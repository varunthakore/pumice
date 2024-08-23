use crate::merkle::hash::Hasher;
use crate::packer_hasher::PackerHasher;
use crate::{CommitmentSchemeProver, CommitmentSchemeVerifier};
use ark_ff::PrimeField;
use channel::fs_prover_channel::FSProverChannel;
use channel::fs_verifier_channel::FSVerifierChannel;
use channel::ProverChannel;
use channel::VerifierChannel;
use randomness::Prng;
use sha3::Digest;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

// These closures are given as an input parameter to packaging commitment scheme prover and verifier
// (correspondingly) to enable creation of inner_commitment_scheme after creating the packer.
type PackagingCommitmentSchemeProverFactory =
    Box<dyn FnOnce(usize) -> Box<dyn CommitmentSchemeProver>>;
type PackagingCommitmentSchemeVerifierFactory =
    Box<dyn FnOnce(usize) -> Box<dyn CommitmentSchemeVerifier>>;

/// Prover of Packaging Commitment Scheme.
pub struct PackagingCommitmentSchemeProver<F: PrimeField, H: Hasher<F>, P: Prng, W: Digest> {
    size_of_element: usize,
    n_elements_in_segment: usize,
    n_segments: usize,
    channel: FSProverChannel<F, P, W>,
    packer: PackerHasher<F, H>,
    inner_commitment_scheme: Box<dyn CommitmentSchemeProver>,
    is_merkle_layer: bool,
    queries: BTreeSet<usize>,
    missing_element_queries: Vec<usize>,
    n_missing_elements_for_inner_layer: usize,
}

#[allow(dead_code)]
impl<F: PrimeField, H: Hasher<F, Output = Vec<u8>>, P: Prng, W: Digest>
    PackagingCommitmentSchemeProver<F, H, P, W>
{
    pub fn new(
        size_of_element: usize,
        n_elements_in_segment: usize,
        n_segments: usize,
        channel: FSProverChannel<F, P, W>,
        inner_commitment_scheme_factory: PackagingCommitmentSchemeProverFactory,
        is_merkle_layer: bool,
    ) -> Self {
        let packer = PackerHasher::new(size_of_element, n_segments * n_elements_in_segment);
        let inner_commitment_scheme = inner_commitment_scheme_factory(packer.n_packages);

        if is_merkle_layer {
            assert!(packer.n_elements_in_package == 2);
        }

        Self {
            size_of_element,
            n_elements_in_segment,
            n_segments,
            channel,
            packer,
            inner_commitment_scheme,
            is_merkle_layer,
            queries: BTreeSet::new(),
            missing_element_queries: vec![],
            n_missing_elements_for_inner_layer: 0,
        }
    }

    pub fn new_with_existing(
        size_of_element: usize,
        n_elements_in_segment: usize,
        n_segments: usize,
        channel: FSProverChannel<F, P, W>,
        inner_commitment_scheme: Box<dyn CommitmentSchemeProver>,
    ) -> Self {
        let commitment_scheme = Self::new(
            size_of_element,
            n_elements_in_segment,
            n_segments,
            channel,
            Box::new(move |_: usize| inner_commitment_scheme),
            true,
        );

        assert_eq!(
            2 * commitment_scheme
                .inner_commitment_scheme
                .segment_length_in_elements(),
            n_elements_in_segment
        );

        commitment_scheme
    }

    fn get_num_of_packages(&self) -> usize {
        self.packer.n_packages
    }

    fn get_is_merkle_layer(&self) -> bool {
        self.is_merkle_layer
    }
}

impl<F: PrimeField, H: Hasher<F, Output = Vec<u8>>, P: Prng, W: Digest> CommitmentSchemeProver
    for PackagingCommitmentSchemeProver<F, H, P, W>
{
    fn element_length_in_bytes(&self) -> usize {
        self.size_of_element
    }

    fn num_segments(&self) -> usize {
        self.n_segments
    }

    fn segment_length_in_elements(&self) -> usize {
        self.n_elements_in_segment
    }

    fn add_segment_for_commitment(&mut self, segment_data: &[u8], segment_index: usize) {
        assert_eq!(
            segment_data.len(),
            self.n_elements_in_segment * self.size_of_element
        );
        assert!(segment_index < self.num_segments());
        let packed = self
            .packer
            .pack_and_hash_internal(segment_data, self.is_merkle_layer);
        self.inner_commitment_scheme
            .add_segment_for_commitment(&packed, segment_index);
    }

    fn commit(&mut self) {
        self.inner_commitment_scheme.commit()
    }

    fn start_decommitment_phase(&mut self, queries: BTreeSet<usize>) -> Vec<usize> {
        self.queries = queries;
        self.missing_element_queries = self
            .packer
            .elements_required_to_compute_hashes(&self.queries);

        let package_queries_to_inner_layer: BTreeSet<usize> = self
            .queries
            .iter()
            .map(|&q| q / self.packer.n_elements_in_package)
            .collect();

        let missing_package_queries_inner_layer = self
            .inner_commitment_scheme
            .start_decommitment_phase(package_queries_to_inner_layer);

        let missing_element_queries_to_inner_layer = self
            .packer
            .get_elements_in_packages(&missing_package_queries_inner_layer);

        self.n_missing_elements_for_inner_layer = missing_element_queries_to_inner_layer.len();

        let mut all_missing_elements = Vec::with_capacity(
            self.missing_element_queries.len() + self.n_missing_elements_for_inner_layer,
        );

        // Add elements from `missing_element_queries` and `missing_element_queries_to_inner_layer`.
        all_missing_elements.extend(self.missing_element_queries.iter().cloned());
        all_missing_elements.extend(missing_element_queries_to_inner_layer);

        all_missing_elements
    }

    fn decommit(&mut self, elements_data: &[u8]) {
        assert_eq!(
            elements_data.len(),
            self.size_of_element
                * (self.missing_element_queries.len() + self.n_missing_elements_for_inner_layer),
        );

        for i in 0..self.missing_element_queries.len() {
            let start = i * self.size_of_element;
            let end = start + self.size_of_element;
            let bytes_to_send = &elements_data[start..end];

            if self.is_merkle_layer {
                let digest = bytes_to_send.to_vec();
                let _ = self.channel.send_decommit_node(digest);
            } else {
                let _ = self.channel.send_data(bytes_to_send);
            }
        }

        let start = self.missing_element_queries.len() * self.size_of_element;
        let end = start + self.n_missing_elements_for_inner_layer * self.size_of_element;
        let data_for_inner_layer = self
            .packer
            .pack_and_hash_internal(&elements_data[start..end], self.is_merkle_layer);

        self.inner_commitment_scheme.decommit(&data_for_inner_layer);
    }

    fn get_proof(&self) -> Vec<u8> {
        self.channel.get_proof()
    }
}

/// Verifier of Packaging Commitment Scheme.
pub struct PackagingCommitmentSchemeVerifier<F: PrimeField, H: Hasher<F>, P: Prng, W: Digest> {
    size_of_element: usize,
    n_elements: usize,
    channel: Rc<RefCell<FSVerifierChannel<F, P, W>>>,
    packer: PackerHasher<F, H>,
    inner_commitment_scheme: Box<dyn CommitmentSchemeVerifier>,
    is_merkle_layer: bool,
}

#[allow(dead_code)]
impl<F: PrimeField, H: Hasher<F, Output = Vec<u8>>, P: Prng, W: Digest>
    PackagingCommitmentSchemeVerifier<F, H, P, W>
{
    /// Constructs a new PackagingCommitmentSchemeVerifier using the commitment scheme factory input.
    ///
    /// # Arguments
    ///
    /// - `size_of_element`: length of element in bytes.
    /// - `n_elements`: number of elements.
    /// - `channel`: Fiat-Shamir verifier channel
    /// - `inner_commitment_scheme_factory`: commitment scheme verifier factory
    /// - `is_merkle_layer`: flag to indicate Merkle layer.
    ///
    /// # Returns
    ///
    /// - `Self`: PackagingCommitmentSchemeVerifier
    pub fn new(
        size_of_element: usize,
        n_elements: usize,
        channel: Rc<RefCell<FSVerifierChannel<F, P, W>>>,
        inner_commitment_scheme_factory: PackagingCommitmentSchemeVerifierFactory,
        is_merkle_layer: bool,
    ) -> Self {
        let packer = PackerHasher::new(size_of_element, n_elements);
        let inner_commitment_scheme = inner_commitment_scheme_factory(packer.n_packages);

        if is_merkle_layer {
            assert_eq!(packer.n_elements_in_package, 2);
        }

        Self {
            size_of_element,
            n_elements,
            channel,
            packer,
            inner_commitment_scheme,
            is_merkle_layer,
        }
    }

    pub fn new_test(
        size_of_element: usize,
        n_elements: usize,
        channel: Rc<RefCell<FSVerifierChannel<F, P, W>>>,
        // inner_commitment_scheme_factory: PackagingCommitmentSchemeVerifierFactory,
        is_merkle_layer: bool,
        packer: PackerHasher<F, H>,
        inner_commitment_scheme: Box<dyn CommitmentSchemeVerifier>,
    ) -> Self {
        // let packer = PackerHasher::new(size_of_element, n_elements);
        // let inner_commitment_scheme = inner_commitment_scheme_factory(packer.n_packages);

        if is_merkle_layer {
            assert_eq!(packer.n_elements_in_package, 2);
        }

        Self {
            size_of_element,
            n_elements,
            channel,
            packer,
            inner_commitment_scheme,
            is_merkle_layer,
        }
    }

    /// Constructs a new PackagingCommitmentSchemeVerifier with the input commitment scheme verifier.
    ///
    /// # Arguments
    ///
    /// - `size_of_element`: length of element in bytes.
    /// - `n_elements`: number of elements.
    /// - `channel`: Fiat-Shamir verifier channel
    /// - `inner_commitment_scheme`: commitment scheme verifier
    /// - `is_merkle_layer`: flag to indicate Merkle layer.
    ///
    /// # Returns
    ///
    /// - `Self`: PackagingCommitmentSchemeVerifier
    pub fn new_with_existing(
        size_of_element: usize,
        n_elements: usize,
        channel: Rc<RefCell<FSVerifierChannel<F, P, W>>>,
        inner_commitment_scheme: Box<dyn CommitmentSchemeVerifier>,
    ) -> Self {
        let commitment_scheme = Self::new(
            size_of_element,
            n_elements,
            channel,
            Box::new(move |_: usize| inner_commitment_scheme),
            true,
        );

        assert_eq!(
            2 * commitment_scheme.inner_commitment_scheme.num_of_elements(),
            n_elements
        );

        commitment_scheme
    }

    fn get_num_of_packages(&self) -> usize {
        self.packer.n_packages
    }

    fn get_is_merkle_layer(&self) -> bool {
        self.is_merkle_layer
    }
}

/// Implement CommitmentSchemeVerifier trait for PackagingCommitmentSchemeVerifier
impl<F: PrimeField, H: Hasher<F, Output = Vec<u8>>, P: Prng, W: Digest> CommitmentSchemeVerifier
    for PackagingCommitmentSchemeVerifier<F, H, P, W>
{
    fn num_of_elements(&self) -> usize {
        self.n_elements
    }

    fn read_commitment(&mut self) -> Result<(), anyhow::Error> {
        self.inner_commitment_scheme.read_commitment()
    }

    fn verify_integrity(&mut self, elements_to_verify: BTreeMap<usize, Vec<u8>>) -> Option<bool> {
        // Get missing elements required to compute hashes
        let keys: BTreeSet<usize> = elements_to_verify.keys().copied().collect();
        let missing_elements_idxs = self.packer.elements_required_to_compute_hashes(&keys);

        let mut full_data_to_verify = elements_to_verify.clone();

        for &missing_element_idx in &missing_elements_idxs {
            if self.is_merkle_layer {
                let mut channel = self.channel.borrow_mut();
                let result_array = channel.recv_decommit_node(H::DIGEST_NUM_BYTES).ok()?;
                full_data_to_verify.insert(missing_element_idx, result_array.to_vec());
            } else {
                let mut channel = self.channel.borrow_mut();
                let data = channel.recv_data(self.size_of_element).ok()?;
                full_data_to_verify.insert(missing_element_idx, data);
            }
        }

        // Convert data to hahses
        let bytes_to_verify = self
            .packer
            .pack_and_hash(&full_data_to_verify, self.is_merkle_layer);

        self.inner_commitment_scheme
            .verify_integrity(bytes_to_verify)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::{make_commitment_scheme_verifier, CommitmentHashes};
    use felt::Felt252;
    use randomness::keccak256::PrngKeccak256;
    use sha3::Sha3_256;

    fn test_packaging_verifier_completeness_with(
        size_of_element: usize,
        n_elements: usize,
        n_verifier_friendly_commitment_layers: usize,
        proof: Vec<u8>,
        elements_to_verify: BTreeMap<usize, Vec<u8>>,
        commitment_hashes: CommitmentHashes,
    ) {
        let channel_prng = PrngKeccak256::new();
        let verifier_channel: Rc<RefCell<FSVerifierChannel<Felt252, PrngKeccak256, Sha3_256>>> =
            Rc::new(RefCell::new(FSVerifierChannel::new(channel_prng, proof)));
        let mut verifier = make_commitment_scheme_verifier(
            size_of_element,
            n_elements,
            verifier_channel,
            n_verifier_friendly_commitment_layers,
            commitment_hashes,
            1,
        );

        let _ = verifier.read_commitment();

        assert!(verifier.verify_integrity(elements_to_verify).unwrap());
    }

    #[test]
    fn test_packaging_completeness() {
        let size_of_element = 1;
        let n_elements = 1;
        let _n_segments = 1;
        let n_verifier_friendly_commitment_layers = 0;
        let _data = vec![218];
        let _queries = vec![0];
        let exp_proof = vec![
            144, 179, 218, 100, 7, 139, 77, 94, 125, 120, 142, 246, 58, 51, 233, 222, 113, 197,
            164, 233, 90, 95, 203, 225, 80, 92, 226, 11, 38, 65, 75, 186,
        ];
        let elements_to_verify = BTreeMap::from([(0, vec![218])]);
        let commitment_hashes = CommitmentHashes::from_single_hash("keccak256".to_string());
        test_packaging_verifier_completeness_with(
            size_of_element,
            n_elements,
            n_verifier_friendly_commitment_layers,
            exp_proof,
            elements_to_verify,
            commitment_hashes,
        );

        // let size_of_element = 32;
        // let n_elements = 4;
        // let _n_segments = 2;
        // let n_verifier_friendly_commitment_layers = 0;
        // let _data: Vec<u8> = vec![
        //     1, 35, 100, 184, 107, 167, 27, 153, 178, 178, 4, 16, 193, 139, 130, 53, 171, 152, 226,
        //     105, 245, 241, 72, 163, 50, 42, 211, 163, 168, 41, 209, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // ];
        // let _queries: Vec<usize> = vec![2, 3];
        // let exp_proof: Vec<u8> = vec![
        //     30, 203, 180, 40, 198, 195, 135, 138, 82, 181, 102, 57, 157, 204, 229, 11, 171, 220,
        //     225, 49, 123, 125, 106, 107, 26, 60, 209, 112, 118, 253, 69, 144, 188, 211, 231, 5,
        //     196, 97, 64, 1, 86, 176, 99, 66, 246, 247, 210, 53, 232, 192, 90, 107, 229, 91, 72, 22,
        //     240, 95, 210, 204, 8, 248, 196, 107,
        // ];
        // let elements_to_verify = BTreeMap::from([
        //     (
        //         2,
        //         vec![
        //             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //             0, 0, 0, 0, 0, 0,
        //         ],
        //     ),
        //     (
        //         3,
        //         vec![
        //             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //             0, 0, 0, 0, 0, 0,
        //         ],
        //     ),
        // ]);
        // let commitment_hashes = CommitmentHashes::from_single_hash("keccak256".to_string());
        // test_packaging_verifier_completeness_with(
        //     size_of_element,
        //     n_elements,
        //     n_verifier_friendly_commitment_layers,
        //     exp_proof,
        //     elements_to_verify,
        //     commitment_hashes,
        // );

        // let size_of_element = 9;
        // let n_elements = 8;
        // let _n_segments = 1;
        // let n_verifier_friendly_commitment_layers = 0;
        // let _data: Vec<u8> = vec![
        //     4, 90, 97, 132, 3, 36, 87, 219, 51, 182, 28, 167, 37, 233, 113, 129, 120, 66, 148, 157,
        //     113, 212, 10, 142, 81, 151, 47, 212, 110, 9, 191, 172, 184, 18, 38, 85, 28, 20, 113,
        //     33, 169, 7, 62, 125, 232, 129, 32, 248, 19, 171, 203, 4, 98, 161, 174, 222, 239, 94,
        //     124, 218, 67, 84, 16, 249, 51, 75, 2, 29, 214, 172, 247, 141,
        // ];
        // let _queries: Vec<usize> = vec![0, 1, 2, 3, 4, 7];
        // let exp_proof: Vec<u8> = vec![
        //     227, 216, 27, 24, 119, 173, 119, 23, 6, 143, 172, 94, 211, 230, 252, 178, 181, 162,
        //     103, 224, 82, 199, 136, 76, 191, 61, 234, 103, 168, 121, 179, 118, 129, 32, 248, 19,
        //     171, 203, 4, 98, 161, 174, 222, 239, 94, 124, 218, 67, 84, 16,
        // ];
        // let elements_to_verify: BTreeMap<usize, Vec<u8>> = BTreeMap::from([
        //     (0, vec![4, 90, 97, 132, 3, 36, 87, 219, 51]),
        //     (1, vec![182, 28, 167, 37, 233, 113, 129, 120, 66]),
        //     (2, vec![148, 157, 113, 212, 10, 142, 81, 151, 47]),
        //     (3, vec![212, 110, 9, 191, 172, 184, 18, 38, 85]),
        //     (4, vec![28, 20, 113, 33, 169, 7, 62, 125, 232]),
        //     (7, vec![249, 51, 75, 2, 29, 214, 172, 247, 141]),
        // ]);
        // let commitment_hashes = CommitmentHashes::from_single_hash("keccak256".to_string());
        // test_packaging_verifier_completeness_with(
        //     size_of_element,
        //     n_elements,
        //     n_verifier_friendly_commitment_layers,
        //     exp_proof,
        //     elements_to_verify,
        //     commitment_hashes,
        // );

        // let size_of_element = 32;
        // let n_elements = 4;
        // let _n_segments = 1;
        // let n_verifier_friendly_commitment_layers = 0;
        // let _data: Vec<u8> = vec![
        //     0, 244, 180, 10, 155, 113, 242, 248, 48, 242, 218, 212, 163, 250, 94, 65, 248, 34, 62,
        //     45, 135, 203, 137, 51, 226, 102, 52, 31, 183, 44, 63, 77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // ];
        // let _queries: Vec<usize> = vec![0, 1, 2, 3];
        // let exp_proof: Vec<u8> = vec![
        //     226, 176, 114, 29, 157, 171, 34, 248, 182, 241, 254, 5, 43, 225, 0, 170, 122, 203, 186,
        //     188, 169, 218, 71, 145, 121, 80, 14, 246, 189, 160, 236, 130,
        // ];
        // let elements_to_verify = BTreeMap::from([
        //     (
        //         0,
        //         vec![
        //             0, 244, 180, 10, 155, 113, 242, 248, 48, 242, 218, 212, 163, 250, 94, 65, 248,
        //             34, 62, 45, 135, 203, 137, 51, 226, 102, 52, 31, 183, 44, 63, 77,
        //         ],
        //     ),
        //     (
        //         1,
        //         vec![
        //             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //             0, 0, 0, 0, 0, 0,
        //         ],
        //     ),
        //     (
        //         2,
        //         vec![
        //             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //             0, 0, 0, 0, 0, 0,
        //         ],
        //     ),
        //     (
        //         3,
        //         vec![
        //             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //             0, 0, 0, 0, 0, 0,
        //         ],
        //     ),
        // ]);
        // let commitment_hashes = CommitmentHashes::from_single_hash("keccak256".to_string());
        // test_packaging_verifier_completeness_with(
        //     size_of_element,
        //     n_elements,
        //     n_verifier_friendly_commitment_layers,
        //     exp_proof,
        //     elements_to_verify,
        //     commitment_hashes,
        // );

        // let size_of_element = 32;
        // let n_elements = 4;
        // let _n_segments = 2;
        // let n_verifier_friendly_commitment_layers = 0;
        // let _data: Vec<u8> = vec![
        //     0, 69, 9, 83, 116, 186, 114, 166, 162, 161, 5, 164, 38, 245, 112, 70, 202, 33, 63, 186,
        //     85, 163, 0, 139, 112, 108, 14, 98, 197, 219, 225, 74, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // ];
        // let _queries: Vec<usize> = vec![0];
        // let exp_proof: Vec<u8> = vec![
        //     146, 94, 184, 144, 65, 226, 66, 73, 78, 28, 77, 148, 37, 124, 200, 171, 209, 97, 100,
        //     196, 17, 231, 9, 39, 122, 53, 63, 43, 193, 194, 119, 162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        //     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 173, 50, 40, 182,
        //     118, 247, 211, 205, 66, 132, 165, 68, 63, 23, 241, 150, 43, 54, 228, 145, 179, 10, 64,
        //     178, 64, 88, 73, 229, 151, 186, 95, 181,
        // ];
        // let elements_to_verify = BTreeMap::from([(
        //     0,
        //     vec![
        //         0, 69, 9, 83, 116, 186, 114, 166, 162, 161, 5, 164, 38, 245, 112, 70, 202, 33, 63,
        //         186, 85, 163, 0, 139, 112, 108, 14, 98, 197, 219, 225, 74,
        //     ],
        // )]);
        // let commitment_hashes = CommitmentHashes::from_single_hash("keccak256".to_string());
        // test_packaging_verifier_completeness_with(
        //     size_of_element,
        //     n_elements,
        //     n_verifier_friendly_commitment_layers,
        //     exp_proof,
        //     elements_to_verify,
        //     commitment_hashes,
        // );
    }
}