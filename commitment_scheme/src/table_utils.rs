use std::cmp::Ordering;

/// Represents a cell in a 2-dimensional array (row and column), with lexicographic comparison.
#[derive(PartialEq, Eq, Clone, Hash)]
pub struct RowCol {
    row: usize,
    col: usize,
}

impl RowCol {
    pub fn new(row: usize, col: usize) -> Self {
        Self { row, col }
    }

    pub fn get_row(&self) -> usize {
        self.row
    }

    pub fn get_col(&self) -> usize {
        self.col
    }

    fn as_pair(&self) -> (usize, usize) {
        (self.row, self.col)
    }
}

impl PartialOrd for RowCol {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.as_pair().cmp(&other.as_pair()))
    }
}

impl Ord for RowCol {
    fn cmp(&self, other: &Self) -> Ordering {
        self.as_pair().cmp(&other.as_pair())
    }
}

pub fn all_query_rows(data_queries: &[RowCol], integrity_queries: &[RowCol]) -> Vec<usize> {
    let mut all_query_rows = Vec::new();
    for query in data_queries.iter() {
        all_query_rows.push(query.get_row());
    }
    for query in integrity_queries.iter() {
        all_query_rows.push(query.get_row());
    }
    all_query_rows
}

pub fn elements_to_be_transmitted(
    n_columns: usize,
    all_query_rows: &[usize],
    integrity_queries: &[RowCol],
) -> Vec<RowCol> {
    let mut to_be_transmitted = Vec::new();

    for &row in all_query_rows.iter() {
        for col in 0..n_columns {
            let query_loc = RowCol::new(row, col);
            // Add the location (row, col) only if it is not part of integrity_queries.
            if !integrity_queries.contains(&query_loc) {
                to_be_transmitted.push(query_loc);
            }
        }
    }

    to_be_transmitted
}