use clam::core::number::Number;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
    None,
}

// New function to compute the distance and direction table using scoring scheme 0; 1; 1
pub fn compute_table<T: Number>(x: &[T], y: &[T]) -> Vec<Vec<(usize, Direction)>> {
    let len_x = x.len();
    let len_y = y.len();

    // Initializing table; subvecs represent rows
    let mut table = vec![vec![(0, Direction::None); len_x + 1]; len_y + 1];

    let gap_penalty = 1;

    // Initialize top row and left column distance values
    for row in 0..(len_y + 1) {
        table[row][0] = (gap_penalty * row, Direction::None);
    }

    for column in 0..(len_x + 1) {
        table[0][column] = (gap_penalty * column, Direction::None);
    }

    // Set values for the body of the table
    for row in 1..(len_y + 1) {
        for col in 1..(len_x + 1) {
            // Check if sequences match at position col-1 in x and row-1 in y
            // Reason for subtraction is that NW considers an artificial gap at the start
            // of each sequence, so the dp tables' indices are 1 higher than that of
            // the actual sequences
            let mismatch_penalty = if x[col - 1] == y[row - 1] { 0 } else { 1 };
            let new_cell = [
                table[row - 1][col - 1].0 + mismatch_penalty,
                table[row - 1][col].0 + gap_penalty,
                table[row][col - 1].0 + gap_penalty,
            ]
            .into_iter()
            .zip([Direction::Diagonal, Direction::Up, Direction::Left].into_iter())
            .min_by(|x, y| x.0.cmp(&y.0))
            .unwrap();

            table[row][col] = new_cell;
        }
    }

    table
}
