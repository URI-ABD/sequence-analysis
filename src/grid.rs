use clam::core::number::Number;

// New function to compute the distance table using scoring scheme 0; 1; 1
fn compute_distance_table<T: Number>(x: &[T], y: &[T]) -> Vec<Vec<usize>> {
    let len_x = x.len();
    let len_y = y.len();

    // Distance_table is vec of vecs (initially populated with zeros)
    // Subvecs represent rows
    let mut distance_table = vec![vec![0; len_x]; len_y];

    let gap_penalty = 1;

    // Initialize top row and left column distance values
    for row in 0..(len_y + 1) {
        distance_table[row][0] = gap_penalty * row;
    }

    for column in 0..(len_x + 1) {
        distance_table[0][column] = gap_penalty * column;
    }

    // Set values for the body of the table
    for row in 1..(len_y + 1) {
        for col in 1..(len_x + 1) {
            // Check if sequences match at position col-1 in x and row-1 in y
            // Reason for subtraction is that NW considers an artificial gap at the start
            // of each sequence, so the dp table's indices are 1 higher than that of
            // the actual sequences
            let mismatch_penalty = if x[col - 1] == y[row - 1] { 0 } else { 1 };
            let new_cell_value = [
                distance_table[row - 1][col - 1] + mismatch_penalty,
                distance_table[row - 1][col] + gap_penalty,
                distance_table[row][col - 1] + gap_penalty,
            ]
            .into_iter()
            .min()
            .unwrap();

            distance_table[row][col] = new_cell_value;
        }
    }

    distance_table
}

// #[derive(Clone, Copy, Debug, PartialEq)]
// pub enum Direction {
//     Diagonal,
//     Up,
//     Left,
//     DiagonalUp,
//     DiagonalLeft,
//     UpLeft,
//     DiagonalUpLeft,
//     None,
// }

// // Creates grid from user entered sequences
// pub fn create_grid<T: Number>(
//     seq_x: &[T],
//     seq_y: &[T],
//     len1: i32,
//     len2: i32,
// ) -> (Vec<i32>, Vec<Direction>) {
//     let (len_x, len_y) = (seq_x.len(), seq_y.len());
//     // Create new vec and set up row and column labels
//     let mut score_grid: Vec<i32> = Vec::new();
//     let total_cells = (len_x as i32 + 1) * (len_y as i32 + 1);
//     let mut directions: Vec<Direction> = Vec::new();
//     let mut i = 0 as i32;
//     while i < total_cells {
//         (score_grid, directions) =
//             max_cell_score(seq_x, seq_y, &mut score_grid, &i, &mut directions);
//         i = i + 1;
//     }
//     return (score_grid, directions);
// }

// // Returns direction of the best score for a given cell
// fn max_cell_score<'a, T: Number>(
//     seq_x: &[T],
//     seq_y: &[T],
//     score_grid: &mut Vec<i32>,
//     cell: &i32,
//     directions: &mut Vec<Direction>, //
// ) -> (Vec<i32>, Vec<Direction>) {
//     // determine if location is a match
//     if cell == &0 {
//         directions.push(Direction::None);
//         score_grid.push(0);
//     } else if cell <= &(seq_y.len() as i32) {
//         directions.push(Direction::Left);
//         score_grid.push(0 - cell);
//     }
//     // Assign every cell in the first column with Up
//     else if cell % &(seq_y.len() as i32 + 1) == 0 {
//         directions.push(Direction::Up);
//         score_grid.push(0 - (cell / (seq_y.len() as i32 + 1)));
//     } else {
//         let seq1_char_index =
//             ((cell - (cell % ((seq_y.len() + 1) as i32))) / (seq_y.len() + 1) as i32 - 1) as usize;
//         let seq2_char_index = (cell % (seq_y.len() + 1) as i32 - 1) as usize;
//         let seq1_char = seq_x.iter().nth(seq1_char_index).unwrap();
//         let seq2_char = seq_y.iter().nth(seq2_char_index).unwrap();
//         let mut match_point = -1;
//         // Get surrounding scores
//         if seq1_char == seq2_char {
//             // Match scoore
//             match_point = 1;
//         }
//         // Gap score
//         let from_above: i32 = score_grid[(cell - (seq_y.len() + 1) as i32) as usize] - 1;
//         let from_left: i32 = score_grid[(cell - 1) as usize] - 1;
//         // Base score
//         let from_diagonal: i32 =
//             score_grid[(cell - (seq_y.len() as i32 + 1) as i32 - 1) as usize] + match_point;
//         // Find best score
//         // Save best directions for each cell
//         if (from_diagonal > from_left) && (from_diagonal > from_above) {
//             directions.push(Direction::Diagonal);
//             score_grid.push(from_diagonal);
//         } else if (from_above > from_left) && (from_above > from_diagonal) {
//             directions.push(Direction::Up);
//             score_grid.push(from_above);
//         } else if (from_left > from_above) && (from_left > from_diagonal) {
//             directions.push(Direction::Left);
//             score_grid.push(from_left);
//         } else if (from_above > from_left) && (from_above == from_diagonal) {
//             directions.push(Direction::DiagonalUp);
//             score_grid.push(from_diagonal);
//         } else if (from_left > from_above) && (from_left == from_diagonal) {
//             directions.push(Direction::DiagonalLeft);
//             score_grid.push(from_diagonal);
//         } else if (from_above == from_left) && (from_above > from_diagonal) {
//             directions.push(Direction::UpLeft);
//             score_grid.push(from_above);
//         } else {
//             directions.push(Direction::DiagonalUpLeft);
//             score_grid.push(from_diagonal);
//         }
//     }
//     return (score_grid.to_vec(), directions.to_vec());
// }
