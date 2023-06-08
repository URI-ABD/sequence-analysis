use clam::core::number::Number;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
    None,
}

// New function to compute the distance and direction tables using scoring scheme 0; 1; 1
// Need for direction table should become obsolete with this scoring scheme and recursive traceback
// function, but I'm leaving it in for now
fn compute_tables<T: Number>(x: &[T], y: &[T]) -> (Vec<Vec<usize>>, Vec<Vec<Direction>>) {
    let len_x = x.len();
    let len_y = y.len();

    // Initializing tables; subvecs represent rows
    let mut distance_table = vec![vec![0; len_x]; len_y];
    let mut direction_table = vec![vec![Direction::None; len_x]; len_y];

    let gap_penalty = 1;

    // Initialize top row and left column distance values
    for row in 0..(len_y + 1) {
        distance_table[row][0] = gap_penalty * row;
        direction_table[row][0] = Direction::Up;
    }

    for column in 0..(len_x + 1) {
        distance_table[0][column] = gap_penalty * column;
        direction_table[0][column] = Direction::Left;
    }

    direction_table[0][0] = Direction::None;

    // Set values for the body of the table
    for row in 1..(len_y + 1) {
        for col in 1..(len_x + 1) {
            // Check if sequences match at position col-1 in x and row-1 in y
            // Reason for subtraction is that NW considers an artificial gap at the start
            // of each sequence, so the dp tables' indices are 1 higher than that of
            // the actual sequences
            let mismatch_penalty = if x[col - 1] == y[row - 1] { 0 } else { 1 };
            let (new_cell_index, new_cell_value) = [
                distance_table[row - 1][col - 1] + mismatch_penalty,
                distance_table[row - 1][col] + gap_penalty,
                distance_table[row][col - 1] + gap_penalty,
            ]
            .into_iter()
            .enumerate()
            .min_by(|x, y| x.1.cmp(&y.1))
            .unwrap();

            distance_table[row][col] = new_cell_value;
            direction_table[row][col] = match new_cell_index {
                0 => Direction::Diagonal,
                1 => Direction::Up,
                2 => Direction::Left,
            }
        }
    }

    (distance_table, direction_table)
}

// Returns a single alignment (disregards ties for best-scoring alignment)
fn get_alignment_recursive<'a, T: Number>(
    distance_table: Vec<Vec<i32>>,
    row_index: usize,
    column_index: usize,
    unaligned_seqs: (&[T], &[T]),
    aligned_seqs: (&'a [T], &'a [T]),
) -> (&'a [T], &'a [T]) {
    if row_index == 0 && column_index == 0 {
        //base case
        return aligned_seqs;
    } else {
        //general case
       
        let (unaligned_x, unaligned_y) = unaligned_seqs;
        let (aligned_x, aligned_y) = aligned_seqs;

        let gap_penalty = 1;
        let mismatch_penalty = if unaligned_x[column_index - 1] == unaligned_y[row_index - 1] {
            0
        } else {
            1
        };
        let (min_cell_index, min_cell_value) = [
            distance_table[row_index - 1][column_index - 1] + mismatch_penalty,
            distance_table[row_index - 1][column_index] + gap_penalty,
            distance_table[row_index][column_index - 1] + gap_penalty,
        ]
        .into_iter()
        .enumerate()
        .min_by(|x, y| x.1.cmp(&y.1))
        .unwrap();

        let (updated_aligned_x, updated_aligned_y) = match min_cell_index {
            0 => (
                Some(unaligned_x[column_index - 1])
                    .iter()
                    .chain(aligned_x.iter())
                    .collect(),
                Some(unaligned_y[row_index - 1])
                    .iter()
                    .chain(aligned_y.iter())
                    .collect()
            ),
            1 => (
                Some(unaligned_x[column_index - 1])
                    .iter()
                    .chain(aligned_x.iter())
                    .collect(),
                Some(T::from('-' as u8).unwrap())
                    .iter()
                    .chain(aligned_y.iter())
                    .collect()
            ),
            2 => (
                Some(T::from('-' as u8).unwrap())
                    .iter()
                    .chain(aligned_x.iter())
                    .collect(),
                Some(unaligned_y[row_index - 1])
                    .iter()
                    .chain(aligned_y.iter())
                    .collect()
            ),
        };

        get_alignment(
            distance_table,
            row_index - 1,
            column_index - 1,
            (unaligned_x, unaligned_y),
            (updated_aligned_x, updated_aligned_y),
        )
    }
}

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
