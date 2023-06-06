use super::grid::Direction;
use clam::core::number::Number;

// Revision of what was originally get_seq_char (now takes index as an argument)
// by taking index as an argument, we remove code that is repeated for both sequence x and y
fn get_next_character<T: Number>(cell: i32, seq: &[T], index: i32) -> T {
    let gap = 0; //TODO: FIGURE THIS OUT
    let seq_char = if index < 0 {
        T::from(gap).unwrap() //TODO: FIGURE THIS OUT
    } else {
        *seq.iter().nth(index as usize).unwrap()
    };
    return seq_char;
}

// Function to compute alignment
// Verify that all lifetimes are necessary
fn first_alignment<'a, T: Number>(
    score_grid: &'a Vec<i32>,
    directions: &'a mut Vec<Direction>,
    seq1: &[T],
    seq2: &[T],
) -> (Vec<T>, Vec<T>, Vec<Direction>) {
    let mut aligned_seq1: Vec<T> = vec![];
    let mut aligned_seq2: Vec<T> = vec![];

    // Save directions of alignment
    let mut alignment_directions: Vec<Direction> = vec![];
    let mut cell = score_grid.len() as i32 - 1;

    // Look for directions if cell is not 0
    while cell > 0 as i32 {
        let seq1_index = ((cell - (cell % (seq2.len() as i32 + 1))) / (seq2.len() as i32 + 1)) - 1;
        let seq2_index = (cell % (seq2.len() as i32 + 1)) - 1;

        let seq1_char = get_next_character(cell, seq1, seq1_index);
        let seq2_char = get_next_character(cell, seq2, seq2_index);

        // If the current direction index includes Diagonal,
        // add the two corresponding characters to the sequence strings
        if (directions[cell as usize] == Direction::Diagonal)
            || (directions[cell as usize] == Direction::DiagonalLeft)
            || (directions[cell as usize] == Direction::DiagonalUp)
            || (directions[cell as usize] == Direction::DiagonalUpLeft)
        {
            aligned_seq1.insert(0, seq1_char);
            aligned_seq2.insert(0, seq2_char);
            alignment_directions.push(Direction::Diagonal);

            // Move to the cell to the diagonally left of the current cell
            if (cell as i32 - seq2.len() as i32 - 2) as i32 >= 0 {
                cell = cell - seq2.len() as i32 - 2;
            }
        }
        // If the current direction index includes Up, add a gap to the first sequence and the corresponding character to the second sequence
        else if (directions[cell as usize] == Direction::Up)
            || (directions[cell as usize] == Direction::UpLeft)
        {
            aligned_seq1.insert(0, seq1_char);
            aligned_seq2.insert(0, '-'); //TODO: fix to use new handling of gaps
            alignment_directions.push(Direction::Up);
            if (cell as i32 - seq2.len() as i32 - 1) as i32 >= 0 && seq1_index >= 0 {
                cell = cell - seq2.len() as i32 - 1;
            }
        }
        // If the current direction index includes Left, add the corresponing character to the second sequence and a gap to the second sequence
        else if directions[cell as usize] == Direction::Left {
            aligned_seq1.insert(0, '-');
            aligned_seq2.insert(0, seq2_char);
            alignment_directions.push(Direction::Left);
            // Move to the cell to the diagonally left of the current cell
            if (cell - 1) as i32 >= 0 && seq2_index >= 0 {
                cell = cell - 1;
            }
        }
    }
    return (aligned_seq1, aligned_seq2, alignment_directions);
}

// Find best alignment
fn priv_best_alignment<'a>(
    score_grid: &'a Vec<i32>,
    directions: &'a mut Vec<Direction>,
    seq1: String,
    seq2: String,
    aligned_seq1: &'a mut Vec<Vec<char>>,
    aligned_seq2: &'a mut Vec<Vec<char>>,
) -> (Vec<String>, Vec<String>) {
    let mut new_path_index: i32 = 0;
    let mut grid_index: i32 = -1;
    let mut prev_aligned_seq1: Vec<char> = vec![];
    let mut prev_aligned_seq2: Vec<char> = vec![];
    let mut prev_aligned_directions: Vec<Direction> = vec![];
    let mut cell: i32 = score_grid.len() as i32 - 1;
    // Check if another sequence can be made
    while new_path_index != -1 {
        // If no alignments have been found, no need to worry which directions have been fully traveled
        // through, can call first alignment
        if aligned_seq1.len() == 0 {
            (
                prev_aligned_seq1,
                prev_aligned_seq2,
                prev_aligned_directions,
            ) = first_alignment(score_grid, directions, seq1.clone(), seq2.clone());
            aligned_seq1.push(prev_aligned_seq1.clone());
            aligned_seq2.push(prev_aligned_seq2.clone());
            (new_path_index, grid_index) =
                is_new_path(directions, prev_aligned_directions.clone(), seq2.clone());
        } else {
            let mut new_aligned_seq1: Vec<char> = vec![];
            let mut new_aligned_seq2: Vec<char> = vec![];
            let mut new_aligned_directions: Vec<Direction> = vec![];
            //copy the beggining of recent alignment into the new alignment
            let mut aligned_index = aligned_seq1[aligned_seq1.len() - 1].len() - 1;
            prev_aligned_seq1 = aligned_seq1[aligned_seq1.len() - 1].clone();
            prev_aligned_seq2 = aligned_seq2[aligned_seq2.len() - 1].clone();
            for i in 0..new_path_index {
                new_aligned_seq1.insert(0, prev_aligned_seq1[aligned_index]);
                new_aligned_seq2.insert(0, prev_aligned_seq2[aligned_index]);
                new_aligned_directions.push(prev_aligned_directions[i as usize]);
                aligned_index = aligned_index - 1;
            }
            // Add new part of sequence
            let seq1_1: String = new_aligned_seq1.to_vec().into_iter().collect();
            let seq2_1: String = new_aligned_seq2.to_vec().into_iter().collect();
            let (seq1_char, seq1_char_index, seq2_char, seq2_char_index) =
                get_seq_char(grid_index, &seq1, &seq2.clone());
            // Choose the new direction
            if (directions[grid_index as usize] == Direction::DiagonalUp
                || directions[grid_index as usize] == Direction::DiagonalUpLeft)
                && prev_aligned_directions[new_path_index as usize] == Direction::Diagonal
            {
                new_aligned_seq2.insert(0, '-');
                new_aligned_seq1.insert(0, seq1_char);
                new_aligned_directions.push(Direction::Up);
                if (grid_index as i32 - seq2.len() as i32 - 1) as i32 >= 0 && seq1_char_index > 0 {
                    grid_index = grid_index - seq2.len() as i32 - 1;
                }
            } else if (directions[grid_index as usize] == Direction::DiagonalLeft
                && prev_aligned_directions[new_path_index as usize] == Direction::Diagonal)
                || ((directions[grid_index as usize] == Direction::DiagonalUpLeft
                    || directions[grid_index as usize] == Direction::UpLeft)
                    && prev_aligned_directions[new_path_index as usize] == Direction::Up)
            {
                new_aligned_seq1.insert(0, '-');
                new_aligned_seq2.insert(0, seq2_char);
                new_aligned_directions.push(Direction::Left);
                if (grid_index - 1) as i32 >= 0 && seq2_char_index > 0 {
                    grid_index = grid_index - 1;
                }
            }
            // Finish alignment
            while new_aligned_seq1.len() < prev_aligned_seq1.len() {
                let (seq1_char, seq1_char_index, seq2_char, seq2_char_index) =
                    get_seq_char(grid_index, &seq1, &seq2);
                // If the current direction index includes Diagonal, add the two corresponding characters to the sequence strings
                if (directions[grid_index as usize] == Direction::Diagonal)
                    || (directions[grid_index as usize] == Direction::DiagonalLeft)
                    || (directions[grid_index as usize] == Direction::DiagonalUp)
                    || (directions[grid_index as usize] == Direction::DiagonalUpLeft)
                {
                    new_aligned_seq1.insert(0, seq1_char);
                    new_aligned_seq2.insert(0, seq2_char);
                    new_aligned_directions.push(Direction::Diagonal);
                    // Move to the cell to the diagonally left of the current cell
                    if (grid_index as i32 - seq2.len() as i32 - 2) as i32 >= 0 {
                        grid_index = grid_index - seq2.len() as i32 - 2;
                    }
                }
                // If the current direction index includes Up, add a gap to the first sequence and the corresponding character to the second sequence
                else if (directions[grid_index as usize] == Direction::Up)
                    || (directions[grid_index as usize] == Direction::UpLeft)
                {
                    new_aligned_seq1.insert(0, seq1_char);
                    new_aligned_seq2.insert(0, '-');
                    new_aligned_directions.push(Direction::Up);
                    if (grid_index as i32 - seq2.len() as i32 - 1) as i32 >= 0
                        && seq1_char_index > 0
                    {
                        grid_index = grid_index - seq2.len() as i32 - 1;
                    }
                }
                // If the current direction index includes Left, add the corresponing character to the second sequence and a gap to the second sequence
                else if directions[grid_index as usize] == Direction::Left {
                    new_aligned_seq1.insert(0, '-');
                    new_aligned_seq2.insert(0, seq2_char);
                    new_aligned_directions.push(Direction::Left);
                    // Move to the cell to the diagonally left of the current cell
                    if (grid_index - 1) as i32 >= 0 && seq2_char_index > 0 {
                        grid_index = grid_index - 1;
                    }
                }
            }
            //Append new alignment and check for more
            if !(alignment_found(
                aligned_seq1,
                aligned_seq2,
                new_aligned_seq1.clone(),
                new_aligned_seq2.clone(),
            )) {
                aligned_seq1.append(&mut vec![new_aligned_seq1.clone()]);
                aligned_seq2.append(&mut vec![new_aligned_seq2.clone()]);
                let seq1_2: String = new_aligned_seq1.to_vec().into_iter().collect();
                let seq2_2: String = new_aligned_seq2.to_vec().into_iter().collect();
                println!("{}", aligned_seq1.len());
                (new_path_index, grid_index) =
                    is_new_path(directions, new_aligned_directions.clone(), seq2.clone());
            }
            prev_aligned_directions = new_aligned_directions.clone();
        }
    }
    // Base case: Turn finished aligned sequences into strings and return them
    let mut str_aligned_seq1: Vec<String> = vec![];
    let mut str_aligned_seq2: Vec<String> = vec![];
    for i in 0..aligned_seq1.len() {
        //let sequence1: String = ;
        str_aligned_seq1.append(&mut vec![aligned_seq1[i].to_vec().into_iter().collect()]);
        str_aligned_seq2.append(&mut vec![aligned_seq2[i].to_vec().into_iter().collect()]);
    }
    return (str_aligned_seq1, str_aligned_seq2);
}

// Find best alignment
pub fn build_best_alignment<'a>(
    score_grid: &'a Vec<i32>,
    directions: &'a mut Vec<Direction>,
    seq1: String,
    seq2: String,
) -> (Vec<String>, Vec<String>) {
    let mut aligned_sequence1: Vec<Vec<char>> = vec![];
    let mut aligned_sequence2: Vec<Vec<char>> = vec![];
    return priv_best_alignment(
        score_grid,
        directions,
        seq1,
        seq2,
        &mut aligned_sequence1,
        &mut aligned_sequence2,
    );
}

// Calculates the edit distance between two sequences.
// THIS FUNCTION IS REDUNDANT-- you don't need to do this calculation if you
// use the scoring scheme (0; 1; 1); you can just get it straight from the DP table
// at an earlier step and avoid this additional computation.
// I'll take it out later
pub fn calculate_edit_distance<T: Number, U: Number>(
    aligned_seq_x: &[T],
    aligned_seq_y: &[T],
) -> U {
    let mut edit_distance: i32 = 0;
    for i in 0..aligned_seq_x.len() {
        // Check for gap or mismatch and add 1 to edit_distance
        if !(aligned_seq_x.iter().nth(i) == aligned_seq_y.iter().nth(i)) {
            edit_distance = edit_distance + 1;
        }
    }
    U::from(edit_distance).unwrap()
}
