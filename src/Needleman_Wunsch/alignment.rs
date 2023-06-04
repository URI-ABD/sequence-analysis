use super::grid::Direction;

// get sequence characters
fn get_seq_char<'a>(cell: i32, seq1: &'a String, seq2: &'a String) -> (char, usize, char, usize) {
    let seq1_char_index = ((cell - (cell % (seq2.len() as i32 + 1))) / (seq2.len() as i32 + 1)) - 1;
    let seq2_char_index = (cell % (seq2.len() as i32 + 1)) - 1;
    let seq1_char: char;
    let seq2_char: char;
    if seq1_char_index < 0 {
        seq1_char = '-';
    } else {
        seq1_char = seq1.chars().nth(seq1_char_index as usize).unwrap();
    }
    if seq2_char_index < 0 {
        seq2_char = '-';
    } else {
        seq2_char = seq2.chars().nth(seq2_char_index as usize).unwrap();
    }
    return (
        seq1_char,
        seq1_char_index as usize,
        seq2_char,
        seq2_char_index as usize,
    );
}

fn first_alignment<'a>(
    score_grid: &'a Vec<i32>,
    directions: &'a mut Vec<Direction>,
    seq1: String,
    seq2: String,
) -> (Vec<char>, Vec<char>, Vec<Direction>) {
    let mut aligned_seq1: Vec<char> = vec![];
    let mut aligned_seq2: Vec<char> = vec![];
    // Save directions of alignment
    let mut alignment_directions: Vec<Direction> = vec![];
    let mut cell = score_grid.len() as i32 - 1;
    // Look for directions if cell is not 0
    while cell > 0 as i32 {
        let (seq1_char, seq1_char_index, seq2_char, seq2_char_index) =
            get_seq_char(cell, &seq1, &seq2);
        // If the current direction index includes Diagonal, add the two corresponding characters to the sequence strings
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
            aligned_seq2.insert(0, '-');
            alignment_directions.push(Direction::Up);
            if (cell as i32 - seq2.len() as i32 - 1) as i32 >= 0 && seq1_char_index >= 0 {
                cell = cell - seq2.len() as i32 - 1;
            }
        }
        // If the current direction index includes Left, add the corresponing character to the second sequence and a gap to the second sequence
        else if directions[cell as usize] == Direction::Left {
            aligned_seq1.insert(0, '-');
            aligned_seq2.insert(0, seq2_char);
            alignment_directions.push(Direction::Left);
            // Move to the cell to the diagonally left of the current cell
            if (cell - 1) as i32 >= 0 && seq2_char_index >= 0 {
                cell = cell - 1;
            }
        }
    }
    return (aligned_seq1, aligned_seq2, alignment_directions);
}

//Checks if another path can be found
fn is_new_path<'a>(
    directions: &'a mut Vec<Direction>,
    previous_alignment_directions: Vec<Direction>,
    seq2: String,
) -> (i32, i32) {
    let mut index: i32 = 0;
    let mut i: i32 = previous_alignment_directions.len() as i32 - 1;
    let mut found_new_path: bool = false;
    // Goes through the path of the most recent alignment to see if a new path can be made
    while i >= 0 as i32 && (found_new_path == false) {
        if previous_alignment_directions[i as usize] == Direction::Diagonal
            && directions.len() > index as usize + seq2.len() + 2
        {
            println!("Index: {}", index);
            if directions[index as usize + seq2.len() + 2] != Direction::Diagonal {
                found_new_path = true;
            }
            index = index + seq2.len() as i32 + 2;
        } else if previous_alignment_directions[i as usize] == Direction::Up
            && directions.len() > index as usize + seq2.len() + 1
        {
            if (directions[index as usize + seq2.len() + 1] == Direction::DiagonalUpLeft)
                || (directions[index as usize + seq2.len() + 1] == Direction::UpLeft)
            {
                found_new_path = true;
            }
            index = index + seq2.len() as i32 + 1;
        } else if previous_alignment_directions[i as usize] == Direction::Left
            && directions.len() > index as usize + 1
        {
            index = index + 1;
        }
        if !(found_new_path) {
            i = i - 1;
        }
    }
    // If there's no new paths, it returns 1
    if !(found_new_path) {
        i = -1;
        index = 0;
    }
    // Else returns i and the index of the grid
    return (i, index);
}

// Check if alignment was already found
fn alignment_found<'a>(
    aligned_seq1: &'a mut Vec<Vec<char>>,
    aligned_seq2: &'a mut Vec<Vec<char>>,
    new_aligned_seq1: Vec<char>,
    new_aligned_seq2: Vec<char>,
) -> bool {
    let mut is_found: bool = false;
    let seq1_1: String = new_aligned_seq1.to_vec().into_iter().collect();
    let seq2_1: String = new_aligned_seq2.to_vec().into_iter().collect();
    if (aligned_seq1.contains(&new_aligned_seq1)) && (aligned_seq2.contains(&new_aligned_seq2)) {
        // Find the index of seq1 in the aligned_seq1 list
        let mut seq_index: usize = 0;
        while aligned_seq1[seq_index] != new_aligned_seq1 {
            seq_index = seq_index + 1;
        }
        // Check if seq2 is in the same location
        if aligned_seq2[seq_index] == new_aligned_seq2 {
            is_found = true;
        }
    }
    return is_found;
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
        // If no alignments have been found, no need to worry which directions have been fully traveled through, can call first alignment
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

// Calculates the score
pub fn score<'a>(str_aligned_seq1: &'a String, str_aligned_seq2: &'a String) -> i32 {
    let mut score: i32 = 0;
    for i in 0..str_aligned_seq1.len() {
        // Check for match and add 1 to score
        if str_aligned_seq1.chars().nth(i) == str_aligned_seq2.chars().nth(i) {
            score = score + 1;
        } else {
            score = score - 1
        }
    }
    return score;
}

// Prints all alignments
pub fn print_alignments<'a>(
    str_aligned_seq1: &'a Vec<String>,
    str_aligned_seq2: &'a Vec<String>,
    score: i32,
) {
    println!("Best Alignments:");
    for i in 0..str_aligned_seq1.len() {
        println!("{}      {}", str_aligned_seq1[i], str_aligned_seq2[i]);
    }
    println!("Score: {}", score)
}
