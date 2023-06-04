#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
    DiagonalUp,
    DiagonalLeft,
    UpLeft,
    DiagonalUpLeft,
    None,
}

// Creates grid from user entered sequences
pub fn create_grid<'a>(
    mut seq1: &'a mut String,
    mut seq2: &'a mut String,
    len1: i32,
    len2: i32,
) -> (Vec<i32>, Vec<Direction>) {
    // Create new vec and set up row and column labels
    let mut score_grid: Vec<i32> = Vec::new();
    let total_cells = (len1 + 1) * (len2 + 1);
    let mut directions: Vec<Direction> = Vec::new();
    let mut i = 0;
    while i < total_cells {
        (score_grid, directions) =
            max_cell_score(&mut seq1, &mut seq2, &mut score_grid, &i, &mut directions);
        i = i + 1;
    }
    return (score_grid, directions);
}

// Returns direction of the best score for a given cell
fn max_cell_score<'a>(
    seq1: &'a mut String,
    seq2: &'a mut String,
    score_grid: &'a mut Vec<i32>,
    cell: &'a i32,
    directions: &'a mut Vec<Direction>,
) -> (Vec<i32>, Vec<Direction>) {
    // determine if location is a match
    if cell == &0 {
        directions.push(Direction::None);
        score_grid.push(0);
    } else if cell <= &(seq2.len() as i32) {
        directions.push(Direction::Left);
        score_grid.push(0 - cell);
    }
    // Assign every cell in the first column with Up
    else if cell % &(seq2.len() as i32 + 1) == 0 {
        directions.push(Direction::Up);
        score_grid.push(0 - (cell / (seq2.len() as i32 + 1)));
    } else {
        let seq1_char_index =
            ((cell - (cell % ((seq2.len() + 1) as i32))) / (seq2.len() + 1) as i32 - 1) as usize;
        let seq2_char_index = (cell % (seq2.len() + 1) as i32 - 1) as usize;
        let seq1_char = seq1.chars().nth(seq1_char_index).unwrap();
        let seq2_char = seq2.chars().nth(seq2_char_index).unwrap();
        let mut match_point = -1;
        // Get surrounding scores
        if seq1_char == seq2_char {
            // Match scoore
            match_point = 1;
        }
        // Gap score
        let from_above: i32 = score_grid[(cell - (seq2.len() + 1) as i32) as usize] - 1;
        let from_left: i32 = score_grid[(cell - 1) as usize] - 1;
        // Base score
        let from_diagonal: i32 =
            score_grid[(cell - (seq2.len() + 1) as i32 - 1) as usize] + match_point;
        // Find best score
        // Save best directions for each cell
        if (from_diagonal > from_left) && (from_diagonal > from_above) {
            directions.push(Direction::Diagonal);
            score_grid.push(from_diagonal);
        } else if (from_above > from_left) && (from_above > from_diagonal) {
            directions.push(Direction::Up);
            score_grid.push(from_above);
        } else if (from_left > from_above) && (from_left > from_diagonal) {
            directions.push(Direction::Left);
            score_grid.push(from_left);
        } else if (from_above > from_left) && (from_above == from_diagonal) {
            directions.push(Direction::DiagonalUp);
            score_grid.push(from_diagonal);
        } else if (from_left > from_above) && (from_left == from_diagonal) {
            directions.push(Direction::DiagonalLeft);
            score_grid.push(from_diagonal);
        } else if (from_above == from_left) && (from_above > from_diagonal) {
            directions.push(Direction::UpLeft);
            score_grid.push(from_above);
        } else {
            directions.push(Direction::DiagonalUpLeft);
            score_grid.push(from_diagonal);
        }
    }
    return (score_grid.to_vec(), directions.to_vec());
}
