use std::ops::RangeBounds;

use clam::number::Number;


#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
    None,
}

// New function to compute the distance and direction table using scoring scheme 0; 1; 1
// Table represents adjacecy matix
pub fn compute_sw_table<T: Number>(x: &[T], y: &[T]) -> (Vec<Vec<usize>>, Vec<(usize, usize)>) {
    let len_x = x.len();
    let len_y = y.len();

    // Initializing table; subvecs represent rows
    let mut table = vec![vec![0; len_x]; len_y];
    let mut matches:Vec<(usize, usize)> = Vec::new();
    // Initialize top row and left column distance values
    #[allow(clippy::needless_range_loop)]

    // Set values for the body of the table
    for row in 0..len_y {
        // Check if it's a duplicate row
        let index = y.iter().position(|&c| c == y[row]).unwrap();
        for col in 0..len_x {
            // Populate table   
            table[row][col] = if x[col] == y[row] { 0 } else { 1 };
            // Store matches
            if table[row][col] == 0 {
                matches.push((col, row));
            }
        }
    }

    (table, matches)
}

pub enum Edit<T: Number> {
    Deletion(usize),
    Insertion(usize, T),
    Substitution(usize, T),
}

// Given an alignement of two sequences, returns the set of edits needed to turn sequence
// x into sequence y
pub fn alignment_to_edits<T: Number>(aligned_x: &[T], aligned_y: &[T]) -> Vec<Edit<T>> {
    aligned_x
        .iter()
        .zip(aligned_y.iter())
        .filter(|(x, y)| x != y)
        .enumerate()
        .map(|(index, (x, y))| {
            if *x == T::from(b'-').unwrap() {
                Edit::Insertion(index, *y)
            } else if *y == T::from(b'-').unwrap() {
                Edit::Deletion(index)
            } else {
                Edit::Substitution(index, *y)
            }
        })
        .collect()
}

// Given a sequence in alignment, determine the indices at which gaps should be inserted
pub fn get_gap_indices<T: Number>(aligned_x: &[T]) -> Vec<usize> {
    aligned_x
        .iter()
        .enumerate()
        .filter(|(_, &value)| value == T::from(b'-').unwrap())
        .map(|(index, _)| index)
        .collect()
}

fn add_third(matches: Vec<(usize, usize)>, curr_position: (usize, usize), consect_gaps: u32, consect_matches: u32, index: usize) -> u8 {
    // The third position should be added if there's only one gap in 2nd ans 3rd and in 3rd and 4th
    // Or if len > 5 and first 2 are 1
    if ((matches.iter().position(|&r| r == (curr_position.0 - index + 2, curr_position.1 - index + 2)) != None) && 
    (matches.iter().position(|&r| r == (curr_position.0 - index + 1, curr_position.1 - index + 1)) == None) &&
    (matches.iter().position(|&r| r == (curr_position.0 - index + 3, curr_position.1 - index + 3)) == None)) ||

    ((matches.iter().position(|&r| r == (curr_position.0 - index + 2, curr_position.1 - index + 2)) == None) &&
    (matches.iter().position(|&r| r == (curr_position.0 - index + 1, curr_position.1 - index + 1)) != None) &&
    (matches.iter().position(|&r| r == (curr_position.0 - index + 3, curr_position.1 - index + 3)) != None))
    {
        return 1;
    }
    return 0;
}

fn add_fourth_fifth(matches: (Vec<usize>, Vec<usize>), curr_position: (usize, usize), consect_gaps: u32, consect_matches: u32, mut index: usize) -> (u8, u8, usize) {
    let mut fourth:u8 = 0;
    let mut fifth:u8 = 0;
    println!("{} ", consect_matches);
    if consect_gaps + consect_matches > 5 || consect_matches == 5
    || (!((matches.0.iter().position(|&r| r == curr_position.0 - index + 2) == matches.1.iter().position(|&r| r == curr_position.1 - index + 2)) && 
    ((matches.0.iter().position(|&r| r == curr_position.0 - index + 3) == matches.1.iter().position(|&r| r == curr_position.1 - index + 3)))) &&
    !((matches.0.iter().position(|&r| r == curr_position.0 - index + 3) == matches.1.iter().position(|&r| r == curr_position.1 - index + 3)) && 
    ((matches.0.iter().position(|&r| r == curr_position.0 - index + 4) == matches.1.iter().position(|&r| r == curr_position.1 - index + 4))))) {
        fourth = 1;
        println!("whyyy {}", curr_position.0 - index + 4);
        if consect_gaps + consect_matches > 5 || matches.0.iter().position(|&r| r == curr_position.0 - index + 4) == matches.1.iter().position(|&r| r == curr_position.1 - index + 4){
            fifth = 1;
            index = index + 1;
        }
    }
    println!("{} {}", fourth, fifth);
    return (fourth, fifth, index);
}

fn cluster_found(alignment_clusters: Vec<Vec<(usize, usize)>>, index: (usize, usize)) -> u8 {
    let mut found:u8 = 0;
    for i in 0..alignment_clusters.len() {
        for j in 0..alignment_clusters[i].len() {
            if alignment_clusters[i][j] == index {
                found = 1;
            }
        }
    }
    return found;
}

// Find all the clusters that represent alignments
fn detect_clusters(matches: (Vec<usize>, Vec<usize>), xlen: usize, ylen: usize) -> Vec<Vec<(usize, usize)>> {
    let mut alignment_clusters: Vec<Vec<(usize, usize)>> = Vec::new();
    for i in 0..(matches.0.len() - 3){
        println!("Start: {} {}", matches.0[i], matches.1[i]);
        if cluster_found(alignment_clusters.clone(), (matches.0[i], matches.1[i])) == 0 {
            println!("Not found");
            let (mut diagonals, mut j) = (1, 1);
            let (mut consect_gaps, mut old_matches, mut consect_matches) = (0, 0, 0);
            let mut curr_position = (matches.0[i], matches.1[i]);
            // Alignment criteria:
            // if consect_matches >= consect_gaps or old_matches >= consect_gaps
            // Cannot go outside the grid
            while (consect_matches >= consect_gaps || old_matches >= consect_gaps) 
            && (Some(&curr_position.0) < matches.0.iter().max() && Some(&curr_position.1) < matches.1.iter().max())  {
                println!("helllooo");
                if matches.0.iter().position(|&r| r == curr_position.0) != None {
                    println!("not none");
                    if matches.0.iter().position(|&r| r == curr_position.0) == Some(23){
                        println!("yay");
                    }
                }
                if matches.0.iter().position(|&r| r == curr_position.0) != None 
                && matches.0.iter().position(|&r| r == curr_position.0) == matches.1.iter().position(|&r| r == curr_position.1) {
                    diagonals = diagonals + 1;
                    consect_matches = consect_matches + 1;
                    println!("Making alignment")
                }
                else {
                    consect_gaps = consect_gaps + 1;
                }
                if consect_gaps + consect_matches == 5 && consect_matches >= 3 || (consect_matches >= 3 && curr_position.0 == xlen - 1 || curr_position.1 == ylen - 1){
                    // Add to cluster if there's more matches
                    let mut ending = 5;
                    while  matches.0.iter().position(|&r| r == curr_position.0 + ending) == matches.1.iter().position(|&r| r == curr_position.0 + ending) {
                        consect_matches = consect_matches + 1;
                        ending = ending + 1;
                    }
                    println!("matches = {}", consect_matches);
                    j = j + ending - 2;
                    let mut curr_alignment:Vec<(usize, usize)> = Vec::new();
                    // The first and second position of cluster should only be added if there's not 2 gaps in a row right after
                    if (matches.0.iter().position(|&r| r == curr_position.0 - j + 1) == matches.1.iter().position(|&r| r == curr_position.1 - j + 1)) || (matches.0.iter().position(|&r| r == curr_position.0 - j + 2) == matches.1.iter().position(|&r| r == curr_position.1 - j + 2)){
                        curr_alignment.push((curr_position.0 - j, curr_position.1 - j));
                        curr_alignment.push((curr_position.0 - j + 1, curr_position.1 - j + 1));
                    }
                    if add_third(matches.clone(), curr_position, consect_gaps, consect_matches, j) == 1 {
                        curr_alignment.push((curr_position.0 - j + 2, curr_position.1 - j + 2));
                    }
                    // The fourth position should be added: always if over 5 length, if not both 3&4 or 4&5 gaps
                    // If cluster is over 5 characters, then it will still have 2 gaps/mismatches
                    // The fifth should be added as long as len not 5 and 5 is 0
                    let (get_fourth, get_fifth, x) = add_fourth_fifth(matches.clone(), curr_position, consect_gaps, consect_matches, j);
                    if get_fourth == 1 {
                        curr_alignment.push((curr_position.0 - j + 3, curr_position.1 - j + 3));
                    }
                    if get_fifth == 1 {
                        curr_alignment.push((curr_position.0 - j + 4, curr_position.1 - j + 4));
                    }
                    alignment_clusters.push(curr_alignment);
                    println!("Alignment added");
                }
                else {
                    curr_position = (curr_position.0 + 1, curr_position.1 + 1);
                }
            }
        }
    }
    return alignment_clusters;
}
// Iterative version of traceback function so we can benchmark both this and recursive option
pub fn traceback_iterative<T: Number>(){

}
#[cfg(test)]
mod tests {
    use crate::alignment_helpers::compute_sw_table;
    //use crate::alignment_helpers::traceback_iterative;
    //use crate::alignment_helpers::traceback_recursive;
    use crate::alignment_helpers::Direction;
    use crate::alignment_helpers::detect_clusters;
    #[test]
    fn test_compute_table() {
        let x_u8 = "NAJIBPEPPERSEATS".as_bytes();
        let y_u8 = "NAJIBEATSPEPPERS".as_bytes();
        let (table, matches) = compute_sw_table(x_u8, y_u8);
        assert_eq!(
            table, vec![
            vec![0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
            vec![1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],  
            vec![1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1],
            vec![1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
            vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
            vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0],
            vec![1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0]]
        );
        assert_eq!(matches, (vec![0, 1, 13, 2, 3, 4, 6, 9, 12, 1, 13, 14, 11,
            15, 5, 7, 8, 6, 9, 12, 5, 7, 8, 5, 7, 8, 6, 9, 12, 10, 11, 15], vec![0,
            1, 1, 2, 3, 4, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 10, 10,
            11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 15, 15]));
        let clusters = detect_clusters(matches, x_u8.len(), y_u8.len());
        println!("{}", clusters.len());
        for i in 0..clusters.len() - 1{
            for j in 0..clusters[0].len() {
                print!("({}, {})\t", clusters[i][j].0, clusters[i][j].1);
            }
            println!("");
        }
        //assert_eq!(clusters, vec![vec![]])
    }
}
