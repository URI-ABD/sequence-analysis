#![feature(test)]
use std::env;
extern crate test;
fn main() {
    println!("hello world");
}

// Levenshtein function for reference
pub fn levenshtein<T: Number, U: Number>(x: &[T], y: &[T]) -> U {
    let (len_x, len_y) = (x.iter().count(), y.iter().count());

    if len_x == 0 {
        // handle special case of 0 length
        U::from(len_x).unwrap()
    } else if len_y == 0 {
        // handle special case of 0 length
        U::from(len_y).unwrap()
    } else if len_x < len_y {
        // require len_a < len_b
        levenshtein(y, x)
    } else {
        let len_y = len_y + 1;

        // initialize DP table for string y
        let mut cur: Vec<usize> = (0..len_y).collect();

        // calculate edit distance
        for (i, cx) in x.iter().enumerate() {
            // get first column for this row
            let mut pre = cur[0];
            cur[0] = i + 1;
            for (j, cy) in y.iter().enumerate() {
                let tmp = cur[j + 1];
                cur[j + 1] = std::cmp::min(
                    // deletion
                    tmp + 1,
                    std::cmp::min(
                        // insertion
                        cur[j] + 1,
                        // match or substitution
                        pre + if cx == cy { 0 } else { 1 },
                    ),
                );
                pre = tmp;
            }
        }
        U::from(cur[len_y - 1]).unwrap()
    }
}

// Wrapper function for calculating distance based on NW-aligned sequences (returns edit distance only)
pub fn needleman_wunsch<T: Number, U: Number>(x: &[T], y: &[T]) -> U {
    (_, _, edit_distance) = needleman_wunsch_with_alignments(x, y); 
    edit_distance
}

// Function for calculating calculating distance based on NW-aligned sequences (returns edit 
// distance AND alignment associated with shortest edit distance)
//
// For now, in cases where there exist ties for the shortest edit distance, we only return 
// one alignment 
//
// When the scoring scheme is 0; 1; 1, the bottom right entry in the DP table for NW IS the edit distance
// For now, I have a separate function which will calculate edit distance separately (with the 0; 1; 1)
// scoring scheme, but that is unnecessary. I will ultimately get rid of it to eliminate redundant computation 
// once the DP table calculations are cleaner
pub fn needleman_wunsch_with_alignments<T: Number, U: Number>(x: &[T], y: &[T]) -> (&[T], &[T], U) {
    // Get sequence lengths
    let len_x = x.len();
    let len_y = y.len();

    // Create the grid
    // TODO: Edit create_grid function to work with references to slices of type T
    // intstead of strings
    let (grid, mut directions) = grid::create_grid(x, y, len_x, len_y);

    // Build and print alignment
    // TODO: Edit build_best_alignment function to work with references to slices of type T
    // instead of strings
    let (aligned_seq_x, aligned_seq_y) =
        alignment::build_best_alignment(&grid, &mut directions, x, y);

    // TODO: Implement edit distance function; this is, under current scoring scheme, inversely
    // related to score
    let edit_distance = calculate_edit_distance(aligned_seq_x, aligned_seq_y);

    // want to return edit distance ultimately
    (aligned_seq_x, aligned_seq_y, edit_distance) 
}


#[cfg(test)]
mod tests {
    mod test_nw {
        mod test_nw_5 {
            use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
            use clam::Metric;
            use test::Bencher;

            #[test]
            fn metric_nw_5to5_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
            fn metric_nw_5to10_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to20_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to30_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to40_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to50_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_10 {
            use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
            use clam::Metric;
            #[test]
            fn metric_nw_10to5_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to10_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to20_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to30_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to40_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to50_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_20 {
            use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
            use clam::Metric;
            #[test]
            fn metric_nw_20to5_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to10_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to20_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to30_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to40_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to50_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_30 {
            use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
            use clam::Metric;
            #[test]
            fn metric_nw_30to5_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to10_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to20_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to30_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to40_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to50_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_40 {
            use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
            use clam::Metric;
            #[test]
            fn metric_nw_40to5_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to10_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to20_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to30_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to40_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to50_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    }

    mod test_nw_50 {
        use algorithms_lib::needleman_wunsch::NeedlemanWunsch;
        use clam::Metric;
        #[test]
        fn metric_nw_50to5_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "CAACC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to10_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TTGCTTTGAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to20_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to30_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to40_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to50_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    }

    mod test_sw {
        mod test_sw_5 {
            use algorithms_lib::smith_waterman::SmithWaterman;
            use clam::Metric;

            #[test]
            fn metric_sw_5to5_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_5to10_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_5to20_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_5to30_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_5to40_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_5to50_test() {
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_sw_10 {
            use algorithms_lib::smith_waterman::SmithWaterman;
            use clam::Metric;
            #[test]
            fn metric_sw_10to5_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_10to10_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_10to20_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_10to30_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_10to40_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_10to50_test() {
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_sw_20 {
            use algorithms_lib::smith_waterman::SmithWaterman;
            use clam::Metric;
            #[test]
            fn metric_sw_20to5_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_20to10_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_20to20_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_20to30_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_20to40_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_20to50_test() {
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_sw_30 {
            use algorithms_lib::smith_waterman::SmithWaterman;
            use clam::Metric;
            #[test]
            fn metric_sw_30to5_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_30to10_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_30to20_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_30to30_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_30to40_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_30to50_test() {
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_sw_40 {
            use algorithms_lib::smith_waterman::SmithWaterman;
            use clam::Metric;
            #[test]
            fn metric_sw_40to5_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_40to10_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_40to20_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_40to30_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_40to40_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_sw_40to50_test() {
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String =
                    "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
                let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    }

    mod test_sw_50 {
        use algorithms_lib::smith_waterman::SmithWaterman;
        use clam::Metric;
        #[test]
        fn metric_sw_50to5_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "CAACC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_sw_50to10_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TTGCTTTGAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_sw_50to20_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_sw_50to30_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let __score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_sw_50to40_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_sw_50to50_test() {
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman { seq1, seq2 };
            let _score: i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    }
}
