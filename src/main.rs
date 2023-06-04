#![feature(test)]
use std::env;
extern crate test;
fn main() {
    let args: Vec<String> = env::args().collect();
    let algorithm = &args[1];
    let seq1 = &args[2];
    let seq2 = &args[3];
    if algorithm == "Needleman-Wunsch" {
        algorithms_lib::needleman_wunsch::align(seq1.to_string(), seq2.to_string());
    } else if algorithm == "Smith-Waterman" {
        algorithms_lib::smith_waterman::align(seq1.to_string(), seq2.to_string());
    } else {
        println!("Algorithm is unknown. Please run program again!");
    }
}


// levenshtein function for reference 
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

// Function for calculating distance based on NW-aligned sequences (returns score only)
// Note to self: This is sort of mimicking what is currently Emily's "align_no_print" fn
// except that it is going to address the (sometimes inverse) relationship between
// score and edit distance
pub fn needleman_wunsch<T: Number, U: Number>(x: &[T], y: &[T]) -> U {
    // Get sequence lengths
    let len_x = x.iter().count();
    let len_y = y.iter().count();

    // Create the grid
    // TODO: Edit create_grid function to work with references to slices of type T
    // intstead of strings
    let (grid, mut directions) = grid::create_grid(x, y, len_x, len_y);

    // Build and print alignment
    // TODO: Edit build_best_alignment function to work with references to slices of type T
    // instead of strings; this should be called something along the lines of
    // "needleman_wunsch_with_alignments"
    let (aligned_seq_x, aligned_seq_y) =
        alignment::build_best_alignment(&grid, &mut directions, x, y);

    // TODO: Implement edit distance function; this is, under current scoring scheme, inversely
    // related to score
    let edit_distance = calculate_edit_distance(aligned_seq_x, aligned_seq_y);

    // want to return edit distance ultimately
    edit_distance
}

// Function for calculating best NW-alignment sequences (returns score only)
// I (Morgan) edited this so that the scoring scheme makes score align with edit 
// distance 
pub fn needleman_wunsch_with_alignments<T: Number, U: Number>(x: &[T], y: &[T]) -> () {
    //Get the aligned sequences (just take 0th one) 

    let mut distance = 0;
    let len_x = x.iter().count();
    for i in 0..len_x {
        // Check for match and add 1 to score
        if !(str_aligned_seq1.chars().nth(i) == str_aligned_seq2.chars().nth(i)) {
            score = score + 1;
        }
    }
}

pub fn smith_waterman<T: Number, U: Number>(x: &[T], y: &[T]) -> U {}

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
