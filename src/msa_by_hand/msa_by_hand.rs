use crate::needleman_wunsch::string;

pub fn msa_by_hand(sequences: &Vec<&[u8]>) -> () {
    let pairwise_distances = distance_many_to_many(sequences, sequences);
    let (center_index, _): (usize, i32) = median(pairwise_distances);
}

pub fn median(pairwise_distances: Vec<Vec<i32>>) -> (usize, i32) {
    pairwise_distances
        .into_iter()
        .map(|v| v.into_iter().sum::<i32>())
        .enumerate()
        .min_by(|(_, l), (_, r)| l.partial_cmp(r).unwrap())
        .unwrap()
}

pub fn distance_many_to_many(left: &Vec<&[u8]>, right: &Vec<&[u8]>) -> Vec<Vec<i32>> {
    left.iter()
        .map(|&l| distance_one_to_many(l, right))
        .collect()
}

pub fn distance_one_to_many(left: &[u8], right: &Vec<&[u8]>) -> Vec<i32> {
    right
        .iter()
        .map(|&r| string::needleman_wunsch(left, r))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::median;
    use crate::msa_by_hand::msa_by_hand::distance_many_to_many;

    #[test]
    fn test_distances_many_to_many() {
        let sequences_0: Vec<&[u8]> = [
            "tgggtactga".as_bytes(),
            "attgttcaga".as_bytes(),
            "gaatgcaagg".as_bytes(),
            "ggccgaatag".as_bytes(),
            "agtgggatgt".as_bytes(),
            "gcacacgtgg".as_bytes(),
            "tcgaaacggc".as_bytes(),
            "cagcggaggg".as_bytes(),
            "agaaagcgca".as_bytes(),
            "gtgacctcat".as_bytes(),
            "ggatgtcaca".as_bytes(),
            "ctttgttacc".as_bytes(),
            "gcagagggcc".as_bytes(),
            "ccttgtgtac".as_bytes(),
            "cttgcctgtt".as_bytes(),
        ]
        .to_vec();

        assert_eq!(
            distance_many_to_many(&sequences_0.clone(), &sequences_0),
            [
                [0, 5, 7, 7, 5, 6, 5, 7, 7, 6, 6, 7, 8, 7, 7],
                [5, 0, 6, 8, 6, 9, 7, 8, 7, 7, 5, 5, 9, 6, 6],
                [7, 6, 0, 6, 6, 5, 7, 5, 6, 8, 5, 7, 6, 8, 7],
                [7, 8, 6, 0, 6, 6, 7, 6, 8, 6, 6, 8, 7, 7, 8],
                [5, 6, 6, 6, 0, 8, 8, 5, 7, 6, 7, 8, 7, 7, 6],
                [6, 9, 5, 6, 8, 0, 5, 6, 7, 7, 8, 9, 5, 7, 7],
                [5, 7, 7, 7, 8, 5, 0, 7, 5, 8, 8, 7, 5, 8, 8],
                [7, 8, 5, 6, 5, 6, 7, 0, 7, 9, 9, 8, 5, 8, 7],
                [7, 7, 6, 8, 7, 7, 5, 7, 0, 6, 5, 9, 5, 8, 9],
                [6, 7, 8, 6, 6, 7, 8, 9, 6, 0, 5, 8, 7, 8, 5],
                [6, 5, 5, 6, 7, 8, 8, 9, 5, 5, 0, 5, 6, 6, 8],
                [7, 5, 7, 8, 8, 9, 7, 8, 9, 8, 5, 0, 7, 3, 6],
                [8, 9, 6, 7, 7, 5, 5, 5, 5, 7, 6, 7, 0, 7, 8],
                [7, 6, 8, 7, 7, 7, 8, 8, 8, 8, 6, 3, 7, 0, 5],
                [7, 6, 7, 8, 6, 7, 8, 7, 9, 5, 8, 6, 8, 5, 0]
            ]
        );

        let sequences_00 = [
            "tgggtactga".as_bytes(),
            "attgttcaga".as_bytes(),
            "agaaagcgca".as_bytes(),
            "gtgacctcat".as_bytes(),
            "ggatgtcaca".as_bytes(),
            "ctttgttacc".as_bytes(),
            "ccttgtgtac".as_bytes(),
            "cttgcctgtt".as_bytes(),
        ]
        .to_vec();

        assert_eq!(
            distance_many_to_many(&sequences_00.clone(), &sequences_00),
            [
                [0, 5, 7, 6, 6, 7, 7, 7],
                [5, 0, 7, 7, 5, 5, 6, 6],
                [7, 7, 0, 6, 5, 9, 8, 9],
                [6, 7, 6, 0, 5, 8, 8, 5],
                [6, 5, 5, 5, 0, 5, 6, 8],
                [7, 5, 9, 8, 5, 0, 3, 6],
                [7, 6, 8, 8, 6, 3, 0, 5],
                [7, 6, 9, 5, 8, 6, 5, 0]
            ]
        );

        let sequences_01 = [
            "gaatgcaagg".as_bytes(),
            "ggccgaatag".as_bytes(),
            "agtgggatgt".as_bytes(),
            "gcacacgtgg".as_bytes(),
            "tcgaaacggc".as_bytes(),
            "cagcggaggg".as_bytes(),
            "gcagagggcc".as_bytes(),
        ]
        .to_vec();

        assert_eq!(
            distance_many_to_many(&sequences_01.clone(), &sequences_01),
            [
                [0, 6, 6, 5, 7, 5, 6],
                [6, 0, 6, 6, 7, 6, 7],
                [6, 6, 0, 8, 8, 5, 7],
                [5, 6, 8, 0, 5, 6, 5],
                [7, 7, 8, 5, 0, 7, 5],
                [5, 6, 5, 6, 7, 0, 5],
                [6, 7, 7, 5, 5, 5, 0]
            ]
        );
    }

    #[test]
    fn test_median() {
        let sequences_0: Vec<&[u8]> = [
            "tgggtactga".as_bytes(),
            "attgttcaga".as_bytes(),
            "gaatgcaagg".as_bytes(),
            "ggccgaatag".as_bytes(),
            "agtgggatgt".as_bytes(),
            "gcacacgtgg".as_bytes(),
            "tcgaaacggc".as_bytes(),
            "cagcggaggg".as_bytes(),
            "agaaagcgca".as_bytes(),
            "gtgacctcat".as_bytes(),
            "ggatgtcaca".as_bytes(),
            "ctttgttacc".as_bytes(),
            "gcagagggcc".as_bytes(),
            "ccttgtgtac".as_bytes(),
            "cttgcctgtt".as_bytes(),
        ]
        .to_vec();

        let pairwise_distances = distance_many_to_many(&sequences_0.clone(), &sequences_0);

        let (center_index, center_distance_sum) = median(pairwise_distances);
        let center_sequence: &str = std::str::from_utf8(&sequences_0[center_index]).unwrap();

        assert_eq!(center_index, 2);
        assert_eq!(center_distance_sum, 89);
        assert_eq!(center_sequence, "gaatgcaagg");

        let sequences_00 = [
            "tgggtactga".as_bytes(),
            "attgttcaga".as_bytes(),
            "agaaagcgca".as_bytes(),
            "gtgacctcat".as_bytes(),
            "ggatgtcaca".as_bytes(),
            "ctttgttacc".as_bytes(),
            "ccttgtgtac".as_bytes(),
            "cttgcctgtt".as_bytes(),
        ]
        .to_vec();

        let pairwise_distances = distance_many_to_many(&sequences_00.clone(), &sequences_00);

        let (center_index, center_distance_sum) = median(pairwise_distances);
        let center_sequence: &str = std::str::from_utf8(&sequences_00[center_index]).unwrap();

        assert_eq!(center_index, 4);
        assert_eq!(center_distance_sum, 40);
        assert_eq!(center_sequence, "ggatgtcaca");

        let sequences_01 = [
            "gaatgcaagg".as_bytes(),
            "ggccgaatag".as_bytes(),
            "agtgggatgt".as_bytes(),
            "gcacacgtgg".as_bytes(),
            "tcgaaacggc".as_bytes(),
            "cagcggaggg".as_bytes(),
            "gcagagggcc".as_bytes(),
        ]
        .to_vec();

        let pairwise_distances = distance_many_to_many(&sequences_01.clone(), &sequences_01);

        let (center_index, center_distance_sum) = median(pairwise_distances);
        let center_sequence: &str = std::str::from_utf8(&sequences_01[center_index]).unwrap();

        assert_eq!(center_index, 5);
        assert_eq!(center_distance_sum, 34);
        assert_eq!(center_sequence, "cagcggaggg");
    }
}
