use abd_clam::cluster;

use crate::needleman_wunsch::{self, string};
use std::{collections::HashMap, vec};

pub struct MSA_Cluster<'a> {
    index: usize,
    name: &'a str,
    center: &'a [u8],
    layer: usize,
    sibling_index: Option<usize>,
}

impl<'a> MSA_Cluster<'a> {
    pub fn new(index: usize, name: &'a str, center: &'a [u8]) -> Self {
        Self {
            index,
            name,
            center,
            layer: name.len(),
            sibling_index: MSA_Cluster::get_sibling_index(index),
        }
    }

    pub fn get_sibling_index(index:usize ) -> Option<usize> {
        if index == 0 {
                None
            }
            else if index % 2 == 0 {
                Some(index -1)
            }
            else {
                Some(index + 1)
            }
    }
}

pub struct MSA_Tree<'a> {
    original_sequences: Vec<&'a [u8]>,
    distance_table: Vec<Vec<i32>>,
    cluster_names: Vec<&'a str>,
    cluster_centers: Vec<&'a [u8]>,
    clusters: Vec<MSA_Cluster<'a>>,
    depth: usize,
}

impl<'a> MSA_Tree<'a> {
    pub fn new(
        original_sequences: Vec<&'a [u8]>,
        cluster_names: Vec<&'a str>,
        cluster_centers: Vec<&'a [u8]>,
    ) -> Self {
        Self {
            original_sequences,
            cluster_names,
            cluster_centers,
            distance_table: Vec::new(),
            clusters: Vec::new(),
            depth: Default::default(),
        }
    }

    fn initialize_clusters(&mut self) -> () {
        self.clusters = self
            .cluster_names
            .iter()
            .zip(self.cluster_centers.iter())
            .enumerate()
            .map(|(index, (&name, &center))| MSA_Cluster::new(index, name, center))
            .collect();


        self.depth = self.clusters.last().unwrap().layer;
    }


    // fn align_layer(&mut self) -> () {

    //     self.clusters.iter_mut().for_each(|mut c| {
    //         if c.layer == self.depth {
    //             c.center = &needleman_wunsch::string::needleman_wunsch_alignment_only::<u8, i32>(c.center, self.clusters[c.sibling_index.unwrap()].center).0;
    //         }
    //     });

    //     self.depth -= 1; 
    // }

    // fn bottom_up_alignment(&mut self) -> () {

    // }
}

pub fn msa_by_hand() -> () {
    let cluster_names = vec![
        "0", "00", "01", "000", "001", "010", "011", "0000", "0001", "0010", "0011", "0100",
        "0101", "0110", "0111", "00010", "00011",
    ];
    let cluster_centers = vec![
        "gaatgcaagg".as_bytes(),
        "ctttgttacc".as_bytes(),
        "ggatgtcaca".as_bytes(),
        "gcacacgtgg".as_bytes(),
        "cagcggaggg".as_bytes(),
        "gtgacctcat".as_bytes(),
        "ctttgttacc".as_bytes(),
        "tgggtactga".as_bytes(),
        "ggatgtcaca".as_bytes(),
        "gcacacgtgg".as_bytes(),
        "gcagagggcc".as_bytes(),
        "ggccgaatag".as_bytes(),
        "cagcggaggg".as_bytes(),
        "attgttcaga".as_bytes(),
        "ctttgttacc".as_bytes(),
    ];

    let sequences = [
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

    let mut msa_tree = MSA_Tree::new(sequences, cluster_names, cluster_centers);
    msa_tree.initialize_clusters();

    // Vector of pairs of centers of sibling clusters based on prebuilt tree
    // in the format (left_child_center, right_child_center)
    let mut sibling_pairs = HashMap::new();
    sibling_pairs.insert("00_01", ("ggatgtcaca".as_bytes(), "cagcggaggg".as_bytes()));
    sibling_pairs.insert(
        "000_001",
        ("ctttgttacc".as_bytes(), "ggatgtcaca".as_bytes()),
    );
    sibling_pairs.insert(
        "010_011",
        ("gcacacgtgg".as_bytes(), "cagcggaggg".as_bytes()),
    );
    sibling_pairs.insert(
        "0000_0001",
        ("gtgacctcat".as_bytes(), "ctttgttacc".as_bytes()),
    );
    sibling_pairs.insert(
        "0010_0011",
        ("tgggtactga".as_bytes(), "ggatgtcaca".as_bytes()),
    );
    sibling_pairs.insert(
        "0100_0101",
        ("gcacacgtgg".as_bytes(), "gcagagggcc".as_bytes()),
    );
    sibling_pairs.insert(
        "0110_0111",
        ("ggccgaatag".as_bytes(), "cagcggaggg".as_bytes()),
    );
    sibling_pairs.insert(
        "00010_00011",
        ("attgttcaga".as_bytes(), "ctttgttacc".as_bytes()),
    );

    let mut alignments = HashMap::new();
    for (key, pair) in sibling_pairs.iter() {
        let (aligned_x, aligned_y) =
            needleman_wunsch::string::needleman_wunsch_alignment_only::<u8, i32>(pair.0, pair.1);

        println!(
            "alignment for {} is {:?}",
            key,
            (
                std::str::from_utf8(&aligned_x).unwrap(),
                std::str::from_utf8(&aligned_y).unwrap()
            )
        );
        alignments.insert(
            key,
            needleman_wunsch::string::needleman_wunsch_alignment_only::<u8, i32>(pair.0, pair.1),
        );
    }
}

pub fn median(pairwise_distances: Vec<Vec<i32>>) -> (usize, i32) {
    pairwise_distances
        .into_iter()
        .map(|v| v.into_iter().sum::<i32>())
        .enumerate()
        .min_by(|(_, l), (_, r)| l.partial_cmp(r).unwrap())
        .unwrap()
}


pub fn anti_median(pairwise_distances: Vec<Vec<i32>>) -> (usize, i32) {
    pairwise_distances
        .into_iter()
        .map(|v| v.into_iter().sum::<i32>())
        .enumerate()
        .max_by(|(_, l), (_, r)| l.partial_cmp(r).unwrap())
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

pub fn get_distances_between_instances(table: &Vec<Vec<i32>>, indices: &[usize], center_index: usize) -> Vec<i32> {
    let mut distances = Vec::new();
    for i in indices {
        distances.push(table[center_index][*i]);

    }  

    distances

}

#[cfg(test)]
mod tests {
    use super::median;
    use crate::msa_by_hand::msa_by_hand::distance_many_to_many;
    use crate::msa_by_hand::msa_by_hand::msa_by_hand;
    use crate::msa_by_hand::msa_by_hand::anti_median;
    use crate::msa_by_hand::msa_by_hand::get_distances_between_instances;

    #[test]
    fn test_distances_many_to_many() {
        let sequences: Vec<&[u8]> = [
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
            distance_many_to_many(&sequences.clone(), &sequences),
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
    }

    #[test]
    fn test_median() {
        let sequences: Vec<&[u8]> = [
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

        let pairwise_distances = distance_many_to_many(&sequences.clone(), &sequences);

        let (center_index, center_distance_sum) = median(pairwise_distances);
        let center_sequence: &str = std::str::from_utf8(&sequences[center_index]).unwrap();

        assert_eq!(center_index, 2);
        assert_eq!(center_distance_sum, 89);
        assert_eq!(center_sequence, "gaatgcaagg");

    }


    #[test]
    fn test_anti_median() {
        // let sequences: Vec<&[u8]> = [
        //     "tgggtactga".as_bytes(),
        //     "attgttcaga".as_bytes(),
        //     "gaatgcaagg".as_bytes(),
        //     "ggccgaatag".as_bytes(),
        //     "agtgggatgt".as_bytes(),
        //     "gcacacgtgg".as_bytes(),
        //     "tcgaaacggc".as_bytes(),
        //     "cagcggaggg".as_bytes(),
        //     "agaaagcgca".as_bytes(),
        //     "gtgacctcat".as_bytes(),
        //     "ggatgtcaca".as_bytes(),
        //     "ctttgttacc".as_bytes(),
        //     "gcagagggcc".as_bytes(),
        //     "ccttgtgtac".as_bytes(),
        //     "cttgcctgtt".as_bytes(),
        // ]
        // .to_vec();

        // let pairwise_distances = distance_many_to_many(&sequences.clone(), &sequences);


        // let (anti_center_index, anti_center_distance_sum) = anti_median(pairwise_distances);
        // let anti_center_sequence: &str = std::str::from_utf8(&sequences[anti_center_index]).unwrap();

        // assert_eq!(anti_center_index, 14);
        // assert_eq!(anti_center_distance_sum, 97);
        // assert_eq!(anti_center_sequence, "cttgcctgtt");

        // let sequences_00 = [
        //     "gaatgcaagg".as_bytes(),
        //     "agtgggatgt".as_bytes(),
        //     "gcacacgtgg".as_bytes(),
        //     "tcgaaacggc".as_bytes(),
        //     "cagcggaggg".as_bytes(),
        //     "agaaagcgca".as_bytes(),
        //     "gtgacctcat".as_bytes(),
        //     "gcagagggcc".as_bytes(),
        // ].to_vec();

        // let pairwise_distances_00 = distance_many_to_many(&sequences_00.clone(), &sequences_00);

        // let (anti_center_index_00, anti_center_distance_sum_00) = anti_median(pairwise_distances_00);
        // let anti_center_sequence_00: &str = std::str::from_utf8(&sequences_00[anti_center_index_00]).unwrap();

        // assert_eq!(anti_center_index_00, 6);
        // assert_eq!(anti_center_distance_sum_00, 51);
        // assert_eq!(anti_center_sequence_00, "gtgacctcat");

        // let sequences_01 = [
        //     "tgggtactga".as_bytes(),
        //     "attgttcaga".as_bytes(),
        //     "ggccgaatag".as_bytes(),
        //     "ggatgtcaca".as_bytes(),
        //     "ctttgttacc".as_bytes(),
        //     "ccttgtgtac".as_bytes(),
        //     "cttgcctgtt".as_bytes(),
        // ].to_vec();

        // let pairwise_distances_01 = distance_many_to_many(&sequences_01.clone(), &sequences_01);
        
        // let (anti_center_index_01, anti_center_distance_sum_01) = anti_median(pairwise_distances_01);
        // let anti_center_sequence_01: &str = std::str::from_utf8(&sequences_01[anti_center_index_01]).unwrap();

        // assert_eq!(anti_center_index_01, 2);
        // assert_eq!(anti_center_distance_sum_01, 44);
        // assert_eq!(anti_center_sequence_01, "ggccgaatag");


        let sequences_000 = [
            "gaatgcaagg".as_bytes(),
            "agtgggatgt".as_bytes(),
            "cagcggaggg".as_bytes(),
            "gcagagggcc".as_bytes(),
        ].to_vec();


        let pairwise_distances_000 = distance_many_to_many(&sequences_000.clone(), &sequences_000); 

        let (anti_center_index_000, anti_center_distance_sum_000) = anti_median(pairwise_distances_000);
        let anti_center_sequence_000: &str = std::str::from_utf8(&sequences_000[anti_center_index_000]).unwrap();

        assert_eq!(anti_center_index_000, 3);
        assert_eq!(anti_center_distance_sum_000, 12);
        assert_eq!(anti_center_sequence_000, "cagcggaggg");

        let sequences_001 = [
            "gcacacgtgg".as_bytes(),
            "tcgaaacggc".as_bytes(),
            "agaaagcgca".as_bytes(),
            "gtgacctcat".as_bytes(),
        ].to_vec();

        let pairwise_distances_001 = distance_many_to_many(&sequences_001.clone(), &sequences_001);

        let (anti_center_index_001, anti_center_distance_sum_001) = anti_median(pairwise_distances_001);
        let anti_center_sequence_001: &str = std::str::from_utf8(&sequences_001[anti_center_index_001]).unwrap();

        assert_eq!(anti_center_index_001, 2);
        assert_eq!(anti_center_distance_sum_001, 12);
        assert_eq!(anti_center_sequence_001, "agaaagcgca");
            


    }

    #[test]
    fn test_msa_by_hand() {
        msa_by_hand()
    }

    #[test]
    fn test_get_distances() {
        let sequences: Vec<&[u8]> = [
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

        let pairwise_distances = distance_many_to_many(&sequences.clone(), &sequences);
    

    // let distances_center = get_distances_between_instances(&pairwise_distances,  &[2, 4, 5, 6, 7, 8, 9, 12], 9);

    // assert_eq!(distances_center, vec![8, 6, 7, 8, 9, 6, 0, 7]);

    // let distances_left_pole = get_distances_between_instances(&pairwise_distances,  &[2, 4, 5, 6, 7, 8, 9, 12], 7);
    
    // assert_eq!(distances_left_pole, vec![5, 5, 6, 7, 0, 7, 9, 5]);

    // let distances_right_pole = get_distances_between_instances(&pairwise_distances,  &[2, 4, 5, 6, 7, 8, 9, 12], 6);
    
    // assert_eq!(distances_right_pole, vec![5, 5, 6, 7, 0, 7, 9, 5]);

    // let distances_center = get_distances_between_instances(&pairwise_distances, &[0, 1, 3, 10, 11, 13, 14], 3);

    // assert_eq!(distances_center, vec![7, 8, 0, 6, 8, 7, 8]);

    // let distances_left_pole = get_distances_between_instances(&pairwise_distances, &[0, 1, 3, 10, 11, 13, 14], 1);

    // assert_eq!(distances_left_pole, vec![5, 0, 8, 5, 5, 6, 6]);

    // let distances_right_pole = get_distances_between_instances(&pairwise_distances, &[0, 1, 3, 10, 11, 13, 14], 13);

    // assert_eq!(distances_right_pole, vec![7, 6, 7, 6, 3, 0, 5]);


}
}
