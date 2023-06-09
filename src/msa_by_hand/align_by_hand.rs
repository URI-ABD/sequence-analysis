use crate::needleman_wunsch;

#[derive(Clone)]
/// Represents a cluster for use in an MSA tree.
pub struct MSACluster<'a> {
    index: usize,
    name: &'a str,
    center: Vec<u8>,
    layer: usize,
    sibling_index: Option<usize>,
    parent_index: Option<usize>,
}

impl<'a> MSACluster<'a> {
    pub fn new(
        index: usize,
        name: &'a str,
        center: Vec<u8>,
        layer: usize,
        parent_index: Option<usize>,
    ) -> Self {
        Self {
            index,
            name,
            center,
            layer,
            parent_index,
            sibling_index: MSACluster::get_sibling_index(index),
        }
    }

    pub fn get_sibling_index(index: usize) -> Option<usize> {
        if index == 0 {
            None
        } else if index % 2 == 0 {
            Some(index - 1)
        } else {
            Some(index + 1)
        }
    }
}

fn get_gap_indices(alignment: Vec<u8>) -> Vec<usize> {
    let mut gap_indices = Vec::new();
    let reversed_alignment = alignment;
    for (i, c) in reversed_alignment.into_iter().enumerate() {
        if c == b'-' {
            gap_indices.push(i);
        }
    }
    gap_indices
}

fn apply_edits(old_center: Vec<u8>, gap_indices: Vec<usize>) -> Vec<u8> {
    let mut new_center = old_center.to_vec();
    let mut gap_indices = gap_indices;
    gap_indices.reverse();
    while !gap_indices.is_empty() {
        let gap_index = gap_indices.pop().unwrap();
        new_center.insert(gap_index, b'-');
    }

    new_center
}

pub struct MSATree<'a> {
    cluster_names: Vec<&'a str>,
    cluster_centers: Vec<Vec<u8>>,
    cluster_layers: Vec<usize>,
    clusters: Vec<MSACluster<'a>>,
    parent_indices: Vec<Option<usize>>,
    depth: usize,
}

impl<'a> MSATree<'a> {
    pub fn new(
        cluster_names: Vec<&'a str>,
        cluster_centers: Vec<Vec<u8>>,
        cluster_layers: Vec<usize>,
        parent_indices: Vec<Option<usize>>,
    ) -> Self {
        Self {
            cluster_names,
            cluster_centers,
            cluster_layers,
            parent_indices,
            clusters: Vec::new(),
            depth: Default::default(),
        }
    }

    pub fn initialize_clusters(&mut self) {
        self.clusters = self
            .cluster_names
            .iter()
            .zip(self.cluster_centers.iter())
            .zip(self.cluster_layers.iter())
            .zip(self.parent_indices.iter())
            .enumerate()
            .map(|(index, (((name, center), &layer), &parent_index))| {
                MSACluster::new(index, name, center.to_vec(), layer, parent_index)
            })
            .collect();

        self.clusters.iter_mut().for_each(|c| {
            c.sibling_index = MSACluster::get_sibling_index(c.index);
        });

        self.depth = self.clusters.last().unwrap().layer;
    }

    fn align_layer(&mut self) {
        let aligned_clusters = self
            .clusters
            .clone()
            .into_iter()
            .map(|c| {
                let new_cluster = if c.layer == self.depth {
                    let center = c.center;
                    let sib_ind = c.sibling_index.unwrap();
                    let sib_center = &self.clusters[sib_ind].center;
                    let aligned_center =
                        needleman_wunsch::string::needleman_wunsch_alignment_only::<u8, i32>(
                            center.as_slice(),
                            sib_center.as_slice(),
                        )
                        .0;
                    MSACluster {
                        index: c.index,
                        name: c.name,
                        center: aligned_center,
                        layer: c.layer,
                        sibling_index: c.sibling_index,
                        parent_index: c.parent_index,
                    }
                } else {
                    c
                };
                new_cluster
            })
            .collect();

        self.clusters = aligned_clusters;
    }

    fn modify_parent_centers(&mut self) {
        let clusters_copy = self.clusters.clone();
        let mut edited_parents = self.clusters.clone();

        self.clusters.iter().for_each(|c| {
            if c.layer == self.depth && c.index % 2 == 0 {
                let index = c.index;
                let center = clusters_copy[index].center.clone();
                let gap_indices = get_gap_indices(center);
                let parent_index = c.parent_index.unwrap();
                let parent_center = clusters_copy[parent_index].center.clone();
                let edited_parent_center = if !gap_indices.is_empty() {
                    apply_edits(parent_center, gap_indices)
                } else {
                    parent_center
                };
                edited_parents[parent_index] = MSACluster {
                    index: parent_index,
                    name: self.clusters[parent_index].name,
                    center: edited_parent_center,
                    layer: self.clusters[parent_index].layer,
                    sibling_index: self.clusters[parent_index].sibling_index,
                    parent_index: self.clusters[parent_index].parent_index,
                };
            }
        });

        self.clusters = edited_parents;
    }

    pub fn align_all(&mut self) {
        while self.depth > 0 {
            self.align_layer();
            self.modify_parent_centers();

            self.depth -= 1;
        }
    }

    pub fn display_aligned_tree(self) {
        for cluster in self.clusters {
            println!(
                "{}: {:?}",
                cluster.name,
                std::str::from_utf8(&cluster.center)
            );
        }
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

pub fn distance_many_to_many(left: &[&[u8]], right: &[&[u8]]) -> Vec<Vec<i32>> {
    left.iter()
        .map(|&l| distance_one_to_many(l, right))
        .collect()
}

pub fn distance_one_to_many(left: &[u8], right: &[&[u8]]) -> Vec<i32> {
    right
        .iter()
        .map(|&r| needleman_wunsch::string::needleman_wunsch(left, r))
        .collect()
}

pub fn get_distances_between_instances(
    table: &[Vec<i32>],
    indices: &[usize],
    center_index: usize,
) -> Vec<i32> {
    let mut distances = Vec::new();
    for i in indices {
        distances.push(table[center_index][*i]);
    }

    distances
}

#[cfg(test)]
mod tests {
    use super::median;
    use crate::msa_by_hand::align_by_hand::distance_many_to_many;
    use crate::msa_by_hand::align_by_hand::get_distances_between_instances;

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
        let sequences_011 = [
            "ctttgttacc".as_bytes(),
            "ccttgtgtac".as_bytes(),
            "cttgcctgtt".as_bytes(),
        ]
        .to_vec();

        let pairwise_distances_011 = distance_many_to_many(&sequences_011.clone(), &sequences_011);
        let (center_index_011, _) = median(pairwise_distances_011);
        let center_sequence_011: &str =
            std::str::from_utf8(&sequences_011[center_index_011]).unwrap();

        assert_eq!(center_index_011, 1);
        assert_eq!(center_sequence_011, "ccttgtgtac");
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

        let distances_center =
            get_distances_between_instances(&pairwise_distances, &[11, 13, 14], 13);
        assert_eq!(distances_center, vec![3, 0, 5]);
    }
}
