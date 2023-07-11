use sequence_analysis::msa_by_hand::align_by_hand::MSATree;

fn main() {
    let _sequences = [
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

    let good_tree_cluster_names = vec![
        "0", "00", "01", "000", "001", "010", "011", "0000", "0001", "0010", "0011", "0100",
        "0101", "0110", "0111", "00010", "00011",
    ];

    let good_tree_parent_indices = vec![0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 8, 8];

    let good_tree_cluster_centers = vec![
        "gaatgcaagg".as_bytes().to_vec(), //0
        "ggatgtcaca".as_bytes().to_vec(), //00
        "cagcggaggg".as_bytes().to_vec(), //01
        "ctttgttacc".as_bytes().to_vec(), //000
        "ggatgtcaca".as_bytes().to_vec(), //001
        "gcacacgtgg".as_bytes().to_vec(), //010
        "cagcggaggg".as_bytes().to_vec(), //011
        "gtgacctcat".as_bytes().to_vec(), //0000
        "ctttgttacc".as_bytes().to_vec(), //0001
        "tgggtactga".as_bytes().to_vec(), //0010
        "ggatgtcaca".as_bytes().to_vec(), //0011
        "gcacacgtgg".as_bytes().to_vec(), //0100
        "gcagagggcc".as_bytes().to_vec(), //0101
        "ggccgaatag".as_bytes().to_vec(), //0110
        "cagcggaggg".as_bytes().to_vec(), //0111
        "attgttcaga".as_bytes().to_vec(), //00010
        "ctttgttacc".as_bytes().to_vec(), //00011
    ];
    let good_tree_cluster_layers = vec![
        0, //0
        1, //00
        1, //01
        2, //000
        2, //001
        3, //010
        3, //011
        3, //0000
        3, //0001
        4, //0010
        4, //0011
        4, //0100
        4, //0101
        4, //0110
        4, //0111
        4, //00010
        4, //00011
    ];

    assert_eq!(
        good_tree_cluster_names.len(),
        good_tree_cluster_centers.len()
    );

    let mut msa_tree = MSATree::new(
        good_tree_cluster_names,
        good_tree_cluster_centers,
        good_tree_cluster_layers,
        good_tree_parent_indices,
    );
    msa_tree.initialize_clusters();
    msa_tree.align_all();
    msa_tree.display_aligned_tree()
}
