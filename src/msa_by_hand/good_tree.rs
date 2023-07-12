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
        "0101", "0110", "0111", "00000", "00001", "00010", "00011", "00110", "00111", "01000",
        "01001", "01100", "01101", "01110", "01111", "000110", "000111",
    ];

    let good_tree_parent_indices: Vec<Option<usize>> = vec![
        None,
        good_tree_cluster_names.iter().position(|&c| c == "0"),
        good_tree_cluster_names.iter().position(|&c| c == "0"),
        good_tree_cluster_names.iter().position(|&c| c == "00"),
        good_tree_cluster_names.iter().position(|&c| c == "00"),
        good_tree_cluster_names.iter().position(|&c| c == "01"),
        good_tree_cluster_names.iter().position(|&c| c == "01"),
        good_tree_cluster_names.iter().position(|&c| c == "000"),
        good_tree_cluster_names.iter().position(|&c| c == "000"),
        good_tree_cluster_names.iter().position(|&c| c == "001"),
        good_tree_cluster_names.iter().position(|&c| c == "001"),
        good_tree_cluster_names.iter().position(|&c| c == "010"),
        good_tree_cluster_names.iter().position(|&c| c == "010"),
        good_tree_cluster_names.iter().position(|&c| c == "011"),
        good_tree_cluster_names.iter().position(|&c| c == "011"),
        good_tree_cluster_names.iter().position(|&c| c == "0000"),
        good_tree_cluster_names.iter().position(|&c| c == "0000"),
        good_tree_cluster_names.iter().position(|&c| c == "0001"),
        good_tree_cluster_names.iter().position(|&c| c == "0001"),
        good_tree_cluster_names.iter().position(|&c| c == "0011"),
        good_tree_cluster_names.iter().position(|&c| c == "0011"),
        good_tree_cluster_names.iter().position(|&c| c == "0100"),
        good_tree_cluster_names.iter().position(|&c| c == "0100"),
        good_tree_cluster_names.iter().position(|&c| c == "0110"),
        good_tree_cluster_names.iter().position(|&c| c == "0110"),
        good_tree_cluster_names.iter().position(|&c| c == "0111"),
        good_tree_cluster_names.iter().position(|&c| c == "0111"),
        good_tree_cluster_names.iter().position(|&c| c == "00011"),
        good_tree_cluster_names.iter().position(|&c| c == "00011"),
    ];

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
        "cttgcctgtt".as_bytes().to_vec(), //00000
        "gtgacctcat".as_bytes().to_vec(), //00001
        "attgttcaga".as_bytes().to_vec(), //00010
        "ctttgttacc".as_bytes().to_vec(), //00011
        "agaaagcgca".as_bytes().to_vec(), //00110
        "ggatgtcaca".as_bytes().to_vec(), //00111
        "tcgaaacggc".as_bytes().to_vec(), //01000
        "gcacacgtgg".as_bytes().to_vec(), //01001
        "agtgggatgt".as_bytes().to_vec(), //01100
        "ggccgaatag".as_bytes().to_vec(), //01101
        "gaatgcaagg".as_bytes().to_vec(), //01110
        "cagcggaggg".as_bytes().to_vec(), //01111
        "ccttgtgtac".as_bytes().to_vec(), //000110
        "ctttgttacc".as_bytes().to_vec(), //000111
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
        5, //00000
        5, //00001
        4, //00010
        4, //00011
        5, //00110
        5, //00111
        5, //01000
        5, //01001
        5, //01100
        5, //01101
        5, //01110
        5, //01111
        5, //000110
        5, //000110
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
