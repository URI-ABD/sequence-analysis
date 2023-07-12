use sequence_analysis::msa_by_hand::align_by_hand::MSATree;

fn main() {
    let _sequences = [
        "tgggtactga",
        "attgttcaga",
        "gaatgcaagg",
        "ggccgaatag",
        "agtgggatgt",
        "gcacacgtgg",
        "tcgaaacggc",
        "cagcggaggg",
        "agaaagcgca",
        "gtgacctcat",
        "ggatgtcaca",
        "ctttgttacc",
        "gcagagggcc",
        "ccttgtgtac",
        "cttgcctgtt",
    ]
    .to_vec();

    let bad_tree_cluster_names = vec![
        "0", "00", "01", "000", "001", "010", "011", "0000", "0001", "0010", "0011", "0100",
        "0101", "0110", "0111", "00000", "00001", "00010", "00011", "00100", "00101", "00110",
        "00111", "01010", "01011", "01110", "01111", "011110", "011111",
    ];

    let bad_tree_cluster_centers = vec![
        "cttgcctgtt".as_bytes().to_vec(), // 0
        "gtgacctcat".as_bytes().to_vec(), // 00
        "ggccgaatag".as_bytes().to_vec(), // 01
        "gcagagggcc".as_bytes().to_vec(), // 000
        "gtgacctcat".as_bytes().to_vec(), // 001
        "ggatgtcaca".as_bytes().to_vec(), // 010
        "ccttgtgtac".as_bytes().to_vec(), // 011
        "agtgggatgt".as_bytes().to_vec(), // 0000
        "gaatgcaagg".as_bytes().to_vec(), // 0001
        "tcgaaacggc".as_bytes().to_vec(), // 0010
        "gcacacgtgg".as_bytes().to_vec(), // 0011
        "tgggtactga".as_bytes().to_vec(), // 0100
        "attgttcaga".as_bytes().to_vec(), // 0101
        "ggccgaatag".as_bytes().to_vec(), // 0110
        "ccttgtgtac".as_bytes().to_vec(), // 0111
        "cagcggaggg".as_bytes().to_vec(), // 00000
        "agtgggatgt".as_bytes().to_vec(), // 00001
        "gcagagggcc".as_bytes().to_vec(), // 00010
        "gaatgcaagg".as_bytes().to_vec(), // 00011
        "agaaagcgca".as_bytes().to_vec(), // 00100
        "tcgaaacggc".as_bytes().to_vec(), // 00101
        "gtgacctcat".as_bytes().to_vec(), // 00110
        "gcacacgtgg".as_bytes().to_vec(), // 00111
        "ggatgtcaca".as_bytes().to_vec(), // 01010
        "attgttcaga".as_bytes().to_vec(), // 01011
        "cttgcctgtt".as_bytes().to_vec(), // 01110
        "ctttgttacc".as_bytes().to_vec(), // 01111
        "ccttgtgtac".as_bytes().to_vec(), // 011110
        "ctttgttacc".as_bytes().to_vec(), // 011111
    ];

    let bad_tree_cluster_layers = vec![
        0, // 0
        1, // 00
        1, // 01
        3, // 000
        3, // 001
        2, // 010
        2, // 011
        4, // 0000
        4, // 0001
        4, // 0010
        4, // 0011
        4, // 0100
        4, // 0101
        3, // 0110
        3, // 0111
        5, // 00000
        5, // 00001
        5, // 00010
        5, // 00011
        5, // 00100
        5, // 00101
        5, // 00110
        5, // 00111
        5, // 01010
        5, // 01011
        4, // 01110
        4, // 01111
        5, // 011110
        5, // 011111
    ];

    let bad_tree_parent_indices = vec![
        None,                                                      // 0
        bad_tree_cluster_names.iter().position(|&c| c == "0"),     // 00
        bad_tree_cluster_names.iter().position(|&c| c == "0"),     // 01
        bad_tree_cluster_names.iter().position(|&c| c == "00"),    // 000
        bad_tree_cluster_names.iter().position(|&c| c == "00"),    // 001
        bad_tree_cluster_names.iter().position(|&c| c == "01"),    // 010
        bad_tree_cluster_names.iter().position(|&c| c == "01"),    // 011
        bad_tree_cluster_names.iter().position(|&c| c == "000"),   // 0000
        bad_tree_cluster_names.iter().position(|&c| c == "000"),   // 0001
        bad_tree_cluster_names.iter().position(|&c| c == "001"),   // 0010
        bad_tree_cluster_names.iter().position(|&c| c == "001"),   // 0011
        bad_tree_cluster_names.iter().position(|&c| c == "010"),   // 0100
        bad_tree_cluster_names.iter().position(|&c| c == "010"),   // 0101
        bad_tree_cluster_names.iter().position(|&c| c == "011"),   // 0110
        bad_tree_cluster_names.iter().position(|&c| c == "011"),   // 0111
        bad_tree_cluster_names.iter().position(|&c| c == "0000"),  // 00000
        bad_tree_cluster_names.iter().position(|&c| c == "0000"),  // 00001
        bad_tree_cluster_names.iter().position(|&c| c == "0001"),  // 00010
        bad_tree_cluster_names.iter().position(|&c| c == "0001"),  // 00011
        bad_tree_cluster_names.iter().position(|&c| c == "0010"),  // 00100
        bad_tree_cluster_names.iter().position(|&c| c == "0010"),  // 00101
        bad_tree_cluster_names.iter().position(|&c| c == "0011"),  // 00110
        bad_tree_cluster_names.iter().position(|&c| c == "0011"),  // 00111
        bad_tree_cluster_names.iter().position(|&c| c == "0101"),  // 01010
        bad_tree_cluster_names.iter().position(|&c| c == "0101"),  // 01011
        bad_tree_cluster_names.iter().position(|&c| c == "0111"),  // 01110
        bad_tree_cluster_names.iter().position(|&c| c == "0111"),  // 01111
        bad_tree_cluster_names.iter().position(|&c| c == "01111"), // 011110
        bad_tree_cluster_names.iter().position(|&c| c == "01111"), // 011111
    ];

    let mut msa_tree = MSATree::new(
        bad_tree_cluster_names,
        bad_tree_cluster_centers,
        bad_tree_cluster_layers,
        bad_tree_parent_indices,
    );
    msa_tree.initialize_clusters();
    msa_tree.align_all();
    msa_tree.display_aligned_tree()
}
