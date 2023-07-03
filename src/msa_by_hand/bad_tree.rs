use sequence_analysis::msa_by_hand::msa_by_hand::MSA_Tree;

fn main() {
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

    let bad_tree_cluster_names = vec![
        "0", "00", "01", "000", "001", "010", "011", "0000", "0001", "0010", "0011", "0100",
        "0101", "0110", "0111", "00000", "00001", "00010", "00011", "00100", "00101", "00110",
        "00111", "01010", "01011", "01110", "01111", "011110", "011111",
    ];

    println!("number of clusters is {}", bad_tree_cluster_names.len());


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

    let mut msa_tree = MSA_Tree::new(sequences, bad_tree_cluster_names, bad_tree_cluster_centers);
    msa_tree.initialize_clusters();
    msa_tree.align_all();
    msa_tree.display_aligned_tree()
}
