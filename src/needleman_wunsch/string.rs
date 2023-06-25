use super::alignment_helpers;
use abd_clam::number::Number;

/// Calculate the edit distance between two strings using Needleman-Wunsch table.
/// This function is only accurate with a scoring scheme for which all penalties
/// are non-negative.
///
/// * [Demo](https://bioboot.github.io/bimm143_W20/class-material/nw/)
/// * [Wikipedia](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
///
/// # Arguments:
///
/// * `x`: unaligned sequence represented as a slice of type T
/// * `y`: unaligned sequence represented as a slice of type T
pub fn needleman_wunsch<T: Number, U: Number>(x: &[T], y: &[T]) -> U {
    let table = alignment_helpers::compute_nw_table(x, y);
    let edit_distance: usize = table[table.len() - 1][table[0].len() - 1].0;

    U::from(edit_distance).unwrap()
}

/// Determine the set of edits needed to turn one unaligned sequence into another,
/// as well as the edit distance between the two sequences.
///
/// Contrast to `needleman_wunsch_with_edits_iterative`, which uses an iterative, as
/// opposed to recursive, traceback function
///
/// For now, in cases where there exist ties for the shortest edit distance, we only
/// return one alignment.
///
/// # Arguments:
///
/// * `x`: unaligned sequence represented as a slice of type T
/// * `y`: unaligned sequence represented as a slice of type T
pub fn needleman_wunsch_with_edits_recursive<T: Number, U: Number>(
    x: &[T],
    y: &[T],
) -> ([Vec<alignment_helpers::Edit<T>>; 2], U) {
    let table = alignment_helpers::compute_nw_table(x, y);
    let (aligned_x, aligned_y) = alignment_helpers::traceback_recursive(&table, (x, y));

    let edit_x_into_y = alignment_helpers::alignment_to_edits(&aligned_x, &aligned_y);
    let edit_y_into_x = alignment_helpers::alignment_to_edits(&aligned_y, &aligned_x);

    let edit_distance: usize = table[table.len() - 1][table[0].len() - 1].0;

    (
        [edit_x_into_y, edit_y_into_x],
        U::from(edit_distance).unwrap(),
    )
}

/// Determine the set of edits needed to turn one unaligned sequence into another,
/// as well as the edit distance between the two sequences.
///
/// Contrast to `needleman_wunsch_with_edits_recursive`, which uses a recursive, as
/// opposed to iterative, traceback function
///
/// For now, in cases where there exist ties for the shortest edit distance, we only
/// return one alignment.
///
/// # Arguments:
///
/// * `x`: unaligned sequence represented as a slice of type T
/// * `y`: unaligned sequence represented as a slice of type T
pub fn needleman_wunsch_with_edits_iterative<T: Number, U: Number>(
    x: &[T],
    y: &[T],
) -> ([Vec<alignment_helpers::Edit<T>>; 2], U) {
    let table = alignment_helpers::compute_nw_table(x, y);
    let (aligned_x, aligned_y) = alignment_helpers::traceback_iterative(&table, (x, y));

    let edit_x_into_y = alignment_helpers::alignment_to_edits(&aligned_x, &aligned_y);
    let edit_y_into_x = alignment_helpers::alignment_to_edits(&aligned_y, &aligned_x);

    let edit_distance: usize = table[table.len() - 1][table[0].len() - 1].0;

    (
        [edit_x_into_y, edit_y_into_x],
        U::from(edit_distance).unwrap(),
    )
}

pub fn needleman_wunsch_alignment_only<T: Number, U: Number>(x: &[T], y: &[T]) -> (Vec<T>, Vec<T>) {
    let table = alignment_helpers::compute_nw_table(x, y);
    let (aligned_x, aligned_y) = alignment_helpers::traceback_recursive(&table, (x, y));

    (aligned_x, aligned_y)
}
