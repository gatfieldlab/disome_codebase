#' Convert Single Letter Amino Acid Sequence into Three Letter Code
#'
#' This function converts a character vector of single letter aminoacid codes
#' into a vector of 3-letter aminoacid codes
#' @param aa_seq Amino acid sequence in uppercase single letter alphabet
#' @keywords translate
#' @export
#' @examples
#' three_letter("DIT")
#'
three_letter <- function(aa_seq) {
  new_codes <- c()
  for (i in 1:nchar(aa_seq)) {
    amino <- substr(aa_seq, i, i)
    new_codes <- c(new_codes, aminoacids_3l[which(aminoacids_1l == amino)])
  }
  new_codes
}