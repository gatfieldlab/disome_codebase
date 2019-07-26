#' Translate Nucleic Acid Sequence
#'
#' This function translates the nucleic acid sequence into amino acid
#' sequence based on the standard genetic code
#' @param nuc_seq Nucleic acid sequence, character.
#' @keywords deduce_nuc_seq
#' @export
#' @examples
#' translate("ATGCGCCGATTAGCA")
#' 
translate <- function(nuc_seq) {
  translatable_len <- nchar(nuc_seq) - nchar(nuc_seq) %% 3
  aa_seq <- ""
  for (i in seq(1, translatable_len, 3)) {
    codon <- substr(nuc_seq, i, i + 2)
    aa <- aminos[which(codons == codon)]
    aa_seq = paste(aa_seq, aa, sep = "")
  }
  aa_seq
}