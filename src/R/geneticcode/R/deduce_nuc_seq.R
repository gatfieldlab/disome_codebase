#' Deduce Nucleic Acid Sequence from Amino Acid Sequence
#'
#' This function produces all possible nucleic acid sequences that
#' would translate into the given amino acid sequence.
#' @param aa_seq Amino acid sequence in uppercase single letter alphabet
#' @keywords translate
#' @export
#' @examples
#' deduce_nuc_seq("DIT")
#'
deduce_nuc_seq <- function(aa_seq) {
  nuc_seqs <- list("")
  for (i in 1:nchar(aa_seq)) {
    amino <- substr(aa_seq, i, i)
    temp_seqs <- list()
    possible_codons <- codons[which(aminos == amino)]
    for (codon in possible_codons) {
      for (base_seq in nuc_seqs) {
        new_seq <- paste(unlist(base_seq), codon, sep = "")
        temp_seqs <- list(temp_seqs, list(new_seq))
      }
    }
    nuc_seqs <- temp_seqs
  }
  unlist(nuc_seqs)
}

# to do: cut it for length input
