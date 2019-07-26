#' Codons
#'
#' Codons in uppercase letters excluding Stop-codons
#'
#' @export
"codons"


codons <- c("ATA", "ATC", "ATT", "ATG", "ACA", "ACC", "ACG", "ACT",
            "AAC", "AAT", "AAA", "AAG", "AGC", "AGT", "AGA", "AGG",
            "CTA", "CTC", "CTG", "CTT", "CCA", "CCC", "CCG", "CCT",
            "CAC", "CAT", "CAA", "CAG", "CGA", "CGC", "CGG", "CGT",
            "GTA", "GTC", "GTG", "GTT", "GCA", "GCC", "GCG", "GCT",
            "GAC", "GAT", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT",
            "TCA", "TCC", "TCG", "TCT", "TTC", "TTT", "TTA", "TTG",
            "TAC", "TAT", "TAA", "TAG", "TGC", "TGT", "TGA", "TGG")

#' Aminos
#'
#' Single letter amino acids corresponding to codons,
#' only for internal usage - use aminoacids or aminoacids_1l
#' instead.
#'
"aminos"
aminos <- c("I", "I", "I", "M", "T", "T", "T", "T",
            "N", "N", "K", "K", "S", "S", "R", "R",
            "L", "L", "L", "L", "P", "P", "P", "P",
            "H", "H", "Q", "Q", "R", "R", "R", "R",
            "V", "V", "V", "V", "A", "A", "A", "A",
            "D", "D", "E", "E", "G", "G", "G", "G",
            "S", "S", "S", "S", "F", "F", "L", "L",
            "Y", "Y", "*", "*", "C", "C", "*", "W")

#' Amino acids
#'
#' Standard amino acids in full names
#'
#' @export
"aminoacids"
aminoacids <- c("Glycine", "Proline", "Alanine", "Valine", "Leucine",
                "Isoleucine", "Methionine", "Cysteine", "Phenylalanine",
                "Tyrosine", "Tryptophan", "Histidine", "Lysine", "Arginine",
                "Glutamine", "Asparagine", "Glutamate", "Aspartate", "Serine",
                "Threonine")

#' Amino acids 1 letter
#'
#' Standard amino acids in single letter format
#'
#' @export
"aminoacids_1l"
aminoacids_1l <- c("G", "P", "A", "V", "L", "I", "M", "C", "F", "Y", "W", "H",
                   "K", "R", "Q", "N", "E", "D", "S", "T")


#' Amino acids 3 letter
#'
#' Standard amino acids in three letter format
#'
#' @export
"aminoacids_3l"
aminoacids_3l <- c("Gly", "Pro", "Ala", "Val", "Leu", "Ile", "Met", "Cys",
                   "Phe", "Tyr", "Trp", "His", "Lys", "Arg", "Gln", "Asn",
                   "Glu", "Asp", "Ser", "Thr")
