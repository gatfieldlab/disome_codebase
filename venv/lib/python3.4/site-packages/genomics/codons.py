#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides codon codes, charges and codon usage related
data structures and tools
"""

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


# Standard genetic code invariables

NTS = ["A", "G", "C", "T"]
CODON_TABLE = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
               'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
               'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
               'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
               'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
               'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
               'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
               'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
               'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
               'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
               'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
STOP_CODONS = set(["TAA", "TGA", "TAG"])


# Charge of side chains by aminoacid

POSITIVE_AA = set(['K', 'R', 'H'])
NEGATIVE_AA = set(['D', 'E'])
NEUTRAL_AA = set(CODON_TABLE.values()) - (POSITIVE_AA | NEGATIVE_AA)
CHARGE_TABLE = {aa: "U" if aa in NEUTRAL_AA else "N" if aa in NEGATIVE_AA else "P" for
                aa in CODON_TABLE.values()}
CHARGE = {'U': 0, 'N': -1, 'P': 1}


# Codon usage tables

CODON_USAGE = {

# This table is obtained from
# http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10090
# Frequencies are per thousand.

    'mouse': {
        'GGG': 15.17,
        'GGA': 16.77,
        'GGT': 11.43,
        'GGC': 21.20,
        'GAG': 39.37,
        'GAA': 26.96,
        'GAT': 20.99,
        'GAC': 26.03,
        'GTG': 28.38,
        'GTA': 7.45,
        'GTT': 10.70,
        'GTC': 15.40,
        'GCG': 6.40,
        'GCA': 15.84,
        'GCT': 20.02,
        'GCC': 26.00,
        'AGG': 12.21,
        'AGA': 12.11,
        'AGT': 12.69,
        'AGC': 19.69,
        'AAG': 33.64,
        'AAA': 21.92,
        'AAT': 15.58,
        'AAC': 20.35,
        'ATG': 22.82,
        'ATA': 7.36,
        'ATT': 15.40,
        'ATC': 22.51,
        'ACG': 5.63,
        'ACA': 15.96,
        'ACT': 13.66,
        'ACC': 18.96,
        'TGG': 12.50,
        'TGA': 1.64,
        'TGT': 11.40,
        'TGC': 12.28,
        'TAG': 0.78,
        'TAA': 0.95,
        'TAT': 12.17,
        'TAC': 16.06,
        'TTG': 13.44,
        'TTA': 6.73,
        'TTT': 17.21,
        'TTC': 21.82,
        'TCG': 4.23,
        'TCA': 11.81,
        'TCT': 16.23,
        'TCC': 18.10,
        'CGG': 10.22,
        'CGA': 6.58,
        'CGT': 4.68,
        'CGC': 9.36,
        'CAG': 34.09,
        'CAA': 11.96,
        'CAT': 10.62,
        'CAC': 15.31,
        'CTG': 39.52,
        'CTA': 8.07,
        'CTT': 13.44,
        'CTC': 20.18,
        'CCG': 6.18,
        'CCA': 17.27,
        'CCT': 18.37,
        'CCC': 18.21}}


# Codon tRNA adaptation indices

CODON_TAI = {

# The values in this table is calculated using tAI.R version 0.2,
# which forms part of codonR, a package for the analysis of codon usage
# The source is slightly modified to include a weight for Met, ATG
# As input, the tRNA gene copy count from GtRNAdb Mus musculus (mm10 Dec 2011)
# {file: mm10_tRNA_gene_count.txt} is used.

    'mouse': {
        "TTT": 0.0680171277997365,
        "TTC": 0.11528326745718,
        "TTA": 0.11528326745718,
        "TTG": 0.102766798418972,
        "TCT": 0.174407114624506,
        "TCC": 0.135046113306983,
        "TCA": 0.0494235836627141,
        "TCG": 0.0652173913043478,
        "TAT": 0.106884057971015,
        "TAC": 0.181159420289855,
        "TGT": 0.599472990777339,
        "TGC": 1,
        "TGG": 0.137022397891963,
        "CTT": 0.0823451910408432,
        "CTC": 0.0592885375494071,
        "CTA": 0.0658843873517787,
        "CTG": 0.185770750988142,
        "CCT": 0.141469038208169,
        "CCC": 0.11133069828722,
        "CCA": 0.131765480895916,
        "CCG": 0.0915678524374176,
        "CAT": 0.113636363636364,
        "CAC": 0.176548089591568,
        "CAA": 0.115284914361001,
        "CAG": 0.234519104084321,
        "CGT": 0.0988142292490119,
        "CGC": 0.0711462450592885,
        "CGA": 0.0823550724637681,
        "CGG": 0.0757575757575758,
        "ATT": 0.207345191040843,
        "ATC": 0.158761528326746,
        "ATA": 0.082364953886693,
        "ATG": 0.312911725955204,
        "ACT": 0.181159420289855,
        "ACC": 0.130434782608696,
        "ACA": 0.0658942687747036,
        "ACG": 0.0869565217391304,
        "AAT": 0.126317523056654,
        "AAC": 0.214097496706192,
        "AAA": 0.230566534914361,
        "AAG": 0.50197628458498,
        "AGT": 0.077733860342556,
        "AGC": 0.131752305665349,
        "AGA": 0.0988142292490119,
        "AGG": 0.113965744400527,
        "GTT": 0.141469038208169,
        "GTC": 0.11133069828722,
        "GTA": 0.0494202898550725,
        "GTG": 0.213438735177866,
        "GCT": 0.417654808959157,
        "GCC": 0.338603425559947,
        "GCA": 0.181197299077734,
        "GCG": 0.22266139657444,
        "GAT": 0.155467720685112,
        "GAC": 0.263504611330698,
        "GAA": 0.131752305665349,
        "GAG": 0.272727272727273,
        "GGT": 0.17868906455863,
        "GGC": 0.270750988142293,
        "GGA": 0.131755599472991,
        "GGG": 0.157444005270092}}


# Codon species-specific tRNA adaptation indices

CODON_STAI = {

# Sij weights were calculated by the stAIcalc utility
# using ALL valid CDS sequences from mouse (mm10)
# then used similarly as above (CODON_TAI) in tAI.R script to calculate the
# normalized weights

    'mouse': {
        "TTT": 0.0958618112175996,
        "TTC": 0.114787138801787,
        "TTA": 0.114787138801787,
        "TTG": 0.180379789545665,
        "TCT": 0.177676171319353,
        "TCC": 0.177500551104233,
        "TCA": 0.173699056599529,
        "TCG": 0.0983889761158174,
        "TAT": 0.150639989056228,
        "TAC": 0.180379789545665,
        "TAA": 0.0163981626859696,
        "TAG": 0.0163981626859696,
        "TGT": 0.838070830265395,
        "TGC": 1,
        "TGA": 0.0288486195401316,
        "TGG": 0.147583464173726,
        "CTT": 0.0819908134298478,
        "CTC": 0.0805511942091317,
        "CTA": 0.127844935014689,
        "CTG": 0.229574277603574,
        "CCT": 0.144879845947414,
        "CCC": 0.14528007342058,
        "CCA": 0.230788956321053,
        "CCG": 0.180379789545665,
        "CAT": 0.15334360728254,
        "CAC": 0.180091865701522,
        "CAA": 0.127237595655949,
        "CAG": 0.311565091033422,
        "CGT": 0.0983889761158174,
        "CGC": 0.096661433050958,
        "CGA": 0.15669355455482,
        "CGG": 0.131185301487756,
        "ATT": 0.210472496691292,
        "ATC": 0.209721028787886,
        "ATA": 0.231396295679793,
        "ATG": 0.311565091033422,
        "ACT": 0.180379789545665,
        "ACC": 0.17721262726009,
        "ACA": 0.202547676139661,
        "ACG": 0.131185301487756,
        "AAT": 0.178029077975542,
        "AAC": 0.213176114917604,
        "AAA": 0.229574277603574,
        "AAG": 0.655926507438782,
        "AGT": 0.109556355677257,
        "AGC": 0.131185301487756,
        "AGA": 0.0983889761158174,
        "AGG": 0.180379789545665,
        "GTT": 0.144879845947414,
        "GTC": 0.14528007342058,
        "GTA": 0.148798142891205,
        "GTG": 0.245972440289543,
        "GCT": 0.431935919615928,
        "GCC": 0.436128144105884,
        "GCA": 0.466740297191393,
        "GCG": 0.344361416405361,
        "GAT": 0.219112711354513,
        "GAC": 0.262370602975513,
        "GAA": 0.131185301487756,
        "GAG": 0.36075957909133,
        "GGT": 0.238214492266795,
        "GGC": 0.278192917973196,
        "GGA": 0.156086215196081,
        "GGG": 0.245972440289543},


# Sij weights were calculated by the stAIcalc utility
# using the highest RPF expressor genes in our liver experiments:
# 'unavailable': 23992, 'too_short': 174, 'untranslatable': 11148,
# 'not_in_list': 34076, 'ok': 6347
# -- note -- this is almost identical if calculated using TR
# then used similarly as above in tAI.R script to calculate the normalized weights

    'mouse-liver': {
        "TTT": 0.0958342182610774,
        "TTC": 0.114754098360656,
        "TTA": 0.114754098360656,
        "TTG": 0.180327868852459,
        "TCT": 0.177625028838234,
        "TCC": 0.180327868852459,
        "TCA": 0.154219793564056,
        "TCG": 0.0983606557377049,
        "TAT": 0.150596628695979,
        "TAC": 0.180327868852459,
        "TAA": 0.0163934426229508,
        "TAG": 0.0163934426229508,
        "TGT": 0.837829599146472,
        "TGC": 1,
        "TGA": 0.0268973891924712,
        "TGG": 0.147540983606557,
        "CTT": 0.0819672131147541,
        "CTC": 0.0819672131147541,
        "CTA": 0.118093503339405,
        "CTG": 0.229508196721312,
        "CCT": 0.144838143592332,
        "CCC": 0.147540983606557,
        "CCA": 0.215179113539769,
        "CCG": 0.180327868852459,
        "CAT": 0.153299468710204,
        "CAC": 0.180327868852459,
        "CAA": 0.125258044930176,
        "CAG": 0.311475409836066,
        "CGT": 0.0983606557377049,
        "CGC": 0.0983606557377049,
        "CGA": 0.144990892531876,
        "CGG": 0.131147540983607,
        "ATT": 0.210411914084135,
        "ATC": 0.213114754098361,
        "ATA": 0.208014571948998,
        "ATG": 0.311475409836066,
        "ACT": 0.180327868852459,
        "ACC": 0.180327868852459,
        "ACA": 0.181117182756527,
        "ACG": 0.131147540983607,
        "AAT": 0.177977833913429,
        "AAC": 0.213114754098361,
        "AAA": 0.229508196721312,
        "AAG": 0.655737704918033,
        "AGT": 0.109524820869803,
        "AGC": 0.131147540983607,
        "AGA": 0.0983606557377049,
        "AGG": 0.180327868852459,
        "GTT": 0.144838143592332,
        "GTC": 0.147540983606557,
        "GTA": 0.133211900425015,
        "GTG": 0.245901639344262,
        "GCT": 0.43181159076277,
        "GCC": 0.442622950819672,
        "GCA": 0.421918639951427,
        "GCG": 0.344262295081967,
        "GAT": 0.219049641739606,
        "GAC": 0.262295081967213,
        "GAA": 0.131147540983607,
        "GAG": 0.360655737704918,
        "GGT": 0.238145924376782,
        "GGC": 0.278688524590164,
        "GGA": 0.152155434122647,
        "GGG": 0.245901639344262}}
