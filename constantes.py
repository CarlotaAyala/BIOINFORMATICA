"""Datos fijos del simulador: colores, complementos y tabla de codones."""

RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"

BASE_COLOR = {
    "A": "\033[42m\033[30m",        # verde
    "T": "\033[41m\033[97m",        # rojo
    "U": "\033[48;5;208m\033[30m",  # naranja
    "G": "\033[44m\033[97m",        # azul
    "C": "\033[43m\033[30m",        # amarillo
}
TITLE_COLOR = "\033[95m"
HEAD_COLOR = "\033[96m"
ENZ_COLOR = "\033[94m"
OK_COLOR = "\033[92m"
STOP_COLOR = "\033[41m\033[97m"

COMP_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}
COMP_RNA = {"A": "U", "T": "A", "G": "C", "C": "G"}

STOP_CODONS_DNA = ["TAA", "TAG", "TGA"]

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}