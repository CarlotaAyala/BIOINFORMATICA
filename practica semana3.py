#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulador del Dogma Central de la Biología Molecular
Bioinformática - ULPGC EII - Práctica 1

Simula, a partir de una cadena de ADN, los procesos de:
  1. Replicación   (ADN -> ADN)   - cadena líder y cadena rezagada (fragmentos de Okazaki)
  2. Transcripción (ADN -> ARN)
  3. Traducción    (ARN -> Proteína)

Uso:
    python3 dogma_central.py
    python3 dogma_central.py --seq ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
    python3 dogma_central.py --random
"""

import argparse
import random
import textwrap
import time

# ---------------------------------------------------------------------------
# Colores ANSI (funcionan en la mayoría de terminales; en Windows usar
# Windows Terminal, PowerShell moderno o WSL para ver el color correctamente)
# ---------------------------------------------------------------------------
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"

BASE_COLOR = {
    "A": "\033[42m\033[30m",   # verde
    "T": "\033[41m\033[97m",   # rojo
    "U": "\033[48;5;208m\033[30m",  # naranja
    "G": "\033[44m\033[97m",   # azul
    "C": "\033[43m\033[30m",   # amarillo
}
TITLE_COLOR = "\033[95m"
HEAD_COLOR = "\033[96m"
ENZ_COLOR = "\033[94m"
OK_COLOR = "\033[92m"
STOP_COLOR = "\033[41m\033[97m"

COMP_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}
COMP_RNA = {"A": "U", "T": "A", "G": "C", "C": "G"}

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

WIDTH = 60  # bases por línea al imprimir


# ---------------------------------------------------------------------------
# Utilidades de impresión
# ---------------------------------------------------------------------------
def title(txt):
    print(f"\n{BOLD}{TITLE_COLOR}{'═' * 70}\n{txt}\n{'═' * 70}{RESET}")


def subtitle(txt):
    print(f"\n{BOLD}{HEAD_COLOR}▶ {txt}{RESET}")


def enzymes(lines):
    print(f"{ENZ_COLOR}", end="")
    for l in lines:
        print(f"   ⚙ {l}")
    print(RESET, end="")


def color_seq(seq):
    """Devuelve la secuencia con cada base coloreada, envuelta en líneas."""
    out_lines = []
    for chunk_start in range(0, len(seq), WIDTH):
        chunk = seq[chunk_start:chunk_start + WIDTH]
        line = "".join(f"{BASE_COLOR.get(b, '')}{b}{RESET}" for b in chunk)
        out_lines.append(line)
    return "\n".join(out_lines)


def pause(anim, secs=0.4):
    if anim:
        time.sleep(secs)


# ---------------------------------------------------------------------------
# Lógica biológica
# ---------------------------------------------------------------------------
def random_dna(min_codons=8, max_codons=14):
    stops = ["TAA", "TAG", "TGA"]
    bases = "ATGC"
    codons = ["ATG"]
    n = random.randint(min_codons, max_codons)
    for _ in range(n):
        c = "".join(random.choice(bases) for _ in range(3))
        while c in stops:
            c = "".join(random.choice(bases) for _ in range(3))
        codons.append(c)
    codons.append(random.choice(stops))
    return "".join(codons)


def replicate(coding, template, anim=False):
    """Simula la cadena líder (continua) y la rezagada (fragmentos de Okazaki)."""
    subtitle("1. Replicación del ADN")
    enzymes([
        "Helicasa: abre la doble hélice rompiendo los puentes de hidrógeno.",
        "Primasa: sintetiza los cebadores (ARN) que inician cada fragmento.",
        "ADN polimerasa: añade nucleótidos leyendo la cadena molde en sentido 5'->3'.",
        "Ligasa: une entre sí los fragmentos de Okazaki de la cadena rezagada.",
    ])

    # Cadena líder: síntesis continua usando 'template' como molde
    leading_new = "".join(COMP_DNA[b] for b in template)
    print(f"\n{DIM}Cadena molde 1 (3'->5'):{RESET}")
    print(color_seq(template))
    pause(anim)
    print(f"\n{OK_COLOR}Cadena líder sintetizada (continua, 5'->3'):{RESET}")
    print(color_seq(leading_new))
    pause(anim)

    # Cadena rezagada: fragmentos de Okazaki usando 'coding' como molde
    lagging_new = "".join(COMP_DNA[b] for b in coding)
    frag_size = 6
    print(f"\n{DIM}Cadena molde 2 (5'->3', dirección opuesta):{RESET}")
    print(color_seq(coding))
    pause(anim)
    print(f"\n{OK_COLOR}Cadena rezagada - fragmentos de Okazaki (cebador + fragmento):{RESET}")
    for i in range(0, len(lagging_new), frag_size):
        frag = lagging_new[i:i + frag_size]
        primer = f"{DIM}[cebador]{RESET}"
        print(f"  {primer} {color_seq(frag)}")
        pause(anim, 0.15)
    print(f"{OK_COLOR}  ↳ La ligasa une todos los fragmentos en una cadena continua.{RESET}")

    return leading_new, lagging_new


def transcribe(template, anim=False):
    subtitle("2. Transcripción (ADN -> ARN)")
    enzymes([
        "ARN polimerasa: lee la cadena molde (3'->5') y sintetiza ARNm "
        "complementario (5'->3'), sustituyendo T por U.",
    ])
    mrna = "".join(COMP_RNA[b] for b in template)
    print(f"\n{DIM}Cadena molde de ADN:{RESET}")
    print(color_seq(template))
    pause(anim)
    print(f"\n{OK_COLOR}ARNm sintetizado:{RESET}")
    print(color_seq(mrna))
    return mrna


def translate(mrna, anim=False):
    subtitle("3. Traducción (ARN -> Proteína)")
    enzymes([
        "Ribosoma: lee el ARNm en codones (grupos de 3 bases).",
        "ARNt: aporta el aminoácido complementario a cada codón, "
        "hasta un codón de STOP.",
    ])
    protein = []
    print()
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]
        aa = CODON_TABLE.get(codon, "?")
        colored_codon = "".join(f"{BASE_COLOR.get(b, '')}{b}{RESET}" for b in codon)
        if aa == "STOP":
            print(f"  {colored_codon}  ->  {STOP_COLOR} STOP {RESET}")
            pause(anim, 0.2)
            break
        print(f"  {colored_codon}  ->  {BOLD}{aa}{RESET}")
        protein.append(aa)
        pause(anim, 0.2)

    print(f"\n{OK_COLOR}Proteína resultante ({len(protein)} aminoácidos):{RESET}")
    print("  " + f" {DIM}-{RESET} ".join(protein))
    return protein


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def run(seq, anim=True):
    seq = seq.upper().strip()
    seq = "".join(ch for ch in seq if ch in "ATGC")
    if not seq or len(seq) < 3:
        print("Secuencia de ADN no válida (usa solo A, T, G, C).")
        return

    coding = seq                                        # cadena codificante 5'->3'
    template = "".join(COMP_DNA[b] for b in coding)      # cadena molde 3'->5'

    title("ADN ORIGINAL (doble hélice)")
    print(f"{DIM}Cadena codificante 5'->3':{RESET}")
    print(color_seq(coding))
    print(f"\n{DIM}Cadena molde 3'->5':{RESET}")
    print(color_seq(template))

    title("FLUJO DE LA INFORMACIÓN GENÉTICA")
    replicate(coding, template, anim)
    mrna = transcribe(template, anim)
    translate(mrna, anim)

    print(f"\n{BOLD}{TITLE_COLOR}{'═' * 70}{RESET}")
    print(f"{BOLD}Resumen del flujo:{RESET}  ADN → ADN (replicación)  |  "
          f"ADN → ARN (transcripción)  |  ARN → Proteína (traducción)")
    print(f"{BOLD}{TITLE_COLOR}{'═' * 70}{RESET}\n")


def main():
    parser = argparse.ArgumentParser(description="Simulador del dogma central de la biología molecular.")
    parser.add_argument("--seq", type=str, help="Secuencia de ADN (cadena codificante, 5'->3').")
    parser.add_argument("--random", action="store_true", help="Genera una secuencia de ADN aleatoria.")
    parser.add_argument("--no-anim", action="store_true", help="Desactiva las pausas de animación.")
    args = parser.parse_args()

    if args.random:
        seq = random_dna()
    elif args.seq:
        seq = args.seq
    else:
        print(f"{BOLD}{TITLE_COLOR}🧬 Simulador del Dogma Central de la Biología Molecular{RESET}")
        entrada = input(
            "\nIntroduce una secuencia de ADN (solo A,T,G,C) "
            "o pulsa Enter para generar una aleatoria: "
        ).strip()
        seq = entrada if entrada else random_dna()

    run(seq, anim=not args.no_anim)


if __name__ == "__main__":
    main()
