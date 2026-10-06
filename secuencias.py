import random

from constantes import COMP_DNA, COMP_RNA, STOP_CODONS_DNA


class SecuenciaADN:

    def __init__(self, texto):
        self.codificante = self.limpiar(texto)
        if len(self.codificante) < 3:
            raise ValueError("Secuencia de ADN no válida (usa solo A, T, G, C).")
        self.molde = self.complementar(self.codificante)

    @staticmethod
    def limpiar(texto):
        texto = texto.upper().strip()
        return "".join(ch for ch in texto if ch in "ATGC")

    @staticmethod
    def complementar(seq, tabla=COMP_DNA):
        return "".join(tabla[b] for b in seq)

    @staticmethod
    def complementar_arn(seq):
        return SecuenciaADN.complementar(seq, COMP_RNA)

    @classmethod
    def aleatoria(cls, min_codones=8, max_codones=14):
        bases = "ATGC"
        codones = ["ATG"]
        for _ in range(random.randint(min_codones, max_codones)):
            c = "".join(random.choice(bases) for _ in range(3))
            while c in STOP_CODONS_DNA:
                c = "".join(random.choice(bases) for _ in range(3))
            codones.append(c)
        codones.append(random.choice(STOP_CODONS_DNA))
        return cls("".join(codones))