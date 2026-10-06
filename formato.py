"""Todo lo relacionado con cómo se ve la salida por pantalla."""
import time

from constantes import (BASE_COLOR, BOLD, DIM, ENZ_COLOR, HEAD_COLOR, OK_COLOR,
                        RESET, STOP_COLOR, TITLE_COLOR)


class Consola:
    """Imprime con colores y controla las pausas de animación."""

    def __init__(self, ancho=60, animar=True):
        self.ancho = ancho
        self.animar = animar

    # --- utilidades ---
    def pausa(self, segundos=0.4):
        if self.animar:
            time.sleep(segundos)

    def colorear(self, seq):
        """Devuelve la secuencia con cada base coloreada, en líneas de `ancho`."""
        lineas = []
        for inicio in range(0, len(seq), self.ancho):
            trozo = seq[inicio:inicio + self.ancho]
            lineas.append("".join(f"{BASE_COLOR.get(b, '')}{b}{RESET}" for b in trozo))
        return "\n".join(lineas)

    # --- elementos de texto ---
    def titulo(self, texto):
        print(f"\n{BOLD}{TITLE_COLOR}{'═' * 70}\n{texto}\n{'═' * 70}{RESET}")

    def subtitulo(self, texto):
        print(f"\n{BOLD}{HEAD_COLOR}▶ {texto}{RESET}")

    def enzimas(self, lineas):
        print(f"{ENZ_COLOR}", end="")
        for linea in lineas:
            print(f"   ⚙ {linea}")
        print(RESET, end="")

    def etiqueta(self, texto, tipo="dim"):
        color = DIM if tipo == "dim" else OK_COLOR
        print(f"\n{color}{texto}{RESET}")

    def secuencia(self, seq):
        print(self.colorear(seq))

    def codon(self, codon, resultado, es_stop=False):
        coloreado = self.colorear_codon(codon)
        if es_stop:
            print(f"  {coloreado}  ->  {STOP_COLOR} STOP {RESET}")
        else:
            print(f"  {coloreado}  ->  {BOLD}{resultado}{RESET}")

    @staticmethod
    def colorear_codon(codon):
        return "".join(f"{BASE_COLOR.get(b, '')}{b}{RESET}" for b in codon)

    def linea_final(self):
        print(f"\n{BOLD}{TITLE_COLOR}{'═' * 70}{RESET}")
