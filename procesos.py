from constantes import BASE_COLOR, CODON_TABLE, DIM, OK_COLOR, RESET
from secuencias import SecuenciaADN

class Replicacion:
    TAM_FRAGMENTO = 6

    def __init__(self, adn, consola):
        self.adn = adn
        self.consola = consola

    def calcular(self):
        lider = SecuenciaADN.complementar(self.adn.molde)
        rezagada = SecuenciaADN.complementar(self.adn.codificante)
        fragmentos = [rezagada[i:i + self.TAM_FRAGMENTO]
                      for i in range(0, len(rezagada), self.TAM_FRAGMENTO)]
        return lider, rezagada, fragmentos

    def mostrar(self):
        c = self.consola
        c.subtitulo("1. Replicación del ADN")
        c.enzimas([
            "Helicasa: abre la doble hélice rompiendo los puentes de hidrógeno.",
            "Primasa: sintetiza los cebadores (ARN) que inician cada fragmento.",
            "ADN polimerasa: añade nucleótidos leyendo la cadena molde en sentido 5'->3'.",
            "Ligasa: une entre sí los fragmentos de Okazaki de la cadena rezagada.",
        ])
        lider, _, fragmentos = self.calcular()

        c.etiqueta("Cadena molde 1 (3'->5'):")
        c.secuencia(self.adn.molde)
        c.pausa()
        c.etiqueta("Cadena líder sintetizada (continua, 5'->3'):", "ok")
        c.secuencia(lider)
        c.pausa()

        c.etiqueta("Cadena molde 2 (5'->3', dirección opuesta):")
        c.secuencia(self.adn.codificante)
        c.pausa()
        c.etiqueta("Cadena rezagada - fragmentos de Okazaki (cebador + fragmento):", "ok")
        for frag in fragmentos:
            print(f"  {DIM}[cebador]{RESET} {c.colorear(frag)}")
            c.pausa(0.15)
        print(f"{OK_COLOR}  ↳ La ligasa une todos los fragmentos en una cadena continua.{RESET}")


class Transcripcion:
    def __init__(self, adn, consola):
        self.adn = adn
        self.consola = consola

    def calcular(self):
        return SecuenciaADN.complementar_arn(self.adn.molde)

    def mostrar(self):
        c = self.consola
        c.subtitulo("2. Transcripción (ADN -> ARN)")
        c.enzimas([
            "ARN polimerasa: lee la cadena molde (3'->5') y sintetiza ARNm "
            "complementario (5'->3'), sustituyendo T por U.",
        ])
        arnm = self.calcular()
        c.etiqueta("Cadena molde de ADN:")
        c.secuencia(self.adn.molde)
        c.pausa()
        c.etiqueta("ARNm sintetizado:", "ok")
        c.secuencia(arnm)
        return arnm


class Traduccion:
    def __init__(self, arnm, consola):
        self.arnm = arnm
        self.consola = consola

    def codones(self):
        for i in range(0, len(self.arnm) - 2, 3):
            codon = self.arnm[i:i + 3]
            aa = CODON_TABLE.get(codon, "?")
            yield codon, aa
            if aa == "STOP":
                return

    def calcular(self):
        return [aa for _, aa in self.codones() if aa != "STOP"]

    def mostrar(self):
        c = self.consola
        c.subtitulo("3. Traducción (ARN -> Proteína)")
        c.enzimas([
            "Ribosoma: lee el ARNm en codones (grupos de 3 bases).",
            "ARNt: aporta el aminoácido complementario a cada codón, "
            "hasta un codón de STOP.",
        ])
        print()
        proteina = []
        for codon, aa in self.codones():
            if aa == "STOP":
                c.codon(codon, aa, es_stop=True)
            else:
                c.codon(codon, aa)
                proteina.append(aa)
            c.pausa(0.2)

        print(f"\n{OK_COLOR}Proteína resultante ({len(proteina)} aminoácidos):{RESET}")
        print("  " + f" {DIM}-{RESET} ".join(proteina))
        return proteina
