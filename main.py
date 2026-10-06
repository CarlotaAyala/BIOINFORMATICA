import argparse

from constantes import BOLD, DIM, RESET, TITLE_COLOR
from formato import Consola
from procesos import Replicacion, Transcripcion, Traduccion
from secuencias import SecuenciaADN


class SimuladorDogma:

    def __init__(self, adn, consola):
        self.adn = adn
        self.consola = consola

    def mostrar_adn_original(self):
        c = self.consola
        c.titulo("ADN ORIGINAL (doble hélice)")
        print(f"{DIM}Cadena codificante 5'->3':{RESET}")
        c.secuencia(self.adn.codificante)
        print(f"\n{DIM}Cadena molde 3'->5':{RESET}")
        c.secuencia(self.adn.molde)

    def ejecutar(self):
        c = self.consola
        self.mostrar_adn_original()
        c.titulo("FLUJO DE LA INFORMACIÓN GENÉTICA")
        Replicacion(self.adn, c).mostrar()
        arnm = Transcripcion(self.adn, c).mostrar()
        Traduccion(arnm, c).mostrar()

        c.linea_final()
        print(f"{BOLD}Resumen del flujo:{RESET}  ADN → ADN (replicación)  |  "
              f"ADN → ARN (transcripción)  |  ARN → Proteína (traducción)")
        c.linea_final()
        print()


def leer_argumentos():
    parser = argparse.ArgumentParser(
        description="Simulador del dogma central de la biología molecular.")
    parser.add_argument("--seq", type=str, help="Secuencia de ADN (cadena codificante, 5'->3').")
    parser.add_argument("--random", action="store_true", help="Genera una secuencia de ADN aleatoria.")
    parser.add_argument("--no-anim", action="store_true", help="Desactiva las pausas de animación.")
    return parser.parse_args()


def obtener_secuencia(args):
    if args.random:
        return SecuenciaADN.aleatoria()
    if args.seq:
        return SecuenciaADN(args.seq)
    print(f"{BOLD}{TITLE_COLOR}🧬 Simulador del Dogma Central de la Biología Molecular{RESET}")
    entrada = input("\nIntroduce una secuencia de ADN (solo A,T,G,C) "
                    "o pulsa Enter para generar una aleatoria: ").strip()
    return SecuenciaADN(entrada) if entrada else SecuenciaADN.aleatoria()


def main():
    args = leer_argumentos()
    try:
        adn = obtener_secuencia(args)
    except ValueError as e:
        print(e)
        return
    SimuladorDogma(adn, Consola(animar=not args.no_anim)).ejecutar()


if __name__ == "__main__":
    main()
