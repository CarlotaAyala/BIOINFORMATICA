# Simulador del Dogma Central de la Biología Molecular

Programa de consola en Python que muestra, paso a paso y con colores, cómo fluye la información genética a partir de una secuencia de ADN:

```
ADN → ADN (replicación)  |  ADN → ARN (transcripción)  |  ARN → Proteína (traducción)
```

Cada base se pinta de un color distinto, y en cada fase se indican las enzimas o moléculas que intervienen.

## Características

- **Replicación**: cadena líder continua y cadena rezagada con fragmentos de Okazaki y cebadores. Explica la función de helicasa, primasa, ADN polimerasa y ligasa.
- **Transcripción**: obtiene el ARNm a partir de la cadena molde (sustituyendo T por U), con la ARN polimerasa.
- **Traducción**: lee el ARNm en codones, muestra el aminoácido de cada uno y se detiene en el codón STOP.
- **Entrada flexible**: secuencia introducida por argumento, por teclado o generada al azar.
- **Animación opcional**: pausas entre pasos que se pueden desactivar.

## Requisitos

- Python 3.8 o superior.
- Un terminal que soporte colores ANSI (la mayoría de terminales de Linux, macOS y Windows Terminal).
- No necesita instalar ninguna librería externa.

## Estructura del proyecto

```
dogma_central/
├── constantes.py   # Colores ANSI, complementos, tabla de codones, codones de parada
├── formato.py      # Clase Consola: impresión con colores y pausas
├── secuencias.py   # Clase SecuenciaADN: validación, complemento y generación aleatoria
├── procesos.py     # Clases Replicacion, Transcripcion y Traduccion
├── main.py         # Clase SimuladorDogma y punto de entrada (argparse)
└── README.md
```

### Clases

| Clase | Archivo | Responsabilidad |
|---|---|---|
| `Consola` | `formato.py` | Todo lo visual: títulos, secuencias coloreadas, codones, pausas de animación. |
| `SecuenciaADN` | `secuencias.py` | Guarda la cadena codificante y calcula la molde. Valida la entrada y genera secuencias aleatorias. |
| `Replicacion` | `procesos.py` | Calcula y muestra la cadena líder y los fragmentos de Okazaki. |
| `Transcripcion` | `procesos.py` | Calcula y muestra el ARNm. |
| `Traduccion` | `procesos.py` | Calcula y muestra la cadena de aminoácidos hasta el codón STOP. |
| `SimuladorDogma` | `main.py` | Orquesta los tres procesos en orden. |

Cada clase de proceso separa el **cálculo** (`calcular()`, devuelve el resultado sin imprimir) de la **presentación** (`mostrar()`, imprime por pantalla). Así la lógica se puede reutilizar y probar sin depender de la salida en consola.

## Uso

Desde la carpeta `dogma_central`:

```bash
# Modo interactivo: pide la secuencia (Enter para una aleatoria)
python main.py

# Con una secuencia concreta (cadena codificante, 5'->3')
python main.py --seq ATGGCCTTTAAATAG

# Con una secuencia aleatoria
python main.py --random

# Sin pausas de animación
python main.py --seq ATGGCCTTTAAATAG --no-anim
```

### Argumentos

| Argumento | Descripción |
|---|---|
| `--seq SECUENCIA` | Secuencia de ADN (cadena codificante, 5'→3'). Solo se aceptan A, T, G y C. |
| `--random` | Genera una secuencia aleatoria que empieza por ATG y termina en un codón de parada. |
| `--no-anim` | Desactiva las pausas de animación. |

## Ejemplo

Para `ATGGCCTTTAAATAG`:

| Etapa | Resultado |
|---|---|
| Cadena codificante | `ATGGCCTTTAAATAG` |
| Cadena molde | `TACCGGAAATTTATC` |
| ARNm | `AUGGCCUUUAAAUAG` |
| Proteína | Met - Ala - Phe - Lys |

## Uso como librería

```python
from secuencias import SecuenciaADN
from procesos import Transcripcion, Traduccion

adn = SecuenciaADN("ATGGCCTTTAAATAG")
arnm = Transcripcion(adn, None).calcular()    
proteina = Traduccion(arnm, None).calcular()  
```

## Notas

- Si la secuencia introducida tiene menos de 3 bases válidas, el programa muestra un mensaje de error y termina.
- Los caracteres que no sean A, T, G o C se descartan al limpiar la secuencia.
- La simulación de la replicación es una representación didáctica simplificada: los fragmentos de Okazaki tienen un tamaño fijo de 6 bases.
- Si un codón no está en la tabla, se muestra `?` como aminoácido.
