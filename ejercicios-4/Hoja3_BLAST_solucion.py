# Hoja de ejercicios 3 - BLAST (versión .py del notebook)

# # Hoja de ejercicios 3 — GBIO · BLAST con Biopython
# 
# Solución de los cuatro ejercicios siguiendo el esquema de `Ejemplos_06_BLAST.ipynb`:
# 
# 1. **Lanzar BLAST** (online con `NCBIWWW.qblast()` o local con `subprocess` + BLAST+).
# 2. **Volcar el resultado a XML** (`-outfmt 5`).
# 3. **Parsear el XML** con la API moderna `Bio.SearchIO` (resultado → hit → hsp) y calcular identidad, coberturas, e-value y bit-score.
# 
# > **Antes de empezar**: ajusta las rutas de la celda de configuración a tu máquina (BLAST+ 2.17, carpeta de bases de datos y carpeta de trabajo).

# %%
import os, sys, platform, subprocess

BLAST_BIN = "C:/Program Files/NCBI/blast-2.17.0+/bin" 
DIR_BD = "C:/Users/lucia/Downloads/Ejercicios3-Practica/BD"
DIR_TRABAJO = os.getcwd()
EMAIL = "lucia.cruz104@alu.ulpgc.es"  
HILOS = "8"                 

EXT = ".exe" if platform.system() == "Windows" else ""
def exe(programa):
    return os.path.join(BLAST_BIN, programa + EXT)

os.makedirs(DIR_TRABAJO, exist_ok=True)
os.makedirs(DIR_BD, exist_ok=True)
def ruta(nombre):
    return os.path.join(DIR_TRABAJO, nombre)

import Bio
print(Bio.__version__)

# %%
# =====================================================================
# FUNCIONES AUXILIARES (reutilizadas en todos los ejercicios)
# =====================================================================
import re, tarfile, urllib.request
from Bio import SearchIO
from Bio.Blast import NCBIWWW

NCBIWWW.email = EMAIL

import concurrent.futures

def blast_online(programa, base_datos, secuencia, archivo_xml, limite_min=15, **opciones):
    """Lanza BLAST en el NCBI; si tarda más de 'limite_min' minutos, se rinde y sigue."""
    print(f"Enviando {programa} contra '{base_datos}' al NCBI (límite {limite_min} min)...")
    ejecutor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    futuro = ejecutor.submit(NCBIWWW.qblast, programa, base_datos, secuencia, **opciones)
    try:
        resultado_web = futuro.result(timeout=limite_min * 60)
    except concurrent.futures.TimeoutError:
        print(f"NCBI no ha respondido en {limite_min} min (cola saturada). Usa el resultado LOCAL.")
        ejecutor.shutdown(wait=False)
        return False
    with open(archivo_xml, "w") as salida:
        salida.write(resultado_web.read())
    resultado_web.close()
    ejecutor.shutdown(wait=False)
    print(f"XML guardado en: {archivo_xml}")
    return True

def blast_local(programa, archivo_fasta, base_datos, archivo_xml, extra=()):
    """Lanza BLAST+ en local (subprocess) con salida XML (-outfmt 5). Si falla, muestra el error de BLAST."""
    comando = [exe(programa),
               "-query", archivo_fasta,
               "-db", base_datos,
               "-out", archivo_xml,
               "-num_threads", HILOS,
               "-outfmt", "5", *extra]
    print("Ejecutando:", " ".join(comando))
    r = subprocess.run(comando, capture_output=True, text=True)
    if r.returncode != 0:
        print("BLAST ha fallado. Mensaje de BLAST:\n", r.stderr)
        return False
    print(f"XML guardado en: {archivo_xml}")
    return True

def descargar_bd_ncbi(nombre_bd, destino=DIR_BD):
    """Descarga y descomprime una BD preformateada de https://ftp.ncbi.nlm.nih.gov/blast/db/
    (equivale a 'update_blastdb.pl --decompress nombre_bd'). Soporta BD multivolumen."""
    base = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
    listado = urllib.request.urlopen(base).read().decode()
    patron = rf'href="({re.escape(nombre_bd)}(?:\.\d+)?\.tar\.gz)"'
    ficheros = sorted(set(re.findall(patron, listado)))
    if not ficheros:
        raise ValueError(f"No encuentro '{nombre_bd}' en {base}")
    carpeta = os.path.join(destino, nombre_bd)
    os.makedirs(carpeta, exist_ok=True)
    for f in ficheros:
        local = os.path.join(carpeta, f)
        if not os.path.exists(local):
            print("Descargando", f, "...")
            urllib.request.urlretrieve(base + f, local)
        with tarfile.open(local) as tar:
            tar.extractall(carpeta)
    print(f"BD '{nombre_bd}' lista en {carpeta} ({len(ficheros)} volumen/es)")
    return os.path.join(carpeta, nombre_bd)   # ruta que se pasa a -db

def resumen_hsp(query_result, hit, hsp):
    """Métricas de la presentación: identidad, coberturas, huecos, e-value, bit-score."""
    return {
        "identidad_pct": 100 * hsp.ident_num / hsp.aln_span,
        "cob_query_pct": 100 * hsp.query_span / query_result.seq_len,
        "cob_subject_pct": 100 * hsp.hit_span / hit.seq_len,
        "huecos_pct": 100 * hsp.gap_num / hsp.aln_span,
        "evalue": hsp.evalue,
        "bitscore": hsp.bitscore,
    }

def tabla_hits(query_result, n=10, ancho=60):
    """Imprime los n mejores hits (mejor HSP de cada uno)."""
    print(f"Query: {query_result.id}  ({query_result.seq_len} posiciones) — {len(query_result)} hits")
    print(f"{'#':>2}  {'Accesión':<16} {'Ident%':>6} {'CobQ%':>6} {'E-value':>9} {'Bits':>7}  Descripción")
    for i, hit in enumerate(query_result[:n], 1):
        hsp = hit[0]
        m = resumen_hsp(query_result, hit, hsp)
        print(f"{i:>2}  {hit.accession or hit.id:<16} {m['identidad_pct']:6.1f} {m['cob_query_pct']:6.1f} "
              f"{m['evalue']:9.2e} {m['bitscore']:7.1f}  {hit.description[:ancho]}")

def ver_alineamiento(hsp):
    print(f"Query:  {hsp.query.seq}")
    print(f"Match:  {hsp.aln_annotation.get('similarity', '')}")
    print(f"Sbjct:  {hsp.hit.seq}")

# ---
# ## Ejercicio 1 · Alerta epidemiológica (16S rRNA, online y local)
# 
# **Planteamiento.** La muestra de LCR es un fragmento del **gen del ARN ribosomal 16S** (ADN) amplificado por PCR. Es nucleótido contra nucleótido → **`blastn`**. La base de datos adecuada no es `nt` (enorme y con mucho "ruido"), sino **`16S_ribosomal_RNA`** (secuencias 16S de cepas tipo, curadas por NCBI), que es la que pide el enunciado.
# 
# - **Online**: `NCBIWWW.qblast("blastn", "16S_ribosomal_RNA", ...)`.
# - **Local**: descargamos `16S_ribosomal_RNA` preformateada del FTP del NCBI y lanzamos `blastn` con `subprocess`.
# 
# Como el *query* es muy corto (56 nt) usamos `-task blastn` (palabra de 11) en local, igual que hace la web con el programa *blastn*, en lugar del *megablast* por defecto.

# %%
# Muestra clínica de líquido cefalorraquídeo (LCR)
secuencia_paciente_lcr = " GACGGCAGCACAGAGAAGCTTGCTTCTCGGGTGGCGAGTGGCGAACGGGTGAGTAA"
secuencia_paciente_lcr = secuencia_paciente_lcr.strip().upper()   # ¡ojo al espacio inicial del enunciado!
print(len(secuencia_paciente_lcr), "nt")

fasta_lcr = ruta("paciente_lcr.fasta")
with open(fasta_lcr, "w") as f:
    f.write(">paciente_LCR_16S\n" + secuencia_paciente_lcr + "\n")

# %%
# ---------- 1A. ONLINE ----------
xml_e1_online = ruta("e1_16S_online.xml")
blast_online("blastn", "16S_ribosomal_RNA", secuencia_paciente_lcr, xml_e1_online,
             hitlist_size=20, expect=1e-5)

# %%
# ---------- 1B. LOCAL ----------
# (Solo la primera vez) descarga de la BD 16S preformateada (~ decenas de MB)
bd_16S = descargar_bd_ncbi("16S_ribosomal_RNA")
# Si ya la tienes: bd_16S = os.path.join(DIR_BD, "16S_ribosomal_RNA", "16S_ribosomal_RNA")

xml_e1_local = ruta("e1_16S_local.xml")
blast_local("blastn", fasta_lcr, bd_16S, xml_e1_local,
            extra=["-task", "blastn", "-evalue", "1e-5", "-max_target_seqs", "20"])

# %%
# ---------- 1C. INTERPRETACIÓN (vale para el XML online y el local) ----------
from collections import Counter

# Rutas de los XML (así esta celda funciona aunque no se haya ejecutado la parte online)
xml_e1_online = ruta("e1_16S_online.xml")
xml_e1_local = ruta("e1_16S_local.xml")

# Principales agentes de meningitis bacteriana
PATOGENOS_MENINGITIS = {
    "Neisseria meningitidis":   "MENINGOCOCO — meningitis/sepsis meningocócica (petequias). EMERGENCIA VITAL + aviso de salud pública (profilaxis a contactos).",
    "Streptococcus pneumoniae": "NEUMOCOCO — meningitis neumocócica. EMERGENCIA VITAL.",
    "Haemophilus influenzae":   "H. influenzae tipo b — meningitis. EMERGENCIA VITAL.",
    "Listeria monocytogenes":   "LISTERIA — meningitis (neonatos/inmunodeprimidos). EMERGENCIA VITAL.",
    "Streptococcus agalactiae": "ESTREPTOCOCO GRUPO B — meningitis neonatal. EMERGENCIA VITAL.",
    "Escherichia coli":         "E. coli K1 — meningitis neonatal. EMERGENCIA VITAL.",
}

def especie(descripcion):
    # 'Neisseria meningitidis strain ATCC 13077 16S ribosomal RNA, partial sequence' -> 'Neisseria meningitidis'
    return " ".join(descripcion.split()[:2])

for etiqueta, xml in [("ONLINE", xml_e1_online), ("LOCAL", xml_e1_local)]:
    if not os.path.exists(xml):
        continue
    print(f"\n================ RESULTADO {etiqueta} ================")
    res = SearchIO.read(xml, "blast-xml")
    if len(res) == 0:
        print("Sin coincidencias: posible contaminación o bacteria no presente en la BD.")
        continue
    tabla_hits(res, n=10)

    mejor = res[0]; hsp = mejor[0]
    print("\nMejor alineamiento:")
    ver_alineamiento(hsp)

    # Todos los hits empatados con la mejor puntuación (fragmento corto => empates)
    empatados = [h for h in res if h[0].bitscore == hsp.bitscore]
    print(f"\nHits empatados en la mejor puntuación: {len(empatados)}")
    print(Counter(especie(h.description) for h in empatados).most_common())
    generos = Counter(especie(h.description).split()[0] for h in empatados)
    print("Género(s) de los mejores hits:", dict(generos))

    sp = especie(mejor.description)
    m = resumen_hsp(res, mejor, hsp)
    print(f"\n>>> Identificación: {sp}  (identidad {m['identidad_pct']:.1f}%, "
          f"cobertura {m['cob_query_pct']:.1f}%, E-value {m['evalue']:.1e})")
    alerta = [v for k, v in PATOGENOS_MENINGITIS.items() if k in {especie(h.description) for h in empatados}]
    if alerta:
        print("🚨 ALERTA:", *alerta, sep="\n   ")
    else:
        print("Ningún patógeno meníngeo clásico entre los mejores hits.")

# ### Resultado y conclusión del ejercicio 1
# 
# **Nota sobre la búsqueda online (1A).** La petición a NCBI (`NCBIWWW.qblast` contra `16S_ribosomal_RNA`) no terminó: superó el límite de 15 minutos que añadí a `blast_online` (Biopython avisó: *"BLAST request GE96TCJN016 is taking longer than 10 minutes"*). La cola de esa base de datos estaba saturada: en otra prueba, la misma consulta llevaba más de 24 h en estado *Searching*. Por eso el resultado se obtuvo con el **BLAST local (1B)** contra la **misma base de datos**, que terminó en segundos.
# 
# **Resultado obtenido (local, `blastn -task blastn`, 20 hits):**
# 
# | # | Accesión | Especie | Identidad | Cobertura | E-value | Bits |
# |---|---|---|---|---|---|---|
# | 1 | NR_104946 | ***Neisseria meningitidis*** M1027 | **100 % (56/56)** | 100 % | 2,1e-22 | 102,3 |
# | 2-10 | NR_026079, NR_121688… | *N. gonorrhoeae, N. sicca, N. cinerea, N. mucosa, N. macacae, N. perflava, N. weaveri, N. dumasiana* | 98,2 % (55/56) | 100 % | 9,0e-21 | 97,8 |
# 
# - **Los 10 primeros hits son todos del género *Neisseria***, y el único con identidad del 100 % es ***N. meningitidis*** (meningococo). No hay empates en la mejor puntuación.
# - Las demás *Neisseria* solo difieren en **1 nucleótido de 56**. El fragmento basta para identificar el género con total seguridad. La especie la marca esa única diferencia y, sobre todo, el cuadro clínico.
# 
# **Conclusión: SÍ, es una emergencia vital.**
# - ADN de *N. meningitidis* en un LCR (que es estéril) y un cuadro de **rigidez de nuca, fiebre alta y petequias** indican una **meningitis o sepsis meningocócica**. Se trata con antibiótico empírico inmediato (cefalosporina de 3ª generación), sin esperar a más pruebas.
# - Es también una **alerta epidemiológica**: es una enfermedad de declaración obligatoria y hay que dar **quimioprofilaxis a los contactos estrechos**.
# - La especie y el serogrupo se confirman con PCR específica (*ctrA*/*sodC*) y cultivo, que también da el antibiograma. BLAST orienta muy rápido, pero el diagnóstico no se basa solo en un fragmento de 56 nt.
# 
# **Online frente a local:** con la misma base de datos, los dos modos deberían dar los mismos hits. Aquí se ve una ventaja práctica del modo local: **no depende de la carga de un servidor externo**. En un hospital tiene además otra: **los datos del paciente no salen de la red del centro**.

# ---
# ## Ejercicio 2 · La caza de enzimas industriales (blastx)
# 
# **¿Qué significa no conocer la dirección ni el marco de lectura?** El ADN es de doble cadena y el código genético se lee en tripletes, así que un fragmento crudo puede traducirse de **6 maneras** (+1, +2, +3 en la cadena directa y −1, −2, −3 en la complementaria inversa). Solo una (como mucho) corresponde a la proteína real; las otras cinco son "basura" de traducción.
# 
# Primero lo comprobamos traduciendo los 6 marcos con Biopython. Después usamos **`blastx`**, que hace exactamente eso de forma automática: **traduce el ADN de consulta en los 6 marcos y compara cada traducción contra una base de datos de proteínas**. Usamos **Swiss-Prot** (curada: los nombres de proteína son fiables) — es la misma BD del ejemplo de clase. El HSP nos dirá además **qué marco** era el bueno (`hsp.query_frame`).

# %%
from Bio.Seq import Seq

# ADN ambiental extraído de una fuente termal (Secuencia codificante cruda)
adn_geiser = "CTCTCCTGGGACCTGGCCAAGGTGCGCACCGACCTGCCCCTGGAGGTGGACTTCGCCAAAAGGCGGGAGC"
adn = Seq(adn_geiser)
gc = 100 * (adn_geiser.count("G") + adn_geiser.count("C")) / len(adn_geiser)
print(f"Longitud: {len(adn)} nt   GC: {gc:.1f}%  (muy alto: típico de termófilos como Thermus)\n")

marcos = {}
for i in range(3):
    fragmento = adn[i:]; fragmento = fragmento[:len(fragmento) // 3 * 3]
    marcos[f"+{i+1}"] = fragmento.translate()
for i in range(3):
    fragmento = adn.reverse_complement()[i:]; fragmento = fragmento[:len(fragmento) // 3 * 3]
    marcos[f"-{i+1}"] = fragmento.translate()

for marco, prot in marcos.items():
    print(f"Marco {marco}: {prot}   (codones STOP: {prot.count('*')})")
print("\nNinguno tiene codones STOP => mirando solo los STOP no sabemos cuál es el bueno. Necesitamos blastx.")

# %%
fasta_geiser = ruta("adn_geiser.fasta")
with open(fasta_geiser, "w") as f:
    f.write(">adn_geiser_Yellowstone\n" + adn_geiser + "\n")

# ---------- 2A. ONLINE ----------
xml_e2_online = ruta("e2_blastx_online.xml")
blast_online("blastx", "swissprot", adn_geiser, xml_e2_online, hitlist_size=10)

# %%
# ---------- 2B. LOCAL (misma BD swissprot del ejemplo de clase) ----------
fasta_geiser = ruta("adn_geiser.fasta")
with open(fasta_geiser, "w") as f:
    f.write(">adn_geiser_Yellowstone\n" + adn_geiser + "\n")

bd_swissprot = descargar_bd_ncbi("swissprot")   # descarga (solo la 1ª vez) y devuelve la ruta para -db
xml_e2_local = ruta("e2_blastx_local.xml")
blast_local("blastx", fasta_geiser, bd_swissprot, xml_e2_local,
            extra=["-evalue", "1e-3", "-max_target_seqs", "10"])

# %%
def nombre_proteina(descripcion):
    """'RecName: Full=DNA polymerase I, thermostable; AltName: ...' -> 'DNA polymerase I, thermostable'"""
    m = re.search(r"Full=([^;\[>]+)", descripcion)
    return m.group(1).strip() if m else descripcion

xml_e2_online = ruta("e2_blastx_online.xml")
xml_e2_local = ruta("e2_blastx_local.xml")

for etiqueta, xml in [("ONLINE", xml_e2_online), ("LOCAL", xml_e2_local)]:
    if not os.path.exists(xml):
        continue
    print(f"\n================ RESULTADO {etiqueta} ================")
    res = SearchIO.read(xml, "blast-xml")
    tabla_hits(res, n=10)
    mejor = res[0]; hsp = mejor[0]
    m = resumen_hsp(res, mejor, hsp)
    print(f"\nMarco de lectura correcto: {hsp.query_frame:+d}")
    ver_alineamiento(hsp)
    print(f"\n>>> PROTEÍNA DESCUBIERTA: {nombre_proteina(mejor.description)}")
    print(f"    Accesión: {mejor.accession} | identidad {m['identidad_pct']:.1f}% | "
          f"E-value {m['evalue']:.1e} | bit-score {m['bitscore']:.1f}")

# ### Resultado y explicación de la búsqueda
# 
# **Nota sobre la búsqueda online (2A).** La petición `blastx` contra Swiss-Prot en NCBI superó el límite de 15 minutos por saturación del servidor. El resultado se obtuvo con **BLAST local contra la misma base de datos (`swissprot`)**.
# 
# **Qué significa no conocer dirección ni marco:** el ADN tiene dos hebras y se lee en tripletes, así que hay **6 marcos de lectura** posibles (+1, +2, +3, −1, −2, −3). Al traducirlos, **ninguno tiene codones STOP**, algo típico en un ADN tan rico en GC (68,6 %). Por tanto, no se puede elegir el marco "a ojo".
# 
# **Búsqueda:** **`blastx`** traduce el ADN en los 6 marcos y compara cada traducción con proteínas conocidas. El marco correcto es el único que da un alineamiento significativo.
# 
# **Resultado obtenido:**
# - Marco correcto: **+1** → `LSWDLAKVRTDLPLEVDFAKRRE`.
# - Único hit con E < 1e-3: **P19821 · `DNA polymerase I, thermostable`**, la **Taq polimerasa** de ***Thermus aquaticus***. Identidad **100 % (23/23 aa)**, cobertura 98,6 % (69 de 70 nt: el último nucleótido no completa un codón), E-value **2,6e-08**, 48,9 bits.
# 
# ```
# >>> PROTEÍNA DESCUBIERTA: DNA polymerase I, thermostable   (Taq DNA polimerasa, Thermus aquaticus)
# ```
# 
# **Comentarios:**
# - El E-value (2,6e-08) es mucho más alto que en los otros ejercicios aunque la identidad sea del 100 %. La razón es que la consulta es muy corta (23 aa): con pocos residuos la puntuación máxima posible es baja. **El E-value depende de la longitud**, no solo de la identidad.
# - Las DNA polimerasas I de otros *Thermus* (*T. thermophilus*, *T. filiformis*…) se parecen un 68-73 % en este tramo, con E-values de 0,003-0,008. Quedan fuera por el umbral `-evalue 1e-3`. Con `-evalue 0.01` aparecerían.
# - *Thermus aquaticus* se aisló precisamente en las fuentes termales de **Yellowstone** (Brock, 1969). Su ADN polimerasa resiste 95 °C y por eso hizo posible la **PCR**. Es una enzima industrial termoestable de enorme valor, aunque no sea propia de detergentes (para eso se buscan proteasas, lipasas o amilasas). El alto contenido en GC encaja con un termófilo.

# ---
# ## Ejercicio 3 · Proteínas célebres contra el genoma humano (local, un solo BLAST)
# 
# ### Las tres proteínas elegidas
# 
# | Proteína | UniProt | Gen | Locus (GRCh38) | Long. | Función |
# |---|---|---|---|---|---|
# | **Hemoglobina** (subunidad β) | P68871 (HBB_HUMAN) | *HBB* | 11p15.4 | 147 aa | Junto a la cadena α forma el tetrámero α2β2 que transporta O₂ en los eritrocitos. Mutaciones: anemia falciforme (Glu6Val) y β-talasemias. |
# | **Mioglobina** | P02144 (MYG_HUMAN) | *MB* | 22q12.3 | 154 aa | Hemoproteína monomérica que almacena O₂ en músculo esquelético y cardíaco. Pariente lejana de las globinas. |
# | **Oxitocina** (precursor oxitocina-neurofisina 1) | P01178 (NEU1_HUMAN) | *OXT* | 20p13 | 125 aa | Se sintetiza como prohormona: péptido señal + oxitocina (9 aa, CYIQNCPLG) + neurofisina 1. Hormona del parto, la lactancia y el vínculo social. |
# 
# Las secuencias FASTA se han obtenido de UniProt (`https://rest.uniprot.org/uniprotkb/<ID>.fasta`).
# 
# ### Estrategia
# - Consulta = **proteína**, base de datos = **genoma (nucleótidos)** → el programa es **`tblastn`** (traduce el genoma en los 6 marcos sobre la marcha).
# - **Un solo lanzamiento de BLAST**: las tres proteínas van juntas en un **multi-FASTA**; BLAST+ procesa todas las consultas en una única ejecución y genera **un único XML** con un *QueryResult* por proteína, que se recorre con `SearchIO.parse()` (como en la última celda del ejemplo de clase).
# - **Genoma humano en local**: BD preformateada del NCBI `human_genome` (GRCh38, cromosomas). Ocupa varios GB: descárgala una vez. Alternativa: descargar el FASTA de GRCh38 y construirla con `makeblastdb` (celda opcional).

# %%
# Multi-FASTA con las tres proteínas (secuencias de UniProt)
proteinas_fasta = """>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH
>sp|P02144|MYG_HUMAN Myoglobin OS=Homo sapiens OX=9606 GN=MB PE=1 SV=2
MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE
DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH
PGDFGADAQGAMNKALELFRKDMASNYKELGFQG
>sp|P01178|NEU1_HUMAN Oxytocin-neurophysin 1 proprotein OS=Homo sapiens OX=9606 GN=OXT PE=1 SV=1
MAGPSLACCLLGLLALTSACYIQNCPLGGKRAAPDLDVRKCLPCGPGGKGRCFGPNICCA
EELGCFVGTAEALRCQEENYLPSPCQSGQKACGSGGRCAVLGLCCSPDGCHADPACDAEA
TFSQR
"""
fasta_celebres = ruta("proteinas_celebres.fasta")
with open(fasta_celebres, "w") as f:
    f.write(proteinas_fasta)

# (Alternativa) descargarlas directamente de UniProt:
# for acc in ["P68871", "P02144", "P01178"]:
#     proteinas_fasta += urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").read().decode()

from Bio import SeqIO
for r in SeqIO.parse(fasta_celebres, "fasta"):
    print(f"{r.id:<24} {len(r.seq):>4} aa")

# %%
# ---------- DESCARGA DEL GENOMA HUMANO (solo la primera vez; varios GB) ----------
# Opción A (recomendada): BD preformateada del NCBI, multivolumen (human_genome.00.tar.gz, .01, ...)
bd_humano = descargar_bd_ncbi("human_genome")
# Ojo: por dentro la BD se llama GCF_000001405.xx_top_level; la celda siguiente busca ese nombre.

# %%
# Opción B (alternativa): FASTA de GRCh38 + makeblastdb (como en la "nota de interés" de la presentación)
# import gzip, shutil
# url = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
#        "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
# gz = os.path.join(DIR_BD, "GRCh38.fna.gz"); fna = os.path.join(DIR_BD, "GRCh38.fna")
# urllib.request.urlretrieve(url, gz)
# with gzip.open(gz, "rb") as src, open(fna, "wb") as dst:
#     shutil.copyfileobj(src, dst)
# bd_humano = os.path.join(DIR_BD, "GRCh38", "GRCh38")
# subprocess.run([exe("makeblastdb"), "-in", fna, "-dbtype", "nucl", "-parse_seqids",
#                 "-title", "GRCh38", "-out", bd_humano], check=True)

# %%
# ---------- UN ÚNICO LANZAMIENTO DE BLAST PARA LAS TRES PROTEÍNAS ----------

# 1) Buscar el nombre REAL de la base de datos dentro de la carpeta human_genome
carpeta_humano = os.path.join(DIR_BD, "human_genome")
archivos = os.listdir(carpeta_humano)
print(len(archivos), "archivos en", carpeta_humano)

nal = [f for f in archivos if f.endswith(".nal")]
print("Archivos .nal:", nal)
if not nal:
    raise FileNotFoundError("No hay archivo .nal: la descarga del genoma no está completa. "
                            "Vuelve a ejecutar descargar_bd_ncbi('human_genome').")

bd_humano = os.path.join(carpeta_humano, nal[0][:-4])   # p. ej. ...\GCF_000001405.39_top_level
print("Base de datos que se usará:", bd_humano)

# 2) Comprobar que BLAST la puede abrir
info = subprocess.run([exe("blastdbcmd"), "-db", bd_humano, "-info"], capture_output=True, text=True)
print(info.stdout[:400] if info.returncode == 0 else "❌ " + info.stderr)

# 3) Lanzar tblastn (puede tardar bastante: déjalo trabajar)
xml_e3 = ruta("e3_tblastn_genoma_humano.xml")
r = subprocess.run([exe("tblastn"),
                    "-query", fasta_celebres,
                    "-db", bd_humano,
                    "-out", xml_e3,
                    "-num_threads", "4",
                    "-outfmt", "5",
                    "-evalue", "1e-5", "-max_target_seqs", "10", "-max_hsps", "50"],
                   capture_output=True, text=True)
if r.returncode == 0:
    print("✅ XML guardado en:", xml_e3)
else:
    print("❌ BLAST ha fallado. Mensaje de BLAST:\n", r.stderr)

# %%
# ---------- ANÁLISIS: un QueryResult por proteína ----------
# En una BD de genoma cada "hit" es un CROMOSOMA entero, y cada HSP es un trozo del gen.
# Agrupamos los HSP en "loci" (modelos de gen aproximados) para reconocer genes y sus exones.

def loci_de_hit(hit, distancia_max=12_000, tolerancia=25):
    """Agrupa los HSP de un cromosoma en 'modelos de gen' aproximados.
    Dos HSP son exones del mismo gen si: están en la misma hebra, a menos de 12 kb,
    y son COLINEALES (el orden en el cromosoma coincide con el orden en la proteína;
    en la hebra - va al revés). Se toleran ~25 aa de solape porque los HSP de tblastn
    suelen 'invadir' unos codones del intrón en los bordes de cada exón."""
    hsps = sorted(hit.hsps, key=lambda h: h.hit_start)
    grupos = []
    for h in hsps:
        mas = h.hit_frame > 0
        destino = None
        for g in reversed(grupos):                      # primero el locus más cercano
            if (g[0].hit_frame > 0) != mas:
                continue
            if h.hit_start - max(x.hit_end for x in g) > distancia_max:
                continue
            if mas:   # hebra +: el nuevo exón debe venir DESPUÉS en la proteína
                colineal = h.query_start >= max(x.query_end for x in g) - tolerancia
            else:     # hebra -: el nuevo exón debe venir ANTES en la proteína
                colineal = h.query_end <= min(x.query_start for x in g) + tolerancia
            if colineal:
                destino = g
                break
        if destino is None:
            grupos.append([h])
        else:
            destino.append(h)
    return grupos

xml_e3 = ruta("e3_tblastn_genoma_humano.xml")

for query_result in SearchIO.parse(xml_e3, "blast-xml"):
    print("\n" + "=" * 75)
    print(f"PROTEÍNA: {query_result.id}  ({query_result.seq_len} aa) — {len(query_result)} cromosomas con hits")
    for hit in query_result[:5]:
        print(f"\n  ▸ {hit.id}  {hit.description[:55]}")
        for locus in sorted(loci_de_hit(hit), key=lambda g: min(h.evalue for h in g)):
            ident = sum(h.ident_num for h in locus); aln = sum(h.aln_span for h in locus)
            q_cubierto = set()
            for h in locus:
                q_cubierto.update(range(h.query_start, h.query_end))
            print(f"     locus {min(h.hit_start for h in locus):>11,}-{max(h.hit_end for h in locus):<11,} "
                  f"HSPs(exones)={len(locus):>2}  identidad={100*ident/aln:5.1f}%  "
                  f"cobertura query={100*len(q_cubierto)/query_result.seq_len:5.1f}%  "
                  f"mejor E={min(h.evalue for h in locus):.1e}  hebra={'+' if locus[0].hit_frame > 0 else '-'}")
            for h in sorted(locus, key=lambda h: h.query_start):
                print(f"        aa {h.query_start+1:>3}-{h.query_end:<3} -> {h.hit_start+1:,}-{h.hit_end:,} "
                      f"(marco {h.hit_frame:+d}, {100*h.ident_num/h.aln_span:.0f}% id)")

# ### Resultados e interpretación
# 
# **Un solo BLAST, tres respuestas.** El XML contiene tres *QueryResult*, uno por proteína, y `SearchIO.parse()` los recorre en orden. Así se cumple la condición del enunciado: **una única ejecución de `tblastn`** para las tres proteínas. BLAST se lanzó contra **GRCh38.p13** (`GCF_000001405.39_top_level`: 639 secuencias, 3.272 Mb).
# 
# **Cómo leer la salida contra un genoma:**
# - Cada **hit** es un **cromosoma entero** (p. ej. `NC_000011.10`), por eso la cobertura del *subject* no tiene sentido aquí.
# - Los genes humanos tienen **intrones**, así que cada proteína aparece como **varios HSP, uno por exón**, separados en el cromosoma y en marcos distintos (cada intrón desplaza la fase). La función `loci_de_hit` agrupa los HSP **colineales** (en el mismo orden en el cromosoma y en la proteína) en "modelos de gen".
# - Los HSP de `tblastn` suelen **invadir unos codones del intrón** en los bordes de cada exón. Por eso algunos HSP solapan unos aminoácidos y tienen identidades algo menores del 100 %.
# - Los **parálogos** (copias duplicadas de un gen ancestral) aparecen como loci adicionales con menor identidad.
# 
# **Resultado obtenido:**
# 
# | Proteína | Cromosoma | Locus principal | Otros loci (parálogos) |
# |---|---|---|---|
# | **Hemoglobina β** | **11** (`NC_000011.10`), hebra − | ***HBB*** 5.225.601-5.227.021: **3 HSP = 3 exones** (aa 1-45 / 31-106 / 105-147), cobertura 100 %, E = 5e-44 | Todo el **clúster β-globina** del mismo cromosoma: *HBD* (δ, 5,23 Mb, 3 exones), pseudogén *HBBP1* (5,24 Mb, solo se detecta 1 exón), *HBG1* y *HBG2* (γ fetales, 5,25 Mb) y *HBE1* (ε embrionaria, 5,27 Mb). En el **cromosoma 16**: el **clúster α-globina** (153-181 kb: *HBA2*, *HBA1*, *HBQ1*…), con un 32-37 % de identidad. |
# | **Mioglobina** | **22** (`NC_000022.11`), hebra − | ***MB*** 35.607.300-35.617.257: **3 HSP = 3 exones codificantes** (aa 1-32 / 32-110 / 93-154), cobertura 100 %, E = 7e-41 | **Ninguno** con E < 1e-5: la mioglobina se separó de las hemoglobinas hace tanto tiempo que su parecido actual (~25 %) no supera el umbral. |
# | **Oxitocina** | **20** (`NC_000020.11`), hebra + | ***OXT*** 3.071.656-3.072.415: **2 HSP**, cobertura 100 %, E = 5e-43. El gen tiene 3 exones, pero sus intrones son muy cortos y, probablemente, `tblastn` une los exones 2 y 3 en un solo HSP atravesando el intrón (por eso su identidad baja al 75 %). | ***AVP*** (vasopresina) a solo **~11 kb** y en la **hebra opuesta** (3.082.981-3.083.181). Alinea el tramo de la **neurofisina** (aa 40-106) con un 90 % de identidad. |
# 
# **Conclusiones:**
# 1. `tblastn` es el programa adecuado cuando la consulta es una proteína y la BD es ADN: **localiza el gen** de una proteína en un genoma sin anotar.
# 2. El número de HSP colineales del locus principal **dibuja la estructura exón-intrón** del gen: 3 exones para *HBB* y *MB*.
# 3. Los loci secundarios muestran la **evolución por duplicación génica**. El clúster β está ordenado en el cromosoma igual que se expresa durante el desarrollo: ε (embrión) → γ (feto) → δ y β (adulto). Oxitocina y vasopresina son una duplicación en tándem de una hormona ancestral.
# 4. Del pseudogén *HBBP1* solo se detecta un exón: es una copia que perdió su función y acumula mutaciones, y el resto ya no se parece lo suficiente.
# 5. Con el genoma en local se pueden lanzar todas las consultas necesarias sin depender de NCBI, a cambio de descargar varios GB.

# ---
# ## Ejercicio 4 · El árbol de la vida (blastp)
# 
# Consulta: **hemoglobina α humana** (UniProt **P69905**, gen *HBA1*, 142 aa). Proteína contra proteína → **`blastp`**.
# - **Online contra `nr`** (*non-redundant*, la mayor BD de proteínas del NCBI), pidiendo los **10 mejores** (`hitlist_size=10`).
# - Como la búsqueda online no terminó (ver nota), la alternativa es **en local contra Swiss-Prot**, que además da **un resultado por especie** y permite ver el gradiente evolutivo. Para obtener el nombre de cada especie se usa la base de taxonomía `taxdb`.

# %%
hba_humana = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
              "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
              "AVHASLDKFLASVSTVLTSKYR")

# 1) blastp contra nr — los 10 mejores
xml_e4 = ruta("e4_blastp_nr_top10.xml")
blast_online("blastp", "nr", hba_humana, xml_e4, hitlist_size=10)

# **Nota sobre la búsqueda online.** La petición `blastp` contra `nr` en NCBI superó el límite de 15 minutos (*"BLAST request GERZV65P014 is taking longer than 10 minutes"*). `nr` es la mayor base de datos de proteínas del NCBI y la cola estaba saturada. El ejercicio se resolvió en **local contra Swiss-Prot** (celda siguiente).

# %%
# Ejercicio 4 en LOCAL: blastp de la HBA humana contra Swiss-Prot
hba_humana = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
              "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
              "AVHASLDKFLASVSTVLTSKYR")
bd_swissprot = os.path.join(DIR_BD, "swissprot", "swissprot")   # descargada en el ejercicio 2

fasta_hba = ruta("hba_humana.fasta")
with open(fasta_hba, "w") as f:
    f.write(">HBA_HUMAN\n" + hba_humana + "\n")

# Base de datos de taxonomía, para ver el nombre de cada especie (pequeña, se descarga una vez)
carpeta_tax = os.path.dirname(descargar_bd_ncbi("taxdb"))
entorno = dict(os.environ, BLASTDB=carpeta_tax)

r = subprocess.run([exe("blastp"), "-query", fasta_hba, "-db", bd_swissprot,
                    "-evalue", "1e-5", "-max_target_seqs", "300", "-num_threads", HILOS,
                    "-outfmt", "6 sacc pident qcovs evalue bitscore sscinames stitle"],
                   capture_output=True, text=True, env=entorno)
if r.returncode != 0:
    print("❌", r.stderr)

filas = [l.split("\t") for l in r.stdout.strip().splitlines()]

print("=== 10 MEJORES HITS ===")
print(f"{'#':>2} {'Accesión':<10} {'Ident%':>6} {'CobQ%':>5} {'E-value':>9} {'Bits':>6}  Especie | Proteína")
for i, (acc, ident, cob, ev, bits, sp, desc) in enumerate(filas[:10], 1):
    print(f"{i:>2} {acc:<10} {float(ident):6.1f} {cob:>5} {float(ev):9.1e} {float(bits):6.1f}  {sp} | {desc[:40]}")

print("\n=== GRADIENTE EVOLUTIVO (un hit por especie) ===")
vistas = set()
for acc, ident, cob, ev, bits, sp, desc in filas:
    if sp not in vistas:
        vistas.add(sp)
        print(f"{float(ident):5.1f}%  {sp:<35} {desc[:45]}")

# ### Comentario de resultados
# 
# **Los 10 mejores hits (Swiss-Prot)** son **todos de primates**, con identidad del 96-100 %, cobertura del 99-100 %, E-values entre 1e-96 y 1e-100 y 276-286 bits:
# 
# | # | Especie | Identidad |
# |---|---|---|
# | 1 | *Homo sapiens*, *Pan troglodytes*, *Pan paniscus* (P69905) | 100 % |
# | 2 | Gorila | 99,3 % |
# | 3 | Gibón (*Hylobates lar*) | 98,6 % |
# | 4 | Orangután | 97,9 % |
# | 5-10 | Langur, macacos, mono araña, tití, cercopiteco | 96-98 % |
# 
# - La HBA de **humano, chimpancé y bonobo es idéntica**: Swiss-Prot las reúne en una sola entrada (P69905).
# - **El top 10 no muestra el árbol de la vida**, solo que la HBA está **casi intacta en primates**. En `nr` pasaría lo mismo, o peor: los 10 primeros serían casi todos copias humanas, estructuras PDB y grandes simios.
# - Con E-values prácticamente nulos en todos, para ordenarlos importan la **identidad** y el **bit-score**.
# 
# **El gradiente evolutivo (un hit por especie, 239 especies con `-max_target_seqs 300`):**
# 
# | Grupo | Ejemplos (identidad con la HBA humana) |
# |---|---|
# | Primates | 92-100 % |
# | Otros mamíferos placentarios | caballo 88 %, vaca 88 %, ratón 85,9 %, perro 83,7 % |
# | Marsupiales | canguro 80,9 %, zarigüeya 71,6 % |
# | Monotremas | equidna 73,8 %, ornitorrinco 72,3 % |
# | Aves | pato 71,8 %, gallina 70,4 %, paloma 69 %, avestruz 68,8 % |
# | Reptiles | cocodrilo 68,8 %, caimán 66,9 %, tortugas y serpientes 57-65 % |
# | Anfibios | *Xenopus laevis* 63,4 %, *X. tropicalis* 59,2 %, tritones 56-60 % |
# | Peces óseos | salmón 57,3 %, carpín 53,2 %, pez cebra 53,1 %, peces antárticos ~49 % |
# 
# **Comentarios:**
# 1. **La identidad baja aproximadamente en el orden en que se separaron los linajes del nuestro**: primates → mamíferos placentarios → marsupiales y monotremas → aves y reptiles → anfibios → peces. Es el **"reloj molecular"**: la hemoglobina conserva su función (los residuos que sujetan el grupo hemo, los contactos α-β) mientras acumula mutaciones neutras con el tiempo.
# 2. El reloj **no es perfecto**: algunas ardillas (~90 %) se parecen más a nosotros que algunos primates lejanos (lémures, ~85-89 %), porque cada linaje evoluciona a su propia velocidad.
# 3. Entre las otras especies aparece una proteína **humana**: la **hemoglobina θ-1** (*HBQ1*, 62 %). Es un **parálogo** de nuestro propio genoma, el mismo que salió en el cromosoma 16 en el ejercicio 3. Una búsqueda por similitud no distingue ortólogos (el mismo gen en otra especie) de parálogos (copias dentro de la misma especie).
# 4. Con 300 hits la lista llega hasta los peces óseos (~48 %). Para ver tiburones, lampreas u otras globinas (β, mioglobina) habría que subir `-max_target_seqs` (p. ej. 1000). El siguiente paso natural sería un alineamiento múltiple de estas secuencias y un árbol filogenético con `Bio.Phylo`.
# 5. **`nr` frente a Swiss-Prot:** `nr` es mucho mayor y redundante, y online tarda mucho. Swiss-Prot está revisada a mano, tiene una entrada por proteína y especie y en local responde en segundos. Para este tipo de análisis evolutivo es la mejor opción.
