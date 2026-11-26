# Hoja de ejercicios 3 · BLAST — Respuestas

Resultados obtenidos al ejecutar el notebook `Hoja3_BLAST_solucion.ipynb` (BLAST+ 2.17, Biopython 1.88).

## Ejercicio 1 · Alerta epidemiológica (blastn, 16S rRNA)

**Nota sobre la búsqueda online (1A).** La petición a NCBI (`NCBIWWW.qblast` contra `16S_ribosomal_RNA`) no terminó: superó el límite de 15 minutos que añadí a `blast_online` (Biopython avisó: *"BLAST request GE96TCJN016 is taking longer than 10 minutes"*). La cola de esa base de datos estaba saturada: en otra prueba, la misma consulta llevaba más de 24 h en estado *Searching*. Por eso el resultado se obtuvo con el **BLAST local (1B)** contra la **misma base de datos**, que terminó en segundos.

**Resultado obtenido (local, `blastn -task blastn`, 20 hits):**

| # | Accesión | Especie | Identidad | Cobertura | E-value | Bits |
|---|---|---|---|---|---|---|
| 1 | NR_104946 | ***Neisseria meningitidis*** M1027 | **100 % (56/56)** | 100 % | 2,1e-22 | 102,3 |
| 2-10 | NR_026079, NR_121688… | *N. gonorrhoeae, N. sicca, N. cinerea, N. mucosa, N. macacae, N. perflava, N. weaveri, N. dumasiana* | 98,2 % (55/56) | 100 % | 9,0e-21 | 97,8 |

- **Los 10 primeros hits son todos del género *Neisseria***, y el único con identidad del 100 % es ***N. meningitidis*** (meningococo). No hay empates en la mejor puntuación.
- Las demás *Neisseria* solo difieren en **1 nucleótido de 56**. El fragmento basta para identificar el género con total seguridad. La especie la marca esa única diferencia y, sobre todo, el cuadro clínico.

**Conclusión: SÍ, es una emergencia vital.**
- ADN de *N. meningitidis* en un LCR (que es estéril) y un cuadro de **rigidez de nuca, fiebre alta y petequias** indican una **meningitis o sepsis meningocócica**. Se trata con antibiótico empírico inmediato (cefalosporina de 3ª generación), sin esperar a más pruebas.
- Es también una **alerta epidemiológica**: es una enfermedad de declaración obligatoria y hay que dar **quimioprofilaxis a los contactos estrechos**.
- La especie y el serogrupo se confirman con PCR específica (*ctrA*/*sodC*) y cultivo, que también da el antibiograma. BLAST orienta muy rápido, pero el diagnóstico no se basa solo en un fragmento de 56 nt.

**Online frente a local:** con la misma base de datos, los dos modos deberían dar los mismos hits. Aquí se ve una ventaja práctica del modo local: **no depende de la carga de un servidor externo**. En un hospital tiene además otra: **los datos del paciente no salen de la red del centro**.

## Ejercicio 2 · La caza de enzimas industriales (blastx)

**Nota sobre la búsqueda online (2A).** La petición `blastx` contra Swiss-Prot en NCBI superó el límite de 15 minutos por saturación del servidor. El resultado se obtuvo con **BLAST local contra la misma base de datos (`swissprot`)**.

**Qué significa no conocer dirección ni marco:** el ADN tiene dos hebras y se lee en tripletes, así que hay **6 marcos de lectura** posibles (+1, +2, +3, −1, −2, −3). Al traducirlos, **ninguno tiene codones STOP**, algo típico en un ADN tan rico en GC (68,6 %). Por tanto, no se puede elegir el marco "a ojo".

**Búsqueda:** **`blastx`** traduce el ADN en los 6 marcos y compara cada traducción con proteínas conocidas. El marco correcto es el único que da un alineamiento significativo.

**Resultado obtenido:**
- Marco correcto: **+1** → `LSWDLAKVRTDLPLEVDFAKRRE`.
- Único hit con E < 1e-3: **P19821 · `DNA polymerase I, thermostable`**, la **Taq polimerasa** de ***Thermus aquaticus***. Identidad **100 % (23/23 aa)**, cobertura 98,6 % (69 de 70 nt: el último nucleótido no completa un codón), E-value **2,6e-08**, 48,9 bits.

```
>>> PROTEÍNA DESCUBIERTA: DNA polymerase I, thermostable   (Taq DNA polimerasa, Thermus aquaticus)
```

**Comentarios:**
- El E-value (2,6e-08) es mucho más alto que en los otros ejercicios aunque la identidad sea del 100 %. La razón es que la consulta es muy corta (23 aa): con pocos residuos la puntuación máxima posible es baja. **El E-value depende de la longitud**, no solo de la identidad.
- Las DNA polimerasas I de otros *Thermus* (*T. thermophilus*, *T. filiformis*…) se parecen un 68-73 % en este tramo, con E-values de 0,003-0,008. Quedan fuera por el umbral `-evalue 1e-3`. Con `-evalue 0.01` aparecerían.
- *Thermus aquaticus* se aisló precisamente en las fuentes termales de **Yellowstone** (Brock, 1969). Su ADN polimerasa resiste 95 °C y por eso hizo posible la **PCR**. Es una enzima industrial termoestable de enorme valor, aunque no sea propia de detergentes (para eso se buscan proteasas, lipasas o amilasas). El alto contenido en GC encaja con un termófilo.

## Ejercicio 3 · Proteínas célebres contra el genoma humano (tblastn local, un solo BLAST)

### Proteínas elegidas

| Proteína | UniProt | Gen | Locus (GRCh38) | Long. | Función |
|---|---|---|---|---|---|
| **Hemoglobina** (subunidad β) | P68871 (HBB_HUMAN) | *HBB* | 11p15.4 | 147 aa | Junto a la cadena α forma el tetrámero α2β2 que transporta O₂ en los eritrocitos. Mutaciones: anemia falciforme (Glu6Val) y β-talasemias. |
| **Mioglobina** | P02144 (MYG_HUMAN) | *MB* | 22q12.3 | 154 aa | Hemoproteína monomérica que almacena O₂ en músculo esquelético y cardíaco. Pariente lejana de las globinas. |
| **Oxitocina** (precursor oxitocina-neurofisina 1) | P01178 (NEU1_HUMAN) | *OXT* | 20p13 | 125 aa | Se sintetiza como prohormona: péptido señal + oxitocina (9 aa, CYIQNCPLG) + neurofisina 1. Hormona del parto, la lactancia y el vínculo social. |

Las secuencias FASTA se han obtenido de UniProt (`https://rest.uniprot.org/uniprotkb/<ID>.fasta`).

**Un solo BLAST, tres respuestas.** El XML contiene tres *QueryResult*, uno por proteína, y `SearchIO.parse()` los recorre en orden. Así se cumple la condición del enunciado: **una única ejecución de `tblastn`** para las tres proteínas. BLAST se lanzó contra **GRCh38.p13** (`GCF_000001405.39_top_level`: 639 secuencias, 3.272 Mb).

**Cómo leer la salida contra un genoma:**
- Cada **hit** es un **cromosoma entero** (p. ej. `NC_000011.10`), por eso la cobertura del *subject* no tiene sentido aquí.
- Los genes humanos tienen **intrones**, así que cada proteína aparece como **varios HSP, uno por exón**, separados en el cromosoma y en marcos distintos (cada intrón desplaza la fase). La función `loci_de_hit` agrupa los HSP **colineales** (en el mismo orden en el cromosoma y en la proteína) en "modelos de gen".
- Los HSP de `tblastn` suelen **invadir unos codones del intrón** en los bordes de cada exón. Por eso algunos HSP solapan unos aminoácidos y tienen identidades algo menores del 100 %.
- Los **parálogos** (copias duplicadas de un gen ancestral) aparecen como loci adicionales con menor identidad.

**Resultado obtenido:**

| Proteína | Cromosoma | Locus principal | Otros loci (parálogos) |
|---|---|---|---|
| **Hemoglobina β** | **11** (`NC_000011.10`), hebra − | ***HBB*** 5.225.601-5.227.021: **3 HSP = 3 exones** (aa 1-45 / 31-106 / 105-147), cobertura 100 %, E = 5e-44 | Todo el **clúster β-globina** del mismo cromosoma: *HBD* (δ, 5,23 Mb, 3 exones), pseudogén *HBBP1* (5,24 Mb, solo se detecta 1 exón), *HBG1* y *HBG2* (γ fetales, 5,25 Mb) y *HBE1* (ε embrionaria, 5,27 Mb). En el **cromosoma 16**: el **clúster α-globina** (153-181 kb: *HBA2*, *HBA1*, *HBQ1*…), con un 32-37 % de identidad. |
| **Mioglobina** | **22** (`NC_000022.11`), hebra − | ***MB*** 35.607.300-35.617.257: **3 HSP = 3 exones codificantes** (aa 1-32 / 32-110 / 93-154), cobertura 100 %, E = 7e-41 | **Ninguno** con E < 1e-5: la mioglobina se separó de las hemoglobinas hace tanto tiempo que su parecido actual (~25 %) no supera el umbral. |
| **Oxitocina** | **20** (`NC_000020.11`), hebra + | ***OXT*** 3.071.656-3.072.415: **2 HSP**, cobertura 100 %, E = 5e-43. El gen tiene 3 exones, pero sus intrones son muy cortos y, probablemente, `tblastn` une los exones 2 y 3 en un solo HSP atravesando el intrón (por eso su identidad baja al 75 %). | ***AVP*** (vasopresina) a solo **~11 kb** y en la **hebra opuesta** (3.082.981-3.083.181). Alinea el tramo de la **neurofisina** (aa 40-106) con un 90 % de identidad. |

**Conclusiones:**
1. `tblastn` es el programa adecuado cuando la consulta es una proteína y la BD es ADN: **localiza el gen** de una proteína en un genoma sin anotar.
2. El número de HSP colineales del locus principal **dibuja la estructura exón-intrón** del gen: 3 exones para *HBB* y *MB*.
3. Los loci secundarios muestran la **evolución por duplicación génica**. El clúster β está ordenado en el cromosoma igual que se expresa durante el desarrollo: ε (embrión) → γ (feto) → δ y β (adulto). Oxitocina y vasopresina son una duplicación en tándem de una hormona ancestral.
4. Del pseudogén *HBBP1* solo se detecta un exón: es una copia que perdió su función y acumula mutaciones, y el resto ya no se parece lo suficiente.
5. Con el genoma en local se pueden lanzar todas las consultas necesarias sin depender de NCBI, a cambio de descargar varios GB.

## Ejercicio 4 · El árbol de la vida (blastp)

**Nota sobre la búsqueda online.** La petición `blastp` contra `nr` en NCBI superó el límite de 15 minutos (*"BLAST request GERZV65P014 is taking longer than 10 minutes"*). `nr` es la mayor base de datos de proteínas del NCBI y la cola estaba saturada. El ejercicio se resolvió en **local contra Swiss-Prot** (celda siguiente).

**Los 10 mejores hits (Swiss-Prot)** son **todos de primates**, con identidad del 96-100 %, cobertura del 99-100 %, E-values entre 1e-96 y 1e-100 y 276-286 bits:

| # | Especie | Identidad |
|---|---|---|
| 1 | *Homo sapiens*, *Pan troglodytes*, *Pan paniscus* (P69905) | 100 % |
| 2 | Gorila | 99,3 % |
| 3 | Gibón (*Hylobates lar*) | 98,6 % |
| 4 | Orangután | 97,9 % |
| 5-10 | Langur, macacos, mono araña, tití, cercopiteco | 96-98 % |

- La HBA de **humano, chimpancé y bonobo es idéntica**: Swiss-Prot las reúne en una sola entrada (P69905).
- **El top 10 no muestra el árbol de la vida**, solo que la HBA está **casi intacta en primates**. En `nr` pasaría lo mismo, o peor: los 10 primeros serían casi todos copias humanas, estructuras PDB y grandes simios.
- Con E-values prácticamente nulos en todos, para ordenarlos importan la **identidad** y el **bit-score**.

**El gradiente evolutivo (un hit por especie, 239 especies con `-max_target_seqs 300`):**

| Grupo | Ejemplos (identidad con la HBA humana) |
|---|---|
| Primates | 92-100 % |
| Otros mamíferos placentarios | caballo 88 %, vaca 88 %, ratón 85,9 %, perro 83,7 % |
| Marsupiales | canguro 80,9 %, zarigüeya 71,6 % |
| Monotremas | equidna 73,8 %, ornitorrinco 72,3 % |
| Aves | pato 71,8 %, gallina 70,4 %, paloma 69 %, avestruz 68,8 % |
| Reptiles | cocodrilo 68,8 %, caimán 66,9 %, tortugas y serpientes 57-65 % |
| Anfibios | *Xenopus laevis* 63,4 %, *X. tropicalis* 59,2 %, tritones 56-60 % |
| Peces óseos | salmón 57,3 %, carpín 53,2 %, pez cebra 53,1 %, peces antárticos ~49 % |

**Comentarios:**
1. **La identidad baja aproximadamente en el orden en que se separaron los linajes del nuestro**: primates → mamíferos placentarios → marsupiales y monotremas → aves y reptiles → anfibios → peces. Es el **"reloj molecular"**: la hemoglobina conserva su función (los residuos que sujetan el grupo hemo, los contactos α-β) mientras acumula mutaciones neutras con el tiempo.
2. El reloj **no es perfecto**: algunas ardillas (~90 %) se parecen más a nosotros que algunos primates lejanos (lémures, ~85-89 %), porque cada linaje evoluciona a su propia velocidad.
3. Entre las otras especies aparece una proteína **humana**: la **hemoglobina θ-1** (*HBQ1*, 62 %). Es un **parálogo** de nuestro propio genoma, el mismo que salió en el cromosoma 16 en el ejercicio 3. Una búsqueda por similitud no distingue ortólogos (el mismo gen en otra especie) de parálogos (copias dentro de la misma especie).
4. Con 300 hits la lista llega hasta los peces óseos (~48 %). Para ver tiburones, lampreas u otras globinas (β, mioglobina) habría que subir `-max_target_seqs` (p. ej. 1000). El siguiente paso natural sería un alineamiento múltiple de estas secuencias y un árbol filogenético con `Bio.Phylo`.
5. **`nr` frente a Swiss-Prot:** `nr` es mucho mayor y redundante, y online tarda mucho. Swiss-Prot está revisada a mano, tiene una entrada por proteína y especie y en local responde en segundos. Para este tipo de análisis evolutivo es la mejor opción.
