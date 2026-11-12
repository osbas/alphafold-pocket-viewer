# Referencias Académicas y Bibliografía de Consulta — AlphaFold y Biología Estructural

Este compendio reúne la literatura fundamental (libros de texto universitarios, artículos científicos con revisión por pares, cursos y manuales interactivos) para comprender a fondo los principios físicos, biológicos y computacionales que sustentan a **AlphaFold** y al microservicio **AlphaFold-Lite Pocket Viewer**.

---

## 1. Libros de Texto Universitarios (Fundamento Bioquímico y Estructural)

### 1.1 Bioquímica y Plegamiento Termodinámico
* **Nelson, D. L., & Cox, M. M. (2017/2021). *Lehninger Principios de Bioquímica* (7ª / 8ª edición). Ediciones Omega / W. H. Freeman.**
  * **ISBN-13:** 978-8428216678 (Ed. en español, Omega) / 978-1319228002 (Macmillan).
  * **Capítulos de estudio prioritario:**
    * **Capítulo 3: *Aminoácidos, péptidos y proteínas*.** Química, polaridad y propiedades ácido-base de los 20 aminoácidos estándar; formación y rigidez planar del enlace peptídico ($C-N$ con carácter parcial de doble enlace).
    * **Capítulo 4: *La estructura tridimensional de las proteínas*.** Ángulos de torsión del esqueleto $\phi$ (phi) y $\psi$ (psi); gráfico de conformaciones estéricamente permitidas de G. N. Ramachandran; estabilidad de $\alpha$-hélices y hojas $\beta$; desnaturalización, renaturalización y el postulado termodinámico de Christian Anfinsen ($\Delta G_{\text{plegamiento}} < 0$).

* **Alberts, B., Johnson, A., Lewis, J., Morgan, D., Raff, M., Roberts, K., & Walter, P. (2015/2022). *Biología Molecular de la Célula* (6ª / 7ª edición). Editorial Omega / Garland Science.**
  * **ISBN-13:** 978-8428216609 / 978-0393884821.
  * **Capítulo de estudio prioritario:**
    * **Capítulo 3: *Proteínas*.** Cómo la cadena polipeptídica adopta una conformación única dictada por interacciones no covalentes (puentes de hidrógeno, atracciones electrostáticas y fuerzas hidrofóbicas que ocultan residuos no polares en el núcleo interior); sitios de unión de ligandos y enzimas como máquinas catalíticas.

### 1.2 Biología Estructural y Bioinformática
* **Brändén, C., & Tooze, J. (1999). *Introduction to Protein Structure* (2ª edición). Garland Science.**
  * **ISBN-13:** 978-0815323051.
  * **Relevancia:** Es el tratado clásico sobre arquitectura molecular. Detalla los motivos supersecundarios ($\beta-\alpha-\beta$, horquillas $\beta$), dominios estructurales (todo-$\alpha$, todo-$\beta$, barriles TIM $\alpha/\beta$) y enseña a interpretar visualmente las representaciones 3D como las generadas en nuestro visor HTML5 Canvas.

* **Gu, J., & Bourne, P. E. (Eds.). (2009). *Structural Bioinformatics* (2ª edición). Wiley-Blackwell.**
  * **ISBN-13:** 978-0470181058.
  * **Capítulos de estudio prioritario:**
    * **Capítulo 4: *Fundamentals of Protein Structure*.**
    * **Capítulo 9: *Coordinate Systems and PDB Format*.** Explica en profundidad la sintaxis del archivo `.pdb` (especificación de líneas `ATOM`, campos de coordenadas cartesianas $x, y, z$ en Ångströms y el factor de temperatura de Debye-Waller o *B-factor*, donde AlphaFold guarda el valor pLDDT).
    * **Capítulo 12: *Geometric and Energetic Properties of Proteins*.** Deducción formal del **Radio de Giro ($R_g$)** como medida de compacidad espacial respecto al centro de masa geométrico:
      $$R_g = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \|\mathbf{r}_i - \mathbf{r}_{\text{cm}}\|^2}$$

---

## 2. Libros de Inteligencia Artificial y Aprendizaje Geométrico 3D

* **Bronstein, M. M., Bruna, J., Cohen, T., & Veličković, P. (2021). *Geometric Deep Learning: Grids, Groups, Graphs, Geodesics, and Gauges*. MIT Press (Preprint abierto).**
  * **Acceso abierto:** [https://geometricdeeplearning.com/](https://geometricdeeplearning.com/)
  * **Relevancia:** Las moléculas tridimensionales son invariantes ante rotaciones y traslaciones (si rotas una proteína en el espacio, sus funciones biológicas son idénticas aunque sus coordenadas $(x, y, z)$ cambien). Este texto explica la teoría de grupos de simetría $SE(3)$ (Special Euclidean Group) que AlphaFold implementa en su módulo de estructura (*Invariant Point Attention / Structure Module*).

* **Goodfellow, I., Bengio, Y., & Courville, A. (2016). *Deep Learning*. MIT Press.**
  * **Acceso abierto:** [https://www.deeplearningbook.org/](https://www.deeplearningbook.org/)
  * **ISBN-13:** 978-0262035613.
  * **Relevancia:** Fundamentos de optimización matemática, retropropagación, representaciones vectoriales (*embeddings*) y mecanismos de autoatención (*Transformers*).

---

## 3. Artículos Científicos Fundacionales (Publicaciones Clave)

1. **El Modelo AlphaFold 2 (Google DeepMind):**
   * Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., Back, T., Petersen, S., Reiman, D., Clancy, E., Zielinski, M., Steinegger, M., Pacholska, M., Berghammer, T., Bodenstein, S., Silver, D., Vinyals, O., Senior, A. W., Kavukcuoglu, K., Kohli, P., & Hassabis, D. (2021).
   * **Título:** *Highly accurate protein structure prediction with AlphaFold*.
   * **Revista:** *Nature*, 596(7873), 583–589.
   * **DOI:** [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
   * **Aspecto clave para el proyecto:** Describe la arquitectura de dos vías *Evoformer*, la función de pérdida *FAPE (Frame Aligned Point Error)* y el cálculo del **pLDDT** discretizado en 50 intervalos de confianza.

2. **La Base de Datos de Estructuras Proteómicas de AlphaFold:**
   * Varadi, M., Anyango, S., Deshpande, M., Nair, S., Natassia, C., Yordanova, G., Yuan, D., Stroe, O., Wood, G., Laydon, A., Žídek, A., Green, T., Tunyasuvunakool, K., Petersen, S., Jumper, J., Clancy, E., Green, R., Vora, A., Lutfi, M., Figurnov, M., Cowie, A., Hobbs, N., Kohli, P., Kleywegt, G., Birney, E., Hassabis, D., & Velankar, S. (2022).
   * **Título:** *AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models*.
   * **Revista:** *Nucleic Acids Research*, 50(D1), D439–D444.
   * **DOI:** [10.1093/nar/gkab1061](https://doi.org/10.1093/nar/gkab1061)
   * **Aspecto clave para el proyecto:** Describe la API pública de consulta en EMBL-EBI consumida por nuestro microservicio (`https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`).

3. **La Métrica Original lDDT (Base del pLDDT):**
   * Mariani, V., Biasini, M., Leaver-Fay, A., & Schwede, T. (2013).
   * **Título:** *lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests*.
   * **Revista:** *Bioinformatics*, 29(21), 2722–2728.
   * **DOI:** [10.1093/bioinformatics/btt473](https://doi.org/10.1093/bioinformatics/btt473)
   * **Aspecto clave para el proyecto:** Demuestra por qué el cálculo local de diferencias de distancia intramolecular en un radio de corte de 15 Å es insensible a deformaciones de bisagra y superior al clásico RMSD global.

4. **El Postulado de Plegamiento Autónomo (Premio Nobel 1972):**
   * Anfinsen, C. B. (1973).
   * **Título:** *Principles that govern the folding of protein chains*.
   * **Revista:** *Science*, 181(4096), 223–230.
   * **DOI:** [10.1126/science.181.4096.223](https://doi.org/10.1126/science.181.4096.223)

5. **La Paradoja de Levinthal (1969):**
   * Levinthal, C. (1969).
   * **Título:** *How to fold graciously*.
   * **Revista:** *Mossbauer Spectroscopy in Biological Systems*, Proceedings of a meeting held at Allerton House, Monticello, Illinois, 22–24.

---

## 4. Cursos Universitarios Abiertos, Manuales y Portales Educativos

1. **EMBL-EBI Training — "AlphaFold: a practical guide to the database and predictions":**
   * **Institución:** European Bioinformatics Institute (Wellcome Genome Campus, Hinxton, Reino Unido) en colaboración con Google DeepMind.
   * **Enlace interactivo gratuito:** [https://www.ebi.ac.uk/training/online/courses/alphafold/](https://www.ebi.ac.uk/training/online/courses/alphafold/)
   * **Módulos recomendados:**
     * *Understanding AlphaFold predictions:* Interpretación visual de las 4 bandas pLDDT.
     * *Confidence scores (pLDDT and PAE):* Cómo distinguir regiones ordenadas de lazos intrínsecamente desordenados (IDRs).

2. **RCSB Protein Data Bank — Educational Portal (PDB-101):**
   * **Enlace:** [https://pdb101.rcsb.org/](https://pdb101.rcsb.org/)
   * **Artículos temáticos recomendados:**
     * *Hemoglobin (Subunidad Alfa, P69905):* Estructura alostérica y transporte de oxígeno.
     * *Insulin (P01308):* Procesamiento proteolítico de proinsulina a insulina madura.
     * *Tumor Suppressor p53 (P04637):* Organización modular en dominios estructurados y regiones intrínsecamente desordenadas.
     * *Guía "Understanding PDB Data":* Desglose pedagógico de los registros `ATOM`, `HETATM` y campos de factor B.

3. **UniProt Knowledgebase (UniProtKB):**
   * **Enlace:** [https://www.uniprot.org/](https://www.uniprot.org/)
   * **Fichas biológicas de las 4 proteínas de la aplicación:**
     * `P69905`: [Hemoglobin subunit alpha (Homo sapiens)](https://www.uniprot.org/uniprotkb/P69905/entry)
     * `P01308`: [Insulin (Homo sapiens)](https://www.uniprot.org/uniprotkb/P01308/entry)
     * `P01542`: [Crambin (Crambe hispanica subsp. abyssinica)](https://www.uniprot.org/uniprotkb/P01542/entry)
     * `P04637`: [Cellular tumor antigen p53 (Homo sapiens)](https://www.uniprot.org/uniprotkb/P04637/entry)

4. **Premio Nobel de Química 2024 — Comunicado y Documento Científico Avanzado:**
   * **Enlace al anuncio:** [https://www.nobelprize.org/prizes/chemistry/2024/press-release/](https://www.nobelprize.org/prizes/chemistry/2024/press-release/)
   * **Scientific Background:** Documento de 18 páginas redactado por el Comité Nobel de Química explicando la evolución desde la predicción por homología clásica hasta AlphaFold 2 y AlphaFold 3.
