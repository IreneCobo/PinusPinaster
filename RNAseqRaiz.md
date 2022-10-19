## RNA-seq de muestras de raíz en *Pinus pinaster*

### 1. Cómo acceder a las raw reads de RNA-seq y descripción de sus nombres

Las raw reads de raíz y tallo están aquí: /home/irene/PinasterRaizRNAseqReads

Yo solo voy a utilizar las de raíz.

```
Lib-Truseq-RNA-173367_18-7-RCO_P_151471_S16_R1_001.fastq
```
**Descripción del nombre**
1. Hasta el guión bajo _ no proporciona información de la muestra

2. Tras el guión bajo nos podemos encontrar o un 1 o un 18 (como en este ejemplo) que habla del portainjertos. El portainjertos contiene la raíz y un trozo grande de tallo

- 1: Portainjertos sensible
- 18: Portainjertos tolerante

3. Después hay un guión -. Lo que viene después del guión es el número del individuo (en este caso es el individuo 7)

4. Luego viene un código de 3 letras (en este caso RCO)

- R: Raíz (o A: acícula)
- C: Habla del tratamiento que puede ser Control (o S sequía)
- O: Habla de la púa, que puede ser Oria, que es tolerante a la sequía (o G Galicia, sensible a la sequía). La púa se colocó muy arriba, luego la púa tiene acículas y un poco de tallo

### 2. Análisis de las raw reads de RNA-seq

Algunos programas son más adecuados que otros para analizar datos de RNA-seq en función de si se tiene un transcriptoma o genoma de referencia, de si el genoma/transcriptoma está completamente anotado etc... Sin embargo, dentro de cada categoría, en muchos casos hay varios programas que hacen esencialmente lo mismo (por ejemplo, Salmon y Kallisto son pseudoaligners para mapear reads contra transcriptomas de referencia parcialmente anotados). Los siguientes artículos son de utilidad para tener una visión más clara sobre cómo tomar decisiones a la hora de analizar datos de RNA-seq (tener en cuenta que esto ha sido escrito en 2022-2023, así que con el tiempo estos artículos pueden quedar desactualizados):

- [**Systematic comparison and assessment of RNA-seq procedures for gene expression quantitative analysis**, Corchete et al. 2020](https://doi.org/10.1038/s41598-020-76881-x)
- [**Approaches to variant discovery for conifer transcriptome sequencing**, Telfer et al. 2018](https://doi.org/10.1371/journal.pone.0205835)
- [**Applications of transcriptome in conifer species**, Wei et al. 2022](https://doi.org/10.1007/s11240-022-02322-4)
- [**A survey of best practices for RNA-seq data analysis**, Conesa et al. 2016](https://doi.org/10.1186/s13059-016-0881-8)

**2.1. Control de calidad con FASTQC/MULTIQC**

Primero, en esa misma carpeta, utilicé el siguiente código para realizar el análisis de FASTQC (utilizo la versión 0.11.7)

```
module load fastqc/0.11.7
fastqc Lib*
```

También se puede utilizar el "path" donde se encuentren las secuencias y correr fastqc en otra carpeta. FASTQ reconoce automáticamente los archivos que contienen las raw reads de la carpeta a la que le dirijas (formato fastq.gz).

```
fastqc /home/FCAM/icobosimon/PineTestRaw/Lib*
```

Como son 24 muestras de raíz, obtienes dos resultados por muestra: un archivo acabado en *.html* y otro en *.zip*

Ahora utilizo [MULTIQC](https://multiqc.info/docs/) v1.12 para obtener un resumen de los resultados de FASTQC de las 24 muestras. MULTIQC te ofrece un report en formato *html* que contiene todos los resultados de FASTQC de todas las muestras en conjunto. Para ello, en la misma carpeta donde están los resultados de FASTC, corres el siguiente comando:

```
multiqc .
```

MULTIQC ya reconoce los archivos que necesita (los acabados en *.zip*). También puedes correr lo siguiente si no te fías

```
multiqc *_fastqc.zip
```

Obtienes como resultado el siguiente report: [multiqc_report](multiqc_report.html). Aquí hay un vídeo de youtube donde te explican cómo interpretar los resultados de MULTIQC: https://www.youtube.com/watch?v=qPbIlO_KWN0

**2.2. Trimming por calidad y longitud**

Aunque en tallo y acícula se utilizó el programa [Reformat de BBMap](https://github.com/BioInfoTools/BBMap/blob/master/sh/reformat.sh) para realizar el trimming por calidad y longitud, yo utilicé Trimmomatic v.0.39 ya que estoy más familiarizada con él y es el más usado en RNA-seq. Los resultados del FASTQC/MULTIQC mostraron ausencia de adaptadores y/o overrepresented sequences, así que no hizo falta remover adaptadores con trimmomatic. Para el trimming por calidad y longitud, utilicé las siguientes flags y probé tres thresholds diferentes: (1) Los usados en acícula y tallo con reformat, (2) los que vienen por defecto en el [manual de trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) y (3) los utilizados en el [manual de la UConn para el análisis de RNA-seq de especies no modelo](https://github.com/CBC-UCONN/RNAseq_nonmodel/blob/master/02_Quality_Control/trimmomatic.sh). Las flags y thresholds por defecto de trimmomatic serían las siguientes: 

- Remove adapters (Flag: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10): No lo usé porque no había adapters
- Remove leading low quality or N bases (below quality 3) (Flag: LEADING:3. El usado con reformat era 20, el manual de la UConn no usa este flag
- Remove trailing low quality or N bases (below quality 3) (Flag: TRAILING:3). EL usado por reformat era 20, el manual de la UConn no usa este flag
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (Flag: SLIDINGWINDOW:4:15). El usado por reformat era 20; en el manual de la UConn, 25. 
- Drop reads below the 36 bases long (Flag: MINLEN:36). El usado por reformat era 30, en el manual de la UConn, 45

El código quedó de la siguiente manera (aquí con los thresholds del programa reformat adaptado para trimmomatic). Como son Single End reads, uso el código de trimmomatic para single ends (SE)

```
module load Trimmomatic/0.39
java -jar $Trimmomatic SE -threads 4 Lib-Truseq-RNA-171792_1-4-RCG_P_171958_S17_R1_001.fastq.gz TrimmedReformat_RaizIndv4PortaSPuaSControl.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30
```
Los resultados de multiqc ([Manual](multiqc_reportDefault.html), [UConn](multiqc_reportUConn.html) y [Reformat](multiqc_reportReformat.html) mostraron que los **parámetros por defecto del manual de Trimmomatic ofrecieron la menor calidad de secuencias (porcentaje de secuencias duplicadas) pero un mayor número de reads**. Los parámetros usados en acículas y tallos con **Reformat, mostraron valores intermedios tanto de calidad (porcentaje de secuencias duplicadas) como de número de reads**. Finalmente, los parámetros del manual de RNA-seq de la **UConn, mostraron los mejores valores de calidad (porcentaje de secuencias duplicadas) pero el menor número de reads, aunque las diferencias no fueron muy significativas entre los tres parámetros**, con la excepción de una de las secuencias que mostró un cambio significativo a mejor en cuanto al número de sequence duplications con los parámetros de la UConn. El resto de parámetros fueron prácticamente iguales entre los tres criterios de trimming. Sin embargo, **decidí usar los parámetros de Reformat**, ya que mostraron valores intermedios y son los usados para acícula y tallo.

**Uno de los individuos mostró valores preocupantes de secuencias duplicadas tras el trimming usando los tres criterios (Indv9PortaSPuaTControl)**. Veremos si es necesario eliminarla en posteriores análisis tras la eliminación del RNA ribosomal y el análisis exploratorio de los datos con PCA. De momento, la dejo. 

**2.3. SortMeRNA: eliminación del RNA ribosomal**

Aunque no hay consenso en la necesidad de eliminar o no el RNA ribosomal o rRNA (yo nunca lo he visto en ninguno de los manuales que he usado para realizar RNA-seq, tanto el de [UCDavis](https://jnmaloof.github.io/BIS180L_web/labs/) como en los de la UConn [with](https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation) and [without reference genome and annotation](https://github.com/CBC-UCONN/RNAseq_nonmodel), muchos argumentan que no es necesario, ya que no se va a cuantificar después (https://www.biostars.org/p/419845/), mientras otros argumentan que remover el rRNA es uno de los pasos más infraestimados en RNA-seq análisis (https://www.qiagen.com/us/knowledge-and-support/knowledge-hub/science-matters/genomics/ribosomal-rna-removal), aunque en la fase de library preparation, no data analisis. 

Dado que se realizó esta eliminación de rRNA en hoja y tallo, voy a realizarlo también en raíz, de este modo aseguro que los resultados sean comparables. 

Cambio de opinión tras leer el [manual de SortMeRNA](https://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.1.pdf). Dice que "The main application of SortMeRNA is filtering rRNA from metatranscriptomic data.". Nosotros no trabajamos con datos meta-transcriptómicos, por eso en ninguno de los manuales de RNA-seq para transcriptómica he visto la eliminación de rRNA de las reads. Sí que he visto, en cambio, eliminar contaminantes en el transcriptoma de referencia antes del mapeo de las reads contra él, y es lo que voy a hacer. Usaré para ello el program enTAP. Mantendré las dos versiones del transcriptoma de referencia (con y sin contaminantes) y realizaré el DE analysis usando ambos, por si acaso. 

Antes de ponerme con el contaminant filtering using EnTAP voy a preguntar cómo se realizó el ensamblaje de novo del transcriptoma (por si acaso este contaminant filtering ya se ha hecho y lo hago otra vez en tonto). Mientras tanto, voy a ir usando Salmon para el mapeo de reads contra el transcriptoma de referencia para ir familiarizándome con el programa y luego haré el DE análisis usando DESeq2. 

**2.4. Salmon: Mapping y conteo**

Aunque he usado siempre kallisto hasta ahora, salmon y kallisto son pseudo-aligners que hacen esencialmente lo mismo (apropiados para transcriptomas de referencia de especies no modelo, que no están completamente anotados). [Aquí](https://gencore.bio.nyu.edu/salmon-kallisto-rapid-transcript-quantification-for-rna-seq-data/) hay una comparación del año 2016 de la velocidad y accuracy de ambos pseudo-aligners. 

Utilizo el [manual de salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

Lo primero que hay que hay que hacer es indexar el transcriptoma de referencia. Para ello uso el siguiente código

```
module load salmon/1.9.0
salmon index -t pinaster.all.cdhit.fasta -i pinaster_index -k 31
```
Flags:
- -t: el nombre del transcriptoma de referencia
- -i el nombre de la carpeta que contendrá los outputs de la indexación
- -k: el tamaño mínimo de los k-mers considerados aceptables (the k size selected here will act as the minimum acceptable length for a valid match). We find that a k of 31 seems to work well for reads of 75bp or longer, but you might consider a smaller k if you plan to deal with shorter reads. (Según el mutiqc report el tamaño de nuestras muestras cae en 75 bp)

Una vez indexado, para obtener las "counts" de las reads mapeadas contra el transcriptoma de referencia (por cada individuo) utilizo el siguiente código. Estas "counts" se usarán para el DE analisis.

```
module load salmon/1.9.0
salmon quant --threads 8 -i pinaster_index -l IU -r TrimmedReformat_RaizIndv10PortaSPuaSSequia.fastq.gz --validateMappings -o TrimmedReformat_RaizIndv10PortaSPuaSSequia_quant
```

Flags:
- --threads: Es el número de cpus que serán usadas por el supercomputador para llevar a cabo el análisis. El propio manual dice que han visto que se obtienen resultados óptimos usando entre 8 y 12, por eso elegí 8
- -i: El nombre de la carpeta con los resultados de la indexación, que dimos en el paso anterior
- -l: library type. Según el manual, para librerías comprimidas en fastq.gz, hay que poner IU
- -r: El nombre de la librería que vamos a mapear (se usa la flag -r para single ends, como es mi caso)
- -o: El nombre de la carpeta que contendrá los outputs del mapeo (con las counts). Le he puesto el mismo nombre que el individuo

```
salmon quant --threads 8 -i pinaster_index -l A -r TrimmedReformat_RaizIndv10PortaSPuaSSequia.fastq.gz --validateMappings --writeUnmappedNames -o TrimmedReformat_RaizIndv10PortaSPuaSSequia_quant_UnmappedNames
```

También usé este código para incluir la flag *--writeUnmappedNames*, que "will tell Salmon to write out the names of reads that do not map to the transcriptome."

Los *.log* files obtenidos en los outputs obtenidos tras correr este código mostraron los siguientes valores de mapeo

**Tabla con los resultados del mapeo contra el transcriptoma de referencia de cada combinación individuo/púa-porta/tratamiento**

| Individuo | Porta | Púa | Tratamiento | Porcentaje de lecturas mapeadas | Número de lecturas mapeadas |
| --------- | ----- | --- | ----------- | ------------------------------- | --------------------------- |
| 1 | Sensible | Oria (Tolerante) | Control | 88.9612% |  18,676,081 |
| 2 | Sensible | Oria (Tolerante) | Sequía | 87.6052% | 9,931,452 |
| 3 | Tolerante | Oria (Tolerante) | Sequía | 85.5583% | 16,921,952 |
| 4 | Sensible | Galicia (Sensible) | Control | 85.604% | 28,207,117 |
| 4 | Tolerante | Galicia (Sensible) | Control | 86.0707% | 15,725,203 |
| 5 | Sensible | Galicia (Sensible) | Control | 85.0118% | 11,937,000 |
| 5 | Tolerante | Galicia (Sensible) | Control | 89.3177% | 14,269,540 |
| 6 | Sensible | Galicia (Sensible) | Control | 86.5574% | 10,716,544 |
| 6 | Tolerante | Galicia (Sensible) | Control | 87.7721% | 27,582,364 |
| 7 | Tolerante | Oria (Tolerante) | Control | 88.8661% | 10,677,096 |
| 8 | Sensible | Galicia (Sensible) | Sequía | 89.7663% | 10,705,819 |
| 9 | Sensible | Oria (Tolerante) | Control | 76.8114% | 9,141,508 |
| 9 | Tolerante | Oria (Tolerante) | Control | 84.3641% | 11,062,264 |
| 10 | Sensible | Galicia (Sensible) | Sequía | 86.0956% | 8,529,399 |
| 10 | Tolerante | Galicia (Sensible) | Sequía | 89.3171% | 12,002,344 |
| 11 | Sensible | Oria (Tolerante) | Control | 85.8028% | 10,052,350 |
| 11 | Tolerante | Oria (Tolerante) | Sequía | 87.9333% | 11,886,579 |










