## RNA-seq de muestras de raíz en Pinus pinaster

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

**2.1. Control de calidad con FASTQC/MULTIQC**

Primero, en esa misma carpeta, utilicé el siguiente código para realizar el análisis de FASTQC

```
fastqc Lib*
```

También se puede utilizar el "path" donde se encuentren las secuencias y correr fastqc en otra carpeta. FASTQ reconoce automáticamente los archivos que contienen las raw reads de la carpeta a la que le dirijas (formato fastq.gz).

```
fastqc /home/FCAM/icobosimon/PineTestRaw/Lib*
```

Como son 24 muestras de raíz, obtienes dos resultados por muestra: un archivo acabado en *.html* y otro en *.zip*

Ahora utilizo [MULTIQC](https://multiqc.info/docs/) para obtener un resumen de los resultados de FASTQC de las 24 muestras. MULTIQC te ofrece un report en formato *html* que contiene todos los resultados de FASTQC de todas las muestras en conjunto. Para ello, en la misma carpeta donde están los resultados de FASTC, corres el siguiente comando:

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

El código quedó de la siguiente manera (aquí con los thresholds del programa reformat adaptado para trimmomatic): 

```
module load Trimmomatic/0.39
java -jar $Trimmomatic SE -threads 4 Lib-Truseq-RNA-171792_1-4-RCG_P_171958_S17_R1_001.fastq.gz TrimmedReformat_RaizIndv4PortaSPuaSControl.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30
```







