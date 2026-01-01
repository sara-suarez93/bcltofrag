# ejecutr en bash
#!/bin/bash

# crear carpeta para reports de fastqc
outdir="../data/fastqc/$pooldir"
mkdir -pv "$outdir"

# for loop para correr fastqc en cada carpeta pool con archivos fastq.gz
for d in ../data/fastq/POOL*; do

# control del proceso
echo "Processing pool: $d"

# guardo el nombre del pool (POOL-numeroN)
pooldir=$(basename "$d")

# guardo el directorio de la carpeta que contiene los archivos .fastq.gz
# (distinto nombre para cada pool. Es el otro directorio que no es ni reports ni stats)
pool_fastq_dir=$(find "$d" -mindepth 1 -maxdepth 1 -type d \
! -name "Stats" \
! -name "Reports"); 

# corro el programa para cada archivo dentro del directorio y guardo en el outdir
fastqc "$pool_fastq_dir"/*.fastq.gz -o "$outdir"

# fin del loop
done

# multiqc genera report multiqc leyendo la carpeta con output de fastqc
multiqc fastqc -o ./data/multiqc