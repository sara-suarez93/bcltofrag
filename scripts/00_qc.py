import os
import re
import glob
import pysam # leer bam en python
import matplotlib.pyplot as plt
import numpy as np
import time

# este script:
# 1) aplica QC descrito en doi: 10.1093/bioadv/vbaf236 
# Only high-quality reads were considered in the analysis (high quality:
# uniquely mapped, no PCR duplicates, both ends are mapped
# with mapping qualities more than 30, and properly paired).

# 2) y genera barplot de lecturas incluidas/excluidas

# output directories
# datos procesddos (QC) - en root (donde data folder)
data_qcdir = '../data/alignments_qc'
os.makedirs(data_qcdir, exist_ok=True)

# calcular tiempo si se ejecuta en notebook/terminal
start = time.time()

# path archivos bam
data_bamsdir = glob.glob("../data/alignments/POOL-*/*.markdup.bam")

# bams a analizar (no vacios) - esto es temporar ya que deberian de haberse ignorado en el alineamiento
bams_lecturas = []

for filepath in data_bamsdir:
    # guardo nombre sample (base del path) sin markdup.bam
    samplename = os.path.basename(filepath).replace('.markdup.bam', '')
    
    # si no es muestra a excluir (vacia)
    if not re.match(r"^ID([1-9]|1[0-9]|2[0-4])_", samplename):
        # anyado a bam con lecturas
        bams_lecturas.append(filepath)

total_files = len(bams_lecturas)

for i, filepath in enumerate(bams_lecturas, start=1):

    # guardo nombre sample (base del path) sin markdup.bam
    samplename = os.path.basename(filepath).replace('.markdup.bam', '')

    # anyado a lista
    #sample_list.append(sample)

    # crear output qc file a escribir
    out_bam_path = os.path.join(data_qcdir, samplename + ".filtered.bam")

    # y escribir cada lectura (los pares que han pasado el filtro)
    print(f'Leyendo y filtrando archivo .markdup.bam #{i}/{total_files} y escribiendo en: {out_bam_path}\n')
        
    # leer bam u escribir filtered.bam correspondiente
    with pysam.AlignmentFile(filepath, "rb") as bam, pysam.AlignmentFile(out_bam_path, "wb", header=bam.header) as out_bam:

        # para cada lectura hasta el final (until_eof=T)
        for read in bam.fetch(until_eof=True):

            # excluir si se cumple cualquiera de estas reglas
            if (
                read.is_unmapped or # lectura no mapeada
                read.mate_is_unmapped or # par de lecturas no alineada ok
                read.is_duplicate or # lectura esta duplicada
                not read.is_proper_pair or # lectura no esta emparejada
                read.mapping_quality <= 30 # la calidad del alineamiento es baja
            ):
            
            out_bam.write(read)
    
    # crear bai (necesario para FinaleToolkit
    pysam.index(out_bam_path)    


end = time.time()

# imprimo tiempo de computacion
print(f"Tiempo total: {end - start/60:.2f} minutos")