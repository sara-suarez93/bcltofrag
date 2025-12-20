import os
import re
import glob
from finaletoolkit.frag import end_motifs

results_endmot = "results/02_endmotifs"
os.makedirs(results_endmot, exist_ok=True)

# ajustar bam directories
# lista de todos los archivos bam - archivos de fragmentos
bam_files = (
    glob.glob("../data/alignments/POOL-*/*.markdup.bam")
    + glob.glob("../data/alignments_qc/*.filtered.bam")
)

klist = [2, 4] # bases de los motivos - igual que colon

# para cada muestra (sort por nombre)
for filename in sorted(bam_files):
    samplename = os.path.basename(filename)

    # ignorar muestras vacias (ID1–ID24)
    if re.match(r"^ID([1-9]|1[0-9]|2[0-4])_", samplename):
        continue

    print(f"Processing {samplename}")

   # quitar .bam para guardar 
    samplename = samplename.replace(".bam", "")

    for ki in klist:

        # nombre archivo para cada motivo
        out_file = os.path.join(results_endmot, f"{samplename}.k{ki}.tsv")
    
        end_motifs(
            input_file = filename,
            # path a 2bit file genoma referencia
            refseq_file = "hg38.2bit",
            k=ki,
            both_strands=True, # defecto, 5´ y 3´
            output_file = out_file)