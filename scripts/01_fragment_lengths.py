import os
import re
import glob
from finaletoolkit.frag import frag_length_bins

results_frag = "results/01_fragment_lengths"
os.makedirs(results_frag, exist_ok=True)

# ajustar bam directories
# lista de todos los archivos bam - archivos de fragmentos
bam_files = (
    glob.glob("../data/alignments/POOL-*/*.markdup.bam")
    + glob.glob("../data/alignments_qc/*.filtered.bam")
)

# para cada muestra (sort por nombre)
for filename in sorted(bam_files):
    samplename = os.path.basename(filename)

    # ignorar muestras vacias (ID1–ID24)
    if re.match(r"^ID([1-9]|1[0-9]|2[0-4])_", samplename):
        continue

    print(f"Processing {samplename}")

   # quitar .bam para guardar 
    samplename = samplename.replace(".bam", "")

    frag_length_bins(
        input_file=filename,
        bin_size=1,
        min_length=50,
        max_length=300,
        quality_threshold=30,
        output_file=os.path.join(results_frag, samplename + ".tsv"),
        # guardar pngs
         histogram_path=os.path.join(results_frag, samplename + ".png")
    )
