import os
import glob


configfile: "config.yaml"

# limitar wildcards al nombre de la muestra
# en las siguientes reglas, el pool y la muestra alojadas en la clave
# del diccionario pool_fastqs se usan como wildcards.
# la constraint siguiente ignora \\. en sample para asegurar que entiende sample
# como el nombre de la muestra (no como los files intermedios i.e, sample.bam) 
wildcard_constraints:
    sample="[^.]+",
    pool="POOL-[0-9A-Za-z]+"

# opcion paralelizar todas las carpetas
#POOLS = config["pools"]

# opcion para inputar la carpeta en el comando y paralelizar files
# (mejor por si falla el pipeline)
POOLS = config.get("pools", [])
if "pool" in config:
    POOLS = [config["pool"]]

FASTQ_DIR = "data/fastq"
ALIGN_DIR = "data/alignments"

os.makedirs(ALIGN_DIR, exist_ok=True)

# ---------------------------------------------------
# funcion para localizar fastq (lo comun al directorio
# data/fastq son las carpetas Stats y Reports. El directorio
# con fastq files tiene nombres distintos segun datos introducidos al secuenciar)
# ---------------------------------------------------
# la function genera:
# samples: lista con nombres de las muestras (L-25000010_S5_R1_001.fastq.gz -> L-25000010_S5_R1_001)
# fastqs: diccionario con directorio de fastqs y los fastq files

def load_fastqs_for_pool(pool):

    pool_path = os.path.join(FASTQ_DIR, pool)

    # find valid FASTQ subfolder (ignore Stats/Reports)
    subfolders = [
        f for f in glob.glob(os.path.join(pool_path, "*"))
        if os.path.isdir(f) and os.path.basename(f) not in ["Stats", "Reports"]
    ]

    if not subfolders:
        raise ValueError(f"No valid FASTQ subfolder found for pool {pool}")

    fastq_dir = subfolders[0]
 
    samples = [] # lista de muestras
    fastqs = {} # diccionario de fastq

    for r1 in glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")):
        samplename = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
        r2 = os.path.join(fastq_dir, samplename + "_R2_001.fastq.gz")

        if not os.path.exists(r2):
            raise ValueError(f"Missing R2 for sample {samplename} in pool {pool}")

        samples.append(samplename)
        fastqs[(pool, samplename)] = (r1, r2)

    return samples, fastqs

# ---------------------------------------------------
# definir pool-muestras-fastq para input del alineamiento
# y comprobar el output
# ---------------------------------------------------
pool_samples = {} # diccionario con lista de muestras por pool
pool_fastqs = {} # diccionario de con key=pool-muestra con values R1.fastq y R2.fastq

for pool in POOLS:
    # genera lista y dicc
    samples, fastqs = load_fastqs_for_pool(pool)
    # diccionario pool_samples
    pool_samples[pool] = samples
    # anyade R1 y R2 fastqs como valores al diccionario pool_fastqs
    pool_fastqs.update(fastqs)

print("Detected samples per pool:")
print(pool_samples)

# ---------------------------------------------------
# output final: .markdup.bam.bai
# ---------------------------------------------------
final_bai = [
    os.path.join(ALIGN_DIR, pool, f"{sample}.markdup.bam.bai")
    for pool, samples in pool_samples.items()
    for sample in samples
]

# el proceso termina cuando esten todos los archivos .bai
# para todas las muestras de un pool
rule all:
    input:
        final_bai

# ---------------------------------------------------
# alineamiento
# ---------------------------------------------------
# bam temporal para marcar duplicados en la siguiente regla
rule align_bwa:
    input:
        # primer value del diccionario pool_fastqs (R1)
        r1 = lambda wc: pool_fastqs[(wc.pool, wc.sample)][0],
        # segundo value del diccionario pool_fastqs (R2)
        r2 = lambda wc: pool_fastqs[(wc.pool, wc.sample)][1],
        # referencia indexada previamiente - ver README
        ref = config["reference"]
    output:
        bam = temp(os.path.join(ALIGN_DIR, "{pool}", "{sample}.bam"))
    threads: config["threads"]["align"]
    shell:
        """
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        mkdir -p {ALIGN_DIR}/{wildcards.pool}

        bwa-mem2 mem -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
            {input.ref} {input.r1} {input.r2} |
        samtools sort -@ {threads} -T $tmpdir/{wildcards.sample} -o {output.bam}
        """

# ---------------------------------------------------
# marcar lecturas duplicadas
# ---------------------------------------------------
rule mark_duplicates:
    input:
        bam = os.path.join(ALIGN_DIR, "{pool}", "{sample}.bam")
    output:
        bam = os.path.join(ALIGN_DIR, "{pool}", "{sample}.markdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{pool}", "{sample}.markdup.bam.bai"),
        metrics = os.path.join(ALIGN_DIR, "{pool}", "{sample}.dup_metrics.txt")
    threads:
        config["threads"]["dup"]
    shell:
        """
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        picard -Xmx2g MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=false \
            TMP_DIR=$tmpdir

        samtools index {output.bam}
        """

