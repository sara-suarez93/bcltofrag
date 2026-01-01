import os

configfile: "config.yaml"

POOLS = config["pools"]
RUNS_DIR = config["runs_dir"]

FASTQ_DIR = config["fastq_dir"]

# creo carpetas para los outputs
os.makedirs(FASTQ_DIR, exist_ok=True)


# ---------------------------------------------------
# output final: done.flag (archivo que marca final regla bcl2fastq)
# ---------------------------------------------------
rule all:
    input:
        expand(
            os.path.join(FASTQ_DIR, "{pool}", "done.flag"),
            pool=POOLS
        )

# ---------------------------------------------------
# bcl2fastq: demultiplex las muestras en cada pool y convertir en FASTQ
# ---------------------------------------------------
rule bcl2fastq:
    input:
        runfolder = lambda wc: os.path.join(RUNS_DIR, wc.pool)
    output:
        # flag file indica que se ha completado el demultiplexing de un POOL
        touch(os.path.join(FASTQ_DIR, "{pool}", "done.flag"))
        
    threads: config["threads"]["bcl2fastq"]
    shell:
        """
        mkdir -p {FASTQ_DIR}

        bcl2fastq \
            --runfolder-dir {input.runfolder} \
            --output-dir    {FASTQ_DIR}/{wildcards.pool} \
            --no-lane-splitting \
            --processing-threads {threads} \
            --use-bases-mask Y*,I8,N10,Y*

        touch {output}
        """
