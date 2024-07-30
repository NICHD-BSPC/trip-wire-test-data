
# See README for details on what data to place in these directories and how to set config options
full_dataset = "full_dataset/"
fastq = "fastq/"
read = ['R1', 'R2']
pool = ['a', 'b']
techrep = ['1', '2']
condition = ['treatment', 'control']
library = ["mapping", "normalization", "expression"]


rule all:
    input:
        r1="results/test_R1.fastq.gz",
        r2="results/test_R2.fastq.gz"


rule getBarcodes:
    """
    Generate all chr19-only barcodes from each pool, extracting them from the
    Chromosome column ('Chr') of the final_TRIP_data_table.txt.
    """
    input:
        expand(full_dataset + "trip_output/pool{{p}}_techrep{t}_{c}/final_TRIP_data_table.txt", t=techrep, c=condition, full_dataset=full_dataset)
    output:
        temporary("results/processing/barcodes_pool{p}.txt")
    resources:
        mem_mb=2 * 1024,
        runtime=60
    shell:
        "cat {input} | grep -w chr19 | cut -f 1 | sort | uniq > {output}"


rule makeFasta:
    """
    Generate a fasta file of chr19-only barcodes.
    """
    input:
        "results/processing/barcodes_pool{p}.txt"
    output:
        temporary("results/processing/barcodes_pool{p}.fasta")
    resources:
        mem_mb=2 * 1024,
        runtime=60
    run:
        with open(str(input), "r") as fin:
            with open(str(output), "w") as fout:
                for i, barcode in enumerate(fin):
                    fout.write(f">{i}\n")
                    fout.write(barcode)


rule cutadaptGetIDMapp:
    """
    Extract the Read ID from the chr19 barcodes fasta for the mapping library
    using cutadapt. Only the R1 read will be extracted.
    """
    input:
        barcodes="results/processing/barcodes_pool{p}.fasta",
        r1=full_dataset+"demux_data/pool{p}_{c}_mapping_R1.fastq.gz",
    output:
        r1o=temporary("results/mapping/r1cutadapt_mapping_pool{p}_{c}.fastq"),
        readID=temporary("results/mapping/readID_mapping_pool{p}_{c}.txt"),
    resources:
        mem_mb=4 * 1024,
        runtime=120
    threads: 16
    log:
        "logs/cutadapt_mapping_pool{p}_{c}.log"
    shell:
        "cutadapt "
        "--action=none "
        "-a file:{input.barcodes} "
        "--discard-untrimmed "
        "-o {output.r1o} "
        "--error-rate 0.0 "
        "--overlap 16 "
        "--cores {threads} "
        "{input.r1} "
        "&> {log}; "
        "awk 'NR%4==1' {output.r1o} | cut -d ' ' -f 1 | sed 's/^@//g' > {output.readID}; "


rule extractReadsWithIDMapp:
    """
    For the mapping library, extract the full read sequence using the ReadID
    from the full dataset. Both R1 and R2 reads are extracted.
    """
    input:
        readID="results/mapping/readID_mapping_pool{p}_{c}.txt",
        fastq=fastq+"Undetermined_R{r}.fastq.gz",
    output:
        temporary("results/mapping/mappingRead_pool{p}_{c}_R{r}.fastq"),
    resources:
        mem_mb=4 * 1024,
        runtime=120
    shell:
        "seqtk subseq {input.fastq} {input.readID} > {output};"


rule cutadaptGetIDNormExpr:
    """
    For the expression and normalization libraries, extract the FastQ
    information from the full dataset, using the Read ID associated with the
    chr19 barcode. The read IDs from only R1 are extracted.
    """
    wildcard_constraints:
        t="1|2",
        lib="expression|normalization"
    input:
        barcodes="results/processing/barcodes_pool{p}.fasta",
        r1=full_dataset+"demux_data/pool{p}_techrep{t}_{c}_{lib}_R1.fastq.gz",
    output:
        r1o=temporary("results/normexpr/r1cutadapt_techrep{t}_{lib}_pool{p}_{c}.fastq"),
        readID=temporary("results/normexpr/readID_techrep{t}_{lib}_pool{p}_{c}.txt"),
    log:
        "logs/cutadapt_techrep{t}_{lib}_pool{p}_{c}.log"
    resources:
        mem_mb=4 * 1024,
        runtime=120
    threads: 16
    shell:
        "cutadapt "
        "--action=none "
        "-a file:{input.barcodes} "
        "--discard-untrimmed "
        "-o {output.r1o} "
        "--error-rate 0.0 "
        "--overlap 16 "
        "--cores {threads} "
        "{input.r1} "
        "&> {log}; "
        "awk 'NR%4==1' {output.r1o} | cut -d ' ' -f 1 | sed 's/^@//g' > {output.readID}; "


rule extractReadsWithIDNormExpr:
    """
    For the expression & normalization libraries, use each chr19 barcode's
    identified read ID to extract the correponding read from the full dataset.
    """
    wildcard_constraints:
        t="1|2",
        lib="normalization|expression"
    input:
        readID="results/normexpr/readID_techrep{t}_{lib}_pool{p}_{c}.txt",
        fastq=fastq + "Undetermined_R{r}.fastq.gz",
    output:
        temporary("results/normexpr/{lib}Read_techrep{t}_pool{p}_{c}_R{r}.fastq"),
    resources:
        mem_mb=4 * 1024,
        runtime=120
    shell:
        "seqtk subseq {input.fastq} {input.readID} > {output};"


rule catAllReads:
    """
    Concatenate all reads from chr19 across all library types, pools,
    conditions, and techreps into a single output file (one for each of R1 and
    R2)"""
    wildcard_constraints:
        r="1|2"
    input:
        expand("results/mapping/{lib}Read_pool{p}_{c}_R{{r}}.fastq", p=pool, c=condition, lib=["mapping"]),
        expand("results/normexpr/{lib}Read_techrep{t}_pool{p}_{c}_R{{r}}.fastq", p=pool, c=condition, lib=["expression", "normalization"], t=techrep),
    output:
        "results/test_R{r}.fastq.gz"
    resources:
        mem_mb=4 * 1024,
        runtime=120
    shell:
        "cat {input} | gzip -c > {output}"
