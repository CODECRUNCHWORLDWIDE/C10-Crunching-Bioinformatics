/*
 * main.nf - C10 Week 12 Exercises (toy RNA-seq quantification, Nextflow)
 *
 * Educational and research use only. Outputs are not clinical software.
 *
 * This Nextflow workflow is the worked reference paired with the
 * Snakefile in the same directory. It runs end to end on the same
 * small didactic dataset: 4 samples, 2 conditions, ~100k reads per
 * sample, a chr21-only transcriptome. Total runtime on a 4-core
 * laptop: 3-6 minutes.
 *
 * Processes (top-down DAG):
 *
 *   TRIM_FASTQ ── (per sample) ──> SALMON_QUANT ──> AGGREGATE_TPM
 *                                                 └> DESEQ2
 *
 *   SALMON_INDEX produces the index once; SALMON_QUANT consumes it.
 *
 * Required tools (declare in nextflow.config or environment.yml):
 *   fastp 0.23.4
 *   salmon 1.10.2
 *   bioconductor-tximport 1.30.0
 *   bioconductor-deseq2 1.42.0
 *
 * Run:
 *   nextflow run main.nf -profile conda
 *
 * Render DAG:
 *   nextflow run main.nf -with-dag dag.svg -preview
 */

nextflow.enable.dsl = 2

params.sample_sheet  = "samples.tsv"
params.fastq_dir     = "data/fastq"
params.transcriptome = "ref/transcripts.fa.gz"
params.tx2gene       = "ref/tx2gene.tsv"
params.outdir        = "results"


workflow {
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id,
                row.condition ?: 'unknown',
                file("${params.fastq_dir}/${row.sample_id}_R1.fastq.gz"),
                file("${params.fastq_dir}/${row.sample_id}_R2.fastq.gz")
            )
        }

    transcriptome_ch = Channel.value(file(params.transcriptome))
    tx2gene_ch       = Channel.value(file(params.tx2gene))

    TRIM_FASTQ(samples_ch)
    SALMON_INDEX(transcriptome_ch)
    SALMON_QUANT(TRIM_FASTQ.out, SALMON_INDEX.out)
    AGGREGATE_TPM(SALMON_QUANT.out.collect())
    DESEQ2(SALMON_QUANT.out.collect(), file(params.sample_sheet), tx2gene_ch)
}


process TRIM_FASTQ {
    tag "${sample_id}"
    cpus 4
    memory '4 GB'
    time '30.min'

    conda 'bioconda::fastp=0.23.4'

    input:
    tuple val(sample_id), val(condition), path(r1), path(r2)

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}_R1.trim.fastq.gz"),
          path("${sample_id}_R2.trim.fastq.gz")

    script:
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${sample_id}_R1.trim.fastq.gz \\
        -O ${sample_id}_R2.trim.fastq.gz \\
        -h ${sample_id}.fastp.html \\
        -j ${sample_id}.fastp.json \\
        --qualified_quality_phred 20 \\
        --length_required 25 \\
        --thread ${task.cpus}
    """
}


process SALMON_INDEX {
    cpus 4
    memory '8 GB'
    time '20.min'

    conda 'bioconda::salmon=1.10.2'

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index -t ${transcriptome} -i salmon_index -k 31 -p ${task.cpus}
    """
}


process SALMON_QUANT {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'
    time '30.min'

    conda 'bioconda::salmon=1.10.2'

    input:
    tuple val(sample_id), val(condition), path(r1), path(r2)
    path salmon_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_quant")

    script:
    """
    salmon quant -i ${salmon_index} -l A \\
        -1 ${r1} -2 ${r2} \\
        -o ${sample_id}_quant \\
        -p ${task.cpus} \\
        --seqBias --gcBias --validateMappings
    """
}


process AGGREGATE_TPM {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 1
    memory '2 GB'
    time '10.min'

    input:
    path quant_tuples

    output:
    path 'tpm_matrix.tsv'

    script:
    """
    python3 <<EOF
    import csv
    import os
    import glob

    dirs = sorted(glob.glob('*_quant'))
    tpm = {}
    samples = []
    for d in dirs:
        sample = d.replace('_quant', '')
        samples.append(sample)
        with open(os.path.join(d, 'quant.sf')) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                tpm.setdefault(row['Name'], {})[sample] = row['TPM']

    with open('tpm_matrix.tsv', 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['transcript_id'] + samples)
        for tx in sorted(tpm):
            writer.writerow([tx] + [tpm[tx].get(s, '0.0') for s in samples])
    EOF
    """
}


process DESEQ2 {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 1
    memory '4 GB'
    time '20.min'

    conda 'bioconda::bioconductor-tximport=1.30.0 bioconda::bioconductor-deseq2=1.42.0'

    input:
    path quant_tuples
    path sample_sheet
    path tx2gene

    output:
    path 'deseq2_results.rds'
    path 'dge_top30.tsv'
    path 'volcano.png'
    path 'pca.png'

    script:
    """
    Rscript - <<'EOF'
    suppressPackageStartupMessages({
        library(tximport)
        library(DESeq2)
    })

    sample_info <- read.table("${sample_sheet}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    tx2gene <- read.table("${tx2gene}", header = FALSE, sep = "\\t", stringsAsFactors = FALSE)

    files <- file.path(paste0(sample_info\$sample_id, "_quant"), "quant.sf")
    names(files) <- sample_info\$sample_id
    stopifnot(all(file.exists(files)))

    txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

    coldata <- data.frame(condition = factor(sample_info\$condition))
    rownames(coldata) <- sample_info\$sample_id

    dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, alpha = 0.05)
    saveRDS(res, file = "deseq2_results.rds")

    top30 <- as.data.frame(res[order(res\$padj), ])[1:30, ]
    top30\$gene_id <- rownames(top30)
    write.table(top30, "dge_top30.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    png("volcano.png", width = 800, height = 800, res = 120)
    with(as.data.frame(res), plot(log2FoldChange, -log10(padj),
        pch = 20, col = ifelse(padj < 0.05, "red", "grey"), main = "Volcano"))
    dev.off()

    vsd <- vst(dds, blind = FALSE)
    png("pca.png", width = 800, height = 800, res = 120)
    plotPCA(vsd, intgroup = "condition")
    dev.off()
    EOF
    """
}


workflow.onComplete {
    println "[main.nf] Workflow complete (status: ${workflow.success ? 'success' : 'failure'})."
    println "[main.nf] Run name:   ${workflow.runName}"
    println "[main.nf] Duration:   ${workflow.duration}"
    println "[main.nf] Work dir:   ${workflow.workDir}"
    println "[main.nf] Disclaimer: Educational and research use only."
}
