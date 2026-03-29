#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────
// Parameters  (override via --param value or nextflow.config)
// ─────────────────────────────────────────────
params.data        = "/home/a-m/fazal2/RNAseq_workflow"
params.type        = "PE"          // "SE" or "PE"
params.fasta       = "/home/a-m/fazal2/xenograft/data/genome/human_mouse.fna"
params.gtf         = "/home/a-m/fazal2/xenograft/data/genome/human_mouse.gtf"
params.junction    = 99
params.strand      = 2
params.trim_length = 30
params.trim_window = "4:18"
params.trim_lead   = 28
params.trim_trail  = 28
params.outdir      = "results"

// ─────────────────────────────────────────────
// STAR genome index
// ─────────────────────────────────────────────
process STAR_INDEX {
    label 'star'

    output:
    path "genome/", emit: genome_index

    script:
    """
    mkdir -p genome
    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir genome/ \\
        --genomeFastaFiles ${params.fasta} \\
        --limitGenomeGenerateRAM 32000000000 \\
        --outTmpDir /scratch/star_index_\${SLURM_JOB_ID:-$$}
    """
}

// ─────────────────────────────────────────────
// Trimmomatic – Single-End
// ─────────────────────────────────────────────
process TRIM_SE {
    label 'trim'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed

    script:
    """
    trimmomatic SE \\
        -threads 2 -phred33 \\
        ${reads} \\
        ${sample_id}.trimmed.fastq.gz \\
        ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-SE.fa:2:30:10 \\
        SLIDINGWINDOW:${params.trim_window} \\
        LEADING:${params.trim_lead} \\
        TRAILING:${params.trim_trail} \\
        MINLEN:${params.trim_length}
    """
}

// ─────────────────────────────────────────────
// Trimmomatic – Paired-End
// ─────────────────────────────────────────────
process TRIM_PE {
    label 'trim'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed

    script:
    """
    trimmomatic PE \\
        -threads 2 -phred33 \\
        ${r1} ${r2} \\
        ${sample_id}_R1.trimmed.fastq.gz  ${sample_id}_R1un.trimmed.fastq.gz \\
        ${sample_id}_R2.trimmed.fastq.gz  ${sample_id}_R2un.trimmed.fastq.gz \\
        ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-PE.fa:2:30:10 \\
        SLIDINGWINDOW:${params.trim_window} \\
        LEADING:${params.trim_lead} \\
        TRAILING:${params.trim_trail} \\
        MINLEN:${params.trim_length}
    """
}

// ─────────────────────────────────────────────
// STAR mapping – Single-End
// ─────────────────────────────────────────────
process MAPPING_SE {
    label 'star'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(trimmed)
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --genomeDir ${genome_index} \\
        --readFilesIn ${trimmed} \\
        --sjdbGTFfile ${params.gtf} \\
        --sjdbOverhang ${params.junction} \\
        --outFileNamePrefix ${sample_id}_ \\
        --limitGenomeGenerateRAM 60000000000 \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --outTmpDir /scratch/star_map_\${SLURM_JOB_ID:-$$}_${sample_id}
    """
}

// ─────────────────────────────────────────────
// STAR mapping – Paired-End
// ─────────────────────────────────────────────
process MAPPING_PE {
    label 'star'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --genomeDir ${genome_index} \\
        --readFilesIn ${r1} ${r2} \\
        --sjdbGTFfile ${params.gtf} \\
        --sjdbOverhang ${params.junction} \\
        --outFileNamePrefix ${sample_id}_ \\
        --limitGenomeGenerateRAM 60000000000 \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --outTmpDir /scratch/star_map_\${SLURM_JOB_ID:-$$}_${sample_id}
    """
}

// ─────────────────────────────────────────────
// featureCounts – Single-End
// ─────────────────────────────────────────────
process COUNTS_SE {
    label 'counts'
    tag "${sample_id}"

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}_featcounts.txt"

    script:
    """
    featureCounts \\
        -T 4 \\
        -s ${params.strand} \\
        -g gene_id \\
        -t exon \\
        -o ${sample_id}_featcounts.txt \\
        -a ${params.gtf} \\
        ${bam}
    """
}

// ─────────────────────────────────────────────
// featureCounts – Paired-End
// ─────────────────────────────────────────────
process COUNTS_PE {
    label 'counts'
    tag "${sample_id}"

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}_featcounts.txt"

    script:
    """
    featureCounts \\
        -T 4 \\
        -p \\
        -s ${params.strand} \\
        -g gene_id \\
        -t exon \\
        -o ${sample_id}_featcounts.txt \\
        -a ${params.gtf} \\
        ${bam}
    """
}

// ─────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────
workflow {

    if (params.type != "SE" && params.type != "PE") {
        error 'Please specify only "SE" or "PE" for the --type parameter.'
    }

    // Build genome index once
    STAR_INDEX()
    genome_index = STAR_INDEX.out.genome_index

    if (params.type == "SE") {

        // Collect SE reads
        reads_ch = Channel
            .fromPath("${params.data}/*.fastq.gz")
            .map { file -> [ file.baseName.replaceAll(/\.fastq$/, ''), file ] }

        TRIM_SE(reads_ch)
        MAPPING_SE(TRIM_SE.out.trimmed, genome_index)
        COUNTS_SE(MAPPING_SE.out.bam)

    } else {

        // Collect PE reads (matched by sample name)
        reads_ch = Channel
            .fromFilePairs("${params.data}/*_R{1,2}.fastq.gz", flat: true)
            .map { sample_id, r1, r2 -> [ sample_id, r1, r2 ] }

        TRIM_PE(reads_ch)
        MAPPING_PE(TRIM_PE.out.trimmed, genome_index)
        COUNTS_PE(MAPPING_PE.out.bam)
    }
}
