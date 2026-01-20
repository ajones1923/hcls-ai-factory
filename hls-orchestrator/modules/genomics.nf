/*
 * Genomics Pipeline Module
 * FASTQ → VCF variant calling
 */

process GENOMICS_FASTQC {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip
    path "versions.yml", emit: versions

    script:
    """
    fastqc -t ${task.cpus} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}

process GENOMICS_BWA_MEM {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${index}/genome.fa \\
        ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

    samtools index ${sample_id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process GENOMICS_GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(bam)
    path reference
    path known_sites

    output:
    tuple val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    gatk HaplotypeCaller \\
        -R ${reference}/genome.fa \\
        -I ${bam} \\
        -O ${sample_id}.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | head -1 | sed 's/.*GATK v//')
    END_VERSIONS
    """
}

process GENOMICS_VCF_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/genomics/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("*.annotated.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # Annotate with VEP or SnpEff
    snpEff -Xmx${task.memory.toGiga()}g \\
        GRCh38.99 \\
        ${vcf} \\
        | bgzip -c > ${sample_id}.annotated.vcf.gz

    tabix -p vcf ${sample_id}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | head -1)
    END_VERSIONS
    """
}

// Main genomics workflow
workflow GENOMICS_PIPELINE {
    take:
    reads  // tuple(sample_id, [read1, read2])

    main:
    ch_versions = Channel.empty()

    // Quality control
    GENOMICS_FASTQC(reads)
    ch_versions = ch_versions.mix(GENOMICS_FASTQC.out.versions)

    // Alignment
    GENOMICS_BWA_MEM(
        reads,
        params.bwa_index ? Channel.fromPath(params.bwa_index) : Channel.empty()
    )
    ch_versions = ch_versions.mix(GENOMICS_BWA_MEM.out.versions)

    // Variant calling
    GENOMICS_GATK_HAPLOTYPECALLER(
        GENOMICS_BWA_MEM.out.bam,
        params.bwa_index ? Channel.fromPath(params.bwa_index) : Channel.empty(),
        params.known_sites ? Channel.fromPath(params.known_sites) : Channel.empty()
    )
    ch_versions = ch_versions.mix(GENOMICS_GATK_HAPLOTYPECALLER.out.versions)

    // Annotation
    GENOMICS_VCF_ANNOTATE(GENOMICS_GATK_HAPLOTYPECALLER.out.vcf)
    ch_versions = ch_versions.mix(GENOMICS_VCF_ANNOTATE.out.versions)

    emit:
    vcf = GENOMICS_VCF_ANNOTATE.out.vcf
    versions = ch_versions
}
