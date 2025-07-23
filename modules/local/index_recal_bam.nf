process INDEX_BAM {
    container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
    publishDir params.outdir_bqsr, mode: 'symlink'

    input:
        tuple val(sample_id), path(bam_file)
    output:
        tuple val(sample_id), path("${bam_file.getBaseName()}.bai"), emit: bam_index
    script:
    """
    samtools index ${bam_file}
    """
}