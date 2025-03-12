process INDEX_BAM {
        container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
        publishDir params.outdir_bamcmp, mode: 'copy'

	input:
	tuple val(sample_id), path(bam_file)

	output:
	tuple val(sample_id), path(bam_file), path("${bam_file}.bai"), emit: indexed_bam

	script:
	"""
	samtools index ${bam_file}
	"""
}
