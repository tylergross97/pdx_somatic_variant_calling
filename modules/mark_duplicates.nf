process MARK_DUPLICATES {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_markduplicates, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file)

	output:
	tuple val(sample_id), path("${sample_id}.marked_duplicates.bam"), emit: marked_dup_bam
	path "${sample_id}.marked_dup_metrics.txt", emit: marked_dup_metrics

	script:
	"""
	gatk MarkDuplicates \
		-I ${bam_file} \
		-O ${sample_id}.marked_duplicates.bam \
		-M ${sample_id}.marked_dup_metrics.txt
	"""
}
