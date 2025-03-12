process APPLY_BQSR_INTERVALS {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file), path(bai), path(recal_table)
	path intervals

	output:
	tuple val(sample_id), path("${sample_id}.recal.bam"), emit: recal_bam

	script:
	"""
	gatk ApplyBQSR \
		-I ${bam_file} \
		--bqsr-recal-file ${recal_table} \
		--intervals ${intervals} \
		--interval-padding 100 \
		-O ${sample_id}.recal.bam
	"""
}
