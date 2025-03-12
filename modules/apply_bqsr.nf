process APPLY_BQSR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file), path(recal_table)

	output:
	tuple val(sample_id), path("${sample_id}.recal.bam"), emit: recal_bam

	script:
	"""
	gatk ApplyBQSR \
		-I ${bam_file} \
		--bqsr-recal-file ${recal_table} \
		-O ${sample_id}.recal.bam
	"""
}
