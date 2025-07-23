process APPLY_BQSR_INTERVALS {
	container "community.wave.seqera.io/library/gatk4_samtools:464a35b5e2f0c13d"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file), path(recal_table)
	path intervals

	output:
	tuple val(sample_id), path("${sample_id}.recal.bam"), emit: recal_bam

	script:
	"""
	samtools index ${bam_file}
	gatk ApplyBQSR \
		-I ${bam_file} \
		--bqsr-recal-file ${recal_table} \
		--intervals ${intervals} \
		--interval-padding 100 \
		-O ${sample_id}.recal.bam
	"""
}
