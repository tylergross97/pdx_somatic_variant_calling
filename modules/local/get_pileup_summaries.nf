process GET_PILEUP_SUMMARIES {
	container "community.wave.seqera.io/library/gatk4_samtools:464a35b5e2f0c13d"
	publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	tuple val(sample_id), path(recal_bam)
	path filtered_vcf
	path filtered_vcf_idx
	
	output:
	tuple val(sample_id), path("${sample_id}.pileups.table"), emit: pileup_table

	script:
	"""
	samtools index ${recal_bam}
	gatk GetPileupSummaries \
		-I ${recal_bam} \
		-V ${filtered_vcf} \
		-L ${filtered_vcf} \
		-O ${sample_id}.pileups.table
	"""
}
