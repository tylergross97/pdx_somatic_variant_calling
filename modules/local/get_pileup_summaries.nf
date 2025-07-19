process GET_PILEUP_SUMMARIES {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	tuple val(sample_id), path(recal_bam)
	path filtered_vcf
	path filtered_vcf_idx
	
	output:
	tuple val(sample_id), path("${sample_id}.pileups.table"), emit: pileup_table

	script:
	"""
	gatk GetPileupSummaries \
		-I ${recal_bam} \
		-V ${filtered_vcf} \
		-L ${filtered_vcf} \
		-O ${sample_id}.pileups.table
	"""
}
