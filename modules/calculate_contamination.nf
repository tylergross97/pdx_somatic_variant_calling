process CALCULATE_CONTAMINATION {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"

	input:
	tuple val(sample_id), path(pileup_table)

	output:
	tuple val(sample_id), path("${sample_id}.contamination.table"), path("${sample_id}.segments.tsv"), emit: contamination_table

	script:
	"""
	gatk CalculateContamination \
		-I ${pileup_table} \
		-O ${sample_id}.contamination.table \
		--tumor-segmentation ${sample_id}.segments.tsv
	"""
}
