process FILTER_MUTECT_CALLS {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	path fasta
	path fasta_idx
	path dict_file
	tuple val(sample_id), path(unfiltered_vcf), path(unfiltered_vcf_stats), path(contamination_table), path(tumor_segmentation_table), path(unfiltered_vcf_idx)

	output:
	tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), emit: filtered_vcf

	script:
	"""
	gatk FilterMutectCalls \
		-R ${fasta} \
		-V ${unfiltered_vcf} \
		--contamination-table ${contamination_table} \
		--tumor-segmentation ${tumor_segmentation_table} \
		-O ${sample_id}.filtered.vcf.gz
	""" 
}
