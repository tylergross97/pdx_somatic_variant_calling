process BASE_RECALIBRATOR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file), path(bai_file)
	path fasta
	path fasta_index
	path dict_file
	path dbsnp_vcf
	path dbsnp_vcf_idx
	path known_indels_vcf
	path known_indels_vcf_idx
	path mills_1000G_vcf
	path mills_1000G_vcf_idx

	output:
	tuple val(sample_id), path("${sample_id}.recal_data.table"), emit: recal_data_table

	script:
	"""
	gatk BaseRecalibrator \
		-I ${bam_file} \
		-R ${fasta} \
		--known-sites ${dbsnp_vcf} \
		--known-sites ${known_indels_vcf} \
		--known-sites ${mills_1000G_vcf} \
		-O ${sample_id}.recal_data.table
	"""
}
