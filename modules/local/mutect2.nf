process MUTECT2 {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_mutect2, mode: "copy"

	input:
	path fasta
	path fasta_index
	path dict_file
	tuple val(sample_id), path(recal_bam), path(recal_bam_idx)
	path germline_resource
	path germline_resource_idx
	path panel_of_normals
	path panel_of_normals_idx

	output:
	tuple val(sample_id), path("${sample_id}.unfiltered.vcf.gz"), path("${sample_id}.unfiltered.vcf.gz.stats"), emit: unfiltered_vcf

	script:
	"""
	gatk Mutect2 \
		-R ${fasta} \
		-I ${recal_bam} \
		--germline-resource ${germline_resource} \
		--panel-of-normals ${panel_of_normals} \
		-O ${sample_id}.unfiltered.vcf.gz
	"""
}
