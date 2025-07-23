process MUTECT2_INTERVALS {
	container "community.wave.seqera.io/library/gatk4_samtools:464a35b5e2f0c13d"
	publishDir params.outdir_mutect2, mode: "copy"

	input:
	path fasta
	path fasta_index
	path dict_file
	tuple val(sample_id), path(recal_bam)
	path germline_resource
	path germline_resource_idx
	path panel_of_normals
	path panel_of_normals_idx
	path intervals

	output:
	tuple val(sample_id), path("${sample_id}.unfiltered.vcf.gz"), path("${sample_id}.unfiltered.vcf.gz.stats"), emit: unfiltered_vcf

	script:
	"""
	samtools index ${recal_bam}
	gatk Mutect2 \
		-R ${fasta} \
		-I ${recal_bam} \
		--germline-resource ${germline_resource} \
		--panel-of-normals ${panel_of_normals} \
		--intervals ${intervals} \
		--interval-padding 100 \
		-O ${sample_id}.unfiltered.vcf.gz
	"""
}
