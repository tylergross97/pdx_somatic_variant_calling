process FUNCOTATOR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"

	input:
	tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_idx)
	path reference_fasta
	path reference_idx
	path reference_dict
	path funcotator_data_sources

	output:
	tuple val(sample_id), path("${sample_id}.filtered.annotated.vcf"), emit: filtered_annotated_vcf

	script:
	"""
	mkdir funcotator_dataSources_extracted
	tar -xzvf ${funcotator_data_sources} -C funcotator_dataSources_extracted
	gatk Funcotator \
                --variant ${filtered_vcf} \
                --reference ${reference_fasta} \
                --ref-version hg38 \
                --data-sources-path funcotator_dataSources_extracted/funcotator_dataSources.v1.8.hg38.20230908s \
                --output ${sample_id}.filtered.annotated.vcf \
                --output-file-format VCF
	"""
}
