process FILTER_HAMA_VARIANTS {
	container "quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
	publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	tuple val(sample_id), path(filtered_vcf)
	path hama_bed
	
	output:
	tuple val(sample_id), path("${sample_id}.filtered.hama_filtered.vcf.gz"), emit: hama_filtered_vcf
	
	script:
	"""
	# Extract HAMA variant IDs from BED file
	cut -f4 ${hama_bed} > hama_variant_ids.txt
	
	# Create a filtered VCF excluding HAMA variants
	bcftools view -h ${filtered_vcf} > header.vcf
	
	# Process the VCF content and create variant IDs
	bcftools view -H ${filtered_vcf} | awk '{OFS="\t"; variant_id=\$1"_"\$2"_"\$4"_"\$5; print variant_id,\$0}' > variants_with_id.txt
	
	# Filter out variants that match HAMA IDs
	grep -v -f hama_variant_ids.txt variants_with_id.txt | cut -f2- > filtered_content.vcf || true
	
	# If no variants were filtered (grep returns non-zero if no matches), create an empty file
	if [ ! -s filtered_content.vcf ]; then
		touch filtered_content.vcf
	fi
	
	# Combine header and filtered content
	cat header.vcf filtered_content.vcf | bgzip > ${sample_id}.filtered.hama_filtered.vcf.gz
	"""
}
