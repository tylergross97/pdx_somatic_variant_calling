process FILTER_HAMA_VARIANTS {
    container "quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
    publishDir params.outdir_mutect2, mode: "copy"
    
    input:
    tuple val(sample_id), path(filtered_vcf)
    path hama_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.hama_annotated.vcf.gz"), emit: hama_filtered_vcf
    
    script:
    """
    # Compress and index the BED file if not already compressed/indexed
    if [ ! -f ${hama_bed}.gz ]; then
      bgzip -c ${hama_bed} > ${hama_bed}.gz
    fi
    
    if [ ! -f ${hama_bed}.gz.tbi ]; then
      tabix -p bed ${hama_bed}.gz
    fi

    # Create custom header for bcftools annotate
    echo '##INFO=<ID=HAMA_ID,Number=1,Type=String,Description="Variant overlaps high-risk HAMA BED region; ID from BED column 4">' > custom_header.txt

    # Run bcftools annotate to add HAMA annotation from BED file
    bcftools annotate \
      -a ${hama_bed}.gz \
      -h custom_header.txt \
      -c CHROM,FROM,TO,INFO/HAMA_ID \
      -o ${sample_id}.filtered.hama_annotated.vcf.gz \
      -O z \
      ${filtered_vcf}

    # Index the resulting annotated VCF
    tabix -p vcf ${sample_id}.filtered.hama_annotated.vcf.gz
    """
}
