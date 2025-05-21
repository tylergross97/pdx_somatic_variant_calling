process INDEX_FILTERED_VCF {
	container "quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
        publishDir params.outdir_mutect2, mode: "symlink"

        input:
        tuple val(sample_id), path(filtered_vcf)

        output:
        tuple val(sample_id), path("${filtered_vcf}.tbi"), emit: filtered_vcf_idx

        script:
        """
        bcftools index -t ${filtered_vcf}
        """
}
