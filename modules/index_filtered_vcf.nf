process INDEX_FILTERED_VCF {
	container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
        publishDir params.outdir_mutect2, mode: "symlink"

        input:
        tuple val(sample_id), path(filtered_vcf)

        output:
        tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"), emit: filtered_vcf_idx

        script:
        """
        bcftools index -t ${filtered_vcf}
        """
}
