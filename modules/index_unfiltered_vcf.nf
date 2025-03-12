process INDEX_UNFILTERED_VCF {
	container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
	publishDir params.outdir_mutect2, mode: "symlink"

	input:
	tuple val(sample_id), path(unfiltered_vcf), path(unfiltered_vcf_stats)

	output:
	tuple val(sample_id), path("${sample_id}.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf_idx

	script:
	"""
	bcftools index -t ${unfiltered_vcf}
	"""
}
