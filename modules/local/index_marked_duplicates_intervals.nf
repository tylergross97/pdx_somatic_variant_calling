process INDEX_MARKED_DUPLICATES_INTERVALS {
	container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"
    	publishDir params.outdir_markduplicates, mode: 'symlink'

    	input:
    	tuple val(sample_id), path(bam_file)

    	output:
    	tuple val(sample_id), path(bam_file), path("${bam_file}.bai"), emit: indexed_bam

    	script:
    	"""
    	samtools index ${bam_file}
    	"""
}
