process SORT_BAM{
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bam, mode: 'copy'

	input:
	tuple val(sample_id), path(bam_file)
	
	output:
	tuple val(sample_id), path("${sample_id}.*.namesorted.bam"), emit: sorted_bam

	script:
	"""
	samtools sort -n -@ 8 -O bam -o ${bam_file.baseName}.namesorted.bam ${bam_file}
	"""
}
