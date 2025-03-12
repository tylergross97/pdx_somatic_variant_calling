process MERGE_BAMS {
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bamcmp, mode: 'copy'
	
	input:
	tuple val(sample_id), path(human_only), path(human_better)

	output:
	tuple val(sample_id), path("${sample_id}.human.bamcmp.coordinatesorted.bam"), emit: human_merged_sorted_bam

	script:
	"""
	samtools merge ${sample_id}.human.bamcmp.bam ${human_only} ${human_better}
	samtools sort -@ 8 -O bam -o ${sample_id}.human.bamcmp.coordinatesorted.bam ${sample_id}.human.bamcmp.bam
	rm ${sample_id}.human.bamcmp.bam
	"""
}
