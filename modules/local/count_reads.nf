process COUNT_READS {
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bamcmp, mode: 'copy'

	input:
	tuple val(sample_id), path(bams)
	
	output:
	path "${sample_id}_read_counts.txt", emit: read_counts

	script:
	"""
	echo "${sample_id}" > ${sample_id}_read_counts.txt
	for bam in ${bams}; do
		count=\$(samtools view -c \$bam)
		echo "\$bam \$count" >> ${sample_id}_read_counts.txt
	done
	"""
}
