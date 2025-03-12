process ALIGN_HUMAN {
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bam, mode: 'copy'
	
	input:
	tuple val(sample_id), path(trimmed_fastq_files)
	tuple val(genome), path(index_files)

	output:
	tuple val(sample_id), path("${sample_id}.human.bam"), emit: bam

	script:
	def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
	"""
	bwa-mem2 mem \
	-t 8 \
	-T 0 \
	-R "${read_group}" \
	human \
	${trimmed_fastq_files[0]} \
	${trimmed_fastq_files[1]} |
	samtools view \
	-Shb \
	-o ${sample_id}.human.bam -
	"""
}
