process INDEX_REFERENCE {
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bqsr, mode: 'copy'

	input:
	path fasta
	
	output:
	path "${fasta}", emit: fasta_file
	path "${fasta}.fai", emit: fasta_index

	script:
	"""
	samtools faidx ${fasta}
	"""
}
