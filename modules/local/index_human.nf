process INDEX_HUMAN {
	container "quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_6"
	publishDir params.outdir_index, mode: 'symlink'
	
	input:
	path fasta

	output:
	tuple val('human'), path("human.*")

	script:
	"""
	bwa-mem2 index -p human ${fasta}
	"""
}
