process CREATE_DICT {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'copy'

	input:
	path fasta

	output:
	path "${fasta.baseName}.dict", emit: dict_file

	script:
	"""
	gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
	"""
}
