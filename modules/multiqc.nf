process MULTIQC {
	container "quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0"
	publishDir params.outdir_fastp, mode: 'copy'

	input:
	path '*'

	output:
	path 'multiqc_report.html', emit: report
	path 'multiqc_data', emit: data

	script:
	"""
	multiqc .
	"""
}
