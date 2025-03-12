process FASTP {
	container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
	publishDir params.outdir_fastp, mode: 'copy'

	input:
	tuple val(sample_id), path(fastq_files)

	output:
	tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz"), emit: trimmed_reads
	path "${sample_id}_fastp.json", emit: json_report
	path "${sample_id}_fastp.html", emit: html_report

	script:
	"""
	fastp \
		-i ${fastq_files[0]} \
		-I ${fastq_files[1]} \
		-o ${sample_id}_trimmed_1.fastq.gz \
		-O ${sample_id}_trimmed_2.fastq.gz \
		-j ${sample_id}_fastp.json \
		-h ${sample_id}_fastp.html
	"""
}
