process BAMCMP {
	container "community.wave.seqera.io/library/bamcmp:2.2--0b3144433037098f"
	publishDir params.outdir_bamcmp, mode: 'copy'

	input:
	tuple val(sample_id), path(bam_human), path(bam_mouse)

	output:
	tuple val(sample_id), path("${sample_id}_humanOnly.bam"), emit: human_only, optional: true
	tuple val(sample_id), path("${sample_id}_humanBetter.bam"), emit: human_better, optional: true
	tuple val(sample_id), path("${sample_id}_mouseOnly.bam"), emit: mouse_only, optional: true
	tuple val(sample_id), path("${sample_id}_mouseBetter.bam"), emit: mouse_better, optional: true
	tuple val(sample_id), path("${sample_id}_humanLoss.bam"), emit: human_loss, optional: true
	tuple val(sample_id), path("${sample_id}_mouseLoss.bam"), emit: mouse_loss, optional: true

	script:
	"""
	bamcmp -n -1 ${bam_human} -2 ${bam_mouse} \
		-a ${sample_id}_humanOnly.bam \
		-A ${sample_id}_humanBetter.bam \
		-b ${sample_id}_mouseOnly.bam \
		-B ${sample_id}_mouseBetter.bam \
		-C ${sample_id}_humanLoss.bam \
		-D ${sample_id}_mouseLoss.bam \
		-s as
	"""
}
