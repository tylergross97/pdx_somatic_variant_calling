#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

ch_fastq = Channel.fromFilePairs("${params.fastq}/*_S*_R{1,2}_001.fastq.gz", checkIfExists: true)

process ACCESSORY_FILES_DOWNLOAD {
    conda "conda-forge::curl=7.86.0"
    publishDir params.outdir_accessory_files, mode: 'copy'

    output:
    // dbSNP files
    path "Homo_sapiens_assembly38.dbsnp138.vcf", emit: dbsnp_vcf
    path "Homo_sapiens_assembly38.dbsnp138.vcf.idx", emit: dbsnp_vcf_idx
    // Known indels
    path "Homo_sapiens_assembly38.known_indels.vcf.gz", emit: known_indels
    path "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi", emit: known_indels_idx
    // Mills indels
    path "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", emit: mills_indels
    path "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", emit: mills_indels_idx
    // gnomAD
    path "af-only-gnomad.hg38.vcf.gz", emit: gnomad
    path "af-only-gnomad.hg38.vcf.gz.tbi", emit: gnomad_idx
    // Filtered VCF (ExAC common SNPs)
    path "small_exac_common_3.hg38.vcf.gz", emit: filtered_vcf
    path "small_exac_common_3.hg38.vcf.gz.tbi", emit: filtered_vcf_idx
    // Panel of Normals (PoN)
    path "1000g_pon.hg38.vcf.gz", emit: pon
    path "1000g_pon.hg38.vcf.gz.tbi", emit: pon_idx

    script:
    """
    # dbSNP files
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

    # Known indels
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

    # Mills indels
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

    # gnomAD
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

    # Filtered VCF (ExAC common SNPs)
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi

    # Panel of Normals (PoN)
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
    """
}

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

process INDEX_MOUSE {
        container "quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_6"
        publishDir params.outdir_index, mode: 'symlink'
        
        input:
        path fasta

        output:
        tuple val('mouse'), path("mouse.*")

        script:
        """
        bwa-mem2 index -p mouse ${fasta}
        """
}

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

process ALIGN_MOUSE {
	container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
	publishDir params.outdir_bam, mode: 'copy'
	
	input:
	tuple val(sample_id), path(trimmed_fastq_files)
	tuple val(genome), path(index_files)

	output:
	tuple val(sample_id), path("${sample_id}.mouse.bam"), emit: bam

	script:
	def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
	"""
	bwa-mem2 mem \
	-t 8 \
	-T 0 \
	-R "${read_group}" \
	mouse \
	${trimmed_fastq_files[0]} \
	${trimmed_fastq_files[1]} |
	samtools view \
	-Shb \
	-o ${sample_id}.mouse.bam -
	"""
}

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

process BAMCMP {
	container "community.wave.seqera.io/library/bamcmp:2.2--0b3144433037098f"
	publishDir params.outdir_bamcmp, mode: 'copy'

	input:
	tuple val(sample_id), path(bam_human), path(bam_mouse)

	output:
	tuple val(sample_id), path("${sample_id}_humanOnly.bam"), emit: human_only
	tuple val(sample_id), path("${sample_id}_humanBetter.bam"), emit: human_better
	tuple val(sample_id), path("${sample_id}_mouseOnly.bam"), emit: mouse_only
	tuple val(sample_id), path("${sample_id}_mouseBetter.bam"), emit: mouse_better
	tuple val(sample_id), path("${sample_id}_humanLoss.bam"), emit: human_loss
	tuple val(sample_id), path("${sample_id}_mouseLoss.bam"), emit: mouse_loss

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

process INDEX_BAM {
        container "community.wave.seqera.io/library/bwa-mem2_samtools:b7ce408fd27b2698"
        publishDir params.outdir_bamcmp, mode: 'copy'

	input:
	tuple val(sample_id), path(bam_file)

	output:
	tuple val(sample_id), path(bam_file), path("${bam_file}.bai"), emit: indexed_bam

	script:
	"""
	samtools index ${bam_file}
	"""
}

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

process MARK_DUPLICATES {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_markduplicates, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file)

	output:
	tuple val(sample_id), path("${sample_id}.marked_duplicates.bam"), emit: marked_dup_bam
	path "${sample_id}.marked_dup_metrics.txt", emit: marked_dup_metrics

	script:
	"""
	gatk MarkDuplicates \
		-I ${bam_file} \
		-O ${sample_id}.marked_duplicates.bam \
		-M ${sample_id}.marked_dup_metrics.txt
	"""
}

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

process BASE_RECALIBRATOR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file)
	path fasta
	path fasta_index
	path dict_file
	path dbsnp_vcf
	path dbsnp_vcf_idx
	path known_indels_vcf
	path known_indels_vcf_idx
	path mills_1000G_vcf
	path mills_1000G_vcf_idx

	output:
	tuple val(sample_id), path("${sample_id}.recal_data.table"), emit: recal_data_table

	script:
	"""
	gatk BaseRecalibrator \
		-I ${bam_file} \
		-R ${fasta} \
		--known-sites ${dbsnp_vcf} \
		--known-sites ${known_indels_vcf} \
		--known-sites ${mills_1000G_vcf} \
		-O ${sample_id}.recal_data.table
	"""
}

process APPLY_BQSR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_bqsr, mode: 'symlink'

	input:
	tuple val(sample_id), path(bam_file), path(recal_table)

	output:
	tuple val(sample_id), path("${sample_id}.recal.bam"), emit: recal_bam

	script:
	"""
	gatk ApplyBQSR \
		-I ${bam_file} \
		--bqsr-recal-file ${recal_table} \
		-O ${sample_id}.recal.bam
	"""
}

process MUTECT2 {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_mutect2, mode: "copy"

	input:
	path fasta
	path fasta_index
	path dict_file
	tuple val(sample_id), path(recal_bam)
	path germline_resource
	path germline_resource_idx
	path panel_of_normals
	path panel_of_normals_idx

	output:
	tuple val(sample_id), path("${sample_id}.unfiltered.vcf.gz"), path("${sample_id}.unfiltered.vcf.gz.stats"), emit: unfiltered_vcf

	script:
	"""
	gatk Mutect2 \
		-R ${fasta} \
		-I ${recal_bam} \
		--germline-resource ${germline_resource} \
		--panel-of-normals ${panel_of_normals} \
		-O ${sample_id}.unfiltered.vcf.gz
	"""
}

process GET_PILEUP_SUMMARIES {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	tuple val(sample_id), path(recal_bam)
	path filtered_vcf
	path filtered_vcf_idx
	
	output:
	tuple val(sample_id), path("${sample_id}.pileups.table"), emit: pileup_table

	script:
	"""
	gatk GetPileupSummaries \
		-I ${recal_bam} \
		-V ${filtered_vcf} \
		-L ${filtered_vcf} \
		-O ${sample_id}.pileups.table
	"""
}

process CALCULATE_CONTAMINATION {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"

	input:
	tuple val(sample_id), path(pileup_table)

	output:
	tuple val(sample_id), path("${sample_id}.contamination.table"), path("${sample_id}.segments.tsv"), emit: contamination_table

	script:
	"""
	gatk CalculateContamination \
		-I ${pileup_table} \
		-O ${sample_id}.contamination.table \
		--tumor-segmentation ${sample_id}.segments.tsv
	"""
}

process INDEX_UNFILTERED_VCF {
	container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
	publishDir params.outdir_mutect2, mode: "symlink"

	input:
	tuple val(sample_id), path(unfiltered_vcf), path(unfiltered_vcf_stats)

	output:
	tuple val(sample_id), path("${sample_id}.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf_idx

	script:
	"""
	bcftools index -t ${unfiltered_vcf}
	"""
}
	
process FILTER_MUTECT_CALLS {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"
	
	input:
	path fasta
	path fasta_idx
	path dict_file
	tuple val(sample_id), path(unfiltered_vcf), path(unfiltered_vcf_stats), path(contamination_table), path(tumor_segmentation_table), path(unfiltered_vcf_idx)

	output:
	tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), emit: filtered_vcf

	script:
	"""
	gatk FilterMutectCalls \
		-R ${fasta} \
		-V ${unfiltered_vcf} \
		--contamination-table ${contamination_table} \
		--tumor-segmentation ${tumor_segmentation_table} \
		-O ${sample_id}.filtered.vcf.gz
	""" 
}

process FUNCOTATOR_DATA_SOURCE_DOWNLOADER {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_resources, mode: "symlink"

	output:
	path "funcotator_dataSources", emit: funcotator_sources

    	script:
    	"""
    	gatk FuncotatorDataSourceDownloader \
        	--somatic \
        	--validate-integrity \
        	--extract-after-download \
		--hg38 \
        	--output funcotator_dataSources
    	"""
}

process INDEX_FILTERED_VCF {
	container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
        publishDir params.outdir_mutect2, mode: "symlink"

        input:
        tuple val(sample_id), path(filtered_vcf)

        output:
        tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"), emit: filtered_vcf_idx

        script:
        """
        bcftools index -t ${filtered_vcf}
        """
}		

process FUNCOTATOR {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_mutect2, mode: "copy"

	input:
	tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_idx)
	path reference_fasta
	path reference_idx
	path reference_dict
	path funcotator_data_sources

	output:
	tuple val(sample_id), path("${sample_id}.filtered.annotated.vcf"), emit: filtered_annotated_vcf

	script:
	"""
	mkdir funcotator_dataSources_extracted
	tar -xzvf ${funcotator_data_sources} -C funcotator_dataSources_extracted
	gatk Funcotator \
                --variant ${filtered_vcf} \
                --reference ${reference_fasta} \
                --ref-version hg38 \
                --data-sources-path funcotator_dataSources_extracted/funcotator_dataSources.v1.8.hg38.20230908s \
                --output ${sample_id}.filtered.annotated.vcf \
                --output-file-format VCF
	"""
}

process FUNCOTATOR_MAF {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
	publishDir params.outdir_mutect2, mode: "copy"

	input:
	tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_idx)
        path reference_fasta
        path reference_idx
        path reference_dict
        path funcotator_data_sources

	output:
	tuple val(sample_id), path("${sample_id}.filtered.annotated.maf.gz"), emit: filtered_annotated_maf
	
	script:
	"""
        mkdir funcotator_dataSources_extracted
        tar -xzvf ${funcotator_data_sources} -C funcotator_dataSources_extracted
	gatk Funcotator \
		--variant ${filtered_vcf} \
		--reference ${reference_fasta} \
		--ref-version hg38 \
		--data-sources-path funcotator_dataSources_extracted/funcotator_dataSources.v1.8.hg38.20230908s \
		--output ${sample_id}.filtered.annotated.maf.gz \
		--output-file-format MAF \
		--remove-filtered-variants
	"""
}

workflow {
	ACCESSORY_FILES_DOWNLOAD()
	FASTP(ch_fastq)
	MULTIQC(FASTP.out.json_report.collect())
	human_index = INDEX_HUMAN(params.hg38_fasta)
	mouse_index = INDEX_MOUSE(params.mm39_fasta)
	ALIGN_HUMAN(FASTP.out.trimmed_reads, human_index)
	ALIGN_MOUSE(FASTP.out.trimmed_reads, mouse_index)
	all_bams = ALIGN_HUMAN.out.bam.mix(ALIGN_MOUSE.out.bam)
	SORT_BAM(all_bams)
	// Separate human and mouse BAMs
    	sorted_bams = SORT_BAM.out.sorted_bam
        	.map { sample_id, bam ->
			def species = bam.name.contains("human") ? "human": "mouse"
            		return tuple(sample_id, species, bam)
        	}
        	.groupTuple(by: [0, 1])
        	.map { sample_id, species, bams -> 
            		return tuple(sample_id, species, bams[0])
        	}

        // Pair human and mouse BAMs
    	paired_bams = sorted_bams
        	.branch {
            		human: it[1] == "human"
            		mouse: it[1] == "mouse"
        	}
    
    	bamcmp_input = paired_bams.human
        	.join(paired_bams.mouse, by: 0)
        	.map { sample_id, human_species, human_bam, mouse_species, mouse_bam ->
            		return tuple(sample_id, human_bam, mouse_bam)
        	}

	BAMCMP(bamcmp_input)
	bams_to_merge = BAMCMP.out.human_only
		.join(BAMCMP.out.human_better)
		.map {sample_id, human_only, human_better ->
			tuple(sample_id, human_only, human_better)
		}
	MERGE_BAMS(bams_to_merge)
	INDEX_BAM(MERGE_BAMS.out.human_merged_sorted_bam)
	bams_for_counting=BAMCMP.out.human_only
		.mix(BAMCMP.out.human_better)
		.mix(BAMCMP.out.mouse_only)
		.mix(BAMCMP.out.mouse_better)
		.mix(BAMCMP.out.human_loss)
		.mix(BAMCMP.out.mouse_loss)
		.groupTuple()
	COUNT_READS(bams_for_counting)
	MARK_DUPLICATES(MERGE_BAMS.out.human_merged_sorted_bam)
	INDEX_REFERENCE(params.hg38_fasta)
	CREATE_DICT(params.hg38_fasta)
	BASE_RECALIBRATOR(
		MARK_DUPLICATES.out.marked_dup_bam,
		params.hg38_fasta,
		INDEX_REFERENCE.out.fasta_index,
		CREATE_DICT.out.dict_file,
		ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf,
		ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf_idx,
		ACCESSORY_FILES_DOWNLOAD.out.known_indels,
		ACCESSORY_FILES_DOWNLOAD.out.known_indels_idx,
		ACCESSORY_FILES_DOWNLOAD.out.mills_indels,
		ACCESSORY_FILES_DOWNLOAD.out.mills_indels_idx

	)
	bqsr_input = MARK_DUPLICATES.out.marked_dup_bam.join(BASE_RECALIBRATOR.out.recal_data_table)
	APPLY_BQSR(bqsr_input)
	MUTECT2(
		params.hg38_fasta,
		INDEX_REFERENCE.out.fasta_index,
		CREATE_DICT.out.dict_file,
		APPLY_BQSR.out.recal_bam,
		ACCESSORY_FILES_DOWNLOAD.out.gnomad,
		ACCESSORY_FILES_DOWNLOAD.out.gnomad_idx,
		ACCESSORY_FILES_DOWNLOAD.out.pon,
		ACCESSORY_FILES_DOWNLOAD.out.pon_idx,
	)
	GET_PILEUP_SUMMARIES(
		APPLY_BQSR.out.recal_bam,
		ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf,
		ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf_idx
	)
	CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup_table)
	INDEX_UNFILTERED_VCF(MUTECT2.out.unfiltered_vcf)
	mutect2_contamination = MUTECT2.out.unfiltered_vcf.join(CALCULATE_CONTAMINATION.out.contamination_table).join(INDEX_UNFILTERED_VCF.out.unfiltered_vcf_idx)
	FILTER_MUTECT_CALLS(
		params.hg38_fasta,
		INDEX_REFERENCE.out.fasta_index,
		CREATE_DICT.out.dict_file,
		mutect2_contamination
	)
	FUNCOTATOR_DATA_SOURCE_DOWNLOADER()
	INDEX_FILTERED_VCF(FILTER_MUTECT_CALLS.out.filtered_vcf)
	FUNCOTATOR(
		FILTER_MUTECT_CALLS.out.filtered_vcf.join(INDEX_FILTERED_VCF.out.filtered_vcf_idx),
		params.hg38_fasta,
		INDEX_REFERENCE.out.fasta_index,
		CREATE_DICT.out.dict_file,
		FUNCOTATOR_DATA_SOURCE_DOWNLOADER.out.funcotator_sources
	)
	FUNCOTATOR_MAF(
                FILTER_MUTECT_CALLS.out.filtered_vcf.join(INDEX_FILTERED_VCF.out.filtered_vcf_idx),
                params.hg38_fasta,
                INDEX_REFERENCE.out.fasta_index,
                CREATE_DICT.out.dict_file,
                FUNCOTATOR_DATA_SOURCE_DOWNLOADER.out.funcotator_sources
        )
}
