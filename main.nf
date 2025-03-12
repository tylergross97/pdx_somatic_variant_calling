#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

ch_fastq = Channel.fromFilePairs("${params.fastq}/*_S*_R{1,2}_001.fastq.gz", checkIfExists: true)

include { ACCESSORY_FILES_DOWNLOAD } from './modules/accessory_files_download.nf'
include { FASTP } from './modules/fastp.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { INDEX_HUMAN } from './modules/index_human.nf'
include { INDEX_MOUSE } from './modules/index_mouse.nf'
include { ALIGN_HUMAN } from './modules/align_human.nf'
include { ALIGN_MOUSE } from './modules/align_mouse.nf'
include { SORT_BAM } from './modules/sort_bam.nf'
include { BAMCMP } from './modules/bamcmp.nf'
include { MERGE_BAMS } from './modules/merge_bams.nf'
include { INDEX_BAM } from './modules/index_bam.nf'
include { COUNT_READS } from './modules/count_reads.nf'
include { MARK_DUPLICATES } from './modules/mark_duplicates.nf'
include { INDEX_REFERENCE } from './modules/index_reference.nf'
include { CREATE_DICT } from './modules/create_dict.nf'
include { BASE_RECALIBRATOR } from './modules/base_recalibrator.nf'
include { APPLY_BQSR } from './modules/apply_bqsr.nf'
include { MUTECT2 } from './modules/mutect2.nf'
include { GET_PILEUP_SUMMARIES } from './modules/get_pileup_summaries.nf'
include { CALCULATE_CONTAMINATION } from './modules/calculate_contamination.nf'
include { INDEX_UNFILTERED_VCF } from './modules/index_unfiltered_vcf.nf'
include { FILTER_MUTECT_CALLS } from './modules/filter_mutect_calls.nf'
include { FUNCOTATOR_DATA_SOURCE_DOWNLOADER } from './modules/funcotator_data_source_downloader.nf'
include { INDEX_FILTERED_VCF } from './modules/index_filtered_vcf.nf'
include { FUNCOTATOR } from './modules/funcotator.nf'
include { FUNCOTATOR_MAF } from './modules/funcotator_maf.nf'

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
