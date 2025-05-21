#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

ch_fastq = Channel.fromFilePairs("${params.fastq}/*_S*_R{1,2}_001.fastq.gz", checkIfExists: true)

// Common modules for both workflows
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
include { GET_PILEUP_SUMMARIES } from './modules/get_pileup_summaries.nf'
include { CALCULATE_CONTAMINATION } from './modules/calculate_contamination.nf'
include { INDEX_UNFILTERED_VCF } from './modules/index_unfiltered_vcf.nf'
include { FILTER_MUTECT_CALLS } from './modules/filter_mutect_calls.nf'
include { INDEX_FILTERED_VCF } from './modules/index_filtered_vcf.nf'
include { HAMA_CSV_TO_BED } from './modules/hama_csv_to_bed.nf'

// Modules for non-interval workflow
include { BASE_RECALIBRATOR } from './modules/base_recalibrator.nf'
include { APPLY_BQSR } from './modules/apply_bqsr.nf'
include { MUTECT2 } from './modules/mutect2.nf'

// Modules for interval workflow
include { INDEX_MARKED_DUPLICATES_INTERVALS } from './modules/index_marked_duplicates_intervals.nf'
include { BASE_RECALIBRATOR_INTERVALS } from './modules/base_recalibrator_intervals.nf'
include { APPLY_BQSR_INTERVALS } from './modules/apply_bqsr_intervals.nf'
include { MUTECT2_INTERVALS } from './modules/mutect2_intervals.nf'

workflow {
    // Common initial steps for both workflows
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
    bams_for_counting = BAMCMP.out.human_only
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
    
    // Generate HAMA BED file
    HAMA_CSV_TO_BED(
        file("${projectDir}/data/high_risk_HAMA_list.csv")
    )

    // Conditional workflow based on whether intervals are provided
    if (params.containsKey('intervals') && params.intervals) {
        // Interval-based workflow
        INDEX_MARKED_DUPLICATES_INTERVALS(MARK_DUPLICATES.out.marked_dup_bam)
        BASE_RECALIBRATOR_INTERVALS(
            INDEX_MARKED_DUPLICATES_INTERVALS.out.indexed_bam,
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf,
            ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf_idx,
            ACCESSORY_FILES_DOWNLOAD.out.known_indels,
            ACCESSORY_FILES_DOWNLOAD.out.known_indels_idx,
            ACCESSORY_FILES_DOWNLOAD.out.mills_indels,
            ACCESSORY_FILES_DOWNLOAD.out.mills_indels_idx,
            params.intervals
        )
        bqsr_input = INDEX_MARKED_DUPLICATES_INTERVALS.out.indexed_bam
            .join(BASE_RECALIBRATOR_INTERVALS.out.recal_data_table)
            .map { sample_id, bam, bai, recal_table ->
                tuple(sample_id, bam, bai, recal_table)
            }
        APPLY_BQSR_INTERVALS(bqsr_input, params.intervals)
        MUTECT2_INTERVALS(
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            APPLY_BQSR_INTERVALS.out.recal_bam,
            ACCESSORY_FILES_DOWNLOAD.out.gnomad,
            ACCESSORY_FILES_DOWNLOAD.out.gnomad_idx,
            ACCESSORY_FILES_DOWNLOAD.out.pon,
            ACCESSORY_FILES_DOWNLOAD.out.pon_idx,
            params.intervals
        )
        GET_PILEUP_SUMMARIES(
            APPLY_BQSR_INTERVALS.out.recal_bam,
            ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf,
            ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf_idx
        )
        CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup_table)
        INDEX_UNFILTERED_VCF(MUTECT2_INTERVALS.out.unfiltered_vcf)
        mutect2_contamination = MUTECT2_INTERVALS.out.unfiltered_vcf
            .join(CALCULATE_CONTAMINATION.out.contamination_table)
            .join(INDEX_UNFILTERED_VCF.out.unfiltered_vcf_idx)
    } else {
        // Non-interval workflow
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
            ACCESSORY_FILES_DOWNLOAD.out.pon_idx
        )
        GET_PILEUP_SUMMARIES(
            APPLY_BQSR.out.recal_bam,
            ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf,
            ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf_idx
        )
        CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup_table)
        INDEX_UNFILTERED_VCF(MUTECT2.out.unfiltered_vcf)
        mutect2_contamination = MUTECT2.out.unfiltered_vcf
            .join(CALCULATE_CONTAMINATION.out.contamination_table)
            .join(INDEX_UNFILTERED_VCF.out.unfiltered_vcf_idx)
    }

    // Common final steps for both workflows
    FILTER_MUTECT_CALLS(
        params.hg38_fasta,
        INDEX_REFERENCE.out.fasta_index,
        CREATE_DICT.out.dict_file,
        mutect2_contamination
    )
    INDEX_FILTERED_VCF(FILTER_MUTECT_CALLS.out.filtered_vcf)
}
