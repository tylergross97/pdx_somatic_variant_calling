#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Common modules for both workflows
include { ACCESSORY_FILES_DOWNLOAD } from './modules/local/accessory_files_download.nf'
include { FASTP } from './modules/local/fastp.nf'
include { MULTIQC } from './modules/local/multiqc.nf'
include { INDEX_HUMAN } from './modules/local/index_human.nf'
include { INDEX_MOUSE } from './modules/local/index_mouse.nf'
include { ALIGN_HUMAN } from './modules/local/align_human.nf'
include { ALIGN_MOUSE } from './modules/local/align_mouse.nf'
include { SORT_BAM } from './modules/local/sort_bam.nf'
include { BAMCMP } from './modules/local/bamcmp.nf'
include { MERGE_BAMS } from './modules/local/merge_bams.nf'
include { INDEX_BAM } from './modules/local/index_bam.nf'
include { MARK_DUPLICATES } from './modules/local/mark_duplicates.nf'
include { INDEX_MARKED_DUPLICATES } from './modules/local/index_marked_duplicates.nf'
include { INDEX_REFERENCE } from './modules/local/index_reference.nf'
include { CREATE_DICT } from './modules/local/create_dict.nf'
include { GET_PILEUP_SUMMARIES } from './modules/local/get_pileup_summaries.nf'
include { CALCULATE_CONTAMINATION } from './modules/local/calculate_contamination.nf'
include { INDEX_UNFILTERED_VCF } from './modules/local/index_unfiltered_vcf.nf'
include { FILTER_MUTECT_CALLS } from './modules/local/filter_mutect_calls.nf'
include { HAMA_CSV_TO_BED } from './modules/local/hama_csv_to_bed.nf'
include { FILTER_HAMA_VARIANTS } from './modules/local/filter_hama_variants.nf'

// Modules for non-interval workflow
include { BASE_RECALIBRATOR } from './modules/local/base_recalibrator.nf'
include { APPLY_BQSR } from './modules/local/apply_bqsr.nf'
include { MUTECT2 } from './modules/local/mutect2.nf'

// Modules for interval workflow

include { BASE_RECALIBRATOR_INTERVALS } from './modules/local/base_recalibrator_intervals.nf'
include { APPLY_BQSR_INTERVALS } from './modules/local/apply_bqsr_intervals.nf'
include { MUTECT2_INTERVALS } from './modules/local/mutect2_intervals.nf'

workflow {
    // Check if we're running in test mode
    def is_test_profile = workflow.profile.contains('test')

    if (is_test_profile) {
        // Use local test reference files - chr22 specific filenames
        ch_dbsnp_vcf = Channel.value(file("${params.test_data_base}/Homo_sapiens_assembly38.dbsnp138.chr22.vcf", checkIfExists: true))
        ch_dbsnp_vcf_idx = Channel.value(file("${params.test_data_base}/Homo_sapiens_assembly38.dbsnp138.chr22.vcf.idx", checkIfExists: true))
        ch_known_indels = Channel.value(file("${params.test_data_base}/Homo_sapiens_assembly38.known_indels.chr22.vcf.gz", checkIfExists: true))
        ch_known_indels_idx = Channel.value(file("${params.test_data_base}/Homo_sapiens_assembly38.known_indels.chr22.vcf.gz.tbi", checkIfExists: true))
        ch_mills_indels = Channel.value(file("${params.test_data_base}/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz", checkIfExists: true))
        ch_mills_indels_idx = Channel.value(file("${params.test_data_base}/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz.tbi", checkIfExists: true))
        ch_gnomad = Channel.value(file("${params.test_data_base}/af-only-gnomad.hg38.chr22.vcf.gz", checkIfExists: true))
        ch_gnomad_idx = Channel.value(file("${params.test_data_base}/af-only-gnomad.hg38.chr22.vcf.gz.tbi", checkIfExists: true))
        ch_filtered_vcf = Channel.value(file("${params.test_data_base}/small_exac_common_3.hg38.chr22.vcf.gz", checkIfExists: true))
        ch_filtered_vcf_idx = Channel.value(file("${params.test_data_base}/small_exac_common_3.hg38.chr22.vcf.gz.tbi", checkIfExists: true))
        ch_pon = Channel.value(file("${params.test_data_base}/1000g_pon.hg38.chr22.vcf.gz", checkIfExists: true))
        ch_pon_idx = Channel.value(file("${params.test_data_base}/1000g_pon.hg38.chr22.vcf.gz.tbi", checkIfExists: true))
        
        // Use test FASTQ input - search recursively in subdirectories
        ch_fastq = Channel.fromFilePairs("${params.fastq}/**/*_S*_R{1,2}_001.fastq.gz", checkIfExists: true)
        
        log.info "Running in TEST mode with chr22 test data"
        log.info "FASTQ input directory: ${params.fastq}"
        log.info "Test accessory files: ${params.test_data_base}"
    } else {
        // Production mode - download accessory files
        ACCESSORY_FILES_DOWNLOAD()
        
        // Map all the process outputs to the same channel variable names
        ch_dbsnp_vcf = ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf
        ch_dbsnp_vcf_idx = ACCESSORY_FILES_DOWNLOAD.out.dbsnp_vcf_idx
        ch_known_indels = ACCESSORY_FILES_DOWNLOAD.out.known_indels
        ch_known_indels_idx = ACCESSORY_FILES_DOWNLOAD.out.known_indels_idx
        ch_mills_indels = ACCESSORY_FILES_DOWNLOAD.out.mills_indels
        ch_mills_indels_idx = ACCESSORY_FILES_DOWNLOAD.out.mills_indels_idx
        ch_gnomad = ACCESSORY_FILES_DOWNLOAD.out.gnomad
        ch_gnomad_idx = ACCESSORY_FILES_DOWNLOAD.out.gnomad_idx
        ch_filtered_vcf = ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf
        ch_filtered_vcf_idx = ACCESSORY_FILES_DOWNLOAD.out.filtered_vcf_idx
        ch_pon = ACCESSORY_FILES_DOWNLOAD.out.pon
        ch_pon_idx = ACCESSORY_FILES_DOWNLOAD.out.pon_idx
        
        // Use production FASTQ input
        ch_fastq = Channel.fromFilePairs("${params.fastq}/*_S*_R{1,2}_001.fastq.gz", checkIfExists: true)
        
        log.info "Running in PRODUCTION mode"
    }

    if (params.containsKey('intervals') && params.intervals) {
        // Convert relative path to absolute path
        def intervals_file = file(params.intervals, checkIfExists: true)
        ch_intervals = Channel.value(intervals_file)
        log.info "Using intervals file: ${intervals_file}"
    } else {
        ch_intervals = Channel.empty()
        log.info "No intervals specified - running whole genome analysis"
    }

    // common initial steps for both workflows
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
            def species = bam.name.contains(".human.") ? "human" : 
                        bam.name.contains(".mouse.") ? "mouse" : "unknown"
            return tuple(sample_id, species, bam)
        }
        .filter { sample_id, species, bam -> species != "unknown" }
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
    MARK_DUPLICATES(MERGE_BAMS.out.human_merged_sorted_bam)
    INDEX_MARKED_DUPLICATES(MARK_DUPLICATES.out.marked_dup_bam)
    INDEX_REFERENCE(params.hg38_fasta)
    CREATE_DICT(params.hg38_fasta)
    
    // Generate HAMA BED file
    HAMA_CSV_TO_BED(
        file("${projectDir}/data/high_risk_HAMA_list.csv")
    )

    // Conditional workflow based on whether intervals are provided
    if (params.containsKey('intervals') && params.intervals) {
        // Interval-based workflow
  
        BASE_RECALIBRATOR_INTERVALS(
            INDEX_MARKED_DUPLICATES.out.indexed_bam,
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            ch_dbsnp_vcf,
            ch_dbsnp_vcf_idx,
            ch_known_indels,
            ch_known_indels_idx,
            ch_mills_indels,
            ch_mills_indels_idx,
            ch_intervals
        )
        bqsr_input = INDEX_MARKED_DUPLICATES.out.indexed_bam
            .join(BASE_RECALIBRATOR_INTERVALS.out.recal_data_table)
            .map { sample_id, bam, bai, recal_table ->
                tuple(sample_id, bam, bai, recal_table)
            }
        APPLY_BQSR_INTERVALS(bqsr_input, ch_intervals)
        MUTECT2_INTERVALS(
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            APPLY_BQSR_INTERVALS.out.recal_bam,
            ch_gnomad,
            ch_gnomad_idx,
            ch_pon,
            ch_pon_idx,
            ch_intervals
        )
        GET_PILEUP_SUMMARIES(
            APPLY_BQSR_INTERVALS.out.recal_bam,
            ch_filtered_vcf,
            ch_filtered_vcf_idx
        )
        CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup_table)
        INDEX_UNFILTERED_VCF(MUTECT2_INTERVALS.out.unfiltered_vcf)
        mutect2_contamination = MUTECT2_INTERVALS.out.unfiltered_vcf
            .join(CALCULATE_CONTAMINATION.out.contamination_table)
            .join(INDEX_UNFILTERED_VCF.out.unfiltered_vcf_idx)
    } else {
        // Non-interval workflow
        BASE_RECALIBRATOR(
            INDEX_MARKED_DUPLICATES.out.indexed_bam,
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            ch_dbsnp_vcf,
            ch_dbsnp_vcf_idx,
            ch_known_indels,
            ch_known_indels_idx,
            ch_mills_indels,
            ch_mills_indels_idx
        )
        bqsr_input = INDEX_MARKED_DUPLICATES.out.indexed_bam.join(BASE_RECALIBRATOR.out.recal_data_table)
        APPLY_BQSR(bqsr_input)
        MUTECT2(
            params.hg38_fasta,
            INDEX_REFERENCE.out.fasta_index,
            CREATE_DICT.out.dict_file,
            APPLY_BQSR.out.recal_bam,
            ch_gnomad,
            ch_gnomad_idx,
            ch_pon,
            ch_pon_idx
        )
        GET_PILEUP_SUMMARIES(
            APPLY_BQSR.out.recal_bam,
            ch_filtered_vcf,
            ch_filtered_vcf_idx
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

    // Filter out HAMA variants
    FILTER_HAMA_VARIANTS(
    	FILTER_MUTECT_CALLS.out.filtered_vcf,
    	HAMA_CSV_TO_BED.out.hama_bed
    )
    
}
