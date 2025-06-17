#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import your processes
include { SIMULATE_PDX_READS } from './modules/testing/simulate_pdx_reads.nf'
include { TEST_FASTQ_FORMAT } from './modules/testing/test_fastq_format.nf'

workflow {    
    // Run simulation of test data
    SIMULATE_PDX_READS()

    // Run the FASTQ format test
    TEST_FASTQ_FORMAT(
        SIMULATE_PDX_READS.out.r1, 
        SIMULATE_PDX_READS.out.r2,
    )
}
