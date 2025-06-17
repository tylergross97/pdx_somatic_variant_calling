#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import your SIMULATE_PDX_READS process from the modules folder
include { SIMULATE_PDX_READS } from './modules/testing/simulate_pdx_reads.nf'

workflow {
    // Run simulation of test data
    SIMULATE_PDX_READS()

    // Access outputs
    SIMULATE_PDX_READS.out.r1.view { println "Generated R1 test FASTQ: $it" }
    SIMULATE_PDX_READS.out.r2.view { println "Generated R2 test FASTQ: $it" }
    SIMULATE_PDX_READS.out.hg38.view { println "Generated human reference: $it" }
    SIMULATE_PDX_READS.out.mm10.view { println "Generated mouse reference: $it" }
}
