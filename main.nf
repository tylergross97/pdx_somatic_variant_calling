#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastp {
    input:
    path inputFile from params.inputFiles  // Get input files

    output:
    path "${inputFile.baseName}_trimmed.fastq.gz"  // Define the output filename based on input

    script:
    """
    fastp -i $inputFile -o ${inputFile.baseName}_trimmed.fastq.gz -h ${inputFile.baseName}_fastp.html
    """
}

workflow {
    fastp()
}
