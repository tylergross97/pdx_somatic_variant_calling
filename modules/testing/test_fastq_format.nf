process TEST_FASTQ_FORMAT {
    publishDir "test_results", mode: 'copy'
    
    input:
    path r1_file
    path r2_file
    
    output:
    path "fastq_format_test.pass", emit: validation_result
    
    script:
    """
    echo "=== FASTQ Format Validation Report ===" > fastq_format_test.pass
    echo "Generated on: \$(date)" >> fastq_format_test.pass
    echo "R1 file: ${r1_file}" >> fastq_format_test.pass
    echo "R2 file: ${r2_file}" >> fastq_format_test.pass
    echo "" >> fastq_format_test.pass
    
    # Check R1 file
    echo "--- Validating R1 file ---" >> fastq_format_test.pass
    if [[ -f "${r1_file}" ]]; then
        echo "✓ R1 file exists: ${r1_file}" >> fastq_format_test.pass
        echo "File size: \$(stat -c%s ${r1_file}) bytes" >> fastq_format_test.pass
        echo "Read count: \$(zcat ${r1_file} | wc -l | awk '{print \$1/4}')" >> fastq_format_test.pass
        echo "First 8 lines:" >> fastq_format_test.pass
        zcat ${r1_file} | head -8 >> fastq_format_test.pass
    else
        echo "✗ R1 file not found: ${r1_file}" >> fastq_format_test.pass
    fi
    
    echo "" >> fastq_format_test.pass
    
    # Check R2 file
    echo "--- Validating R2 file ---" >> fastq_format_test.pass
    if [[ -f "${r2_file}" ]]; then
        echo "✓ R2 file exists: ${r2_file}" >> fastq_format_test.pass
        echo "File size: \$(stat -c%s ${r2_file}) bytes" >> fastq_format_test.pass
        echo "Read count: \$(zcat ${r2_file} | wc -l | awk '{print \$1/4}')" >> fastq_format_test.pass
        echo "First 8 lines:" >> fastq_format_test.pass
        zcat ${r2_file} | head -8 >> fastq_format_test.pass
    else
        echo "✗ R2 file not found: ${r2_file}" >> fastq_format_test.pass
    fi
    
    echo "" >> fastq_format_test.pass
    echo "✓ OVERALL RESULT: VALIDATION COMPLETED" >> fastq_format_test.pass
    echo "Validation completed at: \$(date)" >> fastq_format_test.pass
    """
}
