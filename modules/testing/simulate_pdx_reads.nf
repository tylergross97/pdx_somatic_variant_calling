process SIMULATE_PDX_READS {
    container "community.wave.seqera.io/library/art:2016.06.05--abd910249b0b53f1"
    publishDir params.test_results, mode: 'copy'

    output:
    path "tests/test_fastqs/sample1_S1_R1_001.fastq.gz", emit: r1
    path "tests/test_fastqs/sample1_S1_R2_001.fastq.gz", emit: r2
    path "tests/refs/hg38_chr22.fa", emit: hg38
    path "tests/refs/mm10_chr22.fa", emit: mm10

    script:
    """
    mkdir -p tests/test_fastqs tests/refs synthetic_fastqs

    # Write minimal mock references
    echo '>chr22\\nAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCGAGGCTTGC' > tests/refs/hg38_chr22.fa
    echo '>chr22\\nTCGATCGATCGGATGCTAGCTAGCTAGCTGATCGTAGCTAGCT' > tests/refs/mm10_chr22.fa

    # Simulate reads from each reference
    art_illumina -ss HS25 -sam -i tests/refs/hg38_chr22.fa -l 100 -f 2 -o synthetic_fastqs/human_
    art_illumina -ss HS25 -sam -i tests/refs/mm10_chr22.fa -l 100 -f 2 -o synthetic_fastqs/mouse_

    # Shuffle and mix human and mouse reads
    cat synthetic_fastqs/human_1.fq synthetic_fastqs/mouse_1.fq | shuf | gzip > tests/test_fastqs/sample1_S1_R1_001.fastq.gz
    cat synthetic_fastqs/human_2.fq synthetic_fastqs/mouse_2.fq | shuf | gzip > tests/test_fastqs/sample1_S1_R2_001.fastq.gz
    """
}
