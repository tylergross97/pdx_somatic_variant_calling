process SIMULATE_PDX_READS {
    publishDir "tests/test_fastqs", mode: 'copy'
    
    output:
    path "sample1_S1_R1_001.fastq.gz", emit: r1
    path "sample1_S1_R2_001.fastq.gz", emit: r2
    
    script:
    '''
    python3 << 'EOF'
import gzip
import random

def generate_read(read_id, length=100):
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    quality = ''.join(['I'] * length)  # High quality scores (ASCII 73)
    
    # Proper FASTQ format: header, sequence, plus, quality
    fastq_entry = f"@read_{read_id}\\n{sequence}\\n+\\n{quality}\\n"
    return fastq_entry

print("Generating test FASTQ files...")

# Generate R1 file (in current working directory)
print("Creating R1 file...")
with gzip.open('sample1_S1_R1_001.fastq.gz', 'wt') as f:
    for i in range(100):  # Reduced to 100 reads for faster testing
        entry = generate_read(f"{i:04d}_R1")
        f.write(entry)

# Generate R2 file (in current working directory)
print("Creating R2 file...")
with gzip.open('sample1_S1_R2_001.fastq.gz', 'wt') as f:
    for i in range(100):  # Same number of reads
        entry = generate_read(f"{i:04d}_R2")
        f.write(entry)

print("FASTQ generation completed!")
EOF

    # Verify the generated files
    echo "=== Verification ==="
    echo "Files created:"
    ls -la *.fastq.gz
    
    echo "File sizes:"
    du -h *.fastq.gz
    
    echo "R1 sample (first 12 lines):"
    zcat sample1_S1_R1_001.fastq.gz | head -12
    
    echo "R1 line count: $(zcat sample1_S1_R1_001.fastq.gz | wc -l)"
    echo "R2 line count: $(zcat sample1_S1_R2_001.fastq.gz | wc -l)"
    
    echo "R1 read count: $(zcat sample1_S1_R1_001.fastq.gz | wc -l | awk '{print $1/4}')"
    echo "R2 read count: $(zcat sample1_S1_R2_001.fastq.gz | wc -l | awk '{print $1/4}')"
    '''
}
