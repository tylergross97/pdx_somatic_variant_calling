process HAMA_CSV_TO_BED {
	container "quay.io/biocontainers/python:3.9--1"
	publishDir params.outdir_hama, mode: 'copy'

	input:
	path csv_file

	output:
	path "high_risk_HAMA_filtered.bed", emit: hama_bed

	script:
	"""
	#!/usr/bin/env python3
	import csv
	
	with open("${csv_file}", 'r') as csv_file, open("high_risk_HAMA_filtered.bed", 'w') as bed_file:
		# Skip header
		csv_reader = csv.reader(csv_file)
		header = next(csv_reader)
		
		# Find the index of the Hc column
		hc_index = header.index("Hc")
		chrom_index = header.index("X.chromosome")
		pos_index = header.index("Position")
		ref_index = header.index("ReferenceAllele")
		alt_index = header.index("AlterateAllele")
		
		# Process each row
		count = 0
		for row in csv_reader:
			try:
				# Check if Hc > 1
				hc_value = float(row[hc_index])
				if hc_value > 1:
					chrom = row[chrom_index]
					pos = int(row[pos_index])
					ref = row[ref_index]
					alt = row[alt_index]
					
					# BED format is 0-based, half-open
					start = pos - 1
					end = pos + len(ref) - 1
					
					# Write to BED file (chrom, start, end, name)
					name = f"{chrom}_{pos}_{ref}_{alt}"
					bed_file.write(f"{chrom}\\t{start}\\t{end}\\t{name}\\n")
					count += 1
			except (ValueError, IndexError):
				# Skip rows with invalid data
				continue
		
		print(f"Total rows meeting criteria (Hc > 1): {count}")
	"""
}
