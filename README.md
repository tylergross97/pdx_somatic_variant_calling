# PDX Somatic Variant Calling Nextflow Pipeline

A Nextflow pipeline specifically designed to perform tumor-only SNP and Indel variant calling on Whole Exome Sequencing (WES) data from Patient-Derived Xenograft (PDX) models.

## Table of Contents
1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Pipeline Steps](#pipeline-steps)
6. [Input](#input)
7. [Output](#output)
8. [Configuration](#configuration)
9. [Troubleshooting](#troubleshooting)
10. [Contributing](#contributing)
11. [License](#license)
12. [Citations](#citations)

## Introduction

**Note**: I have provided hyperlinks to helpful learning materials for concepts introduced throughout this README

This pipeline is built using [Nextflow](https://www.nextflow.io/) [3], a workflow management software that utilizes [containerization](https://www.ibm.com/think/topics/containerization#:~:text=Containerization%20is%20the%20packaging%20of,runs%20consistently%20on%20any%20infrastructure.) to allow for portable and reproducible bioinformatics pipelines.

This pipeline is designed to perform [somatic short variant calling](https://www.garvan.org.au/news-resources/science-explained/types-of-variants) (SNPs Indels) from [patient-derived xenograft (PDX) models](https://en.wikipedia.org/wiki/Patient_derived_xenograft). Specifically, the pipeline, was built to handle [whole-exome sequencing (WES)](https://www.illumina.com/techniques/sequencing/dna-sequencing/targeted-resequencing/exome-sequencing.html) data without a matched-normal sample, which is referred to tumor-only variant calling.

Although somatic short variant calling of PDX models without a matched-normal sample is a common task in bioinformatics, particularly in translational oncology research, it introduces a unique set of challenges that many somatic variant calling pipelines do not address. First, it is important to understand that although the tumor is implanted into the mouse, it originated from a human patient, meaning that we are interested in understanding the genome of the human tumor cells. However, during and after implanation of the tumor into the mouse, there is, to a varying extent, some degree of infiltration of mouse cells into the tumor. As discussed and explored throughly in Jo et al., 2019 [2], this can lead to false-positive variant calls that should be minimized through the explicit filtering of the reads f

To achieve this goal, the pipeline begins by 

## Prerequisites

## Installation

## Usage

## Pipeline Steps

## Input

## Output

## Configuration

## Troubleshooting

##  Contributing

## License

## Citations

If you use this pipeline in your work, please cite:
[Tyler Gross] (2025). Tumor-Only and Whole-Exome Sequencing PDX Somatic Variant Calling Nextflow Pipeline [Computer software]. https://github.com/tylergross97/pdx_somatic_variant_calling

This pipeline is based on the following conceptual frameworks and best practices:
1. GATK Best Practices for somatic short variant discovery (SNVs + Indels)
   Broad Institute. (2023). Somatic short variant discovery (SNVs + Indels). Retrieved Jan. 2025, from https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels
2. Jo, S. Y., Kim, E., & Kim, S. (2019). Impact of mouse contamination in genomic profiling of patient-derived models and best practice for robust analysis. Genome Biology, 20, 1-13.

This pipeline uses several tools that should be cited independently:

3. Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature biotechnology, 35(4), 316-319.
