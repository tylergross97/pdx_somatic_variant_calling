# PDX Somatic Variant Calling Nextflow Pipeline

A Nextflow pipeline specifically designed to perform tumor-only SNP and Indel variant calling on Whole Exome Sequencing (WES) data from Patient-Derived Xenograft (PDX) models with built-in functionality for HPC users with Singularity containers.

## Table of Contents
1. [Introduction](#introduction)
2. [Pipeline Workflow](#pipeline-workflow)
3. [Prerequisites](#prerequisites)
4. [Installation](#clone-this-repository)
5. [Usage](#usage)
6. [Pipeline Steps](#pipeline-steps)
7. [Input](#input)
8. [Output](#output)
9. [Configuration](#configuration)
10. [Troubleshooting](#troubleshooting)
11. [Contributing](#contributing)
12. [License](#license)
13. [Citations](#citations)

## Introduction

**Note**: I have provided hyperlinks to helpful learning materials for concepts introduced throughout this README

This pipeline is built using [Nextflow](https://www.nextflow.io/), a workflow management software that utilizes [containerization](https://www.ibm.com/think/topics/containerization#:~:text=Containerization%20is%20the%20packaging%20of,runs%20consistently%20on%20any%20infrastructure.) to allow for portable and reproducible bioinformatics pipelines.

This pipeline is designed to perform [somatic short variant calling](https://www.garvan.org.au/news-resources/science-explained/types-of-variants) (SNPs Indels) from [patient-derived xenograft (PDX) models](https://en.wikipedia.org/wiki/Patient_derived_xenograft). Specifically, the pipeline, was built to handle [whole-exome sequencing (WES)](https://www.illumina.com/techniques/sequencing/dna-sequencing/targeted-resequencing/exome-sequencing.html) data without a matched-normal sample, which is referred to tumor-only variant calling.

Although somatic short variant calling of PDX models without a matched-normal sample is a common task in bioinformatics, particularly in translational oncology research, it introduces a unique set of challenges that many somatic variant calling pipelines do not address.

This pipeline can be conceptually broken down into two main steps:
- Deconvolution (filtering) of mouse reads
- Tumor-only somatic short variant calling of human reads

### Deconvolution of mouse reads
First, it is important to understand that although the tumor is implanted into the mouse, it originated from a human patient, meaning that we are interested in understanding the variants of the human tumor cells. However, during and after implanation of the tumor into the mouse, there is some degree of infiltration of mouse cells into the tumor. As discussed and explored throughly in [Jo et al., 2019](https://link.springer.com/article/10.1186/s13059-019-1849-2), this can lead to false-positive variant calls that should be minimized through the explicit filtering of the reads originating the mouse. This pipeline utilizes the [bamcmp](https://github.com/CRUKMI-ComputationalBiology/bamcmp) tool, although there are others available.

### Tumor-only somatic short variant calling
Performing somatic variant calling without a matched-normal sample also introduces challenges that must be addressed through the use of a database of common germline variants to be filtered out. In the case of a matched-normal sample, germline variants are defined as those present in both the matched-normal sample and the tumor sample. Although the use of a common germline variant database is more prone to rare germline variants showing up as false-positive somatic variant calls, this is often the reality for many researchers working with PDX models. [Mutect2](https://www.biorxiv.org/content/10.1101/861054v1.abstract) is a somatic short variant caller that has a tumor-only mode available and is used in this pipeline, following [GATK's best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels)[2]. The output of the pipeline are called variants in VCF and MAF formats.

### Optional downstream analyses
There is also documentation for downstream analyses of the outputs of the nextflow pipeline (see [Downstream Analyses](#downstream-analyses) section)

## Pipeline Workflow

![Workflow Diagram](images/workflow_diagram.png)

## Prerequisites

Before running this pipeline, ensure you have the following tools and resources installed. This is the hardest part, but I have provided documentation!

1. Nextflow (version 23.10.0 or later)
   - Installation instructions: [Nextflow Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html)
   - Note that if you are using an HPC system, you may be able to load Nextflow using the [module system](https://hpc-wiki.info/hpc/Modules)

2. Singularity (preferred for HPCs) or Docker
   - Singularity: [Singularity Installation Guide](https://sylabs.io/guides/3.0/user-guide/installation.html)
      - Note that this may already be installed on your HPC system
   - Docker: [Docker Installation Guide](https://docs.docker.com/get-docker/)

3. Reference Genomes:
   - Human (hg38)
   - Mouse (mm39)
   - Instructions for obtaining these genomes can be found in the [Reference Genomes](#reference-genomes) section below.
      - Note that as explained in [Zverinova et al., 2021](https://onlinelibrary.wiley.com/doi/10.1002/humu.24311), we recommend using the primary assembly version (see paper for details)

4. GATK Resource Bundle (for hg38)
   - Download from: [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
      - This contains the accessory files needed for the variant calling portion of the pipeline (e.g., database of common germline variants)

5. Input Data:
   - Paired-end FASTQ files from your PDX samples

For optional downstream analysis:

6. R (version 4.0 or later) for downstream analysis with maftools
   - Installation instructions: [R Installation Guide](https://cran.r-project.org/) and [maftools](https://www.bioconductor.org/packages/release/bioc/html/maftools.html)

7. Python (version 3.6 or later) for downstream analysis of contamination
   - Installation instructions: [Python Installation Guide](https://www.python.org/downloads/)

### Reference Genomes

Instructions for obtaining and preparing the reference genomes:

1. Human (hg38):
   - Download from: [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
   - Prepare the genome for BWA-MEM2 indexing:
     ```
     bwa-mem2 index hg38.fa
     ```

2. Mouse (mm39):
   - Download from: [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/)
   - Prepare the genome for BWA-MEM2 indexing:
     ```
     bwa-mem2 index mm39.fa
     ```
Please ensure all prerequisites are properly installed and configured before running the pipeline. For detailed installation instructions of each tool, refer to their respective documentation.

Now all you need to do is define the paths of these downloaded files in your nextflow.config file that you will be creating based on the nextflow.config.template file! Instructions [below](#configuration)

All of the remaining software will be automatically installed via containers specified at the top of each process in the main.nf file!

## Clone this repository

This pipeline is designed to be easy to set up and run. 

## Usage

## Pipeline Steps

## Input

## Output

## Downstream analyses

## Configuration

## Troubleshooting

##  Contributing

## License

## Citations

If you use this pipeline in your work, please cite:
[Tyler Gross] (2025). Tumor-Only and Whole-Exome Sequencing PDX Somatic Variant Calling Nextflow Pipeline [Computer software]. https://github.com/tylergross97/pdx_somatic_variant_calling

This pipeline is based on the following conceptual frameworks and best practices:
1. Jo, S. Y., Kim, E., & Kim, S. (2019). Impact of mouse contamination in genomic profiling of patient-derived models and best practice for robust analysis. Genome Biology, 20, 1-13.
2. GATK Best Practices for somatic short variant discovery (SNVs + Indels)
   Broad Institute. (2023). Somatic short variant discovery (SNVs + Indels). Retrieved Jan. 2025, from https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels

This pipeline uses several tools that should be cited independently:

3. Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature biotechnology, 35(4), 316-319.
4. Garima Khandelwal, Maria Girotti, Christopher Smowton, Sam Taylor, Chris Wirth, Marek Dynowski, Kris Frese, Ged Brady, Deborah Burt, Richard Marais, Crispin Miller. Next-Gen Sequencing Analysis and Algorithms for PDX and CDX Models. Molecular Cancer Research. 2017, 15:8, PMID: 28442585 DOI: 10.1158/1541-7786.MCR-16-0431
5. Benjamin, D., Sato, T., Cibulskis, K., Getz, G., Stewart, C., & Lichtenstein, L. (2019). Calling somatic SNVs and indels with Mutect2. BioRxiv, 861054.
6. Mayakonda A, Lin D, Assenov Y, Plass C, Koeffler PH (2018). “Maftools: efficient and comprehensive analysis of somatic variants in cancer.” Genome Research. doi:10.1101/gr.239244.118.

