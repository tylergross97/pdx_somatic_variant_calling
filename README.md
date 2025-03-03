# PDX Somatic Variant Calling Nextflow Pipeline

A Nextflow pipeline specifically designed to perform tumor-only SNP and Indel variant calling on Whole Exome Sequencing (WES) data from Patient-Derived Xenograft (PDX) models with built-in functionality for HPC users with Singularity containers.

## Table of Contents
1. [Introduction](#introduction)
2. [Pipeline Workflow](#pipeline-workflow)
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

## Installation

### Clone Repository

To get started with this pipeline, clone the repository to your local machine or HPC environment:

```bash
git clone https://github.com/tylergross97/pdx_somatic_variant_calling.git
cd pdx_somatic_variant_calling
```
### Install software and accessory files
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

4. GATK Resource Bundle (for hg38)
   - Download from: [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
      - This contains the accessory files needed for the variant calling portion of the pipeline (e.g., database of common germline variants)

5. Input Data:
   - Paired-end FASTQ files from your PDX samples
   - FASTQ File Naming Convention:
      This pipeline requires a specific naming convention for input FASTQ files. Files should follow this pattern:
      
      *_S*_R{1,2}_001.fastq.gz
      
      Where:
      
      * can be any string (usually sample name or identifier)
      S* represents the sample number (e.g., S1, S2, S3, etc.)
      R{1,2} specifies whether it's the forward (R1) or reverse (R2) read file
      001 is a common suffix in Illumina sequencing output
      Files must be gzipped (.gz extension)

      Examples of correctly named files:
      
      Sample1_S1_R1_001.fastq.gz and Sample1_S1_R2_001.fastq.gz
     
      PDX-tumor_S2_R1_001.fastq.gz and PDX-tumor_S2_R2_001.fastq.gz

      If your files don't match this naming convention, you may need to rename them before running the pipeline.

For optional downstream analysis:

6. R (version 4.0 or later) for downstream analysis with maftools
   - Installation instructions: [R Installation Guide](https://cran.r-project.org/) and [maftools](https://www.bioconductor.org/packages/release/bioc/html/maftools.html)

7. Python (version 3.6 or later) for downstream analysis of contamination
   - Installation instructions: [Python Installation Guide](https://www.python.org/downloads/)

#### Reference Genomes

Instructions for obtaining and preparing the reference genomes
- As explained in [Zverinova et al, 2021](https://onlinelibrary.wiley.com/doi/10.1002/humu.24311), we recommend using primary genome assemblies for references

1. Human (hg38):
   - Download from: [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

2. Mouse (mm39):
   - Download from: [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/)

### Preparing your [nextflow.config](https://www.nextflow.io/docs/latest/config.html) file

In your cloned repository directory, you have a nextflow.config.template file. All you need to do is copy this file as 'nextflow.config' and edit it to reflect the paths of your accessory files you just downloaded and your fastq files

```bash
cp nextflow.config.template nextflow.config
```

Note that the nextflow.config.template file is set up for running Singularity. Adjust as needed.

## Usage

### Running locally

With your nextflow.config and main.nf files in your current working directory and nextflow installed, all you need to do is run the following command:

```bash
nextflow run main.nf
```

### Running on SLURM

If you're using a high-performance computing (HPC) cluster that uses SLURM for job scheduling, you can create a [SLURM script](https://www.arch.jhu.edu/short-tutorial-how-to-create-a-slurm-script/) to run the pipeline. It may look something like this:

```bash
#!/bin/bash
#SBATCH --job-name=pdx_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=pdx_pipeline_%j.out
#SBATCH --error=pdx_pipeline_%j.err

# Load Nextflow module (adjust or remove if Nextflow is in your PATH)
module load nextflow

# Set environment variables
export NXF_WORK=$SCRATCH/pdx_work
export SINGULARITY_CACHEDIR=$HOME/singularity_cache
export NXF_SINGULARITY_CACHEDIR=$HOME/nextflow_singularity_cache

# Run the Nextflow pipeline
nextflow run main.nf
```

## Pipeline Steps

1. [fastp](https://github.com/OpenGene/fastp) for quality control and adapter trimming
2. [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) to align trimmed fastq files to both human and mouse reference genomes
3. [bamcmp](https://github.com/CRUKMI-ComputationalBiology/bamcmp) to perform deconvolution of mouse reads
4. [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) to identify duplicate reads
5. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) and [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR) for base quality score recalibration
6. [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) to call somatic short variants
7. [GATK GetPileupSummaries](https://gatk.broadinstitute.org/hc/en-us/articles/360037593451-GetPileupSummaries), [GATK CalculateContamination](https://gatk.broadinstitute.org/hc/en-us/articles/360036888972-CalculateContamination), and [GATK FilterMutectCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls) to filter variant calls
8. [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator) to annotate variants and generate .vcf and .maf files

## Output

There are many intermediate files generated that will be placed in the results directory you specify in your nextflow.config file. The main files were are interested in are the annotated .vcf and .maf files that can loaded into an R markdown file for analysis with maftools, see [below](#downstream-analysis). These filtered.annotated.vcf.gz and filtered.annotated.maf.gz files were be saved to the ./results/mutect2/directory

If you are looking to analyses the level of contamination of your original samples, you will need to access the files outputted from bamcmp in the ./results/bamcmp directory

## Downstream analyses

### Contamination Analysis (Python)
Provided is a link to a .pdf file with the necessary code and expected outputs to visualize contamination present in your original sample based on  the output of bamcmp. Set the directory variable to ".results/bamcmp/" 
[![Contamination.ipynb](images/contamination_analysis_preview.png)](images/Contamination_bamcmp.pdf)

### maftools analysis (R)
Provided is an html file with the necessary code to leverage maftools to analyze the filtered.annotated.maf.gz files that are the main output from the pipeline. In the provided example, I was looking to identify mutations in genes known to be implicated in Renal Cell Carcinoma (RCC).
[![maftools.Rmd](images/maftools_analysis_preview.png)](images/maftools.pdf)


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
7. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
8. Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997.

