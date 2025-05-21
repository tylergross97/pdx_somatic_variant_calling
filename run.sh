#!/bin/bash
#SBATCH --job-name="pdx"
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --account=rpili
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out        # Standard output file (%j will be replace>
#SBATCH --error=slurm-%j.err         # Standard error file
#SBATCH --mail-user=tgross2@buffalo.edu
#SBATCH --mail-type=ALL

module load nextflow/23.10.0

export NXF_HOME=/projects/academic/rpili/tgross2/tmp/.nextflow
export NXF_WORK=/vscratch/grp-rpili/pdx
export SINGULARITY_LOCALCACHEDIR=/projects/academic/rpili/tgross2/tmp
export SINGULARITY_CACHEDIR=/projects/academic/rpili/tgross2/tmp
export SINGULARITY_TMPDIR=/projects/academic/rpili/tgross2/tmp
export NXF_SINGULARITY_CACHEDIR=/projects/academic/rpili/tgross2/singularity_cache

nextflow run main.nf -resume
