FROM continuumio/miniconda3

# Create a new conda environment and install fastp
RUN conda install -c bioconda fastp -y \
	&& conda clean -a -y

# Verify installation
RUN fastp --version

