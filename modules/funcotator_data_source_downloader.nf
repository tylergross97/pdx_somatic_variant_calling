process FUNCOTATOR_DATA_SOURCE_DOWNLOADER {
	container "community.wave.seqera.io/library/gatk4:4.6.1.0--e3124bcb2431f4a9"
        publishDir params.outdir_resources, mode: "symlink"

	output:
	path "funcotator_dataSources", emit: funcotator_sources

    	script:
    	"""
    	gatk FuncotatorDataSourceDownloader \
        	--somatic \
        	--validate-integrity \
        	--extract-after-download \
		--hg38 \
        	--output funcotator_dataSources
    	"""
}
