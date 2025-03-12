process ACCESSORY_FILES_DOWNLOAD {
    conda "conda-forge::curl=7.86.0"
    publishDir params.outdir_accessory_files, mode: 'copy'

    output:
    // dbSNP files
    path "Homo_sapiens_assembly38.dbsnp138.vcf", emit: dbsnp_vcf
    path "Homo_sapiens_assembly38.dbsnp138.vcf.idx", emit: dbsnp_vcf_idx
    // Known indels
    path "Homo_sapiens_assembly38.known_indels.vcf.gz", emit: known_indels
    path "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi", emit: known_indels_idx
    // Mills indels
    path "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", emit: mills_indels
    path "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", emit: mills_indels_idx
    // gnomAD
    path "af-only-gnomad.hg38.vcf.gz", emit: gnomad
    path "af-only-gnomad.hg38.vcf.gz.tbi", emit: gnomad_idx
    // Filtered VCF (ExAC common SNPs)
    path "small_exac_common_3.hg38.vcf.gz", emit: filtered_vcf
    path "small_exac_common_3.hg38.vcf.gz.tbi", emit: filtered_vcf_idx
    // Panel of Normals (PoN)
    path "1000g_pon.hg38.vcf.gz", emit: pon
    path "1000g_pon.hg38.vcf.gz.tbi", emit: pon_idx

    script:
    """
    # dbSNP files
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

    # Known indels
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

    # Mills indels
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    curl -O https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

    # gnomAD
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

    # Filtered VCF (ExAC common SNPs)
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi

    # Panel of Normals (PoN)
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
    curl -O https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
    """
}
