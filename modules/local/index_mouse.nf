process INDEX_MOUSE {
        container "quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_6"
        publishDir params.outdir_index, mode: 'symlink'
        
        input:
        path fasta

        output:
        tuple val('mouse'), path("mouse.*")

        script:
        """
        bwa-mem2 index -p mouse ${fasta}
        """
}
