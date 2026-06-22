process ANNOTSV_TSV_FILTERING {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
      'nf-core/ubuntu:20.04' }"

    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${prefix}_filtered.tsv") , emit: filtered_tsv
    path "versions.yml"                             , emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: 'print $0'
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk -F'\\t' -v OFS='\\t' '
    NR==1 {
        for(i=1; i<=NF; i++) {
            ${args}
        }
    }
    {
        ${args2}
    }' ${tsv} > ${prefix}_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}