process WHATSHAP_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d837709891c2d98fc0956f6fd0dba18b0f67d96c4db74ccbae7db98fd00afe42/data'
        : 'community.wave.seqera.io/library/whatshap:2.8--7fe530bc624a3e5a' }"

    input:
    tuple val(meta), path(bam), path(bai), path(tsv)

    output:
    tuple val(meta), path("*_h1.bam"),        emit: bam_h1
    tuple val(meta), path("*_h2.bam"),        emit: bam_h2
    tuple val(meta), path("*unassigned.bam"), emit: bam_unassigned
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"), emit: versions_whatshap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    whatshap \\
        split \\
        --output-h1 ${prefix}_h1.bam \\
        --output-h2 ${prefix}_h2.bam \\
        --output-untagged ${prefix}_unassigned.bam \\
        ${args} \\
        ${bam} \\
        ${tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1 | sed 's/whatshap //g')
    END_VERSIONS
    """

}
