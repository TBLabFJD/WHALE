process NANOSTAT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/pip_nanostat:69dc9b3e7f904176"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}_nanostat_report.txt"), emit: report
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    NanoStat \\
        --bam ${bam} \\
        ${args} \\
        --name ${prefix}_nanostat_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(NanoStat -v 2>&1 | sed 's/NanoStat //g')
    END_VERSIONS
    """
}