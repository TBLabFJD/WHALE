process NANOSTAT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/pip_nanostat:69dc9b3e7f904176"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"),                                     emit: report
    tuple val("${task.process}"), val('NanoStat'), eval("NanoStat -v"), emit: versions_nanostat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    NanoStat \\
        --bam ${bam} \\
        ${args} \\
        --name ${prefix}_bam_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NanoStat: \$(NanoStat -v 2>&1 | sed 's/NanoStat //g')
    END_VERSIONS
    """

}
