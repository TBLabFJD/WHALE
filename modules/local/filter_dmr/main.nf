process FILTER_DMR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: dmr_bed
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_asm_dmr"

    """
    awk '\$4 == "different" ${args}' $bed > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}