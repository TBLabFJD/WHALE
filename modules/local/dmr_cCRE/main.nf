process DMR_cCRE {
    tag "$meta.id"
    label 'process_medium'


    container "oras://community.wave.seqera.io/library/bedtools:2.27.1--d5d82fb70c52c058"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(promoters_bed)
    tuple val(meta3), path(enhancers_bed)

    output:
    tuple val(meta), path("${prefix}_overlapped_promoters.bed") , emit: promoters_bed
    tuple val(meta), path("${prefix}_overlapped_enhancers.bed") , emit: enhancers_bed
    path "versions.yml"                                         , emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools intersect \\
        -a ${bed} \\
        -b ${promoters_bed} \\
        -wa \\
        -wb \\
        ${args} \\
        > ${prefix}_overlapped_promoters.bed

    bedtools intersect \\
        -a ${bed} \\
        -b ${enhancers_bed} \\
        -wa \\
        -wb \\
        ${args2} \\
        > ${prefix}_overlapped_enhancers.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}