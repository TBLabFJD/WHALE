process DMR_FILTERING {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)
    path chrom_sizes
    path promoters_bed
    path enhancers_bed

    output:
    tuple val(meta), path("${prefix}_dmr.bed")                  , emit: dmr_bed
    tuple val(meta), path("${prefix}_overlapped_promoters.bed") , emit: promoters_bed
    tuple val(meta), path("${prefix}_overlapped_enhancers.bed") , emit: enhancers_bed
    path "versions.yml"                                         , emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk \\
        -v OFS='\t' \\
        '\$4 == "different" ${args}' \\
        ${bed} \\
        | sort -k1,1 -k2,2n \\
        > ${prefix}_dmr.bed

    bedtools slop \\
        -i ${prefix}_dmr.bed \\
        -g ${chrom_sizes} \\
        -b 50 \\
        ${args2} \\
        > ${prefix}_slopped.bed

    bedtools intersect \\
        -a ${prefix}_slopped.bed \\
        -b ${promoters_bed} \\
        -wa \\
        -wb \\
        ${args3} \\
        > ${prefix}_overlapped_promoters.bed

    bedtools intersect \\
        -a ${prefix}_slopped.bed \\
        -b ${enhancers_bed} \\
        -wa \\
        -wb \\
        ${args4} \\
        > ${prefix}_overlapped_enhancers.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
        sort: \$(sort --version | head -n 1 | awk '{print \$4}')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}