process DMR_FILTERING {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(chrom_sizes)

    output:
    tuple val(meta), path("${prefix}_dmr.bed"), emit: dmr_bed
    path "versions.yml"                       , emit: versions

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk \\
        -v OFS='\t' \\
        '\$4 == "different" ${args} {print \$1, \$2, \$3, \$9, \$10, (\$13 > 0 ? \$13 : -\$13) ${args2}}' \\
        ${bed} \\
        > different_dmr.bed

    bedtools slop \\
        -i different_dmr.bed \\
        -g ${chrom_sizes} \\
        ${args3} \\
        > ${prefix}_dmr.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}