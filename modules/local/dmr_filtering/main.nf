process DMR_FILTERING {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)
    path chrom_sizes
    path promoters_bed
    path enhancers_bed

    output:
    tuple val(meta), path("${prefix}.bed")            , emit: dmr_bed
    tuple val(meta), path("overlapped_promoters.bed") , emit: promoters_bed
    tuple val(meta), path("overlapped_enhancers.bed") , emit: enhancers_bed
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_asm_dmr"

    """
    awk '\$4 == "different" ${args}' $bed | sort -k13,13gr > ${prefix}.bed

    bedtools slop \\
        -i ${prefix}.bed \\
        -g $chrom_sizes \\
        -b 50 \\
        > ${prefix}_slopped.bed

    bedtools intersect \\
        -a  ${prefix}_slopped.bed \\
        -b $promoters_bed \\
        -wa \\
        -wb \\
        > overlapped_promoters.bed

    bedtools intersect \\
        -a ${prefix}_slopped.bed \\
        -b $enhancers_bed \\
        -wa \\
        -wb \\
        > overlapped_enhancers.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
        sort: \$(sort --version | head -n 1 | awk '{print \$4}')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}