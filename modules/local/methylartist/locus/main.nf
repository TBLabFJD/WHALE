process METHYLARTIST_LOCUS {
    tag "$meta.id - $interval"
    label 'process_single'

    conda "bioconda::methylartist=1.2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methylartist:1.2.9--pyh7cba7a3_0':
        'biocontainers/methylartist:1.2.9--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(vcf_tbi)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(gtf), path(gtf_tbi)
    val interval

    output:
    tuple val(meta), path("*.png")    , emit: png, optional: true
    tuple val(meta), path("*.pdf")    , emit: pdf, optional: true
    tuple val(meta), path("*.svg")    , emit: svg, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${interval.replaceAll(/[:-]/, '_')}"

    """
    methylartist locus \\
        -b ${bam} \\
        -i "${interval}" \\
        -g ${gtf} \\
        --ref ${fasta} \\
        --motif CG \\
        --mods m \\
        --phased \\
        --variants ${vcf} \\
        -p ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methylartist: \$(methylartist --version | sed 's/methylartist //g')
    END_VERSIONS
    """
}