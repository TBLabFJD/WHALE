process MODKIT_DMR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta), path(h1_bed), path(h1_tbi), path(h2_bed), path(h2_tbi)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}_every_CpG.bed") , emit: CpG_bed
    tuple val(meta), path("${prefix}_regions.bed")   , emit: regions_bed
    path "versions.yml"                              , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    modkit dmr pair \\
        -a ${h1_bed} \\
        -b ${h2_bed} \\
        -o ${prefix}_every_CpG.bed \\
        --ref ${fasta} \\
        --threads ${task.cpus} \\
        --base C \\
        --segment ${prefix}_regions.bed \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """
}