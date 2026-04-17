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
    tuple val("dmr"), path("*.bed"), emit: dmr_bed
    path "versions.yml"           , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_asm_dmr"

    """
    modkit dmr pair \\
        -a $h1_bed \\
        -b $h2_bed \\
        -o ${prefix}.bed \\
        --ref $fasta \\
        --threads $task.cpus \\
        --base C \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """
}