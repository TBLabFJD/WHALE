process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ontresearch/dorado:0.5.3' :
        'docker.io/ontresearch/dorado:0.5.3' }"

    input:
    tuple val(meta), path(pod5)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(model_dir)

    output:
    tuple val(meta), path("${meta.id}_basecall.bam")  , emit: bam
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    dorado \\
        basecaller \\
        $model_dir \\
        $pod5 \\
        $args \\
        --reference $fasta \\
        > ${meta.id}_basecall.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}