process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(pod5)
    tuple val(meta2), path(fasta), path(fai)
    val dorado_device
    tuple val(meta3), path(model_dir)

    output:
    tuple val(meta), path("${meta.id}_basecall.bam")  , emit: bam
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def device = dorado_device ? "--device ${dorado_device}" : ''
    """
    dorado \\
        basecaller \\
        $model_dir \\
        $pod5 \\
        $device \\
        $args \\
        --reference $fasta \\
        > ${meta.id}_basecall.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}