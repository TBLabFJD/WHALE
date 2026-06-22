process DORADO_DOWNLOAD {
    tag "$dorado_model"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ontresearch/dorado:0.5.3' :
        'docker.io/ontresearch/dorado:0.5.3' }"

    storeDir "${params.dorado_models_dir ?: './dorado_models'}"

    input:
    tuple val(meta), val(dorado_model)

    output:
    tuple val(meta), path("${dorado_model}"), emit: model
    path "versions.yml"                     , emit: versions

    script:
    """
    dorado download --model ${dorado_model}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}