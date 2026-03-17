process DORADO_DOWNLOAD {
    tag "$dorado_model"
    label 'process_single'

    container "docker.io/ontresearch/dorado"

    storeDir "${params.dorado_models_dir ?: './dorado_models'}"

    input:
    tuple val(meta), val(dorado_model)

    output:
    tuple val(meta), path("${dorado_model}"), emit: model

    script:
    """
    dorado download --model ${dorado_model}
    """
}