process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/ontresearch/dorado"

    input:
    tuple val(meta), path(pod5_path)
    tuple val(meta2), path(fasta_path)
    tuple val (meta3), path(fai_path)
    // val dorado_device
    // val dorado_model

    output:
    tuple val(meta), path("${meta.id}_basecall.bam")  , emit: bam_out
    path "versions.yml"                               , emit: versions

    script:
    //def emit_args = (params.dorado_modification == null) ? " --emit-fastq > basecall.fastq && gzip basecall.fastq" : " --modified-bases $params.dorado_modification > basecall.bam"
    """
    dorado download --model dna_r10.4.1_e8.2_400bps_hac@v5.2.0_5mCG_5hmCG@v2
    dorado basecaller hac $pod5_path --device cpu --modified-bases 5mCG_5hmCG --reference $fasta_path > ${meta.id}_basecall.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}