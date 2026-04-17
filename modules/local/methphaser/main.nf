process METHPHASER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/methphaser:0.0.3--90a582334457857e"
    
    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(tbi), path(gtf)
    tuple val(meta_ref), path(fasta), path(fai)

    output:
    tuple val(meta), path("*/*.vcf")       , emit: vcf
    tuple val(meta), path("*/*.{txt,tsv}") , emit: stats
    tuple val(meta), path("*/")            , emit: results_dir
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_methphaser"
    
    """
    meth_phaser_parallel \\
        -t $task.cpus \\
        -b $bam \\
        -r $fasta \\
        -g $gtf \\
        -vc $vcf \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methphaser: \$(meth_phaser_parallel -h | grep "usage" | awk '{print "installed"}')
    END_VERSIONS
    """
}