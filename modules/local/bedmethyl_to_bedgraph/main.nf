process BEDMETHYL_TO_BEDGRAPH {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'ubuntu:22.04' }"

    input:
    tuple val(meta), path(bedgz)

    output:
    tuple val(meta), path("*bedGraph"), emit: bedGraph
    
    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    zcat ${bedgz} | awk -v OFS='\\t' '
    \$4 = "m" {
        pos = \$1"\\t"\$2"\\t"\$3;
        suma[pos] += \$11;
        cuenta[pos]++;
    } 
    END {
        for (p in suma) {
            print p, suma[p]/cuenta[p];
        }
    }' | sort -k1,1 -k2,2n > ${prefix}.bedGraph
    """
}