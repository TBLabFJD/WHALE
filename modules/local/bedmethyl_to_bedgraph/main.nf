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
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # 1. Filtramos primero (drástica reducción de tamaño)
    zcat ${bedgz} \\
        | awk -v OFS='\\t' '\$4 == "m" { print \$1, \$2, \$3, \$11 }' \\
        | sort -k1,1 -k2,2n -T . \\
        | awk -v OFS='\\t' '
            # 2. Leemos la primera línea para inicializar variables
            NR == 1 { chr=\$1; start=\$2; end=\$3; sum=\$4; count=1; next }
            
            # 3. Si la coordenada es la misma que la anterior, sumamos
            \$1 == chr && \$2 == start && \$3 == end { sum += \$4; count++ }
            
            # 4. Si la coordenada cambia, imprimimos el promedio anterior y reseteamos
            \$1 != chr || \$2 != start || \$3 != end {
                print chr, start, end, sum/count;
                chr=\$1; start=\$2; end=\$3; sum=\$4; count=1;
            }
            
            # 5. Imprimimos el último registro al llegar al final del archivo
            END {
                if (NR > 0) print chr, start, end, sum/count;
            }
        ' > ${prefix}.bedGraph
    """
}