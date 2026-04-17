include { DORADO_BASECALLER         } from '../../../modules/local/dorado/basecaller'
include { DORADO_DOWNLOAD           } from '../../../modules/local/dorado/download'
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_MERGE            } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX            } from '../../../modules/nf-core/samtools/index'

workflow BASECALLING {

    take:
    samplesheet
    fasta
    fasta_fai
    fasta_gzi

    main:
    model_ch = channel.of([ 
        [id: 'dorado_model'], 
        params.dorado_model
    ])

    DORADO_DOWNLOAD ( model_ch )

    ch_reference = fasta.map { _meta, file -> file }
        .combine( fasta_fai.map { _meta, file -> file } )
        .map { fa, fai -> [ [id:'genome'], fa, fai ] }

    DORADO_BASECALLER (
        samplesheet,
        ch_reference.first(),
        DORADO_DOWNLOAD.out.model.first()
    )

    SAMTOOLS_SORT (
        DORADO_BASECALLER.out.bam,
        ch_reference.first(),
        "bai"
    )

    ch_bams_for_merge = SAMTOOLS_SORT.out.bam
        .map { _meta, bam -> bam } 
        .collect()
        .map { bams -> [ [id:'merged'], bams ] }

    ch_reference = fasta.map { _meta, file -> file }
        .combine( fasta_fai.map { _meta, file -> file } )
        .combine( fasta_gzi.map { _meta, file -> file } )
        .map { fa, fai, gzi -> [ [id:'genome'], fa, fai, gzi ] }

    SAMTOOLS_MERGE (
        ch_bams_for_merge,
        ch_reference
    )

    SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

    emit:
    bam_bai = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_INDEX.out.bai)

}