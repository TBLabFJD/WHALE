include { DORADO_BASECALLER         } from '../../../modules/local/dorado/basecaller'
include { DORADO_DOWNLOAD           } from '../../../modules/local/dorado/download'
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_MERGE            } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX            } from '../../../modules/nf-core/samtools/index'

workflow BASECALLING {

    take:
    samplesheet
    ch_reference
    fasta_gzi

    main:
    model_ch = channel.of([ 
        [id: 'dorado_model'], 
        params.dorado_model
    ])

    DORADO_DOWNLOAD ( model_ch )

    DORADO_BASECALLER (
        samplesheet,
        ch_reference,
        DORADO_DOWNLOAD.out.model.first()
    )

    SAMTOOLS_SORT (
        DORADO_BASECALLER.out.bam,
        ch_reference,
        "bai"
    )

    ch_bams_for_merge = SAMTOOLS_SORT.out.bam
        .groupTuple() // it merges those BAMs with the same ID

    ch_reference = ch_reference.map {_meta, fa, fai -> [ fa, fai ] }
        .combine( fasta_gzi.map { _meta, gzi -> gzi } )
        .map { fa, fai, gzi -> [ [id:'genome'], fa, fai, gzi ] }
        .first()

    SAMTOOLS_MERGE (
        ch_bams_for_merge,
        ch_reference
    )

    SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

    emit:
    bam_bai = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_INDEX.out.bai)

}