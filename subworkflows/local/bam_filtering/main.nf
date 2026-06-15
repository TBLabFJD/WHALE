include { SAMTOOLS_MARKDUP } from '../../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_VIEW    } from '../../../modules/nf-core/samtools/view'
include { NANOSTAT         } from '../../../modules/local/nanostat'


workflow BAM_FILTERING {

    take:
    not_filtered_bam_bai
    ch_reference
    ch_reads
    ch_intervals
    ch_versions

    main:

    not_filtered_bam = not_filtered_bam_bai.map { meta, bam, _bai -> [ meta, bam ] }

    SAMTOOLS_MARKDUP (      // removal of duplicated reads
        not_filtered_bam,
        ch_reference
    )

    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions_samtools.first())

    SAMTOOLS_INDEX (
        SAMTOOLS_MARKDUP.out.bam
    )

    bam_bai_without_duplicates = SAMTOOLS_MARKDUP.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    SAMTOOLS_VIEW (     // removal of secondary and supplementary alignments
        bam_bai_without_duplicates,
        ch_reference,
        ch_reads,
        ch_intervals,
        "bai"
    )

    filtered_bam_bai = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai)

    emit:
    filtered_bam_bai = filtered_bam_bai
    ch_versions = ch_versions
}