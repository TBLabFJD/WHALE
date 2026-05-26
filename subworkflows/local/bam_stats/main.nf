include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view'
include { NANOSTAT      } from '../../../modules/local/nanostat'

workflow BAM_STATS {

    take:
    bam_bai
    ch_reference
    ch_reads
    ch_intervals

    main:
    ch_versions = channel.empty()

    SAMTOOLS_VIEW (
        bam_bai,
        ch_reference,
        ch_reads,
        ch_intervals,
        ''
    )

    NANOSTAT (
        SAMTOOLS_VIEW.out.bam
    )

    ch_versions = ch_versions.mix(NANOSTAT.out.versions.first())

    emit:
    nanostat_report = NANOSTAT.out.report
    versions        = ch_versions
}