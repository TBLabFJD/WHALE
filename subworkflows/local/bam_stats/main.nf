include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view'
include { NANOSTAT      } from '../../../modules/local/nanostat'


workflow BAM_STATS {

    take:
    bam_bai
    ch_reference
    ch_reads
    ch_intervals

    main:

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

    emit:
    nanostat_report = NANOSTAT.out.report
}