include { SAMTOOLS_VIEW                                     } from '../../../modules/nf-core/samtools/view'
include { NANOSTAT                                          } from '../../../modules/local/nanostat'
include { SAMTOOLS_COVERAGE                                 } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_HAPLOTYPES } from '../../../modules/nf-core/samtools/coverage'
include { WHATSHAP_STATS                                    } from '../../../modules/nf-core/whatshap/stats'    



workflow BAM_STATS {

    take:
    bam_bai
    ch_reference
    ch_reads
    ch_intervals
    phased_vcf
    ch_bam_bai_haplotypes
    ch_versions

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

    SAMTOOLS_COVERAGE (
        bam_bai,
        ch_reference
    )

    SAMTOOLS_COVERAGE_HAPLOTYPES (
        ch_bam_bai_haplotypes,
        ch_reference
    )

    WHATSHAP_STATS (
        phased_vcf,
        true,
        true,
        true
    )

    ch_versions = ch_versions.mix(NANOSTAT.out.versions.first())

    emit:
    nanostat_report            = NANOSTAT.out.report
    samtools_report            = SAMTOOLS_COVERAGE.out.coverage
    samtools_report_haplotypes = SAMTOOLS_COVERAGE_HAPLOTYPES.out.coverage
    whatshap_report            = WHATSHAP_STATS.out.log
    ch_versions                = ch_versions
}