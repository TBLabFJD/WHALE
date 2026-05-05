include { MODKIT_PILEUP                      } from '../../../modules/nf-core/modkit/pileup'
include { TABIX_TABIX                        } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP                        } from '../../../modules/nf-core/tabix/bgzip'
include { MODKIT_DMR                         } from '../../../modules/local/modkit/dmr'
include { DMR_FILTERING                         } from '../../../modules/local/dmr_filtering'

workflow ASM {

    take:
    ch_bam_bai_haplotypes
    ch_reference
    chrom_sizes
    ch_intervals

    main:

    MODKIT_PILEUP ( 
        ch_bam_bai_haplotypes, 
        ch_reference.first(),
        ch_intervals
    )

    TABIX_TABIX (
        MODKIT_PILEUP.out.bedgz
    )

    ch_bed_with_tbi = MODKIT_PILEUP.out.bedgz
        .join(TABIX_TABIX.out.tbi)

    ch_bed_with_tbi.branch {
        h1: it[0].id.endsWith('h1')
        h2: it[0].id.endsWith('h2')
    }.set { ch_haplotypes }

    ch_dmr_input = ch_haplotypes.h1
        .map { meta, bed, tbi -> [ meta.id.replaceAll('_h1$', ''), bed, tbi ] } 
        .join(
            ch_haplotypes.h2.map { meta, bed, tbi -> [ meta.id.replaceAll('_h2$', ''), bed, tbi ] }
        )
        .map { sample_id, h1_bed, h1_tbi, h2_bed, h2_tbi ->
            def new_meta = [ id: sample_id, single_end: false ]
            [ new_meta, h1_bed, h1_tbi, h2_bed, h2_tbi ]
        }

    MODKIT_DMR (
        ch_dmr_input,
        ch_reference.first()
    )

    DMR_FILTERING (
        MODKIT_DMR.out.bed,
        chrom_sizes,
        file(params.promoters_bed),
        file(params.enhancers_bed)
    )

    emit:
    dmr_bed = DMR_FILTERING.out.dmr_bed
    pileup_bedgz = MODKIT_PILEUP.out.bedgz
    pileup_tbi = TABIX_TABIX.out.tbi
}