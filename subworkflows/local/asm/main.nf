include { MODKIT_PILEUP              } from '../../../modules/nf-core/modkit/pileup'
include { TABIX_TABIX                } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP                } from '../../../modules/nf-core/tabix/bgzip'
include { MODKIT_DMR                 } from '../../../modules/local/modkit/dmr'
include { DMR_FILTERING              } from '../../../modules/local/dmr_filtering'
include { DMR_cCRE                   } from '../../../modules/local/dmr_cCRE'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main' 
include { ANNOTSV_TSV_FILTERING      } from '../../../modules/local/annotsv_tsv_filtering' 
include { MERGE_DMR_FILES            } from '../../../modules/local/merge_dmr_files'


workflow ASM {

    take:
    ch_bam_bai_haplotypes
    ch_reference
    chrom_sizes
    ch_intervals

    main:

    MODKIT_PILEUP ( 
        ch_bam_bai_haplotypes, 
        ch_reference,
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
        ch_reference
    )

    DMR_FILTERING (
        MODKIT_DMR.out.differences_bed,
        chrom_sizes
    )

    ch_promoters = channel.value([ [id: params.assembly], file(params.promoters_bed, checkIfExists: true) ])
    ch_enhancers = channel.value([ [id: params.assembly], file(params.enhancers_bed, checkIfExists: true) ])
    
    DMR_cCRE (
        DMR_FILTERING.out.dmr_bed,
        ch_promoters,
        ch_enhancers
    )

    if (params.annotsv_annotations == 'install') {
        ANNOTSV_INSTALLANNOTATIONS()
        
        ch_annotations_dir = ANNOTSV_INSTALLANNOTATIONS.out.annotations.map { it -> [ [id:'annotsv_db'], it ] }
    } else {
        ch_annotations_dir = params.annotsv_annotations ? channel.fromPath(params.annotsv_annotations).map{ it -> [ [id:it.baseName], it ] }.first() : channel.empty()
    }

    ch_transcripts = channel.value( [ [id: 'empty'], [] ] ) // without parameter because it must be always empty

    ch_filtered_dmr = DMR_FILTERING.out.dmr_bed.filter { _meta, bed_file -> 
        bed_file.size() > 0 
    }

    ANNOTSV_ANNOTSV (
        ch_filtered_dmr, // only those which are not empty
        ch_annotations_dir,
        ch_transcripts
    )

    ANNOTSV_TSV_FILTERING (
        ANNOTSV_ANNOTSV.out.tsv
    )

    ch_files_to_merge = ANNOTSV_TSV_FILTERING.out.filtered_tsv
        .join(DMR_FILTERING.out.dmr_bed)
        .join(DMR_cCRE.out.promoters_bed)
        .join(DMR_cCRE.out.enhancers_bed)

    MERGE_DMR_FILES (
        ch_files_to_merge
    )

    emit:
    dmr_bed = DMR_FILTERING.out.dmr_bed
    pileup_bedgz = MODKIT_PILEUP.out.bedgz
    pileup_tbi = TABIX_TABIX.out.tbi
}