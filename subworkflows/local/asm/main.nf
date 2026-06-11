include { MODKIT_PILEUP                } from '../../../modules/nf-core/modkit/pileup'
include { TABIX_TABIX as TABIX_TABIX_1 } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP as TABIX_BGZIP_1 } from '../../../modules/nf-core/tabix/bgzip'
include { MODKIT_DMR                   } from '../../../modules/local/modkit/dmr'
include { DMR_FILTERING                } from '../../../modules/local/dmr_filtering'
include { UCSC_LIFTOVER                } from '../../../modules/nf-core/ucsc/liftover'
include { DMR_cCRE                     } from '../../../modules/local/dmr_cCRE'
include { ANNOTSV_INSTALLANNOTATIONS   } from '../../../modules/nf-core/annotsv/installannotations'
include { ANNOTSV_ANNOTSV              } from '../../../modules/nf-core/annotsv/annotsv' 
include { ANNOTSV_TSV_FILTERING        } from '../../../modules/local/annotsv_tsv_filtering' 
include { MERGE_DMR_FILES              } from '../../../modules/local/merge_dmr_files'
include { METHYLARTIST_LOCUS           } from '../../../modules/local/methylartist/locus'
include { TABIX_TABIX as TABIX_TABIX_2 } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP as TABIX_BGZIP_2 } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_TABIX as TABIX_TABIX_3 } from '../../../modules/nf-core/tabix/tabix'




workflow ASM {

    take:
    ch_bam_bai_haplotypes
    ch_reference
    chrom_sizes
    ch_intervals
    ch_haplotagged_bam_bai
    ch_gtf_tbi
    ch_snv_vcf_gz_tbi
    ch_versions

    main:

    MODKIT_PILEUP ( 
        ch_bam_bai_haplotypes, 
        ch_reference,
        ch_intervals
    )

    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions_modkit.first())

    TABIX_TABIX_1 (
        MODKIT_PILEUP.out.bedgz
    )

    ch_bed_with_tbi = MODKIT_PILEUP.out.bedgz
        .join(TABIX_TABIX_1.out.tbi)

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
        MODKIT_DMR.out.regions_bed,
        chrom_sizes
    )

    ch_versions = ch_versions.mix(DMR_FILTERING.out.versions.first())

    if (params.assembly == 'T2T-CHM13') {
        ch_chain = channel.fromPath(params.chain_T2T_to_hg38, checkIfExists: true)

        UCSC_LIFTOVER (
          DMR_FILTERING.out.dmr_bed,
          ch_chain 
        )

        dmr_bed = UCSC_LIFTOVER.out.lifted
        ch_versions = ch_versions.mix(UCSC_LIFTOVER.out.versions_ucsc.first())
    } else {
        dmr_bed = DMR_FILTERING.out.dmr_bed
    }

    ch_promoters = channel.value([ [id: params.assembly], file(params.promoters_bed, checkIfExists: true) ])
    ch_enhancers = channel.value([ [id: params.assembly], file(params.enhancers_bed, checkIfExists: true) ])
    
    DMR_cCRE (
        dmr_bed,
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

    ch_filtered_dmr = dmr_bed.filter { _meta, bed_file -> 
        bed_file.size() > 0 
    }

    ANNOTSV_ANNOTSV (
        ch_filtered_dmr, // only those which are not empty
        ch_annotations_dir,
        ch_transcripts
    )

    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    ANNOTSV_TSV_FILTERING (
        ANNOTSV_ANNOTSV.out.tsv
    )

    ch_files_to_merge = ANNOTSV_TSV_FILTERING.out.filtered_tsv
        .join(dmr_bed)
        .join(DMR_cCRE.out.promoters_bed)
        .join(DMR_cCRE.out.enhancers_bed)

    MERGE_DMR_FILES (
        ch_files_to_merge
    )

    ch_versions = ch_versions.mix(MERGE_DMR_FILES.out.versions.first())

    if ( params.methylartist_enabled == true ) {

        ch_sample_data = ch_haplotagged_bam_bai
            .join(ch_snv_vcf_gz_tbi)
        
        def intervals_list = params.methylartist_intervals instanceof String 
            ? params.methylartist_intervals.tokenize(',') 
            : params.methylartist_intervals

        ch_plot_intervals = channel.fromList(intervals_list)

        ch_for_methylartist = ch_sample_data.combine(ch_plot_intervals)

        ch_for_methylartist
        .multiMap { meta, bam, bai, vcf, tbi, interval -> 
            sample_data: tuple(meta, bam, bai, vcf, tbi)
            interval:    interval 
        }
        .set { ch_methylartist_inputs }

        ch_methylartist_inputs.sample_data.view { "Data: $it" }

        ch_methylartist_inputs.interval.view { "Interval: $it" }

        METHYLARTIST_LOCUS (
            ch_methylartist_inputs.sample_data,
            ch_reference,
            ch_gtf_tbi,
            ch_methylartist_inputs.interval
        )

        ch_versions = ch_versions.mix(METHYLARTIST_LOCUS.out.versions.first())
    }

    emit:
    dmr_bed = dmr_bed
    pileup_bedgz = MODKIT_PILEUP.out.bedgz
    pileup_tbi = TABIX_TABIX_1.out.tbi
    ch_versions = ch_versions
}