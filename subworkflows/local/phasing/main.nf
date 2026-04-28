include { WHATSHAP_PHASE                     } from '../../../modules/nf-core/whatshap/phase'
include { WHATSHAP_HAPLOTAG                  } from '../../../modules/nf-core/whatshap/haplotag'
include { WHATSHAP_STATS                     } from '../../../modules/nf-core/whatshap/stats'    
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_1 } from '../../../modules/nf-core/samtools/index'
include { WHATSHAP_SPLIT                     } from '../../../modules/local/whatshap/split'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2 } from '../../../modules/nf-core/samtools/index'
include { MODKIT_PILEUP                      } from '../../../modules/nf-core/modkit/pileup'
include { TABIX_TABIX                        } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP                        } from '../../../modules/nf-core/tabix/bgzip'
include { METHPHASER                         } from '../../../modules/local/methphaser'
include { MODKIT_DMR                         } from '../../../modules/local/modkit/dmr'


workflow PHASING {

    take:
    ch_phasing_input
    ch_reference
    fasta
    fasta_fai

    main:

    ch_vcf_tbi = ch_phasing_input.map { meta, vcf, tbi, _bam, _bai -> [ meta, vcf, tbi ] }
    bam_bai    = ch_phasing_input.map { meta, _vcf, _tbi, bam, bai -> [ meta, bam, bai ] }

    WHATSHAP_PHASE (
        ch_vcf_tbi,
        bam_bai,
        ch_reference
    )

    ch_haplotag_input = WHATSHAP_PHASE.out.vcf
        .join(WHATSHAP_PHASE.out.tbi)
        .map { meta, vcf, tbi -> [ meta.id, meta, vcf, tbi ] }
        .join( bam_bai.map { meta, bam, bai -> [ meta.id, bam, bai ] }, by: 0 )
        .map { _id, meta, vcf, tbi, bam, bai -> [ meta, vcf, tbi, bam, bai ] }

    WHATSHAP_HAPLOTAG (
        ch_haplotag_input,
        fasta,
        fasta_fai,
        true
    )

    SAMTOOLS_INDEX_1 (
        WHATSHAP_HAPLOTAG.out.bam
    )

    hap_bam_bai = WHATSHAP_HAPLOTAG.out.bam.join(SAMTOOLS_INDEX_1.out.bai)

    TABIX_BGZIP (
        WHATSHAP_HAPLOTAG.out.tsv
    )

    WHATSHAP_SPLIT (
        hap_bam_bai,
        TABIX_BGZIP.out.output
    )

    ch_bam_h1 = WHATSHAP_SPLIT.out.bam_h1.map { meta, bam ->
        def meta_h1 = meta.clone()
        meta_h1.id  = "${meta.id}_h1"
        [ meta_h1, bam ]
    }

    ch_bam_h2 = WHATSHAP_SPLIT.out.bam_h2.map { meta, bam ->
        def meta_h2 = meta.clone()
        meta_h2.id  = "${meta.id}_h2"
        [ meta_h2, bam ]
    }

    ch_bams_to_index = ch_bam_h1.mix(ch_bam_h2)

    SAMTOOLS_INDEX_2 (
        ch_bams_to_index
    )

    ch_bam_bai_haplotypes = ch_bams_to_index.join(SAMTOOLS_INDEX_2.out.bai)

    WHATSHAP_STATS (
        WHATSHAP_PHASE.out.vcf,
        false,
        true,
        true
    )

    emit:
    phased_vcf  = WHATSHAP_PHASE.out.vcf
    haplotagged_bam_bai = hap_bam_bai
    ch_bam_bai_haplotypes = ch_bam_bai_haplotypes
}