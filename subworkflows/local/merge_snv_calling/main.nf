include { BCFTOOLS_NORM          } from '../../../modules/nf-core/bcftools/norm/main'                     
include { BCFTOOLS_VIEW          } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_MERGE         } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX as TABIX_1 } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_2 } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_3 } from '../../../modules/nf-core/tabix/tabix/main'
include { GET_VCF_FIELDS         } from '../../../modules/local/get_vcf_fields/main'
include { CONSENSUS_GT           } from '../../../modules/local/consensus_gt/main'
include { SOFTWARE_INFO          } from '../../../modules/local/software_info/main'
include { SAMPLE_INFO            } from '../../../modules/local/sample_info/main'
include { FORMAT_VCF             } from '../../../modules/local/format_vcf/main'
include { HEADER_VARIANTS_VCF    } from '../../../modules/local/header_variants_vcf/main'
include { TABIX_BGZIP            } from '../../../modules/nf-core/tabix/bgzip/main'

workflow MERGE_SNV_CALLING {

    take:
    ch_snv_calling_vcfs
    fasta
    fasta_fai

    main:
    BCFTOOLS_NORM (
        ch_snv_calling_vcfs,
        fasta)
    
    BCFTOOLS_VIEW (
        BCFTOOLS_NORM.out.vcf,
        [],
        [],
        [])

    TABIX_1 (
        BCFTOOLS_VIEW.out.vcf)
    
    vcfs = BCFTOOLS_VIEW.out.vcf.groupTuple()
    tbis = TABIX_1.out.tbi.groupTuple()
    vcfs_tbis = vcfs.join(tbis)

    BCFTOOLS_MERGE (
        vcfs_tbis,
        fasta,
        fasta_fai,
        [])

    TABIX_2 (
        BCFTOOLS_MERGE.out.merged_variants)
    
    merged_vcf_tbi = BCFTOOLS_MERGE.out.merged_variants.join(TABIX_2.out.tbi) // channel: [ val(meta), path(merged_vcf), path(tbi) ]

    GET_VCF_FIELDS (
        merged_vcf_tbi)

    CONSENSUS_GT (
        GET_VCF_FIELDS.out.gt,
        GET_VCF_FIELDS.out.caller_order)
    
    SOFTWARE_INFO (
        GET_VCF_FIELDS.out.gt,
        GET_VCF_FIELDS.out.caller_order)
    
    sample_info_input = merged_vcf_tbi.join(CONSENSUS_GT.out.gt_consensus).join(GET_VCF_FIELDS.out.ad_mean).join(GET_VCF_FIELDS.out.dp_mean).join(GET_VCF_FIELDS.out.vaf).join(SOFTWARE_INFO.out.sf).join(CONSENSUS_GT.out.gt_discordances)

    SAMPLE_INFO (
        sample_info_input)
    
    format_vcf_input = merged_vcf_tbi.join(GET_VCF_FIELDS.out.caller_order)

    FORMAT_VCF (
        format_vcf_input)
    
    header_variants_vcf_input = merged_vcf_tbi.join(GET_VCF_FIELDS.out.caller_order).join(FORMAT_VCF.out.format_vcf).join(SAMPLE_INFO.out.sample_info)

    HEADER_VARIANTS_VCF (
        header_variants_vcf_input)

    TABIX_BGZIP (
        HEADER_VARIANTS_VCF.out.final_vcf)
    
    TABIX_3 (
        TABIX_BGZIP.out.output)

    emit:
    final_vcf = HEADER_VARIANTS_VCF.out.final_vcf
    final_tbi = TABIX_3.out.tbi
}