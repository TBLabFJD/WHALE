include { FORMAT2INFO            } from '../../../modules/local/format2info/main'
include { SPLITVCFPVM            } from '../../../modules/local/splitvcfpvm/main'
include { TABIX_BGZIP            } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'
include { AUTOMAP                } from '../../../modules/local/automap/main'
include { ENSEMBLVEP_VEP         } from '../../../modules/nf-core/ensemblvep/vep/main'
include { POSTVEP                } from '../../../modules/local/postvep/main'
include { MERGETSV               } from '../../../modules/local/mergetsv/main'

workflow SNV_ANNOTATION {

    take:
    merged_vcf
    fasta

    main:
    FORMAT2INFO (
        merged_vcf
    )

    SPLITVCFPVM (
        merged_vcf,
        params.n_vcf_variants_split
    )

    TABIX_BGZIP (
        FORMAT2INFO.out.vcf_to_annotate
    )
    
    TABIX_TABIX (
        TABIX_BGZIP.out.output
    )

    sample_info = TABIX_BGZIP.out.output.concat(TABIX_TABIX.out.tbi, FORMAT2INFO.out.fields).groupTuple()

    AUTOMAP (
        merged_vcf,
		params.assembly,
		projectDir
    )

    // Initialize vep_extra_files

    // vep plugins
    vep_plugin_files = []

    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.dbscSNV}", checkIfExists: true))
    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.dbscSNV_tbi}", checkIfExists: true))

    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.dbNSFP}", checkIfExists: true))
    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.dbNSFP_tbi}", checkIfExists: true))

    vep_plugin_files.add(file("${params.vep_annotation_gene_dir}/${params.loFtool}", checkIfExists: true))
    vep_plugin_files.add(file("${params.vep_annotation_gene_dir}/${params.exACpLI}", checkIfExists: true))

    vep_plugin_files.add(file("${params.vep_annotation_gene_dir}/${params.maxEntScan}", checkIfExists: true))

    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.cADD_INDELS}", checkIfExists: true))
    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.cADD_INDELS_tbi}", checkIfExists: true))

    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.cADD_SNVS}", checkIfExists: true))
    vep_plugin_files.add(file("${params.vep_annotation_dir}/${params.cADD_SNVS_tbi}", checkIfExists: true))

    vep_plugin_files_all = vep_plugin_files ? channel.fromPath(vep_plugin_files).collect() : channel.empty()

    // vep custom

    vep_custom_files = []

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.kaviar}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.kaviar_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mAF_FJD_COHORT}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mAF_FJD_COHORT_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cCRS_DB}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cCRS_DB_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cLINVAR}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.cLINVAR_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.dENOVO_DB}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.dENOVO_DB_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_cov}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_cov_tbi}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADe_tbi}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_cov}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_cov_tbi}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.gNOMADg_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mutScore}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.mutScore_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_SNV}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_SNV_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_INDEL}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.spliceAI_INDEL_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.REVEL}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.REVEL_tbi}", checkIfExists: true))

    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.CSVS_dir}/${params.cSVS}", checkIfExists: true))
    vep_custom_files.add(file("${params.vep_annotation_dir}/${params.CSVS_dir}/${params.cSVS_tbi}", checkIfExists: true))

    vep_custom_files_all = vep_custom_files ? channel.fromPath(vep_custom_files).collect() : channel.empty()

    vcf_vep = SPLITVCFPVM.out.vcfs.transpose().combine(sample_info, by: 0).combine(vep_custom_files_all.toList()) // [[meta, vcf], split, [format2info_files], [custom_files]]

    ENSEMBLVEP_VEP (
        vcf_vep,
        params.vep_genome, 
        params.vep_species, 
        params.vep_cache_version, 
        params.vep_cache, 
        fasta,
        vep_plugin_files_all
    )
    
    dbNSFP_gene_path = "${params.vep_annotation_gene_dir}/${params.dbNSFP_gene}"
    dbNSFP_gene = dbNSFP_gene_path ? channel.fromPath(dbNSFP_gene_path).collect() : channel.empty()

    omim_path = "${params.vep_annotation_gene_dir}/${params.omim}"
    omim = omim_path ? channel.fromPath(omim_path).collect() : channel.empty()

    regiondict_path = "${params.vep_annotation_gene_dir}/${params.regiondict}"
    regiondict = regiondict_path ? channel.fromPath(regiondict_path).collect() : channel.empty()

    domino_path = "${params.vep_annotation_gene_dir}/${params.domino}"
    domino = domino_path ? channel.fromPath(domino_path).collect() : channel.empty()

    tissue_expression_path = "${params.vep_annotation_gene_dir}/${params.tissue_expression}"
    tissue_expression = tissue_expression_path ? channel.fromPath(tissue_expression_path).collect() : channel.empty()
    
    postvep_input = ENSEMBLVEP_VEP.out.tsv.groupTuple().join(AUTOMAP.out.roh_automap)

    POSTVEP (
        postvep_input,
        dbNSFP_gene,
        omim,
        regiondict,
        domino,
        tissue_expression,
        params.maf,
        params.genefilter,
        params.glowgenes,
        params.assembly,
        projectDir
    )

    MERGETSV (
        POSTVEP.out.pvm_tsv,
        params.assembly
    )
}
