/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variantpipeline_pipeline'
include { SNV_CALLING            } from '../subworkflows/local/snv_calling'
include { SV_CALLING             } from '../subworkflows/local/sv_calling'
include { MERGE_SNV_CALLING      } from '../subworkflows/local/merge_snv_calling'
include { MERGE_SV_CALLING       } from '../subworkflows/local/merge_sv_calling'
include { SNV_ANNOTATION         } from '../subworkflows/local/snv_annotation'
include { SV_ANNOTATION          } from '../subworkflows/local/sv_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTPIPELINE {

    take:
    samplesheet // channel: samplesheet read in from --input
    fasta
    fasta_fai
    fasta_gzi

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    

 /*    if (params.step == 'mapping') {
        FASTQC (
            samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        MINIMAP2_ALIGN (
            samplesheet,
            fasta,
            true,
            "bai",
            false,
            false)
    
        bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index) // channel: [ val(meta), path(bam), path(bai) ]
    }
    else if (params.step == 'variant_calling') {
        bam_bai = samplesheet
    }
    else if (params.step == 'snv_annotation') {
        merged_vcf = samplesheet
    }
    else if (params.step == 'sv_annotation') {
        merged_final_bed = samplesheet
    }

    //
    // SUBWORKFLOW: Run Deepvariant, Clair3 & NanoCaller
    //

    if (params.snv_calling == true) {
        SNV_CALLING (
            bam_bai,
            fasta,
            fasta_fai,
            fasta_gzi)
        
        snv_calling_vcfs = SNV_CALLING.out.deepvariant_vcf_tbi.concat(SNV_CALLING.out.nanocaller_vcf_tbi, SNV_CALLING.out.clair3_vcf_tbi)
    }

    //
    // SUBWORKFLOW: Run Sniffles, CuteSV & SVIM
    //

    if (params.sv_calling == true) {
        SV_CALLING (
            bam_bai,
            fasta,
            fasta_fai,
            fasta_gzi)
    }

    //
    // SUBWORKFLOW: Run Merge SNV Calling
    //
    
    if (params.merge_snv == true) {
        MERGE_SNV_CALLING (
            snv_calling_vcfs,
            fasta,
            fasta_fai
        )

        merged_vcf = MERGE_SNV_CALLING.out.final_vcf // channel: [ val(meta), path(merged_vcf) ]
    }
    
    //
    // SUBWORKFLOW: Run Merge SV Calling
    //
    
    if (params.merge_sv == true) {
        MERGE_SV_CALLING (
            SV_CALLING.out.sniffles_vcf,
            SV_CALLING.out.cutesv_vcf,
            SV_CALLING.out.svim_vcf
        )

        merged_gt_bed = MERGE_SV_CALLING.out.merged_gt
        merged_final_bed = MERGE_SV_CALLING.out.merged_final
    }

    //
    // SUBWORKFLOW: Run SNV Annotation
    //
    

    if (params.snv_annotation == true) {
        SNV_ANNOTATION (
            merged_vcf, 
            fasta
        )
    }

    //
    // SUBWORKFLOW: Run SV Annotation
    //
    
    if (params.sv_annotation == true) {
        if (params.sv_database == true) {
            SV_ANNOTATION (
                merged_gt_bed
            )
        }
        else {
            SV_ANNOTATION (
                merged_final_bed
            )
        }
    }
 */
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
