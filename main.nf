#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/variantpipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/variantpipeline
    Website: https://nf-co.re/variantpipeline
    Slack  : https://nfcore.slack.com/channels/variantpipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTPIPELINE         } from './workflows/variantpipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_variantpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Rescatamos los valores de igenomes.config automáticamente si se pasa --genome
params.fasta     = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.fasta_gzi = getGenomeAttribute('fasta_gzi')
params.chrom_sizes = getGenomeAttribute('chrom_sizes')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_VARIANTPIPELINE {

    take:
    samplesheet // channel: samplesheet read in from --input
    fasta
    fasta_fai
    fasta_gzi
    chrom_sizes
    
    main:

    //
    // WORKFLOW: Run pipeline
    //
    VARIANTPIPELINE (
        samplesheet,
        fasta,
        fasta_fai,
        fasta_gzi,
        chrom_sizes
    )

    emit:
    multiqc_report = VARIANTPIPELINE.out.multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //

    //
    def ch_fasta     = params.fasta     ? channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()
    def ch_fasta_fai = params.fasta_fai ? channel.fromPath(params.fasta_fai).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()
    def ch_fasta_gzi = params.fasta_gzi ? channel.fromPath(params.fasta_gzi).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()
    def ch_chrom_sizes = params.chrom_sizes ? channel.fromPath(params.chrom_sizes).first() : channel.empty()

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VARIANTPIPELINE (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_gzi,
        ch_chrom_sizes
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_VARIANTPIPELINE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/