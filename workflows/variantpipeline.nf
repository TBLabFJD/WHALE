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
include { BASECALLING            } from '../subworkflows/local/basecalling'  
include { BAM_FILTERING          } from '../subworkflows/local/bam_filtering'
include { PHASING                } from '../subworkflows/local/phasing'
include { ASM                    } from '../subworkflows/local/asm'
include { BAM_STATS              } from '../subworkflows/local/bam_stats'
include { BAM_STATS as HAPLOTYPES_BAM_STATS } from '../subworkflows/local/bam_stats'

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
    chrom_sizes
    gtf_gz
    gtf_tbi

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    def bam_bai = channel.empty()

    ch_reference = fasta.map { meta, fa -> [ meta, fa ] }
        .combine( fasta_fai.map { _meta, fai -> fai } )
        .map { meta, fa, fai -> [ meta, fa, fai ] }
        .first()

    ch_gtf_tbi = gtf_gz.map { meta, gtf -> [ meta, gtf ] }
        .combine( gtf_tbi.map { _meta, tbi -> tbi } )
        .map { meta, gtf, tbi -> [ meta, gtf, tbi ] }
        .first()

    ch_reads = params.reads ?
        channel.value([ [id: 'selected_reads'], params.reads ]) : // in case we want to extract cercain reads
        channel.value([ [id: 'selected_reads'], [] ])
    ch_intervals = params.intervals ?
        // in case we want to extract cercain intervals
        channel.value([ [id: 'selected_intervals'], params.intervals ]) :
        channel.value([ [id: 'selected_intervals'], [] ])

    if (params.step == 'mapping') {
        FASTQC (
            samplesheet
        )

        ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect{ _meta, zip -> zip } )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        MINIMAP2_ALIGN (
            samplesheet,
            fasta,
            true,
            "bai",
            false,
            false
        )

        bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index)
    }
    else if (params.step == 'basecalling') {
        BASECALLING (
            samplesheet,
            ch_reference,
            fasta_gzi,
            ch_versions
        )

        bam_bai = BASECALLING.out.bam_bai
        ch_versions = ch_versions.mix( BASECALLING.out.ch_versions )

    }
    else if (params.step == 'variant_calling') {
            bam_bai = samplesheet
    }
    else if (params.step == 'phasing') {
        def phasing_inputs = samplesheet.multiMap { meta, bam, bai, vcf, tbi ->
            bam_bai_ch: [ meta, bam, bai ]
            vcf_ch:     [ meta, vcf ]
            tbi_ch:     [ meta, tbi ]
        }
        
        bam_bai       = phasing_inputs.bam_bai_ch
        final_snv_vcf = phasing_inputs.vcf_ch
        final_snv_tbi = phasing_inputs.tbi_ch
    }
    else if (params.step == 'snv_annotation') {
        merged_vcf = samplesheet
    }
    else if (params.step == 'sv_annotation') {
        merged_final_bed = samplesheet
    }
    else if (params.step == 'asm') {
        ch_bam_bai_haplotypes = samplesheet
    }

    //
    // SUBWORKFLOW: Run Deepvariant, Clair3 & NanoCaller
    //

    if (params.snv_calling == true) {

        SNV_CALLING (
            bam_bai,
            fasta,
            fasta_fai,
            fasta_gzi
        )

        snv_calling_vcfs = SNV_CALLING.out.deepvariant_vcf_tbi.concat(SNV_CALLING.out.nanocaller_vcf_tbi, SNV_CALLING.out.clair3_vcf_tbi)
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

        merged_vcf = MERGE_SNV_CALLING.out.final_vcf
        final_snv_vcf = MERGE_SNV_CALLING.out.final_vcf
        final_snv_tbi = MERGE_SNV_CALLING.out.final_tbi

        ch_snv_vcf_gz_tbi = MERGE_SNV_CALLING.out.final_vcf_gz
            .join(MERGE_SNV_CALLING.out.final_tbi)

    } else if (params.step == 'snv_annotation') {
        final_snv_vcf = samplesheet.map { meta, vcf, _tbi -> [meta, vcf] }
        final_snv_tbi = samplesheet.map { meta, _vcf, tbi -> [meta, tbi] }
    }

    //
    // PHASING AND DMR
    //

    if (params.phasing == true) {
        ch_snv_vcf_tbi = final_snv_vcf.join(final_snv_tbi)

        ch_phasing_input = ch_snv_vcf_tbi.join(bam_bai) // now the bam and tbi files are grouped by the same meta      

        PHASING (
            ch_phasing_input,
            ch_reference,
            fasta,
            fasta_fai,
            ch_versions
        )

        bam_bai = PHASING.out.haplotagged_bam_bai
        ch_bam_bai_haplotypes = PHASING.out.ch_bam_bai_haplotypes
        ch_multiqc_files = ch_multiqc_files.mix( PHASING.out.whatshap_stats_report.map { _meta, txt -> txt } )
        ch_versions = ch_versions.mix( PHASING.out.ch_versions )
    }

    if (params.asm && !params.phasing) {
        error "ERROR: You must perform phasing (--phasing true) if you want to run ASM (--asm true)"
    } else if (params.asm && params.phasing) {  
        BAM_FILTERING (   
            ch_bam_bai_haplotypes,
            ch_reference,
            ch_reads,
            ch_intervals,
            ch_versions
        )

        asm_input = BAM_FILTERING.out.filtered_bam_bai
        ch_versions = ch_versions.mix( BAM_FILTERING.out.ch_versions )
        // def asm_input = ch_bam_bai_haplotypes ?: samplesheet

        ASM (
            asm_input,
            ch_reference,
            chrom_sizes,
            ch_intervals,
            bam_bai, // haplotagged
            ch_gtf_tbi,
            ch_snv_vcf_gz_tbi,
            ch_versions
        )

        ch_versions = ch_versions.mix( ASM.out.ch_versions )

        HAPLOTYPES_BAM_STATS (
            asm_input,
            ch_reference,
            ch_reads,
            ch_intervals,
            ch_versions
        )

        ch_multiqc_files = ch_multiqc_files.mix( HAPLOTYPES_BAM_STATS.out.samtools_report.map { _meta, txt -> txt } )
    }

    //
    // SUBWORKFLOW: Run Merge SV Calling
    //

    if (params.sv_calling == true) {
        SV_CALLING (
            bam_bai,
            fasta,
            fasta_fai,
            fasta_gzi
        )
    }

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
    // VARIANT ANNOTATION
    //

    if (params.snv_annotation == true) {
        SNV_ANNOTATION (
            merged_vcf, 
            fasta
        )
    }

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

    //
    // BAM STATS
    //

    if (bam_bai) {

        BAM_STATS (
            bam_bai,
            ch_reference,
            ch_reads,
            ch_intervals,
            ch_versions
        )

        ch_multiqc_files = ch_multiqc_files.mix( BAM_STATS.out.nanostat_report.map { _meta, txt -> txt } )
        ch_multiqc_files = ch_multiqc_files.mix( BAM_STATS.out.samtools_report.map { _meta, txt -> txt } )
        
        ch_versions = ch_versions.mix( BAM_STATS.out.versions )
    }

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

    ch_multiqc_config = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    
    MULTIQC (
        ch_multiqc_files.toSortedList(),
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