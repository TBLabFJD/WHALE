include { SV_SAMPLES                 } from '../../../modules/local/sv_samples/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main' 

workflow SV_ANNOTATION {

    take:
    merged_bed

    main:
    if (params.sv_database == true) {
        
        multiinter = params.sv_multiinter ? channel.fromPath(params.sv_multiinter).collect() : channel.empty()

        SV_SAMPLES (
            merged_bed,
            multiinter
        )
    }

    if (params.annotsv_annotations == 'install') {
        ANNOTSV_INSTALLANNOTATIONS()
        
        annotations = ANNOTSV_INSTALLANNOTATIONS.out.annotations
    } else {
        annotations = params.annotsv_annotations ? channel.fromPath(params.annotsv_annotations).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()
    }

    gene_transcripts = params.gene_transcripts ? channel.fromPath(params.gene_transcripts).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()

    annotsv_input = params.sv_database == true ? SV_SAMPLES.out.samples_info_bed : merged_bed

    ANNOTSV_ANNOTSV (
        annotsv_input,
        annotations,
        gene_transcripts
    )

    emit:
    annotated_tsv = ANNOTSV_ANNOTSV.out.tsv
}
