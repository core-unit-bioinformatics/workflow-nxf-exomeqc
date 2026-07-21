/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { FASTQSCREEN_FASTQSCREEN } from '../modules/nf-core/fastqscreen/fastqscreen/main'
include { DRAGEN_CNV_QC           } from '../modules/local/dragen_cnv_qc/main'
include { FINGERPRINTCHECK        } from '../subworkflows/local/fingerprintcheck/main'
include { samplesheetToList       } from 'plugin/nf-schema'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_exomeqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXOMEQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    // initialize empty channels for collecting software versions and multiQC files
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Collect output files for multiQC
    // using the nf-schema samplesheetToList function for adding directories to multiqc via config file
    channel
        .fromList(samplesheetToList(params.multiqc_inputs, "${projectDir}/assets/multiqc_input_schema.json"))
        .set { ch_multiqc_inputs }
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_inputs)

    if (params.enable_fastqscreen) {
        if (!params.fastqscreen_database) {
            error("Please provide --fastqscreen_database when --enable_fastqscreen is set")
        }

        FASTQSCREEN_FASTQSCREEN(
            ch_samplesheet,
            file(params.fastqscreen_database, checkIfExists: true)
        )

        // MultiQC parses FastQ Screen from *_screen.txt output files.
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQSCREEN_FASTQSCREEN.out.txt.map { _meta, txt -> txt }
        )
    }

    if (params.enable_fingerprintcheck) {
        if (!params.fingerprint_haplotype_map) {
            error("Please provide --fingerprint_haplotype_map when --enable_fingerprintcheck is set")
        }
        if (!ch_samplesheet.first().map { meta, _reads -> meta.bam }) {
            error("Fingerprint check requires samplesheet column 'bam' for all samples")
        }

        FINGERPRINTCHECK(
            ch_samplesheet,
            params.fingerprint_haplotype_map
        )

        // MultiQC expects plain file paths; fingerprintcheck outputs are tuple(meta, path).
        ch_multiqc_files       = ch_multiqc_files.mix(
            FINGERPRINTCHECK.out.crosscheck_metrics.map { _meta, crosscheck_metrics -> crosscheck_metrics }
        )
    }

    if (params.enable_cnvqc) {
        if (!params.dragen_dirs) {
            error("Please provide --dragen_dirs when --enable_cnvqc is set")
        }
        
        // add dragen folder to multiqc
        ch_dragen_dirs         = params.dragen_dirs ? Channel.fromPath(params.dragen_dirs).collect() : []
        ch_multiqc_files       = ch_multiqc_files.mix(ch_dragen_dirs)


        DRAGEN_CNV_QC(params.dragen_dirs)

        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.metrics)
        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.json)
        ch_multiqc_files       = ch_multiqc_files.mix(DRAGEN_CNV_QC.out.plots)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'exomeqc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //

    //add DRAGEN dir to multiqc inputs
    if (params.dragen_dirs) {
        ch_multiqc_files = ch_multiqc_files.mix(Channel.fromPath(params.dragen_dirs))
    }
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.mix(ch_multiqc_custom_config).toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
