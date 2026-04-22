/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: FINGERPRINTCHECK
    Run Picard CrosscheckFingerprints on all normal + tumor BAMs from the samplesheet.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_CROSSCHECKFINGERPRINTS } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'

workflow FINGERPRINTCHECK {

    take:
    ch_samplesheet    // channel: [ val(meta), val(fastqs) ]  — each element is one samplesheet row
    haplotype_map     // val/path: path to the haplotype map file

    main:
    ch_versions = channel.empty()

    // -------------------------------------------------------------------------
    // Build per-BAM channels (one element per BAM file) for indexing
    // Each samplesheet row carries meta.normal_bam and meta.tumor_bam
    // -------------------------------------------------------------------------
    ch_bams_to_index = ch_samplesheet
        .filter { meta, _bam -> meta.bam }
        .map { meta, _bam ->
            [ [ id: meta.id ], file(meta.bam.toString(), checkIfExists: true) ]
        }

    SAMTOOLS_INDEX(ch_bams_to_index)

    // -------------------------------------------------------------------------
    // Merge BAM paths with freshly generated indexes
    // -------------------------------------------------------------------------
    ch_bam_with_index = ch_bams_to_index
        .join(SAMTOOLS_INDEX.out.index, by: 0)
        // channel: [ meta_bam, bam, index ]

    // -------------------------------------------------------------------------
    // Build exactly one Picard tuple for the whole cohort:
    // [meta, bam_list, index_list, input2_list, input2_index_list, haplotype_map]
    //
    // collect(flat:false) is critical here; it preserves tuple boundaries so we
    // get rows shaped as [meta, bam, index].
    // -------------------------------------------------------------------------
    ch_fingerprint_input = ch_bam_with_index
        .collect(flat: false)
        .map { rows ->
            if (!rows) {
                error("FINGERPRINTCHECK: no BAM/index rows were collected")
            }

            if (rows.any { row -> row.size() != 3 }) {
                error("FINGERPRINTCHECK: expected [meta, bam, index] rows, observed: ${rows}")
            }

            def bams = rows.collect { row -> row[1] }
            def indexes = rows.collect { row -> row[2] }

            if (bams.size() != indexes.size()) {
                error("FINGERPRINTCHECK: BAM/index list size mismatch (${bams.size()} vs ${indexes.size()})")
            }

            def haplotypeMapFile = file(haplotype_map, checkIfExists: true)

            [
                [ id: 'fingerprintcheck' ],
                bams,
                indexes,
                [],   // input2 (SECOND_INPUT) — unused
                [],   // input2_index           — unused
                haplotypeMapFile
            ]
        }

    // Picard allows the reference input to be empty.
    ch_fingerprint_reference = channel
        .value([[:], []])
        .map { metaRef, fastaRef ->
            [ metaRef, fastaRef ]
        }

    PICARD_CROSSCHECKFINGERPRINTS(
        ch_fingerprint_input,
        ch_fingerprint_reference
    )

    emit:
    crosscheck_metrics = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics
    versions           = ch_versions
}
