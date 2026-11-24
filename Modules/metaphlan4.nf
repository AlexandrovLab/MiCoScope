nextflow.enable.dsl=2

process metaphlan4 {
    label 'per_sample'   // --sample_failure_strategy
    label 'process_high'
    scratch true
    // Intermediate read/alignment files are published only with --save_intermediates true;
    publishDir "${params.metaphlan4_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null }
    conda "${Mpa.env(params)}"  // 4.1.1 or 4.2 by database (lib/Mpa.groovy)

    input:
    tuple val(sampleID), path(reads)
 
    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23).
    val deps
    output:
    // sampleID is carried so downstream steps can join this profile back to its
    // sample instead of looking it up by path in the published output directory.
    tuple val(sampleID),
          path("*.metagenome.bowtie2.bz2"),
          path("*.sam.bz2"),
          path("*.profiled_metagenome.txt")


    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = sampleID
    def read_files = reads instanceof Collection ? reads : [reads]
    def inputs = read_files.collect { "\"${it}\"" }.join(' ')

    """
    # A failing stage inside a pipeline must fail the task. No pipe
    # exists in this script today; the switch protects whatever is added later.
    set -o pipefail
    # The [DEBUG] PATH / which / --version echo that was here is replaced by the tool-version record.
    ${ToolVersions.record(task, params, 'metaphlan bowtie2')}
    out1="${SAMPLE_NAME}.metagenome.bowtie2.bz2"
    out2="${SAMPLE_NAME}.sam.bz2"
    out3="${SAMPLE_NAME}.profiled_metagenome.txt"

    # Actual execution
    gzip -dc ${inputs} > "${SAMPLE_NAME}.trimmed.fastq"

    # MetaPhlAn 4.2 renamed --bowtie2db to --db_dir and --bowtie2out to --mapout (the 4.2.x
    # binary rejects the old names: "unrecognized arguments", vJan26 arm, 22 Nov 2026).
    # --index, -s, -o, --nproc and --input_type are unchanged. Detect the dialect from the
    # installed binary so one module serves 4.1 and 4.2.
    if metaphlan --help 2>&1 | grep -q -- '--db_dir'; then
        DB_FLAG=--db_dir; MAP_FLAG=--mapout
    else
        DB_FLAG=--bowtie2db; MAP_FLAG=--bowtie2out
    fi

    # Run MetaPhlAn4
    metaphlan "${SAMPLE_NAME}.trimmed.fastq" \\
        "\$DB_FLAG" "${Mpa.dbDir(params)}" \\
        --index ${Mpa.index(params)} \\
        "\$MAP_FLAG" "${SAMPLE_NAME}.metagenome.bowtie2.bz2" \\
        -s "${SAMPLE_NAME}.sam.bz2" \\
        --nproc ${task.cpus} \\
        --input_type fastq \\
        -o "${SAMPLE_NAME}.profiled_metagenome.txt"
    """
}
