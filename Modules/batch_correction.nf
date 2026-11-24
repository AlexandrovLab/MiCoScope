// =============================================================================
// BATCH CORRECTION MODULE
// =============================================================================
//
// Batch effect correction and normalization of microbiome data using ConQuR
//
// =============================================================================

nextflow.enable.dsl=2

process PrepareBatchInput {
    tag "${prefix}_${rank}"
    label 'process_low'
    conda "${params.consensus_taxa_env}"
    publishDir path: { "${params.batch_corr_dir}/${prefix}_${rank}" }, mode: 'copy',
        saveAs: { name -> name in ['empty_counts.csv', 'batch_skip_status.tsv'] || name.endsWith('.batch_sample_status.tsv') ? name : null }
    input:
    tuple val(prefix), path(otu_table), path(metadata), val(rank)
    val deps
    output:
    tuple val(prefix), path('ready_counts.tsv'), path(metadata), val(rank), optional: true, emit: ready
    path("${rank}.batch_sample_status.tsv"), optional: true, emit: skipped
    path('empty_counts.csv'), optional: true, emit: empty_table
    path('batch_skip_status.tsv'), optional: true, emit: skip_status
    script:
    """
    set -euo pipefail
    ${ToolVersions.record(task, params, 'python3:python')}
    python3 "${params.scripts}/prepare_batch_input.py" \\
        --input "${otu_table}" --metadata "${metadata}" --rank "${rank}" \\
        --taxon-column "${params.taxon_column}" --sample-column "${params.sample_id_column}" \\
        --batch-column "${params.batch_column}" --type-column "${params.type_column}"
    """
}

process BatchCorrection {

    tag "${prefix}_${rank}"
    label 'process_high'

    conda "${params.batch_corr_env}"
    // Dynamic path: see the note in decontamination.nf. One directory per rank, e.g.
    // decontam_run_species/ and decontam_run_genus/. The QC-only ledger copy is not published.
    publishDir path: { "${params.batch_corr_dir}/${prefix}_${rank}" }, mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith('.batch_sample_status.tsv') ? null : name }

    input:
    // rank: species | genus -- the level of THIS table.
    tuple val(prefix), path(otu_table), path(metadata), val(rank)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23).
    val deps
    output:
    tuple val(prefix),
          path("corrected/*"),
          path("normalized/*.tsv"),
          path("pcoa_plots/*"),
          path("permanova_summary.tsv"),
          path("batch_correction_sample_status.tsv"),
          path("batch_correction_taxon_status.tsv"),
          path("batch_correction_quality.tsv"),
          path("batch_correction_methods.tsv"),
          path("batch_correction_run_settings.tsv"),
          path("batch_correction_effect_sizes.tsv"),
          path("batch_correction_dispersion.tsv"), emit: corrected
    path("${rank}.batch_sample_status.tsv"), emit: qc_status
    // Current-run evidence for the cohort report's batch panel.
    path("${prefix}_${rank}.batch_report"), emit: report_evidence

    script:
    // `--tumor_only false` arrives as the String "false" on Nextflow 25+, which is
    // truthy in Groovy. Compare explicitly rather than testing truthiness.
    def tumor_only_flag = (params.tumor_only.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']) ? "--tumor_only" : ""
    def method_flag = params.batch_corr_method ?: "tune"
    def allow_flag = params.batch_corr_allow_uncorrected.toString() == 'true' ? "--allow_uncorrected" : ""
    // Decontamination controls leave the correction cohort unless the run declares otherwise.
    def include_controls_flag = (params.downstream_include_controls.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']) ? "--include_controls" : ""
    def control_values_arg = params.control_values ? "--control_values \"${params.control_values}\"" : ""
    // Sample type is protected when several types remain and every batch holds two or more.
    // Patient for within-patient permutations: a declared column, else the TCGA barcode with the preset.
    def patient_arg = params.batch_corr_patient_column ? "--patient_column \"${params.batch_corr_patient_column}\"" :
                      (params.cohort_preset == 'tcga' ? "--patient_from_tcga_barcode" : "")
    def protect_flag = (params.batch_corr_protect_sample_type.toString().trim().toLowerCase() in ['false', 'no', 'off', '0']) ? "--no_protect_sample_type" : ""
    // Depth as a covariate: off by default; the cohort report shows depth by batch.
    def depth_flag = (params.batch_corr_depth_covariate.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']) ? "--depth_covariate" : ""
    // Extra before/after effect-size comparison on a named metadata column; none by default.
    def effect_arg = params.batch_corr_effect_column ? "--effect_column \"${params.batch_corr_effect_column}\"" : ""
    """
    ${ToolVersions.record(task, params, 'Rscript:r-base')}
    export R_LIBS="\$CONDA_PREFIX/lib/R/library"
    export R_LIBS_USER="\$CONDA_PREFIX/lib/R/library"
    export R_LIBS_SITE="\$CONDA_PREFIX/lib/R/library"
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
    # ConQuR comes from the repository's local conda channel (conda_envs/local_channel/). A
    # pre-built env passed with --batch_corr_env may lack it: fail here with the fix rather
    # than partway through the R script.
    # Presence is not enough: it must be the pinned fork the script was written against.
    Rscript --vanilla -e 'if (!requireNamespace("ConQuR", quietly = TRUE)) stop("ConQuR is not installed in this environment (env: ", Sys.getenv("CONDA_PREFIX"), "). Use the default environment (conda_envs/lock/batch_correction_env.linux-64.txt, built by Nextflow) or rebuild the prefix with scripts/setup_batch_env.sh"); sha <- utils::packageDescription("ConQuR")\$RemoteSha; if (is.null(sha) || !startsWith(sha, "${params.conqur_sha}")) stop("ConQuR in this environment is not ivartb/ConQuR_par@${params.conqur_sha} (found: ", if (is.null(sha)) "no RemoteSha" else sha, "; env: ", Sys.getenv("CONDA_PREFIX"), "). Use the default environment or rebuild the prefix with scripts/setup_batch_env.sh")'

    # Upstream SHA alone does not identify our compatibility patch.
    Rscript --vanilla -e 'if (!identical(utils::packageDescription("ConQuR")\$CMPipelinePatch, "adonis2-terms-v1")) stop("ConQuR compatibility patch missing; rebuild with scripts/setup_batch_env.sh")'

    # Every value quoted: the TCGA preset's type column is 'Sample Type_x'.
    Rscript --vanilla "${params.batch_corr_script}" \\
        --otu "${otu_table}" \\
        --meta "${metadata}" \\
        --prefix "${prefix}_${rank}" \\
        --batch_column "${params.batch_column}" \\
        --covariates "${params.batch_corr_covariates}" \\
        --type_column "${params.type_column}" \\
        --sample_id_column "${params.sample_id_column}" \\
        --taxon_column "${params.taxon_column}" \\
        --tumor_value "${params.batch_corr_tumor_values}" \\
        ${allow_flag} \\
        ${tumor_only_flag} \\
        ${control_values_arg} \\
        ${include_controls_flag} \\
        ${protect_flag} \\
        ${depth_flag} \\
        ${patient_arg} \\
        ${effect_arg} \\
        --phase ${params.phase} \\
        --r2_threshold ${params.r2_threshold} \\
        --method ${method_flag} \\
        --seed ${params.batch_corr_seed} \\
        --threads ${task.cpus}

    # Rank-named copy of the sample ledger for the cohort QC summary.
    cp batch_correction_sample_status.tsv "${rank}.batch_sample_status.tsv"

    # Evidence bundle for the cohort report: the files the batch panel reads, with their hashes.
    bundle="${prefix}_${rank}.batch_report"
    mkdir -p "\$bundle"
    cp batch_correction_quality.tsv batch_correction_run_settings.tsv batch_correction_methods.tsv \\
       batch_correction_sample_status.tsv permanova_summary.tsv batch_correction_effect_sizes.tsv \\
       batch_correction_dispersion.tsv corrected/BATCH_CORRECTION_STATUS.tsv "\$bundle/"
    {
      printf 'key\\tvalue\\n'
      printf 'report_schema\\t1\\n'
      printf 'rank\\t%s\\n' "${rank}"
      for f in batch_correction_quality.tsv batch_correction_run_settings.tsv batch_correction_methods.tsv \\
               batch_correction_sample_status.tsv permanova_summary.tsv batch_correction_effect_sizes.tsv \\
               batch_correction_dispersion.tsv BATCH_CORRECTION_STATUS.tsv; do
        printf 'report_sha256.%s\\t%s\\n' "\$f" "\$(sha256sum "\$bundle/\$f" | cut -d' ' -f1)"
      done
    } > "\$bundle/manifest.tsv"
    """
}
