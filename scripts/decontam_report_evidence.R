# Reporting only. Sourced by decontamination.R after scientific outputs are final.
# Matrices are read here; neither filtering nor decontam decisions are changed.
report_dir <- sprintf("%s.decontam_report", args$prefix)
if (dir.exists(report_dir)) stopf("Refusing to overwrite existing report evidence: %s", report_dir)
dir.create(report_dir)
report_stages <- list(
  input = report_input_counts,
  declared_removal = count_matrix,
  selected_cohort = selected_counts,
  abundance = filtered_counts,
  prevalence = filtered_counts[evaluated_taxa, , drop = FALSE],
  final = final_counts
)
flow <- do.call(rbind, lapply(names(report_stages), function(stage) {
  m <- report_stages[[stage]]
  k <- if (flag_known) norm_species(rownames(m)) %in% known_species else rep(FALSE, nrow(m))
  data.frame(stage = stage, samples = ncol(m), taxon_rows = nrow(m),
    detected_taxa = sum(rowSums(m) > 0), count_sum = sum(m),
    known_count_sum = if (flag_known) sum(m[k, , drop = FALSE]) else NA_real_,
    known_detected_taxa = if (flag_known) sum(rowSums(m[k, , drop = FALSE]) > 0) else NA_integer_)
}))
flow$lost_count <- c(NA_real_, -diff(flow$count_sum))
flow$lost_detected_taxa <- c(NA_integer_, -diff(flow$detected_taxa))
write_tsv(flow, file.path(report_dir, "stage_flow.tsv"))

# Complete reason ledger, including explicitly removed taxa absent from filter_report.
ledger <- data.frame(taxon_id = rownames(report_input_counts),
  input_all_samples = rowSums(report_input_counts),
  input_selected_samples = rowSums(report_input_counts[, selected_samples, drop = FALSE]),
  after_abundance = 0, final_count = 0, known_contaminant = NA,
  final_action = "remove", terminal_reason = "declared_taxon_removal",
  n_evaluable_batches = NA_integer_, detected_batch_count = NA_integer_,
  stringsAsFactors = FALSE)
if (flag_known) ledger$known_contaminant <- norm_species(ledger$taxon_id) %in% known_species
ix <- match(filter_report$taxon_id, ledger$taxon_id)
ledger$after_abundance[ix] <- filter_report$after_abundance_sum
ledger$final_action[ix] <- filter_report$final_action
ledger$terminal_reason[ix] <- filter_report$terminal_reason
ledger$final_count[match(rownames(final_counts), ledger$taxon_id)] <- rowSums(final_counts)
ix <- match(contaminants_final$taxon_id, ledger$taxon_id)
ledger$n_evaluable_batches[ix] <- contaminants_final$n_evaluable_batches
ledger$detected_batch_count[ix] <- contaminants_final$detected_batch_count
write_tsv(ledger, file.path(report_dir, "taxon_ledger.tsv"))

# Selected cohort only. Denominators are explicit, including zero-count samples.
groups <- unique(data.frame(batch = selected_batches,
  sample_type = validated_metadata[[args$type_column]], role = validated_roles,
  stringsAsFactors = FALSE))
group_rows <- list(); gi <- 0L
for (stage in c("selected_cohort", "abundance", "prevalence", "final")) {
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    ids <- selected_samples[selected_batches == g$batch &
      validated_metadata[[args$type_column]] == g$sample_type & validated_roles == g$role]
    m <- report_stages[[stage]][, ids, drop = FALSE]
    gi <- gi + 1L
    group_rows[[gi]] <- data.frame(stage = stage, batch = g$batch,
      sample_type = g$sample_type, role = g$role, samples = length(ids),
      count_sum = sum(m), detected_taxa = sum(rowSums(m) > 0),
      positive_cells = sum(m > 0), possible_cells = length(m),
      samples_with_counts = sum(colSums(m) > 0),
      eligible = g$batch %in% eligible_batches, stringsAsFactors = FALSE)
  }
}
write_tsv(do.call(rbind, group_rows), file.path(report_dir, "batch_type_summary.tsv"))

# Per-taxon detection frequencies at the actual decontam test input, both roles.
detection_rows <- list(); di <- 0L
for (i in seq_len(nrow(groups))) {
  g <- groups[i, ]
  ids <- selected_samples[selected_batches == g$batch &
    validated_metadata[[args$type_column]] == g$sample_type & validated_roles == g$role]
  m <- report_stages$prevalence[, ids, drop = FALSE]
  if (nrow(m)) {
    di <- di + 1L
    detection_rows[[di]] <- data.frame(taxon_id = rownames(m), batch = g$batch,
      sample_type = g$sample_type, role = g$role, samples = length(ids),
      detected_samples = rowSums(m > 0), count_sum = rowSums(m),
      eligible = g$batch %in% eligible_batches, stringsAsFactors = FALSE)
  }
}
detection <- if (length(detection_rows)) do.call(rbind, detection_rows) else data.frame(
  taxon_id=character(), batch=character(), sample_type=character(), role=character(),
  samples=integer(), detected_samples=integer(), count_sum=numeric(), eligible=logical())
write_tsv(detection, file.path(report_dir, "taxon_detection.tsv"))

for (key in c("sample_validation", "batch_eligibility", "contaminants_by_batch", "contaminants_final", "filter_report", "decontaminated")) {
  suffix <- if (key == "decontaminated") ".csv" else ".tsv"
  if (!file.copy(output_files[[key]], file.path(report_dir, paste0(key, suffix)))) stopf("Report copy failed: %s", key)
}
extra <- c(report_schema = "1", input_source = args$input_source,
  count_unit = if (args$start_mode == "beginning") "Bracken-estimated counts after upstream selection/masking; not unique reads" else "Provided integer count units; upstream read/count semantics not verified",
  input_filename = basename(args$otu_table), remove_taxa = args$remove_taxa,
  blacklist_status = if (args$taxonomic_rank != "species") "not_applicable_genus_rank" else if (is.null(args$known_contaminants) || !nzchar(args$known_contaminants)) "not_provided" else if (!length(known_species)) "empty_list" else "applied_as_flag_only",
  blacklist_filename = if (!is.null(args$known_contaminants) && nzchar(args$known_contaminants)) basename(args$known_contaminants) else "not_provided",
  blacklist_sha256 = if (!is.null(args$known_contaminants) && nzchar(args$known_contaminants)) sha256_file(args$known_contaminants) else "not_provided",
  blacklist_species = as.character(length(known_species)),
  blacklist_matching = "Exact species name after lowercasing and space/underscore normalization; not substring matching",
  sample_blacklist_denominator = "All input samples: per-sample count sum AFTER declared taxon removal, BEFORE abundance/prevalence filtering",
  prevalence_denominator = as.character(length(selected_samples)),
  sensitivity_status = "not_assessed", depth_resampling_status = "not_assessed",
  alternative_mode_status = "not_assessed", false_discovery_calibration = "not_assessed",
  report_helper_sha256 = sha256_file(file.path(dirname(script_path), "decontam_report_evidence.R")))
manifest <- rbind(manifest, data.frame(key = names(extra), value = unname(extra)))
report_files <- list.files(report_dir, full.names = TRUE)
manifest <- rbind(manifest, data.frame(key = paste0("report_sha256.", basename(report_files)),
  value = vapply(report_files, sha256_file, character(1))))
