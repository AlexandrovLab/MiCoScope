# Run production cohort/filter code on independent known-answer fixtures.
# This test requires base R only; it does NOT execute decontam or a complete workflow.
argv <- commandArgs(TRUE); root <- normalizePath(argv[1]); out <- normalizePath(argv[2])
source(file.path(root, "scripts/batch_contract.R"))
expect_error <- function(expr, text) {
  result <- tryCatch({force(expr); NULL}, error=function(e) conditionMessage(e))
  stopifnot(!is.null(result), grepl(text, result, fixed=TRUE))
}
# Two batches, each with ten samples of each TCGA tissue, covariates vary within batch.
ids <- paste0("S",1:60)
types <- rep(rep(c("Primary Tumor","Solid Tissue Normal","Blood Derived Normal"),each=10),2)
meta <- data.frame(sampleid=ids,shipment_batch=rep(c("B1","B2"),each=30),
  type=types,age_diag=rep(40:49,6),sex=rep(c("F","M"),30),row.names=ids)
counts <- data.frame(Taxon=rep(20,60),row.names=ids)
data <- list(otu=counts,meta=meta,unused_metadata=character())
args <- list(sample_id_column="sampleid",batch_column="shipment_batch",type_column="type",
  tumor_only=TRUE,tumor_value="Primary Tumor",prefix=file.path(out,"tumor","run"))
dir.create(dirname(args$prefix),recursive=TRUE,showWarnings=FALSE)
r <- prepare_batch_cohort(data,args,c("age_diag","sex"))
stopifnot(identical(rownames(r$otu),ids[types=="Primary Tumor"]))
ledger <- read.delim(file.path(dirname(args$prefix),"batch_correction_sample_status.tsv"))
stopifnot(sum(ledger$status=="included")==20,
  all(ledger$reason[ledger$status=="excluded"]=="outside requested tumor subset"))
args$tumor_value <- "Tumor"
expect_error(prepare_batch_cohort(data,args,c("age_diag","sex")),"No eligible samples")
args$tumor_only <- FALSE; args$prefix <- file.path(out,"mixed","run")
dir.create(dirname(args$prefix),recursive=TRUE,showWarnings=FALSE)
r <- prepare_batch_cohort(data,args,c("age_diag","sex","type"))
stopifnot(identical(rownames(r$otu),ids))
expect_error(prepare_batch_cohort(data,args,c("age_diag","bmi")),"Missing required columns: bmi")
confounded <- data;confounded$meta$type <- rep(c("Primary Tumor","Blood Derived Normal"),each=30)
expect_error(prepare_batch_cohort(confounded,args,c("age_diag","type")),"Requested covariate does not vary within batch: type")
# Extract the exact production filter and final selection expressions, never duplicate their formula.
code <- readLines(file.path(root,"scripts/decontamination.R"))
block <- function(start,end) {
 a <- which(code==start);b <- which(code==end)
 stopifnot(length(a)==1,length(b)==1,b>a)
 code[a:(b-1)]
}
filter <- block("selected_counts <- count_matrix[, selected_samples, drop = FALSE]",
 "if (length(evaluated_taxa) > 0L && length(eligible_batches) < args$min_batches) {")
final <- block("final_taxa <- evaluated_taxa[!final_contaminant_flag]",
 "# Known-contaminant flag (June-Young's list). Species rank only: the list names species, and a genus")
writeLines(filter,file.path(out,"production_filter_fragment.R"))
writeLines(final,file.path(out,"production_final_fragment.R"))
selected_samples <- paste0("S",1:40)
taxon_ids <- c("Always","Count4","Count5","Below","At","Above")
count_matrix <- matrix(0,6,41,dimnames=list(taxon_ids,paste0("S",1:41)))
count_matrix["Always",] <- 20
count_matrix["Count4",] <- 4;count_matrix["Count5",] <- 5
count_matrix["Below",1] <- 5;count_matrix["At",c(1,21)] <- 5
count_matrix["Above",c(1,2,21)] <- 5
count_matrix[,41] <- 999  # Unselected sample must not enter denominator/count sums.
args <- list(min_abundance=5,min_prevalence=.05)
eval(parse(text=filter))
stopifnot(ncol(filtered_counts)==40,all(filtered_counts["Count4",]==0),
 all(filtered_counts["Count5",]==5),filtered_prevalence["Below"]==.025,
 filtered_prevalence["At"]==.05,filtered_prevalence["Above"]==.075,
 identical(evaluated_taxa,c("Always","Count5","At","Above")))
final_contaminant_flag <- rep(FALSE,length(evaluated_taxa));eval(parse(text=final))
stopifnot(identical(rownames(final_counts),c("Always","Count5","At","Above")),
 sum(final_counts["At",])==10,!"Below" %in% rownames(final_counts))
for (threshold in c(.0499999,.0500001)) {
 args$min_prevalence <- threshold;eval(parse(text=filter))
 stopifnot(("At" %in% evaluated_taxa)==(threshold<.05),"Above" %in% evaluated_taxa,
  !"Below" %in% evaluated_taxa)
}
writeLines(c("PASS: tumor exact label and mixed inclusion", "PASS: missing BMI and complete tissue/batch confounding rejected",
 "PASS: production count4/5, prevalence below/at/above5%, selected denominator40, final removal expressions",
 "LIMIT: decontam package execution and whole workflow not tested here"),file.path(out,"base_R_checks.txt"))
cat("Base-R production contract checks PASS\n")
