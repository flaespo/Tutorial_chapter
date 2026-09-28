##########
#Numerical examples and complete workflow reported in the chapter:
#Nonnegative Matrix Factorization for RNA-seq Data Analysis: From Foundations to Benchmark Applications
#authored by A. Boccarelli, C. Castiello, N. Del Buono, F. Esposito, L.Selicato
#Corrected and simplified version: fixes SEQC row indexing, removes unused or
#redundant commands, and consolidates repeated NMF scaling operations.
##########
#


#
# -----------------------------------------------------------------------------
# 0. Install if necessary. The packages are required for a first run without 
#    reinstalling packages every time.
# -----------------------------------------------------------------------------
#
#######
# Uncomment the following line the first time you run the R-code
#######
#install.packages(c("NMF", "pheatmap"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "airway", "SummarizedExperiment", "seqc"))

#######
#NMF supplies factorization and rank functions; pheatmap supplies heatmaps. 
#Suppressing messages keeps output compact
######
suppressPackageStartupMessages({
  library(NMF)
  library(pheatmap)
})

####
#Set the global pseudorandom-number seed
#It makes the simulated W, H, depth factors, and Poisson counts reproducible
###


set.seed(20260904)

# Utility functions----------------------------------------------------------

####
#It computes library sizes, converts each sample column to counts per million, 
#adds a pseudocount, and applies log base 2.
#The transformation removes gross library-depth differences, 
#compresses the count range, avoids log(0), and preserves nonnegativity.
####

log_cpm <- function(counts, prior_count = 1) {
  lib_size <- colSums(counts)
  cpm <- t(t(counts) / lib_size) * 1e6
  log2(cpm + prior_count)
}

####
#It provides label-free feature selection and safely handles 
#matrices with fewer than n rows
####

top_variable_rows <- function(x, n = 500) {
  rv <- apply(x, 1, var)
  order(rv, decreasing = TRUE)[seq_len(min(n, nrow(x)))]
}

#####
#Return the row index of the largest H value in every sample column
#The result is the dominant NMF component for each sample
###
dominant_component <- function(H) {
  apply(H, 2, which.max)
}

# Resolve the NMF scaling ambiguity in a reusable way. Columns of W are given
# unit L1 norm; the corresponding scale is transferred to the rows of H. The
# final H columns are normalized to sum to one and can be read as relative
# component exposures.
normalize_nmf_factors <- function(W, H) {
  w_scale <- colSums(W)
  stopifnot(all(is.finite(w_scale)), all(w_scale > 0))

  W_l1 <- sweep(W, 2, w_scale, "/")
  H_adjusted <- sweep(H, 1, w_scale, "*")
  h_total <- colSums(H_adjusted)
  stopifnot(all(is.finite(h_total)), all(h_total > 0))

  H_fraction <- sweep(H_adjusted, 2, h_total, "/")
  list(W_l1 = W_l1, H_adjusted = H_adjusted, H_fraction = H_fraction)
}

# -----------------------------------------------------------------------------
# 1. Simulated data: we know that the true rank is 3
# -----------------------------------------------------------------------------

#Stores the numbers of genes, samples, and true latent components.
n_genes <- 300
n_samples <- 18
k_true <- 3

#Create three six-sample groups and systematic sample and gene identifiers 
#with zero padding
sample_group <- factor(rep(c("A", "B", "C"), each = 6))
sample_names <- paste0("S", sprintf("%02d", seq_len(n_samples)))
gene_names <- paste0("gene_", sprintf("%03d", seq_len(n_genes)))

#####
#Initialize W_true with a small background loading, names its dimensions, 
#and give each component a block of 40 high-loading marker genes.
#This creates three known gene signatures against which estimated factors 
#can be interpreted.
#####

# W_true: gene signatures; each component has 40 main marker genes.
W_true <- matrix(0.05, nrow = n_genes, ncol = k_true,
                 dimnames = list(gene_names, paste0("C", 1:k_true)))
for (j in seq_len(k_true)) {
  idx <- ((j - 1) * 40 + 1):(j * 40)
  W_true[idx, j] <- runif(length(idx), 2, 5)
}

#####
#Initialize H_true with low background activity
#Samples are mixed rather than perfectly pure, 
#which makes the recovery problem more realistic.
#####

# H_true: activity of the components in the samples.
H_true <- matrix(0.2, nrow = k_true, ncol = n_samples,
                 dimnames = list(paste0("C", 1:k_true), sample_names))
for (j in seq_len(k_true)) {
  H_true[j, sample_group == levels(sample_group)[j]] <- 3
}
H_true <- H_true + matrix(runif(k_true * n_samples, 0, 0.4), k_true)

#####
# Multiply W_true and H_true to obtain expected expression, 
#applies sample-specific depth factors, and draws Poisson counts.
#This introduces two RNA-seq features: variable library depth 
#and count sampling noise.
#####

# Simplified RNA-seq counts: Poisson distribution with library depth factors
mu <- W_true %*% H_true
depth <- runif(n_samples, 0.8, 1.3)
mu <- sweep(mu, 2, depth * 30, "*")
counts_sim <- matrix(rpois(length(mu), lambda = as.vector(mu)),
                     nrow = n_genes, ncol = n_samples,
                     dimnames = list(gene_names, sample_names))

X_sim <- log_cpm(counts_sim)
stopifnot(all(is.finite(X_sim)), min(X_sim) >= 0)

# Exploring the matrix.
dim(X_sim)
summary(as.vector(X_sim))

ann_sim <- data.frame(group = sample_group, row.names = sample_names)
pheatmap(X_sim[top_variable_rows(X_sim, 60), ],
         scale = "row", annotation_col = ann_sim,
         show_rownames = FALSE,
         main = "Simulated data: the 60 most variable genes")

####
#Fits repeated NMF models for ranks 2 through 5 and plots rank diagnostics
#The survey illustrates how stability and consensus can inform rank selection
####

# Choice of rank: use a small number of runs in the classroom; increase nrun for a real-world analysis.
rank_sim <- nmfEstimateRank(X_sim, range = 2:5, method = "brunet", nrun = 10, seed = 20260904)
plot(rank_sim)

#Rank 3 matches the known generating model
#Fit the final rank-3 Brunet NMF with 30 starts and extracts W and H.
fit_sim <- nmf(X_sim, rank = 3, method = "brunet",
               nrun = 30, seed = 20260904)
W_sim <- basis(fit_sim)  # genes x components
H_sim <- coef(fit_sim)   # components x samples

dim(W_sim)
dim(H_sim)

pheatmap(H_sim, scale = "row", annotation_col = ann_sim,
         main = "Activity of the components in the simulated samples")

pred_sim <- factor(dominant_component(H_sim), levels = 1:3)
table(componente_NMF = pred_sim, group_known = sample_group)

# Genes with the greatest weighting in each component.
top_genes_sim <- lapply(seq_len(ncol(W_sim)), function(j) {
  head(sort(W_sim[, j], decreasing = TRUE), 10)
})
names(top_genes_sim) <- colnames(W_sim)
top_genes_sim

# Relative reconstruction error.
Xhat_sim <- W_sim %*% H_sim
rel_error_sim <- sqrt(sum((X_sim - Xhat_sim)^2)) / sqrt(sum(X_sim^2))
rel_error_sim

# EXERCISE 1
# Repeat fit_sim with rank 2 and rank 4. Compare:
# (a) relative error; (b) separation of clusters; (c) interpretability of the signatures.


# -----------------------------------------------------------------------------
# 2. Real-world airway dataset
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(SummarizedExperiment)
})

data(airway)
se <- airway
colData(se)$dex <- droplevels(colData(se)$dex)

####
#Extract counts, retains genes with at least 10 counts in at least three samples,
#and subset the complete container
####
counts_airway <- assay(se)
keep <- rowSums(counts_airway >= 10) >= 3
se_filt <- se[keep, ]

####
#Construct a DESeq2 object with cell and dex in the design, 
#estimate size factors, and extracts normalized counts.
####
dds <- DESeqDataSet(se_filt, design = ~ cell + dex)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

####
#Log-transforms normalized counts, selects 1,000 variable genes, 
#subsets the matrix, and checks NMF validity.
#
# log2(normalised counts + 1) preserves the non-negativity required by NMF.
####

X_airway_all <- log2(norm_counts + 1)
sel <- top_variable_rows(X_airway_all, n = 1000)
X_airway <- X_airway_all[sel, ]
stopifnot(all(is.finite(X_airway)), min(X_airway) >= 0)

meta <- as.data.frame(colData(se_filt)[, c("cell", "dex")])
meta[] <- lapply(meta, as.factor)

dim(X_airway)
table(meta$dex)

####
#Plot the 60 most variable airway genes with sample annotations.
####
pheatmap(X_airway[top_variable_rows(X_airway, 60), ],
         scale = "row", annotation_col = meta,
         show_rownames = FALSE,
         main = "airway: 60 of the most variable genes")

# With only 8 samples, the rank-based analysis is illustrative, not conclusive.
rank_airway <- nmfEstimateRank(X_airway, range = 2:4,method = "brunet", nrun = 15,seed = 20260904)
plot(rank_airway)

fit_airway <- nmf(X_airway, rank = 2, method = "brunet",nrun = 50, seed = 20260904)
W_airway <- basis(fit_airway)
H_airway <- coef(fit_airway)

pheatmap(H_airway, scale = "row", annotation_col = meta,
         main = "airway: activity of NMF components")

# Descriptive table: dominant component and experimental metadata.
sample_summary <- data.frame(
  campione = colnames(X_airway),
  cellula = meta$cell,
  trattamento = meta$dex,
  componente_dominante = dominant_component(H_airway),
  t(H_airway),
  row.names = NULL,
  check.names = FALSE
)
sample_summary

top_genes_airway <- lapply(seq_len(ncol(W_airway)), function(j) {
  head(sort(W_airway[, j], decreasing = TRUE), 20)
})
names(top_genes_airway) <- colnames(W_airway)
top_genes_airway

# To what extent does each component distinguish between treatments and controls?
old_par <- par(mfrow = c(1, nrow(H_airway)))
for (j in seq_len(nrow(H_airway))) {
  boxplot(H_airway[j, ] ~ meta$dex,
          xlab = "Dexamethasone", ylab = "Weight",
          main = paste("Component", j))
}
par(old_par)

# EXERCISE
# (a) Repeat the analysis using 500 and 2000 variable genes.
# (b) Check whether the dominant component follows ‘dex’ or ‘cell’.
# (c) Examine the first 20 genes of each component: what further biological analyses
#     would you suggest to assign a functional significance?

# -----------------------------------------------------------------------------
# 3. SEQC/MAQC-III benchmark: recovery of the known mixtures with k = 2
# -----------------------------------------------------------------------------

# The SEQC/MAQC-III experiment contains two reference RNAs and two mixtures:
#   A = Universal Human Reference RNA (UHRR)
#   B = Human Brain Reference RNA (HBRR)
#   C = 75% A + 25% B
#   D = 25% A + 75% B
#
# Following the manuscript, we begin with a single Illumina sequencing centre
# (AGR). This avoids asking NMF to explain mixture composition and cross-centre
# technical variation simultaneously. The physical construction suggests k = 2.

# data(..., package = "seqc") loads the object without attaching the package.
data("ILM_refseq_gene_AGR", package = "seqc")
seqc_df <- ILM_refseq_gene_AGR

# The seqc package documents the first four columns as feature metadata:
# EntrezID, Symbol, GeneLength and IsERCC. All remaining columns are counts.
stopifnot(is.data.frame(seqc_df), ncol(seqc_df) > 4)
feature_info <- seqc_df[, 1:4, drop = FALSE]

# ERCC spike-ins are deliberately set aside. They can be used for an external
# validation, but they are not included in the matrix factorized here.
stopifnot(!anyNA(feature_info[[4]]))
is_ercc <- tolower(trimws(as.character(feature_info[[4]]))) %in%
  c("true", "t", "1")

counts_seqc <- as.matrix(seqc_df[!is_ercc, -(1:4), drop = FALSE])
storage.mode(counts_seqc) <- "integer"

# Use the RefSeq/Entrez identifier in the first feature column. Missing or
# duplicated identifiers are made unique so that every matrix row has a stable
# name; no rows are reordered by this operation.
gene_id_seqc <- make.unique(as.character(feature_info[[1]][!is_ercc]))
rownames(counts_seqc) <- gene_id_seqc

stopifnot(
  nrow(counts_seqc) > 0,
  ncol(counts_seqc) > 0,
  !anyNA(counts_seqc),
  all(counts_seqc >= 0)
)

# Illumina SEQC library names have the form
# (Sample)_(Replicate)_(Lane)_(FlowCell), for example A_1_L01_FlowCellA.
# The explicit checks make the code fail visibly if the package naming scheme
# changes, rather than silently assigning incorrect A/B/C/D labels.
sample_class <- factor(
  sub("_.*", "", colnames(counts_seqc)),
  levels = c("A", "B", "C", "D")
)
if (anyNA(sample_class)) {
  stop(
    "Could not infer the A/B/C/D classes from the SEQC column names. ",
    "Inspect colnames(counts_seqc) and update the parsing rule."
  )
}
if (!all(levels(sample_class) %in% unique(as.character(sample_class)))) {
  stop("The AGR table does not contain all four expected SEQC samples A-D.")
}

table(sample_class)

# Filter genes with very low information: retain genes with at least 10 counts
# in at least three libraries. This rule uses no knowledge of A/B/C/D labels.
keep_seqc <- rowSums(counts_seqc >= 10) >= 3
counts_seqc <- counts_seqc[keep_seqc, , drop = FALSE]

# DESeq2 size factors correct between-library scale and composition effects.
# design = ~ 1 is sufficient because the labels are used for validation only;
# they are not used to create the matrix supplied to NMF.
dds_seqc <- DESeqDataSetFromMatrix(
  countData = counts_seqc,
  colData = data.frame(
    sample_class = sample_class,
    row.names = colnames(counts_seqc)
  ),
  design = ~ 1
)
dds_seqc <- estimateSizeFactors(dds_seqc)
norm_seqc <- counts(dds_seqc, normalized = TRUE)

# Select the 2,000 genes with the largest variance on the log-normalized scale.
# Selection is unsupervised: it does not use the known mixture labels, avoiding
# leakage of the benchmark answer into the feature-selection step.
log_seqc <- log2(norm_seqc + 1)
gene_var_seqc <- apply(log_seqc, 1, var)
selected_seqc <- order(gene_var_seqc, decreasing = TRUE)[
  seq_len(min(2000, length(gene_var_seqc)))
]

# Two complementary representations are retained, as discussed in the paper.
# X_seqc_log compresses the dynamic range and is useful for exploration.
X_seqc_log <- log_seqc[selected_seqc, , drop = FALSE]

# Physical RNA mixing is additive on the linear scale. Therefore, the primary
# benchmark is the normalized-count matrix. Dividing each row by its mean keeps
# highly abundant genes from dominating the NMF objective while preserving the
# A/B mixture pattern of each selected gene.
X_seqc_linear <- norm_seqc[selected_seqc, , drop = FALSE]
row_mean_seqc <- rowMeans(X_seqc_linear)
positive_rows_seqc <- is.finite(row_mean_seqc) & row_mean_seqc > 0
X_seqc_linear <- X_seqc_linear[positive_rows_seqc, , drop = FALSE]
X_seqc_linear <- sweep(
  X_seqc_linear,
  1,
  row_mean_seqc[positive_rows_seqc],
  "/"
)

stopifnot(
  all(is.finite(X_seqc_log)),
  all(is.finite(X_seqc_linear)),
  min(X_seqc_log) >= 0,
  min(X_seqc_linear) >= 0,
  identical(colnames(X_seqc_log), colnames(X_seqc_linear))
)

# A heatmap of the most variable genes provides an initial view of the sample
# structure. The A/B/C/D annotation is displayed but was not used for selection.
meta_seqc <- data.frame(
  sample = sample_class,
  row.names = colnames(X_seqc_log)
)
pheatmap(
  X_seqc_log[top_variable_rows(X_seqc_log, 60), , drop = FALSE],
  scale = "row",
  annotation_col = meta_seqc,
  show_rownames = FALSE,
  main = "SEQC AGR: 60 most variable genes"
)

# Primary analysis: rank 2 corresponds to the two source RNAs. Multiple random
# starts are used because NMF is non-convex and can converge to local solutions.
fit_seqc <- nmf(X_seqc_linear,rank = 2,method = "brunet",nrun = 50,seed = 20260904)
W_seqc <- basis(fit_seqc)
H_seqc <- coef(fit_seqc)

stopifnot(ncol(W_seqc) == 2,nrow(H_seqc) == 2,all(is.finite(W_seqc)),
          all(is.finite(H_seqc)),min(W_seqc) >= 0,min(H_seqc) >= 0)

# Resolve the scaling ambiguity W H = (W D)(D^-1 H) with the helper defined
# above. The resulting H columns are relative component exposures.
scaled_seqc <- normalize_nmf_factors(W_seqc, H_seqc)
W_seqc_l1 <- scaled_seqc$W_l1
H_seqc_fraction <- scaled_seqc$H_fraction

# Component numbers are arbitrary. Orient them after fitting: the A-like
# component is whichever component has the larger mean exposure in pure A.
mean_in_A <- rowMeans(H_seqc_fraction[, sample_class == "A", drop = FALSE])
component_A <- which.max(mean_in_A)
component_B <- setdiff(seq_len(nrow(H_seqc_fraction)), component_A)
stopifnot(length(component_A) == 1, length(component_B) == 1)

estimated_A <- as.numeric(H_seqc_fraction[component_A, ])
estimated_B <- as.numeric(H_seqc_fraction[component_B, ])

# Known physical fraction of source A in each SEQC sample. The principal test
# is recovery of the order A > C > D > B, rather than exact calibration.
expected_A_by_class <- c(A = 1.00, B = 0.00, C = 0.75, D = 0.25)
expected_A <- unname(expected_A_by_class[as.character(sample_class)])
stopifnot(!anyNA(expected_A))

benchmark_seqc <- data.frame(library = colnames(X_seqc_linear),
  sample = sample_class,
  expected_A = expected_A, estimated_A = estimated_A, estimated_B = estimated_B,
  row.names = NULL
)

# Mean exposures should reproduce the intended A-C-D-B progression. Spearman
# correlation evaluates this ordering; mean absolute error gives a transparent
# calibration summary but is not expected to be zero because normalization,
# gene selection and NMF scaling affect the numerical exposures.
benchmark_summary_seqc <- aggregate(
  cbind(expected_A, estimated_A, estimated_B) ~ sample,
  data = benchmark_seqc,
  FUN = mean
)

# IMPORTANT: the comma selects rows. Without it, a data frame would be indexed
# as a list and the command would reorder columns instead of sample-class rows.
benchmark_summary_seqc <- benchmark_summary_seqc[
  match(c("A", "C", "D", "B"), benchmark_summary_seqc$sample),
  ,
  drop = FALSE
]

spearman_seqc <- cor(
  benchmark_seqc$expected_A,
  benchmark_seqc$estimated_A,
  method = "spearman"
)
mae_seqc <- mean(abs(benchmark_seqc$expected_A - benchmark_seqc$estimated_A))

benchmark_summary_seqc
spearman_seqc
mae_seqc

# Visual comparison with the identity line. The colours distinguish the four
# known sample types while every point remains an individual RNA-seq library.
seqc_colours <- c(A = "#1B9E77", B = "#D95F02", C = "#7570B3", D = "#E7298A")
plot(benchmark_seqc$expected_A,benchmark_seqc$estimated_A,
  xlab = "Known fraction of sample A",
  ylab = "NMF-estimated A-like exposure",
  pch = 19,
  col = unname(seqc_colours[as.character(benchmark_seqc$sample)]),
  xlim = c(0, 1),
  ylim = c(0, 1)
)
abline(0, 1, col = "steelblue", lwd = 2)
legend(
  "topleft",
  legend = c("A", "B", "C", "D"),
  col = seqc_colours[c("A", "B", "C", "D")],
  pch = 19,
  title = "SEQC sample"
)

# The exposure heatmap shows the same result at library level. No row scaling
# is used because H_seqc_fraction already contains comparable fractions.
rownames(H_seqc_fraction) <- c("Component 1", "Component 2")
rownames(H_seqc_fraction)[component_A] <- "A-like component"
rownames(H_seqc_fraction)[component_B] <- "B-like component"
pheatmap(
  H_seqc_fraction,
  scale = "none",
  annotation_col = meta_seqc,
  main = "SEQC AGR: normalized component exposures"
)
###
#Sorts each column of W_sim and retains the ten largest gene weights, then prints the list.
###
# Genes with the largest weights are candidates for biological interpretation;
# they are not, by themselves, formal differentially expressed genes.
top_genes_seqc <- list(
  A_like = head(sort(W_seqc_l1[, component_A], decreasing = TRUE), 20),
  B_like = head(sort(W_seqc_l1[, component_B], decreasing = TRUE), 20)
)
top_genes_seqc

# Relative reconstruction error is reported for completeness, but rank should
# never be selected from reconstruction error alone.
Xhat_seqc <- W_seqc %*% H_seqc
rel_error_seqc <- sqrt(sum((X_seqc_linear - Xhat_seqc)^2)) /
  sqrt(sum(X_seqc_linear^2))
rel_error_seqc

# Optional sensitivity analysis described in the manuscript: fit k = 2 to the
# log-scale representation. Agreement of the A-C-D-B ordering across scales
# strengthens the conclusion; disagreement identifies preprocessing sensitivity.
fit_seqc_log <- nmf(
  X_seqc_log,
  rank = 2,
  method = "brunet",
  nrun = 50,
  seed = 20260904
)
scaled_seqc_log <- normalize_nmf_factors(
  basis(fit_seqc_log),
  coef(fit_seqc_log)
)
H_seqc_log_fraction <- scaled_seqc_log$H_fraction
component_A_log <- which.max(rowMeans(
  H_seqc_log_fraction[, sample_class == "A", drop = FALSE]
))
estimated_A_log <- as.numeric(H_seqc_log_fraction[component_A_log, ])

scale_comparison_seqc <- data.frame(
  library = colnames(X_seqc_linear),
  sample = sample_class,
  expected_A = expected_A,
  linear_estimated_A = estimated_A,
  log_estimated_A = estimated_A_log,
  row.names = NULL
)

# Display and quantify the optional log-scale sensitivity analysis.
spearman_log_seqc <- cor(
  scale_comparison_seqc$expected_A,
  scale_comparison_seqc$log_estimated_A,
  method = "spearman"
)
scale_comparison_seqc
spearman_log_seqc

# EXERCISE
# (a) Repeat the analysis with 500, 1000 and 5000 variable genes.
# (b) Survey ranks 2:5 with nmfEstimateRank; rank 2 is the biological working
#     hypothesis, but a reproducible extra component may represent lane/flowcell.
# (c) Repeat the analysis at another Illumina centre in the seqc package and
#     compare the A-like ordering and the overlap of the top-loading genes.

# -----------------------------------------------------------------------------
# 4. Reproducibility
# -----------------------------------------------------------------------------

sessionInfo()

