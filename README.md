# Nonnegative Matrix Factorization for RNA-seq Data Analysis

R tutorial for the chapter **Nonnegative Matrix Factorization for RNA-seq Data Analysis: From Foundations to Benchmark Applications**, by **A. Boccarelli, C. Castiello, N. Del Buono, F. Esposito, and L. Selicato**.

This repository provides the R workflow underlying the chapter's worked examples. It introduces nonnegative matrix factorization (NMF) for exploring gene-expression data, extracting gene signatures (metagenes), and interpreting their activity across samples. The tutorial progresses from a controlled simulation to two publicly available RNA-seq datasets, connecting methodological concepts with practical analysis.

The workflow creates fitted NMF objects, gene-loading and component-activity matrices, lists of top-loading genes, reconstruction-error summaries, and diagnostic plots. The SEQC analysis also computes relative component exposures, Spearman correlation with known A fractions, and mean absolute error. Results remain in the R session; the script does not explicitly export tables or figures.

The examples are intended for teaching and methodological exploration. The eight-sample `airway` analysis provides an illustration of component interpretation. In the SEQC benchmark, normalized component exposures are assessed primarily for recovery of the expected **A > C > D > B** ordering; they need not equal the physical RNA mixing fractions. Top-loading genes are candidates for biological follow-up, rather than formal differential-expression results.


## Worked examples

| Example | Data | Main purpose |
| --- | --- | --- |
| Controlled simulation | 300 genes, 18 samples, and three latent signatures | Explore recovery of known structure and the effect of rank choice. |
| `airway` | Eight RNA-seq samples from four human airway smooth-muscle cell lines, with and without dexamethasone treatment | Prepare real expression data and interpret components against treatment and cell-line metadata. |
| SEQC/MAQC-III | Reference RNAs and known mixtures from the AGR Illumina centre, provided by the `seqc` package | Assess whether rank-two NMF captures the expected mixture ordering and compare linear and log-scale representations. |

The SEQC example uses two reference RNAs: A (Universal Human Reference RNA) and B (Human Brain Reference RNA), together with mixtures C (75% A + 25% B) and D (25% A + 75% B). The primary analysis uses row-scaled normalized counts; a second fit explores sensitivity to log transformation.

## Requirements

The tutorial uses R and the following packages:

- **CRAN:** `NMF`, `pheatmap`, `RColorBrewer`.
- **Bioconductor:** `DESeq2`, `airway`, `SummarizedExperiment`, `seqc`.

Install the dependencies in R before running the examples:

```r
install.packages(c("NMF", "pheatmap", "RColorBrewer"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "airway", "SummarizedExperiment", "seqc"))
```

The real datasets are loaded from the Bioconductor packages; no separate count-matrix download is needed.

## Running the tutorial

Download the R script from this repository and open it in R or RStudio. Run the sections in order, starting with the package-loading block and utility functions. Running the code interactively makes it easier to inspect each plot, table, and fitted object.

The script is organized into four sections:

1. Simulated data.
2. The `airway` dataset.
3. The SEQC/MAQC-III benchmark.
4. Computational environment information (`sessionInfo()`).

All fits use `method = "brunet"` and an explicit random seed (`20260904`). Rank surveys and repeated NMF fits may take time; their `nrun` settings are specified in the script.



## Reproducibility

The script fixes random seeds and reports the computational environment with `sessionInfo()`. Retain this output when reporting results, together with the chosen ranks, number of runs, and preprocessing settings. Package versions and computational environments can affect numerical results.

## Associated chapter

Boccarelli, A., Castiello, C., Del Buono, N., Esposito, F., and Selicato, L. *Nonnegative Matrix Factorization for RNA-seq Data Analysis: From Foundations to Benchmark Applications*.

Please acknowledge the associated chapter when using this tutorial in research or teaching, and cite the relevant software and original datasets where appropriate.
