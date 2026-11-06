# Wilcoxon test for CNV effect on fitness

###### SET UP ######
# Install the Bioconductor manager if you don't have it yet
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install the missing S4Vectors package
BiocManager::install("S4Vectors")

install.packages("dplyr")
library(dplyr)

#Google drive location
setwd("/Users/siyuanmelaniewu/siyuanwu@umass.edu - Google Drive/My Drive/Research/Ma lab/Projects/Dry lab projects/STEE mouse/8996 STEE paper/Data & Analysis/Genotype_Phenotype Association/260913_geno-pheno-cnv")

#local location
setwd("/Users/siyuanmelaniewu/Desktop/STEE related/Genotype_phenotype analysis/260913_geno-pheno-cnv")

# ── Load data ──────────────────────────────────────────────────────────────
X <- read.csv("cnv_rf_input.csv", row.names = 1)
Y <- read.csv("phenotype_matrix.csv", row.names = 1)
Y <- t(Y)  # Transpose to samples x conditions

RF_imp <- read.csv("SHAP(new cnv calc)/shap_mean_abs.csv", row.names = 1)

cat("Genotype matrix:", nrow(X), "samples x", ncol(X), "genes\n")
cat("Phenotype matrix:", nrow(Y), "samples x", ncol(Y), "conditions\n")

# ── Parameters ─────────────────────────────────────────────────────────────
CNV_threshold <- 0.5  # 
min_mutants <- 2      # Minimum samples with CNV for test
fdr_cutoff <- 0.10    # FDR < 10% = significant

# ── Binarize genotype: strong CNV vs WT ────────────────────────────────────
# 1 = has strong CNV (|scaled_value| > 1.0), 0 = wild-type
X_binary <- abs(X) > CNV_threshold
X_binary <- X_binary * 1

cat("\nBinarization threshold: |scaled_value| >", CNV_threshold, "\n")

# ── Wilcoxon test function ─────────────────────────────────────────────────
run_wilcoxon_test <- function(fitness, has_cnv) {
  fitness_with_cnv <- fitness[has_cnv == 1]
  fitness_without_cnv <- fitness[has_cnv == 0]
  
  if (length(fitness_with_cnv) < 2 | length(fitness_without_cnv) < 2) {
    return(NA)
  }
  
  test <- wilcox.test(fitness_with_cnv, fitness_without_cnv)
  return(test$p.value)
}

# ── Filter by RF importance & run Wilcoxon tests ────────────────────────────
results_list <- list()

# ── Replace your Wilcoxon test loop section with this: ────────────────
for (condition in colnames(Y)) {
  
  rf_imp_values <- RF_imp[, condition]
  names(rf_imp_values) <- rownames(RF_imp)
  
  mean_imp <- mean(rf_imp_values, na.rm = TRUE)
  sd_imp <- sd(rf_imp_values, na.rm = TRUE)
  threshold_imp <- mean_imp + 0.5 * sd_imp
  
  important_genes <- names(rf_imp_values)[rf_imp_values > threshold_imp]
  
  for (gene in important_genes) {
    
    # Check if gene exists in original X (rows)
    if (!(gene %in% rownames(X))) {
      cat("Warning: gene", gene, "not found in X\n")
      next
    }
    
    # Get CNV values for this gene across all samples
    cnv_values <- X[gene, ]
    
    # Binarize: 1 if |CNV| > threshold, 0 otherwise
    has_cnv <- (abs(cnv_values) > CNV_threshold) * 1
    n_with_cnv <- sum(has_cnv)
    
    # Filter: minimum samples with CNV
    if (n_with_cnv >= min_mutants) {
      fitness <- Y[, condition]
      pval <- run_wilcoxon_test(fitness, has_cnv)
      
      if (!is.na(pval)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          gene = gene,
          condition = condition,
          n_with_cnv = n_with_cnv,
          n_without_cnv = sum(has_cnv == 0),
          p_value = pval
        )
      }
    }
  }
}









# # ── ADD THIS DIAGNOSTIC SECTION ────────────────────────────────────────
# for (condition in colnames(Y)) {
#   rf_imp_values <- RF_imp[, condition]
#   names(rf_imp_values) <- rownames(RF_imp)
#   
#   mean_imp <- mean(rf_imp_values, na.rm = TRUE)
#   sd_imp <- sd(rf_imp_values, na.rm = TRUE)
#   threshold_imp <- mean_imp + 0.5 * sd_imp
#   
#   important_genes <- names(rf_imp_values)[rf_imp_values > threshold_imp]
#   
#   cat("\n=== Condition:", condition, "===\n")
#   cat("RF importance - Mean:", mean_imp, "SD:", sd_imp, "\n")
#   cat("Threshold:", threshold_imp, "\n")
#   cat("Genes passing importance filter:", length(important_genes), "\n")
#   
#   # Check CNV distribution for these genes
#   for (gene in important_genes[1:min(3, length(important_genes))]) {  # Show first 3
#     has_cnv <- X_binary[, gene]
#     n_with <- sum(has_cnv)
#     cat("  Gene", gene, "- samples with CNV:", n_with, "\n")
#   }
# }
# 
# for (condition in colnames(Y)) {
#   
#   # Extract RF importance values and add gene names
#   rf_imp_values <- RF_imp[, condition]
#   names(rf_imp_values) <- rownames(RF_imp)
#   
#   mean_imp <- mean(rf_imp_values, na.rm = TRUE)
#   sd_imp <- sd(rf_imp_values, na.rm = TRUE)
#   threshold_imp <- mean_imp + 0.25 * sd_imp
#   
#   # Get genes above importance threshold
#   important_genes <- names(rf_imp_values)[rf_imp_values > threshold_imp]
#   
#   for (gene in important_genes) {
#     
#     # Check if gene exists in X_binary
#     if (!(gene %in% colnames(X_binary))) {
#       cat("Warning: gene", gene, "not found in X_binary\n")
#       next
#     }
#     
#     has_cnv <- X_binary[gene, ]  # ← Use [, ] not [[]]
#     n_with_cnv <- sum(has_cnv)
#     
#     # Filter: minimum 3 samples with CNV
#     if (n_with_cnv >= min_mutants) {
#       fitness <- Y[, condition]
#       pval <- run_wilcoxon_test(fitness, has_cnv)
#       
#       if (!is.na(pval)) {
#         results_list[[length(results_list) + 1]] <- data.frame(
#           gene = gene,
#           condition = condition,
#           n_with_cnv = n_with_cnv,
#           n_without_cnv = sum(has_cnv == 0),
#           p_value = pval
#         )
#       }
#     }
#   }
# }

# ── Combine results and apply FDR correction ────────────────────────────────
if (length(results_list) > 0) {
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # FDR correction (Benjamini-Hochberg)
  results_df$fdr <- p.adjust(results_df$p_value, method = "BH")
  
  # Sort by FDR
  results_sorted <- results_df[order(results_df$fdr), ]
  
  significant <- results_sorted$fdr < fdr_cutoff
  
  cat("\n=== RESULTS ===\n")
  cat("Pairs tested:", nrow(results_df), "\n")
  cat("Significant (FDR <", fdr_cutoff, "):", sum(significant), "\n\n")
  
  print(results_sorted)
  
  # Save all results
  write.csv(results_sorted, "wilcoxon_results.csv", row.names = FALSE)
  
  # Save significant results
  if (sum(significant) > 0) {
    write.csv(results_sorted[significant, ], "wilcoxon_significant.csv", row.names = FALSE)
    cat("\n=== SIGNIFICANT FINDINGS (FDR <", fdr_cutoff, ") ===\n")
    print(results_sorted[significant, ])
  } else {
    cat("\nNo significant findings at FDR <", fdr_cutoff, "\n")
  }
  
} else {
  cat("No pairs passed filtering!\n")
}

cat("\nSaved: wilcoxon_results.csv and wilcoxon_significant.csv\n")




###############################################################################

## making wilcoxon wide

w <- read.csv("wilcoxon_significant.csv", header = TRUE, sep = ",")
a <- read.csv("cnv_cutoff0.5_wide.csv", header = TRUE, sep = ",")

wa <- w %>%
  mutate(GENE = gene) %>%
  inner_join(a, by = "GENE")
write.csv(wa, "wilcoxon_sig_wide.csv", row.names = FALSE)
