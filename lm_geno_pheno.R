library(tidyverse)

geno  <- read.csv("genotype_filtered_matrix.csv", row.names=1)  # samples x genes
pheno <- read.csv("phenotype_matrix.csv",         row.names=1)  # conditions x samples

geno  <- as.data.frame(geno)
pheno <- as.data.frame(t(pheno))  # samples x conditions

genes      <- colnames(geno)
conditions <- colnames(pheno)
results    <- list()

for (gene in genes) {
  for (cond in conditions) {
    af      <- geno[[gene]]
    fitness <- pheno[[cond]]
    
    # Simple linear regression: fitness ~ allele_frequency
    fit <- lm(fitness ~ af)
    s   <- summary(fit)
    
    # Extract coefficient for allele frequency
    if ("af" %in% rownames(s$coefficients)) {
      coef_row <- s$coefficients["af", ]
      results[[paste(gene, cond, sep = "_")]] <- data.frame(
        gene      = gene,
        condition = cond,
        estimate  = coef_row["Estimate"],       # effect size: change in log2FC per unit AF
        se        = coef_row["Std. Error"],
        t_value   = coef_row["t value"],
        p_raw     = coef_row["Pr(>|t|)"],
        r_squared = s$r.squared                 # variance explained by this gene in this condition
      )
    }
  }
}

results_df       <- bind_rows(results)
results_df$p_adj <- p.adjust(results_df$p_raw, method = "BH")
results_df       <- results_df %>% arrange(p_adj, p_raw)

# Filter meaningful hits
hits <- results_df %>%
  filter(
    p_adj        < 0.1,
    abs(estimate) > 0.2,
    r_squared    > 0.3   # gene explains at least 10% of fitness variance in this condition
  )

print(hits)
write.csv(results_df, "lm_all_results.csv",    row.names = FALSE)
write.csv(hits,       "lm_hits_filtered.csv",  row.names = FALSE)
