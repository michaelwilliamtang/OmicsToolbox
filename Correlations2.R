#BiocManager::install("Hmisc")
correlations <- function(corr_mat) {
  require(Hmisc)
  # cor <- rcorr(format(corr_mat, digits=20), type="spearman")
  cor <- rcorr(format(corr_mat, digits=20), type="pearson")
  cor.data <- cor$r
  cor.data[upper.tri(cor.data, diag = T)] <- 0
  pval.data <- cor$P
  pval.data[upper.tri(pval.data, diag = T)]<- NA
  # FDR.data <- apply(pval.data, 2, p.adjust, method="BH", n = length(pval.data))
  FDR.data <- apply(pval.data, 2, p.adjust, method="BH")
  cor.data[FDR.data > 0.05]=0
  return(cor.data)
}
