#BiocManager::install("Mfuzz")
standardise_data2 <- function(n) { # accepts numeric matrix
  require(Mfuzz)
  eset1 <- ExpressionSet(n)
  eset1 <- standardise(eset1) #Running standarise 
  o <- exprs(eset1)
  o <- na.omit(o)
  return(o)
}
