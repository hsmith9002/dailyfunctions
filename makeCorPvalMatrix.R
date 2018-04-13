#makeCorPvalMatrix is a function that calculates correlations (pearson or spearman) and p-values. 
#It then takes these values and creates a matrix where the upper triangle
#contains pair-wise correlations, and the lower tirangle contains p-values.
#This function also allows for multiple testing corrections if the user
#deems this necessary.

#########################################################################################
#Arguments:
#dataset: a dataframe containing variables to be correlated. Must contain numbers only.
#type: the type of correlation to perform (pearson or spearman)
#adjust: the type of multiple testing complarison to apply ("holm", "hochberg", "hommel",
#"bonferroni", "BH", "BY", "fdr", if none use NULL. This is the default)
#n: number of samples (used for multiple testing) 
#nround: nuber of digits to round values included in final matrix
#alpha: signifigance level.
#############################################################################################

###########################################
#Function version of creating a correlation
#matrix that contains correlations on top
#and p-values on bottom
###########################################

makeCorPvalMatrix <- function(dataset, type = "pearson", adjust = NULL, n = NULL, nround = 3, alpha = 0.05) {
  library(Hmisc)
  library(psych)
  if (is.null(adjust)) {
 
  #calculate correlation matrix
  CorMatrix <- round(rcorr(as.matrix(dataset), type = type)$r, nround)
  CorMatrix.DF <- as.data.frame(CorMatrix)
  colnames(CorMatrix.DF) <- colnames(dataset)
  rownames(CorMatrix.DF) <- colnames(dataset)
  #calculate p-values
  pvalMatrix <- round(rcorr(as.matrix(dataset), type = type)$P, nround)
  pvalMatrix.DF <- as.data.frame(pvalMatrix)
  colnames(pvalMatrix.DF) <- colnames(dataset)
  rownames(pvalMatrix.DF) <- colnames(dataset)
  # create a new matrix that is the pvalues
  FinalCorMatrix <- pvalMatrix
  # make diaginal 1.00
  diag(FinalCorMatrix) <- 1.00
  # replace the upper triangle of this new matrix with the 
  # upper triangle of the correlation matrix
  FinalCorMatrix[upper.tri(FinalCorMatrix)] <- CorMatrix[upper.tri(CorMatrix)]
  FinalCorMatrix.DF <- as.data.frame(FinalCorMatrix)
  colnames(FinalCorMatrix.DF) <- colnames(dataset)
  rownames(FinalCorMatrix.DF) <- colnames(dataset)
  return(FinalCorMatrix.DF)
  }
  ## Perform pairwise correlations, and report p-values adjusted for multiple comparisons
  else {
  
  #calculate correlation matrix
  CorMatrix <- round(rcorr(as.matrix(dataset), type = type)$r, nround)
  CorMatrix.DF <- as.data.frame(CorMatrix)
  colnames(CorMatrix.DF) <- colnames(dataset)
  rownames(CorMatrix.DF) <- colnames(dataset)
  #calculate p-values
  pvalMatrix <- round(rcorr(as.matrix(dataset), type = type)$P, nround)
  pvalMatrix.DF <- as.data.frame(pvalMatrix)
  colnames(pvalMatrix.DF) <- colnames(dataset)
  rownames(pvalMatrix.DF) <- colnames(dataset)
  #Perform multiple testing corrections on p-values
  pvalMatrix.adjusted <- corr.p(CorMatrix, n, adjust = adjust, alpha = alpha)
  # create a new matrix that is the pvalues
  FinalCorMatrix <- pvalMatrix
  # make diaginal 1.00
  diag(FinalCorMatrix) <- 1.00
  # replace the bottom triangle of this new matrix with its own upper triangle
  FinalCorMatrix[lower.tri(FinalCorMatrix)] <- FinalCorMatrix[upper.tri(FinalCorMatrix)]
  # replace the upper triangle of this new matrix with the 
  # upper triangle of the correlation matrix
  FinalCorMatrix[upper.tri(FinalCorMatrix)] <- CorMatrix[upper.tri(CorMatrix)]
  FinalCorMatrix.DF <- as.data.frame(FinalCorMatrix)
  colnames(FinalCorMatrix.DF) <- colnames(dataset)
  rownames(FinalCorMatrix.DF) <- colnames(dataset)
  return(FinalCorMatrix.DF)
  }
}




