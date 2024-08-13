
# R version of GUIDE model given an M x T matrix of summary statistics (`betas'), where M = number of genetic variants
# and T is the number of traits

# Dependencies:
  library(fastICA)
  library(Matrix)

guide <- function(betas, L = 100) {


  M <- max(dim(betas))
  
  betas_m <- scale(betas, center = TRUE, scale = FALSE)
  betas_m <- scale(t(betas_m), center = TRUE, scale = FALSE)
  
  svd_res <- svd(betas_m, nu = L, nv = L)
  Uc <- svd_res$u[, 1:L]
  Sc <- svd_res$d[1:L]
  Vc <- t(svd_res$v[, 1:L])
  UVc <- rbind(Uc, t(Vc)) / sqrt(2)
  
  ica_res <- fastICA(UVc, n.comp = L, maxit = 10000, tol = 1e-6)
  
  W_XL <- ica_res$S[1:M, ] *sqrt(2)
  W_LT <- t(ica_res$S[-(1:M), ]) * sqrt(2)
  mix <- ica_res$A
  
  return(list(W_XL = W_XL, W_LT = W_LT, Sc = Sc, mix = mix))
}


# The key differences between the Python and R code are:

# 1. In R, we use the `fastICA` package to perform the Independent Component Analysis (ICA) instead of the `sklearn.decomposition.FastICA` module used in the Python code.
# 2. The `scale()` function in R is used for mean-centering the `betas_m` matrix, instead of the manual subtraction of means used in the Python code.
# 3. The singular value decomposition (SVD) is performed using the `svd()` function in R, with the `nu` and `nv` parameters to specify the number of left and right singular vectors to compute.
# 4. The output of the function is a list containing the required components, instead of returning them as separate variables.

# The functionality of the two implementations should be the same, with the R version using the appropriate R packages and syntax.
