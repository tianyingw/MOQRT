source("MOQRT_main.R")

#### Test for X,Y,Z ####
# generate X,Y,Z
set.seed(123)
N_sample = 5000
sample_beta = 0.1
sample_gamma = 0.05
MAF = 0.2
tau_list = c(0.1,0.25,0.5,0.75,0.9)

X = rnorm(N_sample)
X = as.numeric(cut(X,breaks = c(-Inf,qnorm((1-MAF)^2),qnorm((1-MAF)^2+2*MAF*(1-MAF)),Inf),
                   labels = c(0,1,2)))-1
X = matrix(X,N_sample,1)
Z = cbind(1,rnorm(N_sample))

epsilon = rcauchy(N_sample)
Y = Z %*% c(1,1)+sample_beta*X+ (1+sample_gamma*X)*epsilon

# run MOQRT
res_INT = MOQRT_main(X, Y, Z, tau_list,INT = TRUE)
cat("Single-SNP MOQRT p-value:\n", signif(res_INT$pval_MOQRT))

#### Test for a group of X and (Y,Z) ####
set.seed(123)
N_sample = 5000
MAF = 0.2
tau_list = c(0.1,0.25,0.5,0.75,0.9)

X_all = sapply(1:100, function(i){
  X = rnorm(N_sample)
  X = as.numeric(cut(X,breaks = c(-Inf,qnorm((1-MAF)^2),qnorm((1-MAF)^2+2*MAF*(1-MAF)),Inf),
                     labels = c(0,1,2)))-1
  X = matrix(X,N_sample,1)
})

Z = cbind(1,rnorm(N_sample))
 
epsilon = rcauchy(N_sample)
Y = Z %*% c(1,1) + epsilon

# INT 
Y_resid_irnt = residual_INT(Y,Z)

# Fit the null quantile regressions once for all SNPs.
rqfit0_list = MOQRT_rqfit0(Y_resid_irnt, Z, tau_list,INT = FALSE)

# Test all 100 SNPs individually; retain one result per SNP.
results = data.frame(SNP = seq_len(ncol(X_all)), pval_MOQRT = NA_real_)
for(i in seq_len(ncol(X_all))){
  X = as.matrix(X_all[,i])

  res_INT = MOQRT_main(X, Y_resid_irnt, Z, tau_list = tau_list,
             INT = FALSE,
             rqfit0_list = rqfit0_list)
  results$pval_MOQRT[i] = as.numeric(res_INT$pval_MOQRT)
}
cat("\nMOQRT p-values for 100 simulated SNPs under the null:\n")
head(results, row.names = FALSE)