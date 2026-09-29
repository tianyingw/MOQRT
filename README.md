# MOQRT implementation and examples

R implementation of the minimax-optimal quantile rank-score test (MOQRT). 

## Dependencies

The example uses `quantreg` to fit the null quantile regressions and `MORST::Get_Q_pval()` to compute tail probabilities for the weighted chi-square distribution.

```r
install.packages("quantreg")
install.packages("remotes")  # needed only to install MORST from GitHub
remotes::install_github("yaowuliu/MORST")
```

Validation environment: R 4.5.3, quantreg 6.1, MORST 0.9.0.

## Single SNP

```r
source("MOQRT_main.R")
res_INT = MOQRT_main(X, Y, Z, tau_list,INT = TRUE)
res_INT$pval_MOQRT
```

- `X`: genotype dosages, one per subject.
- `Y`: responses in the same subject order as `X`.
- `Z`: covariate matrix, one row per subject. If it has no intercept column, one is added with a warning.
- `tau_list`: increasing quantile levels; default `c(0.1,0.25,0.5,0.75,0.9)`.
- `INT`: whether to apply residual INT before testing; default `TRUE`.
- `opt`: reserved for the two-parameter minimax search. Currently, `TRUE` gives a warning and still returns only the closed-form result.

The result contains `Sn` (rank scores), `Vn` (estimated score covariance, `Qn * V` in manuscript notation), `phi_dens_cut` (estimated density values), `orthogonal_score` (the squared null correlation of the two projected scores), and `pval_MOQRT`.

## Multiple SNPs with shared null fits

```r
Y_INT <- residual_INT(Y, Z)
rqfit0_list <- MOQRT_rqfit0(Y_INT, Z, tau_list, INT = FALSE)
pvals <- vapply(seq_len(ncol(X_all)), function(j) {
  res <- MOQRT_main(X_all[,j], Y_INT, Z, tau_list,
                    INT = FALSE, rqfit0_list = rqfit0_list)
  as.numeric(res$pval_MOQRT)
}, numeric(1))
```

For multiple SNPs, let `X_all` be a genotype matrix with subjects in rows and SNPs in columns, ordered to match `Y` and `Z`. We apply residual INT once and fit the null quantile regressions for the chosen quantile grid. The resulting `rqfit0_list` is reused for each SNP. 
