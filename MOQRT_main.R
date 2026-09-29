library(quantreg)
library(MORST)

add_intercept <- function(Z){
  Z = as.matrix(Z)
  if(!any(apply(Z, 2, function(z) all(z == 1)))){
    warning("Z has no intercept column; adding a column of ones.")
    Z = cbind(1, Z)
  }
  return(Z)
}

get_tau_mat <- function(tau_list){
  L = length(tau_list)
  tau_mat = matrix(0,L,L)
  for(t in 1:L){
    tau = tau_list[t]
    tau_mat[t,t] = tau - tau^2
  }
  for(t1 in 1:(L-1)){
    for(t2 in (t1+1):(L)){
      tau1 = tau_list[t1]
      tau2 = tau_list[t2]
      tau_mat[t1,t2] =tau_mat[t2,t1] = min(tau1,tau2) - tau1* tau2
    }
  }
  return(tau_mat)
}

Get_chisq1 <- function(q,lambda){
  if( 0 %in% lambda){
    lambda = lambda[-which(lambda==0)]
  }
  pval = MORST::Get_Q_pval(Q=q,w = lambda)
  return(pval)
}

theta_minmax <- function(V,D_phi,tau_list= c(0.1,0.25,0.5,0.75,0.9)){
  L = length(tau_list)
  d1 = rep(1,L) / sqrt(sum(rep(1,L)^2))
  d2 = qnorm(tau_list)
  d2 = d2 / sqrt(sum(d2^2))
  
  d1 = D_phi %*% d1
  d2 = D_phi %*% d2
  
  lambda1 = c(t(d1) %*%solve(V)%*%solve(V) %*% d1) / c(t(d1) %*%solve(V) %*% d1)
  lambda2 = c(t(d2) %*%solve(V)%*%solve(V) %*% d2) / c(t(d2) %*%solve(V) %*% d2)
  return(  lambda1 / (lambda1 + lambda2) )
}

center.V.theta<- function(theta,V,D_phi,tau_list = c(0.1,0.25,0.5,0.75,0.9)) {
  L = length(tau_list)
  
  d1 = rep(1,L) / sqrt(sum(rep(1,L)^2))
  d2 = qnorm(tau_list)
  d2 = d2 / sqrt(sum(d2^2))
  
  d1 = D_phi %*% d1
  d2 = D_phi %*% d2
  
  u1 = solve(V) %*% d1
  u2 = solve(V) %*% d2
  u1 = u1 / norm(u1,'2')
  u2 = u2 / norm(u2,'2')
  
  mat = theta*u1 %*% t(u1) + (1-theta)* u2 %*% t(u2)
  return(mat)
}

MOQRT_sym <- function(Sk,Vk,thetac = NULL,D_phi_est,tau_list){
  if(is.null(thetac)){
    thetac = theta_minmax(V = Vk,D_phi = D_phi_est,tau_list)
  }
  eigen_Vk = eigen(Vk)
  Vk_sqrt = eigen_Vk$vectors %*% diag(sqrt(eigen_Vk$values)) %*% t(eigen_Vk$vectors)
  
  L = length(tau_list)
  
  center_V = center.V.theta(thetac,Vk,D_phi_est,tau_list)
  Test_T = t(Sk) %*% center_V %*% Sk
  center_eigenvalues = eigen(Vk_sqrt %*% center_V %*% Vk_sqrt)$values
  if(max(Im(center_eigenvalues))>1e-12){
    warning('The imaginary part of the eigenvalue > 1e-12')
  }
  if(min(Re(center_eigenvalues))< -1e-12){
    warning('The minimum eigenvalue < -1e-12')
  }
  center_eigenvalues[which(abs(center_eigenvalues/center_eigenvalues[1])<1e-12)] = 0
  center_eigenvalues = Re(center_eigenvalues)
  if(max(center_eigenvalues)<1e-12){
    warning('max(center_eigenvalues)<1e-12')
    pval = NA
  }else{
    pval  = Get_chisq1(Test_T,center_eigenvalues)
  }
  return(pval)
}



# Squared null correlation of the two projected score directions.
Get.orthogonal_score <- function(D_phi_est,tau_list){
  L = length(tau_list)
  V_inv = solve(get_tau_mat(tau_list))
  d1 = D_phi_est %*% rep(1,L)
  d2 = D_phi_est %*% qnorm(tau_list)
  Gamma = c(t(d1) %*% V_inv %*% d2)
  r1 = c(t(d1) %*% V_inv %*% d1)
  r2 = c(t(d2) %*% V_inv %*% d2)
  orthogonal_score = Gamma^2/(r1*r2)
}

# Apply residual INT after adjusting Y for Z.
# Y: responses; Z: covariate matrix (an intercept is added if absent).
# Returns: inverse-normal-transformed residuals.
residual_INT<- function(Y,Z){
  Z = add_intercept(Z)
  Y_ori = Y
  fit_LR =lm(Y_ori~Z-1)
  Y_resid = fit_LR$residuals
  percentile = ecdf(Y_resid)
  Y_resid_irnt = qnorm(percentile(Y_resid)-min(percentile(Y_resid))/2)
  return(Y_resid_irnt)
}

# Compute the MOQRT test for one SNP.
# X: genotype dosages; Y: responses; Z: covariate matrix.
# tau_list: quantile levels; INT: apply residual INT internally.
# rqfit0_list: optional matching null fits from MOQRT_rqfit0().
# opt: reserved for the two-parameter minimax search.
# Returns: Sn, Vn, phi_dens_cut, orthogonal_score, and pval_MOQRT.
MOQRT_main <- function(X, Y, Z, tau_list = c(0.1,0.25,0.5,0.75,0.9),
                       INT = TRUE, opt = FALSE, rqfit0_list = NULL){
  Z = add_intercept(Z)
  if(INT == TRUE){
    Y = residual_INT(Y,Z) 
  }
  Y = as.matrix(Y)
  Z = as.matrix(Z)
  X = as.matrix(X)
  L = length(tau_list)
  
  N_sample = length(X)
  xstar = as.matrix(qr.resid(qr(Z), X))
  Sn = c()
  # Reuse MOQRT_rqfit0() when testing multiple SNPs against the same Y and Z.
  if(is.null(rqfit0_list)){
    phi_dens_cut = c()
    for(t_tau in 1:L){
      tau = tau_list[t_tau]
      rqfit0 = quantreg::rq.fit.br(Z,Y,tau)
      an0 = rqfit0$dual
      # Centred dual rank score at this quantile.
      Sn[t_tau] = 1/sqrt(N_sample)*sum(xstar*(an0- (1-tau)))
      
      # error density
      uhat0 = rqfit0$residuals
      d1 = density(uhat0, from=quantile(uhat0,0.001), to=quantile(uhat0,0.999),
                   cut=0.01, kernel="gaussian",n = 2^10)
      min_id = which(abs(d1$x)==min(abs(d1$x)))
      if(length(min_id)>1){min_id= min_id[1]}
      phi_dens_cut[t_tau] = d1$y[min_id]
    }
  }else{
    for(t_tau in 1:L){
      tau = tau_list[t_tau]
      an0 = rqfit0_list$an0_list[[t_tau]]
      Sn[t_tau] = 1/sqrt(N_sample)*sum(xstar*(an0- (1-tau)))
    }
    phi_dens_cut = rqfit0_list$phi_dens_cut
  }

  Qn = 1/N_sample*sum(xstar^2)
  Vn = Qn * get_tau_mat(tau_list)
  D_phi_est = diag(phi_dens_cut)
  orthogonal_score = Get.orthogonal_score(D_phi_est,tau_list)
  pval_MOQRT = MOQRT_sym(Sn,Vn,thetac = NULL,D_phi_est,tau_list)

  res = list(Sn = Sn, Vn = Vn, phi_dens_cut = phi_dens_cut,
             orthogonal_score = orthogonal_score,
             pval_MOQRT = pval_MOQRT)
  if(opt == TRUE){
    warning("The two-parameter minimax search is computationally intensive; returning only the closed-form MOQRT result.")
  }
  return(res)
}

# Fit null quantile regressions once for reuse across SNPs.
# Y, Z, tau_list, and INT: as in MOQRT_main().
# Returns: an0_list (null dual scores) and phi_dens_cut (densities) in a list.
MOQRT_rqfit0 <- function(Y, Z, tau_list = c(0.1,0.25,0.5,0.75,0.9),INT = TRUE){
  Z = add_intercept(Z)
  if(INT == TRUE){
    Y = residual_INT(Y,Z) 
  }
  Y = as.matrix(Y)
  Z = as.matrix(Z)
  
  L = length(tau_list)
  rqfit0_list= list()
  
  phi_dens_cut = c()
  for(t_tau in 1:L){
    tau = tau_list[t_tau]
    rqfit0 = quantreg::rq.fit.br(Z,Y,tau)
    rqfit0_list$an0_list[[t_tau]] = rqfit0$dual
    
    # density
    uhat0 = rqfit0$residuals
    d1 = density(uhat0, from=quantile(uhat0,0.001), to=quantile(uhat0,0.999),
                 cut=0.01, kernel="gaussian",n = 2^10)
    min_id = which(abs(d1$x)==min(abs(d1$x)))
    if(length(min_id)>1){min_id= min_id[1]}
    phi_dens_cut[t_tau] = d1$y[min_id]
  }
  rqfit0_list$phi_dens_cut = phi_dens_cut
  return(rqfit0_list)
}
