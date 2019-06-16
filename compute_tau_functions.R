compute_tau_mean<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
  N<- N_eps*gamma
  return(N)
}


compute_tau_mean_large_dim<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  log_val<- lgamma(1+theta) + lgamma(1+ theta/alpha + 1/(1-alpha) ) - lgamma(1+theta/alpha) - lgamma(1+ alpha/(1-alpha) + theta)
  #gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
  gamma_val<- exp(log_val)
  N<- N_eps*gamma_val
  return(N)
}

##############  compute variance
compute_tau_var<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma2<- (gamma(1+theta)*gamma(1+ theta/alpha + 2/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (2*alpha)/(1-alpha) + theta))
  gamma1<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (1*alpha)/(1-alpha) + theta))
  gamma<-(gamma2-gamma1*gamma1)
  N<- (N_eps^2)*gamma
  return(sqrt(N))
}


compute_tau_var_large_dim<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  log_gamma_2<- lgamma(1+theta) + lgamma(1+ theta/alpha + 2/(1-alpha)) - lgamma(1+theta/alpha) - lgamma(1+ (2*alpha)/(1-alpha) + theta)
  gammaval_2<- exp(log_gamma_2)
  # gamma2<- (gamma(1+theta)*gamma(1+ theta/alpha + 2/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (2*alpha)/(1-alpha) + theta))
  log_gamma_1<- lgamma(1+theta) + lgamma(1+ theta/alpha + 1/(1-alpha)) - lgamma(1+theta/alpha) - lgamma(1+ (1*alpha)/(1-alpha) + theta)
  # gamma1<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (1*alpha)/(1-alpha) + theta))
  gammaval_1<- exp(log_gamma_1)
  gamma_val<-(gammaval_2 -gammaval_1*gammaval_1)
  N<- (N_eps^2)*gamma_val
  return(sqrt(N))
}

