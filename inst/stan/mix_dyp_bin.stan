data {
  int<lower = 1> M; // number of years 
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower=1> Nt; // number of obs / block
  real<lower=0> y[M,Nt]; // all rainfall values
  int<lower=0> N[M];  // number of events - sample
    // prior parameters for hyper-distributions
    // each shape and scale (a, b) of a gamma for C and w
    real<lower = 0> C10prior[2];
    real<lower = 0> w10prior[2];
    real<lower = 0> C20prior[2];
    real<lower = 0> w20prior[2];
    real<lower = 0> srho0prior[2];
    real<lower = 0> mrho0prior[2];
    real<lower = 0> pn0prior[2];
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
parameters {
    real<lower=0> C2;
    real<lower=0> w2;
    real<lower=0> C1;
    real<lower=0> w1;
    // real<lower=0> DC;
    // real<lower=0> Dw;
  real<lower=0, upper = 1> pn;
  real<lower=0, upper = 1> rho[M];
  // real<lower=0, upper = 1> alphas;
    real<lower=0, upper=1> mrho;
  // real<lower=0, upper=1> srho;
  real<lower=0, upper=sqrt(mrho*(1-mrho))> srho;
}
transformed parameters{
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  
  alpha_rho = mrho^2*(1-mrho)/srho^2 - mrho;
  beta_rho = alpha_rho*(1-mrho)/mrho;
  
}
model {
  // Priors
  // alphap ~ beta(2, 2);
  rho ~ beta(alpha_rho, beta_rho);
  // mean_rho ~ beta(60, 60);
  // stdv_rho ~ beta(10, 60);
  mrho ~ beta(mrho0prior[1], mrho0prior[2]);
  srho ~ beta(srho0prior[1], srho0prior[2]);
  // stdv_rho ~ beta(10, 60);
  pn ~ beta(pn0prior[1], pn0prior[2]);
  // C1 ~ gamma(C0prior[1], 0.8*C0prior[2]);
  // w1 ~ gamma(w0prior[1], w0prior[2]);
  // C2 ~ gamma(C0prior[1], 1.2*C0prior[2]);
  // w2 ~ gamma(w0prior[1], 1.2*w0prior[2]);
  C1 ~ gamma(C10prior[1], C10prior[2]);
  w1 ~ gamma(w10prior[1], w10prior[2]);
  C2 ~ gamma(C20prior[1], C20prior[2]);
  w2 ~ gamma(w20prior[1], w20prior[2]);
  
  // DC ~ inv_gamma(10*1, 10);
  // Dw ~ inv_gamma(10*0.1, 10);
  // Model specification
    for (m in 1:M) {
    target += binomial_lpmf(N[m] | Nt, pn);
      for (j in 1:Nt) {
        if (y[m,j] >1e-6) {
                    target += log_mix(rho[m], 
                            weibull_lpdf(y[m,j] | w1, C1),
                            weibull_lpdf(y[m,j] | w2, C2));
      // target += weibull_lpdf(y[m,j] | w, C);
      }
    }
  }
}
generated quantities {
  // real log_lik[nwets]; // log like for the accumulations only (not for the N model)
  // real log_lik_N[M]; // log like for the N model only
  // real<lower=0> ygen[Mgen,Nt]; // all observed values
  // real<lower=0>maxgen[Mgen];
  int<lower=0> Ngen[Mgen]; // all observed values
  real<lower=0>Cgen1[Mgen]; // generated yearly parameters
  real<lower=0>wgen1[Mgen]; // generated yearly parameters
  real<lower=0>Cgen2[Mgen]; // generated yearly parameters
  real<lower=0>wgen2[Mgen]; // generated yearly parameters
  // real<lower=0, upper = 1>alphap[Mgen]; // generated yearly parameters
  // int <lower=0> count;
  
  
  // count = 0;
  // for (m in 1:M){
  //   log_lik_N[m] = binomial_lpmf(N[m] | Nt, pn);
  //   for (j in 1:Nt){
  //     if (y[m,j] > 1e-6){
  //       count = count + 1;
  //       log_lik[count] = weibull_lpdf(y[m,j] | w, C);
  //       }
  //     }
  //   }
    
  for (m in 1:Mgen) {
    Cgen1[m] = C1;
    wgen1[m] = w1;
    Cgen2[m] = C2;
    wgen2[m] = w2;
    Ngen[m] = binomial_rng(Nt, pn);
    // alphap[m] = alphas;
    // for (j in 1:Nt) {
    //   if (j <= Ngen[m]) {
    //      ygen[m,j] = weibull_rng( wgen[m], Cgen[m]);
    //   } else {
    //     ygen[m,j] = 0.0;
    //   }
    // }
    // maxgen[m] = max(ygen[m]);
  }

    
    
    
}












