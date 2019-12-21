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
    real<lower = 0> pn0prior[2];
    real<lower = 0> rhoj0prior[2];
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
  // real<lower=0, upper = 1> alphap[M];
  real<lower=0, upper = 1> rhoj;
}
// transformed parameters{
//   real<lower=0> C2;
//   real<lower=0> w2;
//   C2 = C1 + DC;
//   w2 = w1 + Dw;
// }
model {
  // Priors
  // alphas ~ beta(60, 60);
  rhoj ~ beta(rhoj0prior[1], rhoj0prior[2]);
  pn ~ beta(pn0prior[1], pn0prior[2]);
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
                    target += log_mix(rhoj, 
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
  real<lower=0, upper = 1>rho[Mgen]; // generated yearly parameters
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
    rho[m] = rhoj;
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












