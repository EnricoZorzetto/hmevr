data {
  int<lower = 1> M; // number of years 
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower=1> Nt; // number of obs / block
  real<lower=0> y[M,Nt]; // all rainfall values
  int<lower=0> N[M];  // number of events - sample
  // prior parameters for hyper-distributions
  // each shape and scale (a, b) of a gamma for C and w
  real<lower = 0> C0prior[2];
  real<lower = 0> w0prior[2];
  // real<lower = 0> pn0[2];
  real<lower=0> mun0prior[2];
  real<lower=0> omegan0prior[2]; 
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
parameters {
    real<lower=0> C;
    real<lower=0> w;
  real<lower=0, upper = 1> pn[M];
  real<lower=0> mun;    // mean
  real<lower=0> omegan; // overdispersion
}
transformed parameters{
 real<lower=0> an; // beta shape1
 real<lower=0> bn; // beta shape2
 // real K;
  // K = (omega + 1)/(Nt-mu);
  // an = mu*(1-K)/(K*Nt-1);
  an = mun*(1- ((omegan + 1)/(Nt-mun)) )/( ((omegan + 1)/(Nt-mun)) *Nt-1);
  bn = an/mun*(Nt-mun); 
}
model {
    // Priors
      mun  ~ gamma(mun0prior[1], mun0prior[2]); 
    omegan ~ gamma(omegan0prior[1], omegan0prior[2]); 
    // Model specification
    pn ~ beta(an, bn);
    for (m in 1:M) {
    target += binomial_lpmf(N[m] | Nt, pn[m]);
       }
  // Priors
  C ~ gamma(C0prior[1], C0prior[2]);
  w ~ gamma(w0prior[1], w0prior[2]);
  // Model specification
    for (m in 1:M) {
      for (j in 1:Nt) {
        if (y[m,j] >1e-6) {
      target += weibull_lpdf(y[m,j] | w, C);
      }
    }
  }
}
generated quantities {
  // real log_lik[nwets]; // log like for the accumulations only (not for the N model)
  // real log_lik_N[M]; // log like for the N model only
  // real<lower=0> ygen[Mgen,Nt]; // all observed values
  // real<lower=0>maxgen[Mgen];
  real<lower=0> varn;
  int<lower=0> Ngen[Mgen]; // all observed values
  real<lower=0>Cgen[Mgen]; // generated yearly parameters
  real<lower=0>wgen[Mgen]; // generated yearly parameters
  // int <lower=0> count;
  
  varn = mun*(omegan + 1);
  
  // count = 0;
  // for (m in 1:M){
  //   log_lik_N[m] = beta_binomial_lpmf(N[m] | Nt, an, bn);
  //   for (j in 1:Nt){
  //     if (y[m,j] > 1e-6){
  //       count = count + 1;
  //       log_lik[count] = weibull_lpdf(y[m,j] | w, C);
  //       }
  //     }
  //   }
    
  for (m in 1:Mgen) {
    Cgen[m] = C;
    wgen[m] = w;
    Ngen[m] = beta_binomial_rng(Nt, an, bn);
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












