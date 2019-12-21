// The hbevm model - weibull dynamic lognormal
functions {
  // real fr_lpdf(real x, real loc, real scale, real shape) {
  // real frechet3_lpdf(real x, real m, real s, real a) {
  //   // real lpdf;
  //   real norm;
  //   norm = (x-m)*inv(s);
  //   return log(a)-log(s)-(1+a)*log(norm)-exp(-a*log(norm));
  // }
  //   real frechet3_lpdf(real x, real m, real s, real a) {
  //   real xmm;
  //   xmm = x-m;
  //   return frechet_lpdf(xmm | a, s);
  // }
    real exploc_lpdf(real x, real m, real s) {
    real xmm;
    xmm = x-m;
    // return frechet_lpdf(xmm | a, s);
    return exponential_lpdf(xmm | inv(s));
  }
}
data {
  int<lower = 1> M; // number of years 
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of observations /year
  real<lower=0> y[M,Nt]; // all observed values
  int<lower=0> N[M];  // number of events - sample
  real<lower=0> cexploc0prior[2];
  real<lower=0> wexploc0prior[2];
  real<lower=0> cexpsc0prior[2];
  real<lower=0> wexpsc0prior[2];
  real<lower=0> pn0prior[2];
}
parameters {
  real <lower=0> cloc; // average of c 
  real <lower=0> wloc; 
  real<lower=cloc> C[M]; 
  real<lower=wloc> w[M]; 
  real <lower=0> csc; // standard deviation of c
  real <lower=0> wsc; 
  real<lower=0, upper = 1> pn;
}
model {
  // Priors
  pn ~ beta(pn0prior[1], pn0prior[2]);

  // cloc ~ gamma(10*10, 10);
  // wloc ~ gamma(0.7*10, 10);
  // csc ~ gamma(10*0.2*200, 200);
  // wsc ~ gamma(0.7*0.1*200, 200);
  
  cloc ~ gamma(cexploc0prior[1], cexploc0prior[2]);
  wloc ~ gamma(wexploc0prior[1], wexploc0prior[2]);
  csc ~ gamma(cexpsc0prior[1], cexpsc0prior[2]);
  wsc ~ gamma(wexpsc0prior[1], wexpsc0prior[2]);

  
  for (m in 1:M) {
      C[m] ~ exploc(cloc, csc); // frechet shape and scale (loc = 0)
      w[m] ~ exploc(wloc, wsc); // frechet shape and scale (loc = 0)
  target += binomial_lpmf(N[m] | Nt, pn);
    for (j in 1:Nt) {
      if (y[m,j] >1e-6) { 
         y[m,j] ~ weibull( w[m], C[m]);
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
  real<lower=0>Cgen[Mgen]; // generated yearly parameters
  real<lower=0>wgen[Mgen]; // generated yearly parameters
  // int <lower=0> count;
  
  // count = 0;
  // for (m in 1:M){
  //   log_lik_N[m] = binomial_lpmf(N[m] | Nt, pn);
  //   for (j in 1:Nt){
  //     if (y[m,j] > 1e-6){
  //       count = count + 1;
  //       log_lik[count] = weibull_lpdf(y[m,j] | w[m], C[m]);
  //       }
  //     }
  //   }
    
for (m in 1:Mgen) {
    // Cgen[m] = lognormal_rng(muc, sigmac);
    // wgen[m] = lognormal_rng(muw, sigmaw);
    // Cgen[m] = cloc + frechet_rng(csh, csc);
    // wgen[m] = wloc + frechet_rng(wsh, wsc);
    Cgen[m] = cloc + exponential_rng(inv(csc));
    wgen[m] = wloc + exponential_rng(inv(wsc));
    Ngen[m] = binomial_rng(Nt, pn);
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
