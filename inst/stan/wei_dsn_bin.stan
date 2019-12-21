// The hbevm model - weibull dynamic lognormal
// functions {
//   real skew_normal_lb_rng(real loc, real sc, real sh, real lb) {
//     real p = skew_normal_cdf(lb, loc, sc, sh);  // cdf for bounds
//     real u = uniform_rng(p, 1);
//     return sigma * (-log1m(u))^(1 / alpha);  // inverse cdf for value
//   }
//   // real fr_lpdf(real x, real loc, real scale, real shape) {
//   // real frechet3_lpdf(real x, real m, real s, real a) {
//   //   // real lpdf;
//   //   real norm;
//   //   norm = (x-m)*inv(s);
//   //   return log(a)-log(s)-(1+a)*log(norm)-exp(-a*log(norm));
//   // }
//     real frechet3_lpdf(real x, real m, real s, real a) {
//     real xmm;
//     xmm = x-m;
//     return frechet_lpdf(xmm | a, s);
//   }
//     real exploc_lpdf(real x, real m, real s) {
//     real xmm;
//     xmm = x-m;
//     // return frechet_lpdf(xmm | a, s);
//     return exponential_lpdf(xmm | inv(s));
//   }
// }
data {
  int<lower = 1> M; // number of years 
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of observations /year
  real<lower=0> y[M,Nt]; // all observed values
  int<lower=0> N[M];  // number of events - sample
  // real<lower=0> sigmac0prior[2]; 
  // real<lower=0> sigmaw0prior[2];
  // real muc0prior[2]; 
  // real muw0prior[2]; 
  // real<lower=0> sc0prior[2]; 
  // real<lower=0> sw0prior[2];
  // real mc0prior[2]; 
  // real mw0prior[2]; 
  // real<lower = 0> pn0prior[2];
  real<lower=0> csnloc0prior[2];
  real<lower=0> wsnloc0prior[2];
  real<lower=0> csnsc0prior[2];
  real<lower=0> wsnsc0prior[2];
  real<lower=0> csnsh0prior[2];
  real<lower=0> wsnsh0prior[2];
  real<lower=0> pn0prior[2];
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
parameters {

real<lower=0> C[M]; 
real<lower=0> w[M]; 
 real <lower=0> cloc; // average of c 
 real <lower=0> wloc; 
 real <lower=0> csc; // standard deviation of c
 real <lower=0> wsc; 
 real <lower=0> csh;
 real <lower=0> wsh;
  // real muc;
  // real muw;
  // real<lower = 0> sigmac;
  // real<lower = 0> sigmaw;
  real<lower=0, upper = 1> pn;
}
// transformed parameters{
// 
//  // real <lower=0> mc; // average of c 
//  // real <lower=0> mw; 
//  // real <lower=0> sc; // standard deviation of c
//  // real <lower=0> sw; 
//   //  real muc;
//   // real muw;
//   // real<lower = 0> sigmac;
//   // real<lower = 0> sigmaw;
//   // muc = log(mc)-0.5*log1p( pow(sc/mc, 2));
//   // muw = log(mw)-0.5*log1p( pow(sw/mw, 2));
//   // sigmac = sqrt(log1p(pow(sc/mc, 2)));
//   // sigmaw = sqrt(log1p(pow(sw/mw, 2)));
// 
//  C = C0 + cloc;
//  w = w0 + wloc;
// }
model {
  // Priors
  pn ~ beta(pn0prior[1], pn0prior[2]);

  
  cloc ~ gamma(csnloc0prior[1], csnloc0prior[2]);
  wloc ~ gamma(wsnloc0prior[1], wsnloc0prior[2]);
  csc ~ gamma(csnsc0prior[1], csnsc0prior[2]);
  wsc ~ gamma(wsnsc0prior[1], wsnsc0prior[2]);
  csh ~ gamma(csnsh0prior[1], csnsh0prior[2]);
  wsh ~ gamma(wsnsh0prior[1], wsnsh0prior[2]);
  
  // cloc ~ gamma(10*10*0.9, 10);
  // wloc ~ gamma(10*0.7*0.9, 10); 
  // csc ~ gamma(10*0.2*200, 200);
  // wsc ~ gamma(0.7*0.1*200, 200);
  // csh ~gamma(10*200, 200);
  // wsh ~gamma(10*200, 200);
  
  for (m in 1:M) {
      // C[m] ~ frechet3(cloc, csc, csh); // frechet shape and scale (loc = 0)
      // w[m] ~ frechet3(wloc, wsc, wsh); // frechet shape and scale (loc = 0)
      // C[m] ~ skew_normal(cloc, csc, csh)T[0.01, positive_infinity()]; 
      // w[m] ~ skew_normal(wloc, wsc, wsh)T[0.01, positive_infinity()]; 
      C[m] ~ skew_normal(cloc, csc, csh); 
      w[m] ~ skew_normal(wloc, wsc, wsh); 
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
      int count_w;
    real genw;
    int count_c;
    real genc;

    count_c = 0;
    while (count_c < Mgen) {
      genc =  skew_normal_rng(cloc, csc, csh);
      if (genc > 0){
        count_c = count_c + 1;
        Cgen[count_c] = genc;
      }
    }


    count_w = 0;
    while (count_w < Mgen){
      genw =  skew_normal_rng(wloc, wsc, wsh);
      if (genw > 0){
        count_w = count_w + 1;
        wgen[count_w] = genw;
      }
    }
    
for (m in 1:Mgen) {
    // Cgen[m] = lognormal_rng(muc, sigmac);
    // wgen[m] = lognormal_rng(muw, sigmaw);
    // Cgen[m] = cloc + frechet_rng(csh, csc);
    // wgen[m] = wloc + frechet_rng(wsh, wsc);
    // Cgen[m] = skew_normal_rng(cloc, csc, csh);
    // wgen[m] = skew_normal_rng(wloc, wsc, wsh);
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
