  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of years
  real<lower=0> y[M,Nt]; // all rainfall values
  int<lower=0> N[M];  // number of events - sample
    // prior parameters for hyper-distributions (shape and rate)
    real<lower = 0> mc0prior[2]; // a, b of an inverse gamma
    real<lower = 0> sc0prior[2]; // a, b of an inverse gamma
    real<lower = 0> mw0prior[2]; // a, b of an inverse gamma
    real<lower = 0> sw0prior[2]; // a, b of an inverse gamma
    real<lower = 0> pn0prior[2];
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
parameters {
real<lower=0> C[M]; // vectorized
real<lower=0> w[M]; // vectorized
  real<lower=0> mc;
  real<lower=0> sc;
  real<lower=0> mw;
  real<lower=0> sw;
  real<lower=0, upper = 1> pn;
}
// transformed parameters{
//  real<lower=0> ac;
//  real<lower=0> bc;
//  real<lower=0> aw;
//  real<lower=0> bw;
// //  // For inverse gamma::
//  // ac = 2 + (mc/sc)^2;
//  // aw = 2 + (mw/sw)^2;
//  // bc = mc*(1 + (mc/sc)^2);
//  // bw = mw*(1 + (mw/sw)^2);
// //  // For Gamma::
//  // ac = (mc/sc)^2;
//  // aw = (mw/sw)^2;
//  // bc = mc/sc^2;
//  // bw = mw/sw^2;
// }
model {
  // Priors
  pn ~ beta(pn0prior[1], pn0prior[2]);
  mw ~ inv_gamma(mw0prior[1], mw0prior[2]);
  sw ~ inv_gamma(sw0prior[1], sw0prior[2]);
  mc ~ inv_gamma(mc0prior[1], mc0prior[2]);
  sc ~ inv_gamma(sc0prior[1], sc0prior[2]);
    // C ~ inv_gamma(ac, bc); // vectorized
    // w ~ inv_gamma(aw, bw); // vectorized
    C ~ gumbel(mc, sc); // vectorized
    w ~ gumbel(mw, sw); // vectorized
//model for the number of events/year
    for (m in 1:M) {
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
      genc =  gumbel_rng(mc, sc);
      if (genc > 0){
        count_c = count_c + 1;
        Cgen[count_c] = genc;
      }
    }


    count_w = 0;
    while (count_w < Mgen){
      genw =  gumbel_rng(mw, sw);
      // genw =  skew_normal_rng(wloc, wsc, wsh);
      if (genw > 0){
        count_w = count_w + 1;
        wgen[count_w] = genw;
      }
    }
    
for (m in 1:Mgen) {
    // Cgen[m] = gumbel_rng(mc, mc);
    // wgen[m] = gumbel_rng(mw, mw);
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





