  // Hierarchical model for extreme events
data {
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of years
  real<lower=0> y[M,Nt]; // all rainfall values
  int<lower=0> N[M];  // number of events - sample
    // prior parameters for hyper-distributions (shape and rate)
    real<lower = 0> mc10prior[2]; // a, b of an inverse gamma
    real<lower = 0> sc10prior[2]; // a, b of an inverse gamma
    real<lower = 0> mw10prior[2]; // a, b of an inverse gamma
    real<lower = 0> sw10prior[2]; // a, b of an inverse gamma
    real<lower = 0> mc20prior[2]; // a, b of an inverse gamma
    real<lower = 0> sc20prior[2]; // a, b of an inverse gamma
    real<lower = 0> mw20prior[2]; // a, b of an inverse gamma
    real<lower = 0> sw20prior[2]; // a, b of an inverse gamma
    
    real<lower = 0> rhoj0prior[2];
    real<lower = 0> pn0prior[2];
}
// transformed data{
//   int<lower=0> nwets;
//   nwets = sum(N);
// }
parameters {
real<lower=0> C2[M]; // vectorized
real<lower=0> w2[M]; // vectorized
  real<lower=0> mc2;
  real<lower=0> sc2;
  real<lower=0> mw2;
  real<lower=0> sw2;
  
  // real<lower=0> DC[M]; // vectorized
  // real<lower=0> Dw[M]; // vectorized
  // real<lower=0> Dmc;
  // real<lower=0> Dsc;
  // real<lower=0> Dmw;
  // real<lower=0> Dsw;
  
  real<lower=0> C1[M]; // vectorized
  real<lower=0> w1[M]; // vectorized
  real<lower=0> mc1;
  real<lower=0> sc1;
  real<lower=0> mw1;
  real<lower=0> sw1;
  real<lower=0, upper = 1> pn;
  real<lower=0, upper = 1> rhoj; // partition for the mixture
}
transformed parameters{
 real<lower=0> ac1;
 real<lower=0> bc1;
 real<lower=0> aw1;
 real<lower=0> bw1;
 real<lower=0> ac2;
 real<lower=0> bc2;
 real<lower=0> aw2;
 real<lower=0> bw2;
 // real Dac;
 // real Dbc;
 // real Daw;
 // real Dbw;
 // real<lower=0> C2[M];
 // real<lower=0> w2[M];

    // for (m in 1:M) {
    //    // C2[m] = (C1[m]+DC[m]); // difference > 0: C2 > C1
    //    w2[m] = (w1[m]+Dw[m]); // difference > 0: w2 > w1
    // }
    // 
 ac1 = (mc1/sc1)^2;
 aw1 = (mw1/sw1)^2;
 bc1 =  mc1/sc1^2;
 bw1 =  mw1/sw1^2;
 
 // Dac = (Dmc/Dsc)^2;
 // Daw = (Dmw/Dsw)^2;
 // Dbc =  Dmc/Dsc^2;
 // Dbw =  Dmw/Dsw^2;
 
   ac2 = (mc2/sc2)^2;
   bc2 =  mc2/sc2^2;
   aw2 = (mw2/sw2)^2;
   bw2 =  mw2/sw2^2;
}
model {
  // Priors
  // alphaps ~ beta(60, 60);
  rhoj ~ beta(rhoj0prior[1], rhoj0prior[2]);
  pn ~ beta(pn0prior[1], pn0prior[2]);
  
  // mw1 ~ inv_gamma( mw0prior[1], 0.9*mw0prior[2]);
  // sw1 ~ inv_gamma( sw0prior[1], sw0prior[2]);
  // mc1 ~ inv_gamma( mc0prior[1], 0.8*mc0prior[2]);
  // sc1 ~ inv_gamma( sc0prior[1], sc0prior[2]);
  // 
  // mw2 ~ inv_gamma( mw0prior[1], 1.2*mw0prior[2]);
  // sw2 ~ inv_gamma( sw0prior[1], 1*sw0prior[2]);
  // mc2 ~ inv_gamma( mc0prior[1], 1.2*mc0prior[2]);
  // sc2 ~ inv_gamma( sc0prior[1], sc0prior[2]);
  
  mw1 ~ inv_gamma( mw10prior[1], mw10prior[2]);
  sw1 ~ inv_gamma( sw10prior[1], sw10prior[2]);
  mc1 ~ inv_gamma( mc10prior[1], mc10prior[2]);
  sc1 ~ inv_gamma( sc10prior[1], sc10prior[2]);
  
  mw2 ~ inv_gamma( mw20prior[1], mw20prior[2]);
  sw2 ~ inv_gamma( sw20prior[1], sw20prior[2]);
  mc2 ~ inv_gamma( mc20prior[1], mc20prior[2]);
  sc2 ~ inv_gamma( sc20prior[1], sc20prior[2]);
  
  // Dmw ~ normal(0, 0.2);
  // Dsw ~ normal(0, 0.1);
  // Dmc ~ normal(0, 2);
  // Dsc ~ normal(0, 0.1);
  
  // Dmw ~ inv_gamma(0.01*0.1, 0.01);
  // Dsw ~ inv_gamma(0.01*0.1, 0.01);
  // Dmc ~ inv_gamma(10*2  , 10);
  // Dsc ~ inv_gamma(10*0.1, 10);
  
  
    // DC ~ gamma(Dac, Dbc); // vectorized
    // Dw ~ gamma(Daw, Dbw); // vectorized
    C2 ~ gamma(ac2, bc2); // vectorized
    w2 ~ gamma(aw2, bw2); // vectorized
    C1 ~ gamma(ac1, bc1); // vectorized
    w1 ~ gamma(aw1, bw1); // vectorized
//model for the number of events/year
    for (m in 1:M) {
    target += binomial_lpmf(N[m] | Nt, pn);
      for (j in 1:Nt) {
        if (y[m,j] >1e-6) {
          target += log_mix(rhoj, 
                            weibull_lpdf(y[m,j] | w1[m], C1[m]),
                            weibull_lpdf(y[m,j] | w2[m], C2[m]));
           // y[m,j] ~ weibull( w[m], C[m]);
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
  real<lower=0>Cgen2[Mgen]; // generated yearly parameters
  real<lower=0>Cgen1[Mgen]; // generated yearly parameters
  real<lower=0>wgen2[Mgen]; // generated yearly parameters
  real<lower=0>wgen1[Mgen]; // generated yearly parameters
  real<lower=0, upper=1>rho[Mgen]; // generated yearly parameters
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
    Cgen1[m] = gamma_rng(ac1, bc1);
    wgen1[m] = gamma_rng(aw1, bw1);
    Cgen2[m] = gamma_rng(ac2, bc2);
    wgen2[m] = gamma_rng(aw2, bw2);
    // Cgen2[m] = Cgen1[m] + gamma_rng(Dac, Dbc);
    // wgen2[m] = wgen1[m] + gamma_rng(Daw, Dbw);
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





