// old version of the GEV stan code
// without bounds for xi
functions {
    real gev_lpdf(vector y, real mu, real k, real psi) {
    // gev log pdf 
    // int N = rows(y);
    real N = rows(y);
    real inv_k = inv(k);
    if (k<0 && max(y-mu)/psi > -inv_k)
      reject("k<0 and max(y-mu)/psi > -1/k; found k, psi =", k, psi)
    if (psi<=0)
      reject("psi<=0; found psi =", psi)
    if (fabs(k) > 1e-15)
      // return -(1+inv_k)*sum(log1p((y-mu) * (k/sigma))) -N*log(sigma);
      return (-N*log(psi) -(1+inv_k)*sum(log1p( (y-mu)*(k/psi) ))
               -sum( exp( -inv_k*log1p((y-mu)*(k/psi)) )));
    else
      // return -sum(y-mu)/sigma -N*log(sigma); // limit k->0
      return -N*log(psi) - sum( (y - mu)/psi + exp(-(y - mu)/psi) );
  }
    real gev_lcdf(vector y, real mu, real k, real psi) {
    // gev log cdf
    real inv_k = inv(k);
    if (k<0 && max(y-mu)/psi > -inv_k)
      reject("k<0 and max(y-mu)/psi > -1/k; found k, psi =", k, psi)
    if (psi<=0)
      reject("psi<=0; found psi =", psi)
    if (fabs(k) > 1e-15)
      return sum(  -exp( -inv_k*log( (1+(y-mu)*(k/psi)) ) ));
    else
      return  exp( sum( log( -exp(-(y-mu)/psi)   )));
  }
    real gev_lccdf(vector y, real mu, real k, real psi) {
    // gev log ccdf
    real inv_k = inv(k);
    if (k<0 && max(y-mu)/psi > -inv_k)
      reject("k<0 and max(y-mu)/psi > -1/k; found k, psi =", k, psi)
    if (psi<=0)
      reject("psi<=0; found psi =", psi)
    if (fabs(k) > 1e-15)
      return sum( log1m_exp( -exp( -inv_k*log( (1+(y-mu)*(k/psi)) ) )));
    else
      return sum(log1m_exp(- exp( -(y-mu)/psi))); // limit k->0
  }
    real gev_rng(real mu, real k, real psi) {
    if (psi<=0)
      reject("[psi<=0; found psi =", psi)
    if (fabs(k) > 1e-15)
      return mu +  psi/k*(  -1 + exp( -k*log(-log( uniform_rng(0, 1)))  ) );
    else
      return mu - psi*log( -log(uniform_rng(0,1))); // limit k->0
  }
}
data {
  int<lower=0> N;
  int<lower=0> Mgen;
  vector[N] y;
  real pr_psi[2]; // Normal parameters mu, sigma
  real pr_mu[2];    // Normal parameters mu, sigma
  real pr_k[2];   // Beta  parameters a, b
}
parameters {
    real k; 
    real<lower=0> psi;
    real mu; 
}
model {
  mu ~ normal(pr_mu[1], pr_mu[2]);
  psi ~ lognormal(pr_psi[1], pr_psi[2]);
  k ~ normal(pr_k[1], pr_k[2]);
  y ~ gev(mu, k, psi);
}
generated quantities {
  real yrep[Mgen];
  real log_lik[N];
  for (i in 1:Mgen){
    yrep[i] = gev_rng( mu, k, psi);
  }
    for (i in 1:N) {
    log_lik[i] = gev_lpdf(rep_vector(y[i],1) | mu, k, psi);
  }
}



  
