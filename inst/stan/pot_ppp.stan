// partially based on the code by Aki Vehtari
// https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
//
functions {
  real gpareto_lpdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log pdf 
    int N = rows(y);
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return -(1+inv_k)*sum(log1p((y-ymin) * (k/sigma))) -N*log(sigma);
    else
      return -sum(y-ymin)/sigma -N*log(sigma); // limit k->0
  }
  real gpareto_lcdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log cdf
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return sum(log1m_exp((-inv_k)*(log1p((y-ymin) * (k/sigma)))));
    else
      return sum(log1m_exp(-(y-ymin)/sigma)); // limit k->0
  }
  real gpareto_lccdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log ccdf
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return (-inv_k)*sum(log1p((y-ymin) * (k/sigma)));
    else
      return -sum(y-ymin)/sigma; // limit k->0
  }
  real gpareto_rng(real ymin, real k, real sigma) {
    // generalised Pareto rng
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return ymin + (uniform_rng(0,1)^-k -1) * sigma / k;
    else
      return ymin - sigma*log(uniform_rng(0,1)); // limit k->0
  }
  // For computing log like for annual maxima:
      real gev_lpdf(vector y, real mu, real k, real psi) {
    // gev log pdf 
    // int N = rows(y);
    real N = rows(y);
    real inv_k = inv(k);
    if (k<0 && max(y-mu)/psi > -inv_k)
        return (-9*log(10)); // very small probability: set to 10^-9
    if (k>0 && min(y-mu)/psi < -inv_k)
        return (-9*log(10)); // very small probability: set to 10^-9
    if (psi<=0)
      reject("psi<=0; found psi =", psi);
    if (fabs(k) > 1e-15)
      // return -(1+inv_k)*sum(log1p((y-mu) * (k/sigma))) -N*log(sigma);
      return (-N*log(psi) -(1+inv_k)*sum(log1p( (y-mu)*(k/psi) ))
               -sum( exp( -inv_k*log1p((y-mu)*(k/psi)) )));
    else
      // return -sum(y-mu)/sigma -N*log(sigma); // limit k->0
      return -N*log(psi) - sum( (y - mu)/psi + exp(-(y - mu)/psi) );
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
  real ymin; // threshold
  int<lower=0> M; // sample size (exceedances)
  int<lower = 1> My; // number of years of observations
  vector<lower=ymin>[M] y;
  int<lower=0> N[My];      // sample (number of events)
  int Mgen;     //
  real maxima[My]; // annual maxima for computing likelihood
  real pr_sigma[2]; // normal mu and sigma
  real pr_k[2];   //  normal mu and sigma
  real pr_lambda[2];   // prior, gamma a, b
}
transformed data {
  real ymax;
  ymax = max(y);
}
parameters {
  real<lower=0> sigma; 
  real<lower=-sigma/(ymax-ymin)> k;
  real<lower=0> lambda;
}
transformed parameters{
  real<lower=0> psi; // GEV scale
  real mu;  // GEV location
  mu = ymin + sigma/k*(lambda^(k) - 1);
  psi  = sigma*lambda^k;
}
model {
  // Model for N
  lambda ~ gamma(pr_lambda[1], pr_lambda[2]);
  for (m in 1:My) {
    target += poisson_lpmf(N[m] | lambda);
  }
  // Priors
  // k ~ normal(0.1, 0.1)T[-sigma/(ymax-ymin), positive_infinity()];
  // k ~ normal(pr_k[1], pr_k[2])T[-sigma/(ymax-ymin), positive_infinity()];
  k ~ normal(pr_k[1], pr_k[2]);
  // sigma ~ normal(pr_sigma[1], pr_sigma[2]);
  sigma ~ gamma(pr_sigma[1], pr_sigma[2]);
  // Model for exceedances
  y ~ gpareto(ymin, k, sigma);
}
generated quantities {
  // real<lower=0> psi; // GEV scale
  // real mu;  // GEV location
  real yrep[Mgen];
  real log_lik[My];
  // mu = ymin + sigma/k*(lambda^(k) - 1);
  // psi  = sigma*lambda^k;

  for (i in 1:Mgen){
    yrep[i] = gev_rng( mu, k, psi);
  }
    for (i in 1:My) {
    log_lik[i] = gev_lpdf(rep_vector(maxima[i],1) | mu, k, psi);
    // log_lik[i] = gev_lpdf(rep_vector(maxima[i],1) | mu, k, psi);
  }
}


