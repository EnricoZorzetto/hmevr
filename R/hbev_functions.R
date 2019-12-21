

reshuffle_events <- function(df, Nt = 366){
"----------------------------
  from a dataframe df produce a new
  dataframe dfr obtained reshuffling
  the number of events/ year
  and all the events
  --------------------------"
  table = table_max(df, Nt = Nt, reshuffle_days = TRUE)
  # now rebuild a dataframe from it
  # it's ok if years have all 366 events now
  PRCP = as.vector(t( table$data)) # by rows
  YEAR = rep(table$years, each = Nt)
  dfr = as.data.frame(cbind(YEAR, PRCP))
  return(dfr)
}


# hbev_cdf_fun <- function(x, C = C, W = W, N = N){
# "------------------------------------------
# cumulative distribution function
# for the block maxima in the hbev model
# x = quantile of interest
# C = vector with yearly Weibull scale
# W = vector with yearly Weibull shape
# N = vector with yearly number of events
# ------------------------------------------"
# cdfi = (1 - exp( -(x/C)^W))^N
# cdf = mean(cdfi)
# return(cdf)
# }
#
#
# hbev_pdf_fun <- function(x, C = C, W = W, N = N){
# "------------------------------------------
# probability density function
# for the block maxima in the hbev model
# x = value of interest
# C = vector with yearly Weibull scale
# W = vector with yearly Weibull shape
# N = vector with yearly number of events
# ------------------------------------------"
# M_gen = length(C)
# pdfi = N*(1 - exp( -(x/C)^W))^(N-1)*W/C*(x/C)^(W-1)*exp(-(x/C)^W)
# pdf = mean(pdfi)
# return(pdf)
# }

hbev_mix_cdf <- function(x, C1=Cs1, W1=ws1, C2=Cs2, W2=ws2, alp = alp,  N=Ns){
  # mixture of Weibull
  cdf = mean((alp*pweibull(x, W1, C1) + (1-alp)*pweibull(x, W2, C2))^N)
}

hbev_mix_pdf <- function(x, C1=Cs1, W1=ws1, C2=Cs2, W2=ws2, alp = alp,  N=Ns){
  # mixture of Weibull
  cdf = mean(  N*(alp*pweibull(x, W1, C1) + (1-alp)*pweibull(x, W2, C2))^(N-1)*
                 (alp*dweibull(x, W1, C1) + (1-alp)*dweibull(x, W2, C2)) )
}

hbev_mgw_cdf <- function(x, C=Cs, W=ws, a=as, b=bs, alp = alp,  N=Ns){
  # mixture of Gamma and Weibull
  cdf = mean((alp*pgamma(x, a, b) + (1-alp)*pweibull(x, W, C))^N)
}

hbev_mgw_pdf <- function(x, C=Cs, W=ws, a=as, b=bs, alp = alp,  N=Ns){
  # mixture of Gamma and Weibull
  cdf = mean(  N*(alp*pgamma(x, a, b) + (1-alp)*pweibull(x, W, C))^(N-1)*
                 (alp*dgamma(x, a, b) + (1-alp)*dweibull(x, W, C)) )
}


rmixweibull <- function(Ngen, C1, w1, C2, w2, alp){
  'draw Ngen values from a Weibull mixture'
  flipcoin = runif(Ngen)
  # print(flipcoin)
  R = rep(0, Ngen)
  # print(R)
  for (i in 1:Ngen){
    if (flipcoin[i] < alp) {
      R[i] = rweibull(1, w1, C1)
    } else {
      R[i] = rweibull(1, w2, C2)
    }
  }
  return(R)
}


hbev_wei_cdf <- function(x, C = C, W = W, N = N){
  # only for scalar x!
  cdf = mean(   pweibull(x, W, C )^N  )
return(cdf)
}


hbev_wei_pdf <- function(x, C = C, W = W, N = N){
  # only for scalar x!
  pdf = mean( N*pweibull(x, W, C)^(N-1)*dweibull(x, W, C))
return(pdf)
}


hbev_gam_cdf <- function(x, a = a, b = b, N = N){
  # only for scalar x!
  cdf = mean(   pgamma(x, a, b )^N  )
  return(cdf)
}


hbev_gam_pdf <- function(x, a = a, b = b, N = N){
  # only for scalar x!
  pdf = mean( N*pgamma(x, a,b)^(N-1)*dgamma(x, a, b))
  return(pdf)
}

# generalized gamma
hbev_gga_cdf <- function(x, C = C, W = W, A = A, N = N){
  nu = W
  sigma = A/nu
  mu = C*sigma
  # only for scalar x!
  cdf = mean(   rmutil::pggamma(x, sigma, mu, nu)^N)
return(cdf)
}




hbev_gga_pdf <- function(x, C = C, W = W, A = A, N = N){
  nu = W
  sigma = A/nu
  mu = C*sigma
  # only for scalar x!
  pdf = mean(rmutil::pggamma(x, sigma, mu, nu)^(N-1)*rmutil::dggamma(x, sigma, mu, nu))
return(pdf)
}


LOGN_pars<-function(m, s){
  # return lognormal param mu, sigma
  # given its mean and stdv
  s2 = s^2
  sigma = sqrt(log(s2/m^2 + 1))
  mu = log(m) - 0.5*sigma^2
  return(list(mu = mu, sigma = sigma))
}
#
#
LOGN_moms<-function(mu, sigma){
  # return lognormal moments (mean and stdv)
  # given its parameters mu and sigma
  sigma2 = sigma^2
  m = exp(mu + sigma2/2)
  s2 = (exp(sigma2)-1)*exp(2*mu + sigma2)
  stdv = sqrt(s2)
  return(list(m = m, stdv = stdv))
 }
#
 IG_mom<-function(mean, var){
  a = mean^2/var + 2
  b = mean*(mean^2/var + 1)
  return(list(a = a, b = b))
 }


fit_betabin <-function(sample, Nt = 365){
  # fit binomial distribution
  # by matching first two sample raw moments
  raw1 <- mean(sample)
  raw2 <- mean(sample^2)
  alpha <- (raw1*Nt - raw2)/( Nt*(raw2/raw1 - raw1- 1) + raw1)
  beta <- (Nt - raw1)*(Nt-raw2/raw1)/( Nt*(raw2/raw1 - raw1- 1) + raw1)
  # these parameter estimates can be negative when the sample is underdispersed.
  # In this case we assume a binomial distribution as prior
  if (alpha < 0){
  print("fit_betabin WARNING: negative parameters found underdispersion. Assuming binomial instead")
    # assume alpha + beta = theirsum
    # large so small variance, and given mean = binomial arrival rate
    # theirsum = 500
  # alpha = raw1*theirsum
  # beta = alpha - theirsum
    alpha = 150
    beta = 150
  }
  parhat = list(alpha = alpha, beta = beta)
  return(parhat)
}


BB_mom<-function(mean, var, Nt = 365){
  raw1 = mean
  raw2 = var + raw1^2
  a <- (raw1*Nt - raw2)/( Nt*(raw2/raw1 - raw1- 1) + raw1)
  b <- (Nt - raw1)*(Nt-raw2/raw1)/( Nt*(raw2/raw1 - raw1- 1) + raw1)
  if (a < 0){
    print("BB_mom WARNING: betabinomial model works only
                             for overdisperse variables")
  }
  return(list(an=a, bn=b))
}



BB_meanvar <-function(a, b, Nt = 365){
  bbmean = Nt*a/(a + b)
  bbvar = Nt*a*b*(a + b + Nt)/(a + b)^2/(a + b + 1)
  return(list(mean = bbmean, var = bbvar))
}


# decluster <- function(datamat, up_to = 50, first_indep = 50, maxlag = 365,
#                        #upper_quantile = 0.9,
#                       tau_max = 15,
#                       signif_lim = 0.1){
#   if (is.matrix(datamat)){
#   shape = dim(datamat)
#   M = shape[1]
#   Nt = shape[2]
#   dimN = Nt*M
#   my_vec = as.vector(datamat) # original dataset
#   } else if (is.data.frame(datamat)) {
#     my_vec = as.vector(datamat$PRCP)
#   dimN = length(my_vec)
#   } else {
#     print('decluster ERROR: input datamat must be either
#                             a matrix or a dataframe!')
#   }
#   dec_vec = my_vec  # initialize declustered dataset
#
#   # compute correlation:
#   myacf = acf(my_vec, lag.max = maxlag)
#   short_range = myacf$acf[1:up_to]
#   long_range  = myacf$acf[first_indep:maxlag]
#   # signif_lim = quantile(long_range, upper_quantile)
#   # signif_lim = qnorm(upper_quantile, mean = mean(long_range), sd = sd(long_range))
#   # signif_lim = mean(long_range) + 5*sd(long_range)
#   signif_lim = signif_lim
#   Tau <- which(myacf$acf <= signif_lim)[1]-1 # last point above limit
#   # Tau = 20
#
#   # check maximum reasonable value for Tau:
#   print(sprintf("decorrelation Tau = %s", Tau))
#   if (Tau > up_to){
#     Tau = up_to
#     print('decluster WARNING: time series seems too much correlated. Increase up_to')
#   }
#
#   if (Tau > tau_max){
#     Tau = tau_max
#     print('Decluster WARNING: reached tau_max')
#   }
#
#   # now slide a running window of size Tau,
#   # and in each window keep only the maximum
#
#   count = 0
#   end = 0
#   while(end < dimN){
#     start = count + 1
#     end = count + Tau
#     window = my_vec[start:end]
#     # posmax = which(window == max(window)) + count
#     for (i in start:end){
#       if (my_vec[i] < max(window)){
#         dec_vec[i] = 0.0
#       }
#     }
#     count = count + 1
#   }
#   if (is.matrix(datamat)){
#   decluster_datamat <- matrix(dec_vec, nrow = M, byrow = TRUE)
#   } else {
#     df$PRCP = dec_vec
#     decluster_datamat = df
#   }
#   return(list( decdata = decluster_datamat,
#                tau = Tau, acf = myacf, upquant = signif_lim))
# }


#
# decluster <- function(datamat, up_to = 50, first_indep = 50, maxlag = 365,
#                       #upper_quantile = 0.9,
#                       tau_max = 5,
#                       signif_lim = 0.1){
#   '-----------------------------------------------------------------------------
#   decluster a time series of daily rainfall
#   INPUT:
#   datamat -> data in matrix nyears*ndays, or data frame with field PRCP
#
#
#
#   -----------------------------------------------------------------------------'
#   if (is.matrix(datamat)){
#     shape = dim(datamat)
#     M = shape[1]
#     Nt = shape[2]
#     dimN = Nt*M
#     my_vec = as.vector(datamat) # original dataset
#   } else if (is.data.frame(datamat)) {
#     my_vec = as.vector(datamat$PRCP)
#     dimN = length(my_vec)
#   } else {
#     print('decluster ERROR: input datamat must be either
#                             a matrix or a dataframe!')
#   }
#   dec_vec = my_vec  # initialize declustered dataset
#
#   # compute correlation:
#   myacf = acf(my_vec, lag.max = maxlag, plot = FALSE)
#   short_range = myacf$acf[1:up_to]
#   # long_range  = myacf$acf[first_indep:maxlag]
#   # signif_lim = quantile(long_range, upper_quantile)
#   # signif_lim = qnorm(upper_quantile, mean = mean(long_range), sd = sd(long_range))
#   # signif_lim = mean(long_range) + 5*sd(long_range)
#   # signif_lim = signif_lim
#   Tau <- which(myacf$acf <= signif_lim)[1]-1 # last point above limit
#   # print(Tau)
#   # Tau = 20
#
#   # check maximum reasonable value for Tau:
#   print(sprintf("decorrelation Tau = %s", Tau))
#   if (Tau > up_to){
#     Tau = up_to
#     print('decluster WARNING: time series seems too much correlated. Increase up_to')
#   }
#
#   if (Tau > tau_max){
#     Tau = tau_max
#     print('Decluster WARNING: reached tau_max')
#   }
#
#   # now slide a running window of size Tau,
#   # and in each window keep only the maximum
#
#   count = 0
#   end = 0
#   while(end < dimN){
#     start = count + 1
#     end = count + Tau
#     window = my_vec[start:end]
#     print(length())
#     # posmax = which(window == max(window)) + count
#     for (i in start:end){
#       if (my_vec[i] < max(window)){
#         dec_vec[i] = 0.0
#       }
#     }
#     count = count + 1
#   }
#   if (is.matrix(datamat)){
#     decluster_datamat <- matrix(dec_vec, nrow = M, byrow = TRUE)
#   } else {
#     df$PRCP = dec_vec
#     decluster_datamat = df
#   }
#   return(list( decdata = decluster_datamat,
#                tau = Tau, acf = myacf, upquant = signif_lim))
# }



decluster <- function(datamat, up_to = 50, first_indep = 50, maxlag = 365,
                      #upper_quantile = 0.9,
                      tau_max = 5,
                      signif_lim = 0.1){
  '-----------------------------------------------------------------------------
  decluster a time series of daily rainfall
  INPUT:
  datamat -> data in matrix nyears*ndays, or data frame with field PRCP



  -----------------------------------------------------------------------------'
  if (is.matrix(datamat)){
    shape = dim(datamat)
    M = shape[1]
    Nt = shape[2]
    dimN = Nt*M
    my_vec = as.vector(datamat) # original dataset
  }else if (is.vector(datamat)) {
    dimN = length(datamat)
    my_vec = datamat
  } else if (is.data.frame(datamat)) {
    my_vec = as.vector(datamat$PRCP)
    dimN = length(my_vec)
  } else {
    print('decluster ERROR: input datamat must be either
                            a matrix or a dataframe!')
  }
  dec_vec = my_vec  # initialize declustered dataset

  # compute correlation:
  myacf = acf(my_vec, lag.max = maxlag, plot = FALSE)
  short_range = myacf$acf[1:up_to]
  # long_range  = myacf$acf[first_indep:maxlag]
  # signif_lim = quantile(long_range, upper_quantile)
  # signif_lim = qnorm(upper_quantile, mean = mean(long_range), sd = sd(long_range))
  # signif_lim = mean(long_range) + 5*sd(long_range)
  # signif_lim = signif_lim
  Tau <- which(myacf$acf <= signif_lim)[1]-1 # last point above limit
  # print(Tau)
  # Tau = 20

  # check maximum reasonable value for Tau:
  # print(sprintf("decorrelation Tau = %s", Tau))
  if (Tau > up_to){
    Tau = up_to
    print('decluster WARNING: time series seems too much correlated. Increase up_to')
  }

  if (Tau > tau_max){
    Tau = tau_max
    print('Decluster WARNING: reached tau_max')
  }

  # now slide a running window of size Tau,
  # and in each window keep only the maximum

  count = 0
  end = 0
  while(end < dimN){
    start = count + 1
    end = count + Tau
    window = my_vec[start:end]
    # print(length(window))
    posmax = count + which.max(window)
    for (i in start:end){
      if (i != posmax){
      # if (my_vec[i] < max(window)){
        dec_vec[i] = 0.0
      }
    }
    count = count + 1
  }

  if (is.matrix(datamat)){
    decluster_datamat <- matrix(dec_vec, nrow = M, byrow = TRUE)
  } else if (is.data.frame(datamat)){
    df$PRCP <- dec_vec
    decluster_datamat <- df
  } else if (is.vector(datamat)){
   decluster_datamat <- dec_vec
  }

  return(list( decdata = decluster_datamat,
               tau = Tau, acf = myacf, upquant = signif_lim))
}





  # use the dgev from the Stan mode::
  # expose_stan_functions("gev.stan")
# check that integral of exp(gpareto_lpdf) (from ymin to ymax) matches with gpareto_cdf
# mu = 10
# k = 0.1
# psi = 5
# y = rgev(1000, mu, k, psi)
# gev_pdf <- function(y, mu, psi, k) {
#   exp(sapply(y, FUN = gev_lpdf, mu = mu, k = k, psi = psi))
# }

# gev_pdf(c(1, 2, 3, 4, 5), mu, k, psi)


df_2_mat <- function(df, Nt=366, decluster = FALSE, signif_lim = 0.05){
    if (decluster){
    decres = decluster(df, signif_lim = signif_lim)
    df = decres$decdata
    # print(sprintf('Declustering time series: interval = %s days', decres$tau))
    # tau
    # acf -> also available
  }
years = unique(df$YEAR)
nyears = length(years)
datamat=matrix(0, nrow = nyears, ncol = Nt)
for (i in 1:nyears){
  df_ii = subset(df, df$YEAR == years[i])
  data_ii = df_ii$PRCP
  nobs_ii = length(data_ii)
  datamat[i, 1:nobs_ii] = data_ii
}
return(datamat)
}


reshapemat <- function(datamat, Nt){
  " reshape an array of elements of size nyears0*Nt0 (ordered row-wise)
  to another elements of size nyears*Nt with Nt specified.
  If not divisible throw away the eccess values
  * ALL ROW - WISE *"
datavec = as.vector(t(datamat))
# print(datavec)
nelem_old = length(datavec)
nelem = floor(nelem_old/Nt)*Nt
myvec = datavec[1:nelem]
mat = matrix(myvec, ncol = Nt, byrow = TRUE)
return(mat)
}

info_crit <- function(pdfmat, psis = TRUE){
  "-----------------------------------------------------------------------------
  given a matrix of p(yi | \theta^s) compute various information criteria
  pdfmat -> matrix with shape S*M, where
  S -> number of draw from the posterior = number of rows
  M -> sample size used for evaluation = number of columns
  if psis = TRUE: compute also pareto smoothed loo

  NB: In current version, lppd-lpml-elpd are computed with a -1/M factor,
  where M is the sample size, to make them positive (smallest = better)
  and more narrowly distributed
  -----------------------------------------------------------------------------"
  S = dim(pdfmat)[1]
  M = dim(pdfmat)[2]
  loglik = log(pdfmat)

  cpoi = rep(0, M)
  mean_pdfi = rep(0, M)
  mean_logl = rep(0, M)
  var_logl = rep(0, M)
  # for (i in 1:M){
  #     cpoi[i] = 1/mean(1/pdfmat[,i])
  #     mean_pdfi[i] = mean(pdfmat[,i])
  #     mean_logl[i] = mean(loglik[,i])
  #     var_logl[i] = var(loglik[,i])
  # }
    for (i in 1:M){
      cpoi[i] = 1/mean(1/pdfmat[,i])
      mean_pdfi[i] = mean(pdfmat[,i])
      mean_logl[i] = mean(loglik[,i])
      var_logl[i] = var(loglik[,i])
  }
  lpml_sum = sum(log(cpoi))
  lppd_sum = sum(log(mean_pdfi))
  lpml = -mean(log(cpoi))       # should be +sums
  lppd = -mean(log(mean_pdfi))  # should be +sums

  p_waic1 = 2*(lppd_sum - sum(mean_logl)) # effective number of parameters, waic 1
  p_waic2 = sum(var_logl) # effective number of parameters, waic 2

  # elpd_waic1 = lppd_sum - p_waic1
  # elpd_waic2 = lppd_sum - p_waic2
  elpd_waic1 = -1/M*(lppd_sum - p_waic1) # should be +sums
  elpd_waic2 = -1/M*(lppd_sum - p_waic2) # should be +sums

  if (psis == TRUE){
    psis_loo = loo(loglik)
    elpd_loo = psis_loo$estimates[1] # elpd, psis-loo
    p_loo = psis_loo$estimates[2] # eff number of parameters, psis-loo
  } else {
    elpd_loo = 0
    p_loo = 0
  }

 infocrits = list(lpml = lpml, lppd = lppd, elpd_waic1 = elpd_waic1,
                  elpd_waic2 = elpd_waic2,
                  p_waic1 = p_waic1, p_waic2 = p_waic2,
                  elpd_loo = elpd_loo, p_loo = p_loo)
 return(infocrits)
 }


skewness <- function(x){
  m <- mean(x)
  s = sd(x)
  skew = mean((x-m)^3)/(s^3)
}
