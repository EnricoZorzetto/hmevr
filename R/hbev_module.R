

#' Compute block (yearly) rainfall statistics for a rainfall dataset.
#'
#' @param df A data frame with the precipitation data. It must be a data frame
#'           with fields PRCP (daily precipitation) and
#'           YEAR (integer in format YYYY)
#' @param Nt The number of events/year (default is Nt = 366). If record length in an
#'           year is shorter than Nt, zero values are added.
#'           If longer, excess values are dropped.
#' @param reshuffle_days if TRUE, the number of events/block
#'                    and the intensities of all events
#'                    in the time series are
#'                    resampled without resubstitution
#'                    (default is reshuffle_days = FALSE)
#' @return a list with the following quantities:
#' \describe{
#'   \item{data}{matrix of size nyears*Nt with the daily values.
#'     Each row is an yearly sample (with non-zero values at the beginning)}
#'   \item{max}{block maxima values (array with length = nyears)}
#'   \item{Xi}{annual maxima sorted in ascending order (array of size nyears)}
#'   \item{Fi}{ Empirical non exceedance frequency for the annual maxima Xi
#'         (array of size nyears)
#'         computed by means of the Weibull plotting position formula}
#'   \item{Tr}{ Empirical return time for the annual maxima Xi
#'         (array of size nyears)
#'         computed by means of the Weibull plotting position formula}
#'  \item{N}{number of events in each year (array of size nyears)}
#'  \item{years}{years in the dataset (array of size nyears)}
#'  \item{nyears}{number of years in the dataset}
#'  \item{min}{block minima of non-zero values in each year of record
#'        (array of size nyears)}
#'  \item{totals}{yearly precipitation totals (array of size nyears)}
#'  \item{sdwets}{standard deviation of non-zero precipitation values
#'        for each year (array of size nyears)}
#'  \item{skews}{skewness of non-zero precipitation values
#'        for each year (array of size nyears)}
#' }
#' @examples
#' \dontrun{
#' table_max(nycp)
#' }
#' @export
table_max <- function(df, Nt = 366, reshuffle_days = FALSE){
  years  = unique(df$YEAR)
  nyears = length(years)
  maxima = array(0, nyears)
  minima = array(0, nyears)
  totals = array(0, nyears)
  sdwets = array(0, nyears)
  skews = array(0, nyears)
  N = array(0, nyears)
  datamat = matrix(0, nyears, Nt)
  count = 0
  for (i in 1:nyears){
    # print(years[i])
    all_events = df$PRCP[df$YEAR == years[i]]
    if (length(all_events)>Nt){
      print('table_max ERROR: some years of record have more than Nt obs.')
    }
    wets = all_events[all_events > 0]
    N[i] = length(wets)
    if (N[i] > 0){
      maxima[i] = max(wets)
      minima[i] = min(wets)
      totals[i] = sum(wets)
      sdwets[i] = sd(wets)
      skews[i] = skewness(wets)
      datamat[i,1:length(wets)] = wets
    }
    else{
      count = count + 1
      sprintf('table_max WARNING: some %s years have no non-zero observations!', count)
    }
  }
  # reshuffle all days and number of events in the series if requested
  if (reshuffle_days){
      N2 = sample(N, nyears, replace = FALSE)
      all_wets = datamat[datamat > 0]
      numwets = length(all_wets)
      all_wets2 = sample(all_wets, numwets, replace = FALSE)
      datamat2 = matrix(0, nyears, Nt)
      maxima2 = array(0, nyears)
      minima2 = array(0, nyears)
      totals2 = array(0, nyears)
      sdwets2 = array(0, nyears)
      skews2  = array(0, nyears)
      count2 = 0
      for (i in 1:nyears){
        if (N2[i] > 0){
          datamat2[i, 1:N2[i]] = all_wets2[(count2 + 1):(count2 + N2[i])]
          maxima2[i] = max(all_wets2[(count2 + 1):(count2 + N2[i])])
          minima2[i] = min(all_wets2[(count2 + 1):(count2 + N2[i])])
          totals2[i] =  sum(all_wets2[(count2 + 1):(count2 + N2[i])])
          sdwets2[i] =  sd(all_wets2[(count2 + 1):(count2 + N2[i])])
          skews2[i] =  skewness(all_wets2[(count2 + 1):(count2 + N2[i])])
          count2 = count2 + N2[i]
        }
      }
  } else{
    datamat2 = datamat
    N2 = N
    maxima2 = maxima
    minima2 = minima
    totals2 = totals
    sdwets2 = sdwets
    skews2 = skews
  }
  Fi = (1:nyears)/(1+nyears)
  Xi = sort(maxima)
  Tr = 1/(1-Fi)
  res = list(data = datamat2, max = maxima2, Fi = Fi, Xi = Xi, Tr = Tr,
             N = N2, years = years, nyears = nyears, min = minima2,
             totals = totals, sdwets = sdwets, skews = skews)
  return(res)
}

#'   Generate syntetic data according to a given specification.
#'   generates a matrix xij of synthetic data with dimension S*Nt
#'   where S = number of blocks (e.g., years), and
#'         Nt = number of observations/block
#'   with nj the number of non-zero events in each block
#'   Data are i.i.d., and non-zero events are stacked at the 'beginning' of each
#'   block (e.g., their random inter-arrival times are not reproduced here!)
#'   The distribution for the xij magnitudes and number of events nj
#'   must be specified from the possible available.
#'
#'   ----------------------------------------------------------------------------
#'      'bin'      -> nj ~ binomial(Nt, pn)
#'      'betabin'  -> nj ~ betabinomial(Nt, an, bn) or [mn, varn]
#'      'constant' -> nj ~ diracDelta(nc)
#'  -----------------------------------------------------------------------------

#'
#'  The magnitudes of the events xij are drawn from one of the following models::
#'  -----------------------------------------------------------------------------
#'      Model specification:                         | Parameters to pass:
#'  -----------------------------------------------------------------------------
#'      'gam'     -> xij ~ gamma(a, b)               |  [a, b]
#'      'gpd'     -> xij ~ gpd(xi, sigma)            |  [xi, sigma]
#'      'wei'     -> xij ~ Weibull(C, w)             |  [C, w]
#'      'wei_dyl' -> xij ~ Weibull(Cj, wj) with      |  [muc, sigmac, muw, sigmaw]
#'                                                   |  or:   [mc, sc, mw, sw]
#'                   Cj  ~ lognormal(muc, sigmac)    |
#'                   wj  ~ lognormal(muw, sigmaw)    |
#'      'wei_dsl' -> xij ~ Weibull(Cj, w) with       |  [muc, sigmac, w]
#'                   Cj  ~ lognormal(muc, sigmac)    |
#'      'wei_dyn' -> xij ~ Weibull(Cj, wj) with      |  [ac, bc, aw, bw]
#'                                                   |  or: [mc, sc, mw, sw]
#'                    Cj   ~ gamma(ac, bc)           |
#'                    wj   ~ gamma(aw, bw)           |
#'      'gam_dyn' ->  xij ~ gamma(aj, bj) with       |  [aga, bga, agb, bgb]
#'                                                   |  or: [ma, sa, mb, sb]
#'                    aj  ~ gamma(aga, bga)          |
#'                    bj  ~ gamma(agb, bgb)          |
#'      'gan_dsc' ->  xij ~ gamma(a, bj) with        |  [a, agb, bgb]
#'                    bj   ~ gamma(agb, bgb)         |
#'  -----------------------------------------------------------------------------
#'
#'  NB: for n~betabinomial, I can pass one of the following:
#'  ntrue = list(an=, bn=), or ntrue = list(mn=, varn=)
#'  i.e., mean and variance prescribed instead of the beta-binomial parameters.
#'
#'  NB: for Ptrue, If *-> you can pass vector of means and stdvs instead
#'  as follows: list(ma=, sa=, mb=, sb=) or list(mc=, sc=, mw=, sw=)
#'
#'
#'   ptrue and ntrue must be named lists
#'   with parameter names as mentioned above for each distribution
#'   e.g., pass ptrue = list(pn = 0.4) in case of a binomial distr.
#'   ptrue = list(muw = 0, muc = 2, sigmac = 0.05, sigmaw = 0.05)
#'   in case of 'wei_dyl' distribution.
#'
#'   Returns a named list containing::
#'   data   = the matrix with the data, of size S*Nt
#'   N      = array with number of events / block (length S)
#'   maxima = array of block maxima (of length S)
#'   Xi     = sorted annual maxima (length S)
#'   Fi     = non exceedance frequency of the Xi (length S)
#'   Tr     = empirical return time of the Xi (length S)
#'  ------------------------------------------------------------
#'
#' @param S length of the time series to be generated (in years / blocks)
#' @param ptrue list of parameters describing the distribution of event magnitudes
#' @param ntrue list of parameters describing the distribution of the number of events/block
#' @param Nt the length of each block (default Nt = 366)
#' @param ndist distribution of the number of events/block
#'   (default is ndist='betabin' for the betabinomial distribution)
#' @param dist distribution of the event magnitudes
#'   (default is dist='gam' for the gamma distribution)
#'
#' @return a list with the following quantities:
#' \describe{
#'   \item{data}{matrix of size nyears*Nt with the daily values.
#'     Each row is an yearly sample}
#'   \item{maxima}{block maxima values (array with length = nyears)}
#'   \item{Xi}{annual maxima sorted in ascending order (array of size nyears)}
#'   \item{Fi}{ Empirical non exceedance frequency for the annual maxima Xi
#'         (array of size nyears)
#'         computed by means of the Weibull plotting position formula}
#'   \item{Tr}{ Empirical return time for the annual maxima Xi
#'         (array of size nyears)
#'         computed by means of the Weibull plotting position formula}
#'  \item{N}{number of events in each year (array of size nyears)}
#'
#' @export
load_synth_data <- function(S, ptrue, ntrue, Nt = 366,
                            ndist = 'betabin', dist = 'gam'){


  if (ndist == 'bin'){
    nj = rbinom(S, Nt, ntrue$pn)
  } else if (ndist == 'bbn') {
        if (is.null(ntrue[["an"]]) ){
          ntrue = BB_mom(ntrue$mn, ntrue$varn, Nt=Nt)
        }
    # nj = rbetabinom.ab(S, Nt, ntrue$an, ntrue$bn)
    nj = extraDistr::rbbinom(S, Nt, alpha = ntrue$an, beta = ntrue$bn)
  } else if (ndist == 'constant'){
    nj = rep(ntrue$nc, S)
  } else {
    print('load_synth_data ERROR: invalid ndist')
  }
data = matrix(0, S, Nt)
for (j in 1:S){
  if (dist == 'gam'){
      data[j,1:nj[j]] = rgamma(nj[j], ptrue$a, ptrue$b)
  } else if (dist == 'wei'){
      data[j,1:nj[j]] = rweibull(nj[j], ptrue$w, ptrue$C)
  } else if (dist == 'gpd'){
      # data[j,1:nj[j]] = evd::rgpd(nj[j], loc = 0,
      #               scale = ptrue$sigma, shape = ptrue$xi)
            data[j,1:nj[j]] = extraDistr::rgpd(nj[j], mu = 0,
                    sigma = ptrue$sigma, xi = ptrue$xi)
  } else if (dist == 'wei_dyl'){

    if (is.null(ptrue[["muc"]]) ){
      cpars = LOGN_pars(ptrue$mc, ptrue$sc)
      ptrue$muc =  cpars$mu
      ptrue$sigmac =  cpars$sigma
      wpars = LOGN_pars(ptrue$mw, ptrue$sw)
      ptrue$muw = wpars$mu
      ptrue$sigmaw = wpars$sigma
    }
      Cj = rlnorm(1, ptrue$muc, ptrue$sigmac)
      wj = rlnorm(1, ptrue$muw, ptrue$sigmaw)
      data[j,1:nj[j]] = rweibull(nj[j], wj, Cj)
  } else if (dist == 'wei_dsl'){
      Cj = rlnorm(1, ptrue$muc, ptrue$sigmac)
    data[j,1:nj[j]] = rweibull(nj[j], ptrue$w, Cj)
      } else if (dist == 'wei_dsc'){
      Cj = rgamma(1, ptrue$ac, ptrue$bc)
    data[j,1:nj[j]] = rweibull(nj[j], ptrue$w, Cj)
  } else if (dist == 'wei_dyn'){

    if (is.null(ptrue[["ac"]]) ){
      ptrue$ac = (ptrue$mc/ptrue$sc)^2
      ptrue$bc = ptrue$mc/ptrue$sc^2
      ptrue$aw = (ptrue$mw/ptrue$sw)^2
      ptrue$bw = ptrue$mw/ptrue$sw^2
    }
      Cj = rgamma(1, ptrue$ac, ptrue$bc)
      wj = rgamma(1, ptrue$aw, ptrue$bw)
      data[j,1:nj[j]] = rweibull(nj[j], wj, Cj)
  } else if (dist == 'wei_dgu'){
      # print('Hello World!')
      Cj = extraDistr::rgumbel(1, mu = ptrue$mc, sigma = ptrue$sc)
      wj = extraDistr::rgumbel(1, mu = ptrue$mw, sigma = ptrue$sw)
      data[j,1:nj[j]] = rweibull(nj[j], wj, Cj)
  } else if (dist == 'gam_dyn'){

    if (is.null(ptrue[["aga"]]) ){
      ptrue$aga = (ptrue$ma/ptrue$sa)^2
      ptrue$bga = ptrue$ma/ptrue$sa^2
      ptrue$agb = (ptrue$mb/ptrue$sb)^2
      ptrue$bgb = ptrue$mb/ptrue$sb^2
    }
    aj = rgamma(1, ptrue$aga, ptrue$bga)
    bj = rgamma(1, ptrue$agb, ptrue$bgb)
    data[j,1:nj[j]] = rgamma(nj[j], aj, bj)
  } else if (dist == 'gam_dsc'){
    bj = rgamma(1, ptrue$abg, ptrue$bgb)
    data[j,1:nj[j]] = rgamma(nj[j], ptrue$a, bj)
  }
}
  # compute annual maxima and their frequency
  maxima = apply(data, 1, max) # max in each row
  Fi = (1:S)/(1+S)
  Xi = sort(maxima)
  Tr = 1/(1-Fi)
  # return list with values of interest
  return( list(data = data, N = nj, maxima = maxima,
          Xi = Xi, Fi = Fi, Tr = Tr))
}

#' Load a dataset of rainfall observations
#'
#' @param filepath
#' @param readdata
#' @param maxmiss
#' @param min_nevents
#' @param dividebyten
#' @param Nt
#' @return df, dataframe with the precipitation data
#' @export
load_obs_data <- function(filepath, readdata=TRUE,
                          maxmiss = 30, min_nevents = 0,
                          dividebyten = TRUE , Nt = 366){
  "------------------------------------------------------------
  loads an observed time series from a csv file
  from the path + filename 'filepath', (if readdata =TRUE)

  else if readdata = FALSE, provide already a dataframe as filepath instead
  of the actual file path - must have fields YEAR and PRCP.


  csv file must include the variable of interest
  in a $PRCP field (e.g., rainfall accumulations)
  and a $DATE field, integer in format YYYYMMDD,
  OR a $YEAR field, integer in format YYYY.
  Remove from the record years with more than 'maxmiss' obs.
  (either missing data, negative values of NANs are removed)

  If dividebyten -> divide by 10 the values in $PRCP
  (original GHCN data are stored in [tenths of mm] -> GOTO [mm])
  Nt = expected number of events / year (default 366 = daily scale)

  Remove from the record all years with less than
  'min_nevents', minimum number of non zero events
  (default = 0).
  Some records might have many missing events in a few years.

  Returns a dataframe df with the dataset.
  -------------------------------------------------------------"

  # load file as a data frame or not
  if (readdata){
    df <-read.csv(file = filepath, header = TRUE, sep = ',')
  } else {
    df <- filepath
  }


  # read years if provided - else
  # compute years from date in format YYYYMMDD
  if(!("YEAR" %in% colnames(df))){
      myfun <- function(x) {floor(x/10000)}
      df$YEAR = mapply(myfun, df$DATE)
  }


  # change unit if requested
  if (dividebyten){
      df$PRCP = df$PRCP/10
  }

  # remove NaNs and -9999 if any
  df  = df[df$PRCP >= 0.0,]
  df = df[!is.nan(df$PRCP), ]

  # remove years with more than maxmiss missing observations
  # or with less than 'min_events' obs. greater than zero.
  years0 = unique(df$YEAR)
  nyears0= length(years0)
  for (i in 1:nyears0){
      year = years0[i]
      dfii = df[df$YEAR == year,]
      dfiiwets = dfii[dfii$PRCP > 0, ] # or some other detection threshold
      all_obs = dim(dfii)[1]
      all_events = dim(dfiiwets)[1]
      # two condition to remove current year from the record::
      cond1 = all_obs < (Nt - maxmiss)
      cond2 = all_events < min_nevents

      if ( cond1 || cond2) {
          df = df[df$YEAR != year, ]
      }
  }
  return(df)
}


#' @export
split_obs_data <- function(df, M_cal = 20,
                           M_val = 200,
                           Nt = 366,
                           cross_val = TRUE,
                           reshuffle = FALSE,
                           flip_time = TRUE,
                           reshuffle_days = FALSE,
                           decluster = FALSE,
                           signif_lim = 0.05){
  "---------------------------------------------------------
  split a dataframe df in two datasets.

  Extract two datasets for validation and calibration,
  with sizes M_val and M_cal years respectively.
  If M_cal and M_val are longer than the available
  time series, all available years are used instead.

  If cross_val -> the two dataset are independent.
  If not       -> the calibration sample is usually
                  contained in the validation one.

  If M_cal or M_val are longer than the sample size,
  only the available data is used, thus producing
  shorter samples, except in the case
  of cross_validation when is not possible to produce an
  independent validation sample as expected
  (e.g., if M_cal > sample size an independent sample is not
  created and an error is raised.)

  Id decluster = TRUE, before fitting the models we apply a
  declustering technique (default is FALSE):
  The time lag Tau at which the serial correlation decays below
  a given value 'signif_min' (default 0.05) is computed.
  Then the time series is declustered with a running window of
  length tau, in which only the largest observations is kept.

  signif_lim = Value of autocorrelation below which
  the time series is assumed to be uncorretated (used to determine the
  window size for declustering if decluster = TRUE)

  If reshuffle -> the order of the years in the record
  is reshuffled before splitting the series in the two datasets.
  (default id FALSE)

  If reshuffle_days -> the number of events and the intensities
                       of all events in the time series are
                       resampled without resubstitution
                       (default is FALSE)

  If flip_time -> If not reshuffle, start extracting years
                  from the end of the record,
                  flipping the order of observed yerars
                  (e.g., to simulate longer time series that
                  go progressively back longer in time)
                  extract as calibration the last M_cal years,
                  (beacuse more recent data may have better resolution),
                  and then extract the M_val years before (if cross_val),
                  or again extract the last M_val years (with overlapping)
                  if cross_val is FALSE.

  If neither reshuffle or flip time, start extracting
  samples from the beginning of the time series.

  Nt = number of events year (default 366, exclude leap years)

  Return a list of two lists (datacal, dataval), each containing:
  $data   = matrix of size nyears*Nt with the data
  $N      = yearly number of events
  $max    = block (annual) maxima
  $Xi     = sorted maxima
  $Fi     = empirical non exceedanc frequency of the maxima Xi
  $Tr     = empirical return time of the maxima Xi
  $years  = years from the record used for the dataset
  $nyears = length of the dataset in years
  ---------------------------------------------------------"

  if (decluster){
    decres = decluster(df, signif_lim = signif_lim)
    df = decres$decdata
    lagtau = decres$tau
    # print(sprintf('Declustering time series: interval = %s days', decres$tau))
    # tau
    # acf -> also available
  } else {
    lagtau = 0
  }


  # resuffle order of the years in the record if requested
  # else flip them to start sampling calibration years
  # from most recent observations
  years = unique(df$YEAR)
  nyears= length(years)
  if (reshuffle){
      perm_years = sample(years, nyears, replace = FALSE)
  } else if (flip_time){
      perm_years = sort(years, decreasing = TRUE)
  } else {
      perm_years = years
  }

  # if reshuffle days, do it here
  if (reshuffle_days){
   df = reshuffle_events(df, Nt = Nt)
  }

  # extract a datasets for calibration
  # from shuffled or flipped time series
  if (M_cal > nyears){
    writeLines('split_obs_data WARNING:
                requested calibration sample
                is longer than available time series,
                using all available years instead.')
  }
  upper_cal = min(M_cal, nyears)
  cal_years = perm_years[1:upper_cal] # last M_cal elements

  # now select years for validation
  if (cross_val){

      if ((M_cal < nyears) && ( (M_cal + M_val) > nyears)){
          writeLines("split_obs_data WARNING:
              requested valibration sample is longer
              than remaining independent time series,
              using all reamining years instead.")
          # upper_val = min( (M_cal+M_val), nyears)
          val_years = perm_years[(M_cal+1):nyears]
      }
      else if ( (M_cal + M_val) > nyears){
          writeLines('split_obs_data ERROR: not enough data
                      for cross validation. Independent
                      cross validation is not feasible.
                      Reduce requested sample size "M_val".')
          upper_val = NaN
          val_years = NaN
      }
      else {
          # upper_val = min( (M_cal+M_val), nyears)
          val_years = perm_years[(M_cal+1):(M_cal+M_val)]
      }
  }
  else{ # case of no cross validation
          if (M_cal > nyears){
              writeLines("split_obs_data WARNING:
                          requested valibration sample is longer
                          than available time series,
                          using all available years instead.")
          }
          upper_val = min(M_val, nyears)
          val_years = perm_years[1:upper_val]
  }

  # extract the two datasets
  dfcal = df[df$YEAR %in% cal_years,]
  dfval = df[df$YEAR %in% val_years,]

  # now compute maxima for each dataset
  nyears_val = length(val_years)
  M_val = nyears_val
  cal1 = table_max(dfcal, Nt = Nt)
  val1 = table_max(dfval, Nt = Nt)

  # save lists with values of interest
  datacal = list(data = cal1$data, N = cal1$N, max = cal1$max,
                    Fi  = cal1$Fi, Xi = cal1$Xi, Tr = cal1$Tr,
                    years = cal1$years, nyears = cal1$nyears,
                    lagtau = lagtau)
  dataval = list(data = val1$data, N = val1$N, max = val1$max,
                    Fi  = val1$Fi, Xi = val1$Xi, Tr = val1$Tr,
                    years = val1$years, nyears = val1$nyears,
                    lagtau = lagtau)
  return( list(datacal = datacal, dataval = dataval))
}

#' Fit an extreme value model to data
#'
#' @param data matrix with daily precipitation data, with size
#' (number of blocks)*(number of observations in each block)
#' so that each row is the sample within one block.
#' Alternatively, it can be a vector with the block maxima if one sets
#' the parameter onlymax_gev = TRUE (only for the GEV model).
#' @param model which extreme value model.
#'              Default is model='wei_dyn_bin'
#' @param Mgen Number of samples to draw for the latent level variable models
#' ( for models of the HBEV family only)
#' @param thresh_pot threshold for the Peak-Over-Threshold / Point Poisson Process
#' (default thresh_pot = 10 units) not used by default!
#' @param thresh_hbev threshold for the hbev family models
#' (default thresh_pot = 0 units)
#' @param pot_frac if true, use frac_exc to determine POT model threshold
#' instead of the given thresh_pot (default TRUE)
#' @param frac_exc fraction of non-zero observations above threshold
#' (default is 0.05, for POT model only)
#' @param iter number of iterations in each chains for HMC.
#' By default 50% burn it period.
#' @param chains number of parallel chains for HMC sampler
#' @param onlymax_gev if TRUE, allow for data = vector of annual maxima
#' (default FALSE, only for GEV model)
#' @param refresh To show output or not (default refresh = 0, no output shown)
#' @param empirical_prior If TRUE, use some empirical priors
#' (default empirical_prior = FALSE)
#' @param adapt_delta for HMC sampler (default adapt_delta = 0.8)
#' @param priorpar to provide specific prior parameter values (default NULL)
#' @param draw_priors if TRUE produce density for the priors
#' (mostly for plotting)

#' @return a list with the following quantities:
#' \describe{
#'   \item{model}{model name}
#'   \item{Nt}{number of observations / block}
#'   \item{prior}{prior distributions used (parameters, type and possibly rng)}
#'   \item{Mgen}{number of samples for latent level paramaters}
#'   \item{thresh_pot}{threshold used for POT model}
#'  \item{thresh_hbev}{threshold used for hbev models}
#'  \item{model_fit}{Stan object with the model fit}
#' }
#' @export
fit_ev_model <- function(data, model = 'wei_dyn_bin', Mgen = 50,
                         thresh_pot = 10, thresh_hbev = 0.0,
                         pot_frac = TRUE, frac_exc = 0.05,
                         iter = 1000, chains = 4,
                         onlymax_gev = FALSE, refresh = 0,
                         empirical_prior = FALSE,
                         adapt_delta = 0.8,
                         max_treedepth = 10,
                         priorpar = NULL,
                         draw_priors=FALSE)
{
  "-------------------------------------------------------------------
  fit one of the following extreme value models to data:
  model = 'gev' for the Generalized Extreme Value Distribution
        = 'pot' for the Poisson Process / Peak Over Threshold
        = 'hbevs' for the static Weibull model
        = 'hbevm' for the complete hierarchical model (latent lognormal)
        = 'hbevd' for the complete hierarchical model, with latent invgamma

  *** NEW MODELS *****************************************************
  'wdm' -> weibull dynamic model [complete] (latent lognormal)
  'wds' -> weibull dynamic scale [only]
  'wst' -> weibull static
  'gdm' -> gamma dynamic model [complete]
  'gds' -> gamma dynamic scale [only]
  'gst' -> gamma static
  *******************************************************************

  For the models of the hbev* family only, uose one of the
  following models for the number of events/block:
  (not applicabile for GEV; POT uses Poisson model only)
  N_model = 'bin' for Binomial
          = 'betabin' for Betabinomial

  thresh_pot = threshold for the POT model (default = 10 units)
  thresh_hbev = threshold for the hbev* models (default = 0 units)
  niter = 1000 number of iterations per chain
          (default sampling after 500 warmup iterations)
  nchains = 4

  data = matrix of observations, with size nblocks*nobs oer block
                      (e.g., nyears*ndays)
  If 'onlymax_gev = TRUE', provide a sample of annual maxima instead of
  the daily matrix data.

  Generate Mgen years of data from the posterior distribution

  If pot_frac = TRUE, instead of using thresh_pot, select a threshold
  such that only a fraction 'frac exc' of the nonzero values is above
  threshold (by default is 0.05, i.e., keeps 5% of the non zero
  events as values above threshold)

  refresh=0 -> no output displayed

  Returns a list with the following named quantitites:
  model       ->  (model name   - one of the above)
  N_model     ->  (N model name - one of the above)
  Nt          ->  (number of observations / block)
  thresh_pot  ->  (as in input)
  thresh_hbev ->  (as in input)
  model_fit   ->  (fit object as returned by Stan)
  N_model_fit ->  (fit object as returned by Stan)
  prior       ->  (information on the prior distributions used)
  N_prior     ->  (information on the prior distributions used for N)
  ---------------------------------------------------------------------"
  nyears = dim(data)[1]
  Nt     = dim(data)[2]
 if (model == 'gev' | model == 'gev_nobounds'){
    print(sprintf('Fitting the %s model', model))
    if (onlymax_gev){
      maxima = data
    } else{
      maxima = apply(data, 1, max)
    }
    # prior distributions for each parameters
    mean_y = mean(maxima)
    sd_y   = sd(maxima)
    pr_mu = c(mean_y, 0.3*mean_y)   # normal mu, sigma
    # pr_psi = c(log(sd_y), 0.3)    # lognormal mu, sigma
    pr_psi = c(0.5*sd_y, 0.5)    # gamma mu, sigma
    pr_k = c(0.114, 0.125)           # normal mu, sigma
    priorlist = list(pr_mu=pr_mu, pr_psi=pr_psi, pr_k=pr_k)
    model_data0 <- list(N=length(maxima), y=maxima, Mgen=Mgen)
    model_data = append(model_data0, priorlist)
    modelfile = sprintf('%s.stan', model)

    model_fit  <- rstan::sampling(stanmodels[[model]], data = model_data,
                       chains = chains, iter = iter,
                       init_r = 0.02, refresh = refresh,
                    control = list(adapt_delta = adapt_delta,
                                   max_treedepth = max_treedepth))
  }
  else if (model == 'pot_ppp'){
    print(sprintf('Fitting the %s model', model))
    # compute threshold based on given exceedance prob:
    if (pot_frac){
      wets = data[data > 0]
      thresh_pot = quantile(wets, 1-frac_exc)
    }
    print(sprintf('pot-ppp: threshold value = %s', thresh_pot))
    # print("Hello World!")
    maxima = apply(data, 1, max)
    # print(length(maxima))
    Nex = array(0, nyears)
    for (i in 1:nyears){
      rain = data[i,]
      Nex[i] = length(rain[rain > thresh_pot])
    }
    datavec = array(data)
    exceedances = datavec[datavec > thresh_pot]
    # not excesses, the threshold is passed as location parameter.

    # prior distributions for each parameters:
    sd_y   = sd(exceedances) # no need to subtract threshold here
    # pr_sigma = c(0.4*sd_y, 0.4)  # gamma a, b
    # pr_xi01 = c(9,6)             # beta a, b
    pr_sigma = c(0.5*sd_y, 0.5)    # gamma mu, sigma
    pr_k = c(0.114, 0.125)             # normal mu, sigma [truncated]
    # pr_lambda = c(1.2, 0.1)      # gamma a, b for Poisson
    pr_lambda = c(2, 0.5)      # gamma a, b for Poisson
    print(sprintf('pot-ppp: %s excesses in %s years', length(exceedances), nyears))

    priorlist = list(pr_sigma = pr_sigma, pr_k  = pr_k, pr_lambda = pr_lambda)
    model_data0 = list(ymin = thresh_pot, M = length(exceedances),
                      My = nyears, y = exceedances, N = Nex, Mgen=Mgen,
                       maxima = maxima)
    model_data = append(model_data0, priorlist)
    modelfile = sprintf('%s.stan', model)

    # model_fit = stan( file = modelfile, data = model_data,
    #                   chains = chains, iter=iter,
    #                   refresh = refresh,
    #                 control = list(adapt_delta = adapt_delta))

        model_fit = rstan::sampling(stanmodels[[model]], data = model_data,
                      chains = chains, iter=iter,
                      refresh = refresh,
                    control = list(adapt_delta = adapt_delta,
                                   max_treedepth = max_treedepth))
  }  else {
    # FIT HBEV MODEL
  print(sprintf('Fitting the %s model', model))
  priorlist = hbev_priors(model, empirical=empirical_prior,  data = data,
                          thresh_hbev=thresh_hbev,
                          draw_rng=draw_priors, ndraws = iter*chains/2,
                          priorpar = priorpar)
  exceedances = pmax(data-thresh_hbev, 0.0)
    N = rep(0, nyears)
    for (i in 1:nyears){
      samplei = exceedances[i,]
      N[i] = length(samplei[samplei > 0])
    }
  model_data0 <- list(M = nyears, Mgen=Mgen, Nt = Nt, y = exceedances, N=N)
  model_data = append(model_data0, priorlist$prior_par)
  modelfile = sprintf('%s.stan',model)
  model_fit <- rstan::sampling(stanmodels[[model]], data = model_data,
                    iter = iter, chains = chains, refresh = refresh,
                    control = list(adapt_delta = adapt_delta,
                                   max_treedepth = max_treedepth))
  }
    return( list(model = model, Nt = Nt, prior = priorlist, Mgen=Mgen,
               thresh_pot = thresh_pot, thresh_hbev = thresh_hbev,
               model_fit = model_fit))
}




#' Compute quantiles and goodness of fit measures for a fitted model
#'
#' @param model_fit list as returned from the fit_ev_model function
#' @param maxval array of annual maxima to be used for
#' model validation and for computing quantiles
#' @param trmin minimum return time for computing quantiles,
#' and to be included in the fractional square error / mean bias measures
#' (default trmin = 2)
#' @return a list with the following named quantities
#' \describe{
#'   \item{qmean}{expected value for the computed quantiles (for Tr >= trmin)}
#'   \item{qupper}{upper credidility interval for quantiles (prob = 0.95)
#'   (quantiles only for Tr >= trmin)}
#'   \item{qlower}{lower credibility interval for quantiles (prob = 0.05)
#'   (quantiles only for Tr >= trmin)}
#'   \item{qmedian}{median value of the computed quantiles  (prob = 0.50)
#'   (quantiles only for Tr >= trmin)}
#'   \item{fmean}
#'   \item{fupper}
#'   \item{flower}
#'   \item{fmedian}
#'   \item{dmean}
#'   \item{dupper}
#'   \item{dlower}
#'   \item{dmedian}
#'   \item{lpml}
#'   \item{lppd}
#'   \item{fse}
#'   \item{mbias}
#'   \item{mwidth}
#'   \item{fse_Tr}
#'   \item{bias_Tr}
#'   \item{wisth_Tr}
#'   \item{elpd_loo}
#'   \item{p_loo}
#'   \item{elpd_waic1}
#'   \item{p_waic1}
#'   \item{elpd_waic2}
#'   \item{p_waic2}
#'   \item{quants}
#'   \item{cdfs}
#'   \item{pdfs}
#'   \item{epsi}
#'   \item{Fival}
#'   \item{Trval}
#'   \item{Xival}
#'   \item{FivalQ}
#'   \item{TrvalQ}
#'   \item{XivalQ}
#'   \item{log_lik}{log likelihood (nyears*ndraws)}
#'   \item{model}{name of the model}
#'   \item{Mgen}{number of draws for latent variable level for hbev models}
#'   \item{thresh_hbev}{theshold used for the hbev-family models}
#'   \item{nwarning_x0}{number of warnings for optimization convergence problems}
#'   }
#' @export
comp_quant <- function(model_fit, maxval, trmin = 2){
  Mgen = model_fit$Mgen
  nwarning_x0 = 0 # if true there are numerical problems with optimization
  Np = length(maxval)     # number of years for the validation
  Fival = (1:Np)/(1+Np)   # Weibull plotting position non exceedance frequency
  Xival = sort(maxval)    # sorted maxima
  Trval = 1/(1-Fival)     # return time empirical values

  # for quantiles only, compute only for Tr > trmin
  TrvalQ = Trval[Trval >= trmin]
  FivalQ = Fival[Trval >= trmin]
  XivalQ = Xival[Trval >= trmin]
  NpQ = length(XivalQ)

  model = model_fit$model
  Nt = model_fit$Nt
  model_sim = rstan::extract(model_fit$model_fit)
  S = dim(model_sim$lp__) # number of draws from the posterior

  # quantities to compute:
  quants = matrix(0, nrow = S, ncol =  NpQ)
  cdfs   = matrix(0, nrow = S, ncol =  Np )
  pdfs   = matrix(0, nrow = S, ncol =  Np )
  epsi   = matrix(0, nrow = S, ncol =  NpQ)

  if (model == 'gev' | model == 'pot_ppp' | model == 'gev_nobounds'){
      print(sprintf('computing quantiles for the %s model', model))
      xi     = model_sim$k
      psi  = model_sim$psi
      mu     = model_sim$mu
      for (s in 1:S){
        for (i in 1:NpQ){
          # from evd
          # quants[s, i] = qgev(FivalQ[i], loc = mu[s],
          #               scale = psi[s], shape = xi[s])
          # from extraDistr
          quants[s, i] = extraDistr::qgev(FivalQ[i], mu = mu[s],
                        sigma = psi[s], xi = xi[s])
          # epsi[s, i] = (XivalQ[i] - quants[s, i])/XivalQ[i]
          epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
        }
        for (i in 1:Np){
          # from evd
          # cdfs[s, i] = pgev(Xival[i], loc = mu[s],
          #               scale = psi[s], shape = xi[s])
          # pdfs[s, i] = dgev(Xival[i], loc = mu[s],
          #               scale = psi[s], shape = xi[s])
          # from extraDistr
          cdfs[s, i] = extraDistr::pgev(Xival[i], mu = mu[s],
                        sigma = psi[s], xi = xi[s])
          pdfs[s, i] = extraDistr::dgev(Xival[i], mu = mu[s],
                        sigma = psi[s], xi = xi[s])
                    # pdfs[s, i] = gev_pdf(Xival[i],mu[s],psi[s],xi[s])
          if (pdfs[s, i] < 1e-9){
            pdfs[s, i] <- 1e-9
          }
          }
        }
      }
      else if (substr(model, 1, 3) == 'mix') { # mixture model
        thresh = model_fit$thresh_hbev
        for (s in 1:S){
            Cs1 = model_sim$Cgen1[s,]
            ws1 = model_sim$wgen1[s,]
            Cs2 = model_sim$Cgen2[s,]
            ws2 = model_sim$wgen2[s,]
            alp = model_sim$rho[s, ]
            Ns = model_sim$Ngen[s,]
            for (i in 1:NpQ){
                myfun<-function(x) hbev_mix_cdf(x, C1=Cs1, W1=ws1, C2=Cs2,
                                          W2=ws2, alp = alp,  N=Ns) - FivalQ[i]
                # x0 = 8*mean(Cs)
                # x0 = 6*mean(Cs)/mean(ws)*(1 + log10(TrvalQ[i])/3)
                # F0 <- 0.1 # Tr 10 anni to start :: hip constant params::
                # x0 <-  mean(Cs1)*(log(mean(Ns)/(1-F0)))^(1/mean(ws1))
                # x0 = 80
                Cs = mean(c(mean(Cs1), mean(Cs2)))
                ws = mean(c(mean(ws1), mean(ws2)))
                # F0 <- 0.8*FivalQ[] # Tr 10 anni to start :: hip constant params::
                F0 <- 0.9 # Tr 10 anni to start :: hip constant params::
                x0 <-  mean(Cs)*(log(mean(Ns)/(1-F0)))^(1/mean(ws))
                # x0 <-  mean(Cs2)*(log(mean(Ns)/(1-F0)))^(1/mean(ws2))
                optim = nleqslv(x0, myfun)
                # optim = nleqslv(x0, myfun)
                # print(optim$x)
                # print(dim(quants))
                quants[s, i] = optim$x + thresh
                termcd = optim$termcd
                if (termcd != 1){
                  print("comp_quant WARNING: nleqslv might not be converging")
                  print("If 6 Jacobian is singular, new guess")
                  print("termcd = ")
                  print(termcd)
                  nwarning_x0 = nwarning_x0 + 1
                }
                # epsi[s, i]   = (XivalQ[i] - quants[s, i])/XivalQ[i]
                epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
            }
            for (i in 1:Np){
                qfval = Xival[i] - thresh
                if (qfval <= 0){ # below threshold
                  # print("comp_quant ERROR: quantile below threshold")
                  cdfs[s, i] = 0.0
                  pdfs[s, i] = 0.0
                } else {
                cdfs[s, i] = hbev_mix_cdf(qfval, C1=Cs1, W1=ws1, C2=Cs2,
                                          W2=ws2, alp = alp,  N=Ns)
                pdfs[s, i] = hbev_mix_pdf(qfval, C1=Cs1, W1=ws1, C2=Cs2,
                                          W2=ws2, alp = alp,  N=Ns)
                }
            }
        }

         }    else if (substr(model, 1, 3) == 'mgw') { # mixture gamma / weibull
        thresh = model_fit$thresh_hbev
        for (s in 1:S){
            Cs = model_sim$Cgen[s,]
            ws = model_sim$wgen[s,]
            as = model_sim$agen[s,]
            bs = model_sim$bgen[s,]
            alps = model_sim$rho[s, ]
            Ns = model_sim$Ngen[s,]
            for (i in 1:NpQ){
                myfun<-function(x) hbev_mgw_cdf(x, C=Cs, W=ws, a=as,
                                          b=bs, alp = alps,  N=Ns) - FivalQ[i]
                # x0 = 8*mean(Cs)
                # x0 = 6*mean(Cs)/mean(ws)*(1 + log10(TrvalQ[i])/3)
                # F0 <- 0.1 # Tr 10 anni to start :: hip constant params::
                # x0 <-  mean(Cs1)*(log(mean(Ns)/(1-F0)))^(1/mean(ws1))
                # x0 = 80
                Cs = mean(c(mean(Cs1), mean(Cs2)))
                ws = mean(c(mean(ws1), mean(ws2)))
                # F0 <- 0.8*FivalQ[] # Tr 10 anni to start :: hip constant params::
                F0 <- 0.9 # Tr 10 anni to start :: hip constant params::
                x0 <-  mean(Cs)*(log(mean(Ns)/(1-F0)))^(1/mean(ws))
                # x0 <-  mean(Cs2)*(log(mean(Ns)/(1-F0)))^(1/mean(ws2))
                optim = nleqslv(x0, myfun)
                # optim = nleqslv(x0, myfun)
                # print(optim$x)
                # print(dim(quants))
                quants[s, i] = optim$x + thresh
                termcd = optim$termcd
                if (termcd != 1){
                  print("comp_quant WARNING: nleqslv might not be converging")
                  print("If 6 Jacobian is singular, new guess")
                  print("termcd = ")
                  print(termcd)
                  nwarning_x0 = nwarning_x0 + 1
                }
                # epsi[s, i]   = (XivalQ[i] - quants[s, i])/XivalQ[i]
                epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
            }
            for (i in 1:Np){
                qfval = Xival[i] - thresh
                if (qfval <= 0){ # below threshold
                  # print("comp_quant ERROR: quantile below threshold")
                  cdfs[s, i] = 0.0
                  pdfs[s, i] = 0.0
                } else {
                cdfs[s, i] = hbev_mgw_cdf(qfval, C=Cs, W=ws, a=as,
                                          b=bs, alp = alps,  N=Ns)
                pdfs[s, i] = hbev_mgw_pdf(qfval, C=Cs, W=ws, a=as,
                                          b=bs, alp = alps,  N=Ns)
                }
            }
        }

  }      else if (substr(model, 1, 3) == 'wei') {
        thresh = model_fit$thresh_hbev
        for (s in 1:S){
            Cs = model_sim$Cgen[s,]
            ws = model_sim$wgen[s,]
            Ns = model_sim$Ngen[s,]
            for (i in 1:NpQ){
                myfun<-function(x) hbev_wei_cdf(x, C=Cs, W=ws, N=Ns) - FivalQ[i]
                # x0 = 8*mean(Cs)
                # x0 = 6*mean(Cs)/mean(ws)*(1 + log10(TrvalQ[i])/3)
                F0 <- 0.1 # Tr 10 anni to start :: hip constant params::
                x0 <-  mean(Cs)*(log(mean(Ns)/(1-F0)))^(1/mean(ws))
                optim = nleqslv(x0, myfun)
                # print(optim$x)
                # print(dim(quants))
                quants[s, i] = optim$x + thresh
                termcd = optim$termcd
                if (termcd != 1){
                  print("comp_quant WARNING: nleqslv might not be converging")
                  print("If 6 Jacobian is singular, new guess")
                  print("termcd = ")
                  print(termcd)
                  nwarning_x0 = nwarning_x0 + 1
                }
                # epsi[s, i]   = (XivalQ[i] - quants[s, i])/XivalQ[i]
                epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
            }
            for (i in 1:Np){
                qfval = Xival[i] - thresh
                if (qfval <= 0){ # below threshold
                  # print("comp_quant ERROR: quantile below threshold")
                  cdfs[s, i] = 0.0
                  pdfs[s, i] = 0.0
                } else {
                cdfs[s, i] = hbev_wei_cdf(qfval, C=Cs, W=ws, N=Ns)
                pdfs[s, i] = hbev_wei_pdf(qfval, C=Cs, W=ws, N=Ns)
                }
            }
      }
  }

  else if (substr(model, 1, 3) == 'gam') {
      thresh = model_fit$thresh_hbev
          for (s in 1:S){
            as = model_sim$agen[s,]
            bs = model_sim$bgen[s,]
            Ns = model_sim$Ngen[s,]
            for (i in 1:NpQ){
                myfun<-function(x) hbev_gam_cdf(x, a=as, b=bs, N=Ns) - FivalQ[i]
                x0 = 6/mean(bs)
                optim = nleqslv(x0, myfun)
                quants[s, i] = optim$x + thresh
                termcd = optim$termcd
                if (termcd != 1){
                  print("comp_quant WARNING: nleqslv might not be converging")
                  print("If 6 Jacobian is singular, new guess")
                  print("termcd = ")
                  print(termcd)
                  nwarning_x0 = nwarning_x0 + 1
                }
                # epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
                epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
            }
            for (i in 1:Np){
                qfval = Xival[i] - thresh
                if (qfval <= 0){ # below threshold
                  # print("comp_quant ERROR: quantile below threshold")
                  cdfs[s, i] = 0.0
                  pdfs[s, i] = 0.0
                } else {
                cdfs[s, i] = hbev_gam_cdf(qfval, a=as, b=bs, N=Ns)
                pdfs[s, i] = hbev_gam_pdf(qfval, a=as, b=bs, N=Ns)
                }
            }
          }
  }
    else if (substr(model, 1, 3) == 'gga') {
      # generalized gamma
      thresh = model_fit$thresh_hbev
          for (s in 1:S){
            as = model_sim$agen[s,]
            ws = model_sim$wgen[s,]
            Cs = model_sim$Cgen[s,]
            Ns = model_sim$Ngen[s,]
            for (i in 1:NpQ){
                myfun<-function(x) hbev_gga_cdf(x, C=Cs, W=ws, A=as, N=Ns) - FivalQ[i]
                # x0 = 6/mean(bs)
                #   nucc = mean(ws)
                #   sigmacc = mean(as/ws)
                #    mucc = mean(Cs*as/ws)
                # x0 = rmutil::qggamma(FivalQ[i], sigmacc, mucc, nucc)
                x0 = 100
                optim = nleqslv(x0, myfun)
                quants[s, i] = optim$x + thresh
                termcd = optim$termcd
                if (termcd != 1){
                  print("comp_quant WARNING: nleqslv might not be converging")
                  print("If 6 Jacobian is singular, new guess")
                  print("termcd = ")
                  print(termcd)
                  nwarning_x0 = nwarning_x0 + 1
                }
                # epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
                epsi[s, i]   = (quants[s, i] - XivalQ[i])/XivalQ[i]
            }
            for (i in 1:Np){
                qfval = Xival[i] - thresh
                if (qfval <= 0){ # below threshold
                  # print("comp_quant ERROR: quantile below threshold")
                  cdfs[s, i] = 0.0
                  pdfs[s, i] = 0.0
                } else {
                cdfs[s, i] = hbev_gga_cdf(qfval, C=Cs, W=ws, A=as, N=Ns)
                pdfs[s, i] = hbev_gga_pdf(qfval, C=Cs, W=ws, A=as, N=Ns)
                }
            }
          }
    }
  # compute uncertainty bands and average / mean values
  qmean =  apply(quants, 2, mean)
  qquant = apply(quants, 2, quantile, prob = c(0.05, 0.5, 0.95))
  fmean = apply(cdfs, 2, mean)
  fquant = apply(cdfs, 2, quantile, prob = c(0.05, 0.5, 0.95))
  dmean = apply(pdfs, 2, mean)
  dquant = apply(pdfs, 2, quantile, prob = c(0.05, 0.5, 0.95))

  # now compute some goodness of fit measures
    # such as CPOi and LPML, LPPD, and FSE

  # fse = sqrt(mean(epsi^2))

  # COMPUTE THE FRACTIONAL SQUARE ERROR ABOVE THE GIVEN TRMIN:
  fse_Tr = rep(NA, NpQ) # fractional seuare error for a given non exceedance frequency
  bias_Tr = rep(NA, NpQ) # bias
  width_Tr = rep(NA, NpQ) # width of credibility interval
  for (i in 1:NpQ){
    fse_Tr[i]  = sqrt( mean(epsi[,i]^2)) # averaged over the S draws from the posterior
    bias_Tr[i] = mean(epsi[, i])
    # width_Tr[i] = qupper[i] - qlower[i]
    width_Tr[i] = qquant[3, i] - qquant[1, i]
  }
  fse = mean(fse_Tr)
  mbias = mean(bias_Tr)
  mwidth = mean(width_Tr)

  log_lik = log(pdfs)

  infc = info_crit(pdfs, psis = FALSE) # compute information criteria

  # cpoi = rep(0, Np)
  # mean_pdfi = rep(0, Np)
  # for (i in 1:Np){
  #     cpoi[i] = 1/mean(1/pdfs[,i])
  #     mean_pdfi[i] = mean(pdfs[,i])
  # }
  # # lpml = mean(log(cpoi))       # should be sums
  # # lppd = mean(log(mean_pdfi))  # should be sums
  # lpml = mean(log(cpoi))       # should be sums
  # lppd = mean(log(mean_pdfi))  # should be sums
  #
  # # fse = sqrt(apply((epsi^2), 1, mean))
  #
  # # fse = sqrt(apply((epsi^2), mean))
  #
  # # compute log likelihood for maxima and PSIS - LOO:
  # log_lik = log(pdfs)
  #
  # # psis_loo = loo(log_lik)
  # # elpd_loo = psis_loo$estimates[1]/Np # divide by the sample size of validation
  # # p_loo = psis_loo$estimates[2]
  # psis_loo = 0
  # p_loo = 0
  # elpd_loo = 0



  return( list(qmean = qmean, qupper = qquant[3,],
               qlower = qquant[1, ], qmedian = qquant[2, ],
               fmean = fmean, fupper = fquant[3,],
               flower = fquant[1, ], fmedian = fquant[2, ],
               dmean = dmean, dupper = dquant[3,],
               dlower = dquant[1, ], dmedian = dquant[2, ],
               lpml = infc$lpml, lppd = infc$lppd,
               fse = fse, fse_Tr = fse_Tr,
               mbias = mbias, bias_Tr = bias_Tr,
               mwidth = mwidth, width_Tr = width_Tr,
               elpd_loo = infc$elpd_loo, p_loo = infc$p_loo,
               elpd_waic1 = infc$elpd_waic1, elpd_waic2 = infc$elpd_waic2,
               p_waic1 = infc$p_waic1, p_waic2 = infc$p_waic2,
               quants = quants, cdfs = cdfs, pdfs = pdfs, epsi=epsi,
               Fival=Fival, Xival=Xival, Trval=Trval,
               FivalQ=FivalQ, XivalQ=XivalQ, TrvalQ=TrvalQ,
               log_lik = log_lik,
               model = model, Mgen = Mgen,
               nwarning_x0 = nwarning_x0,
               thresh_hbev = model_fit$thresh_hbev
               ))
}




#' @export
hbev_max_ppd <- function(model_fit){
  "---------------------------------------------------------------------
   comments -> compute ppd for yearly maxima

   arguments:
   model_fit -> result from fit_ev_model

   output: max -> array of simulated annual maxima ndraws*nyears
  (the nyears is the number of years generated in the model)
  ----------------------------------------------------------------------"
  x_mod = unlist(strsplit(model_fit$model, "_"))[[1]]
  sim = rstan::extract(model_fit$model_fit)
  ndraws = dim(sim$Ngen)[1]
  nyears = dim(sim$Ngen)[2]
  max = matrix(0, nrow = ndraws, ncol = nyears)
  for (s in 1:ndraws){
    for (j in 1:nyears){
      if (x_mod == 'gam'){
        xij = rgamma(sim$Ngen[s, j], sim$agen[s, j], sim$bgen[s, j])
      } else if (x_mod == 'wei'){
        xij = rweibull(sim$Ngen[s, j], sim$wgen[s, j], sim$Cgen[s, j])
      }
      max[s, j] <- max(xij) + model_fit$thresh_hbev
    }
  }
  return(max)
}


#' @export
hbev_xij_ppd <- function(model_fit){
  "---------------------------------------------------------------------
   comments -> compute ppd for yearly maxima

   arguments:
   model_fit -> result from fit_ev_model

   output: list with array of simulated annual maxima ndraws*nyears
   and of daily data xij (ndraes*(nyears*Nt)) - including zeros
  (the nyears is the number of years generated in the model)
  ----------------------------------------------------------------------"
  x_mod = unlist(strsplit(model_fit$model, "_"))[[1]]
  sim = rstan::extract(model_fit$model_fit)
  ndraws = dim(sim$Ngen)[1]
  nyears = dim(sim$Ngen)[2]
  max = matrix(0, nrow = ndraws, ncol = nyears)
  daily_data = matrix(0, nrow = ndraws, ncol = nyears*model_fit$Nt)
  for (s in 1:ndraws){
    daily_data_s = matrix(0, nrow = nyears, ncol = model_fit$Nt)
    for (j in 1:nyears){

      if (sim$Ngen[s, j] > 0){ ### ADDED
      if (x_mod == 'gam'){
        xij = rgamma(sim$Ngen[s, j], sim$agen[s, j], sim$bgen[s, j])
      } else if (x_mod == 'wei'){
        xij = rweibull(sim$Ngen[s, j], sim$wgen[s, j], sim$Cgen[s, j])
      } else if (x_mod == 'mix'){
        xij <- rmixweibull(sim$Ngen[s, j], sim$Cgen1[s, j], sim$wgen1[s, j], sim$Cgen2[s, j],
         sim$wgen2[s, j], sim$alphap[s, j])
        # xij = rmixwei(sim$Ngen[s, j], sim$wgen1[s, j], sim$Cgen1[s, j],
        # sim$wgen1[s, j], sim$Cgen1[s, j], sim$alphas[s, j])
      }
      daily_data_s[j, 1:sim$Ngen[s, j]] = xij
      max[s, j] <- max(xij) + model_fit$thresh_hbev
    }
    daily_data[s, ] = as.vector(daily_data_s)
  } ### ADDED
  }
  return(list(max = max, xij = daily_data))
}


#' @export
plot_quants <- function(listq){
  "listq => list of outputs from compute_quants function"
  colors = c("black", "red", "blue", "green", "orange", "cyan")
  nmods = length(listq)
  # colors = colors[1:nmods]
  modelcolors = list()
  # print(nmods)
p1 <- ggplot() +
  geom_point(aes( listq[[1]]$TrvalQ, listq[[1]]$XivalQ),  size = 1.5, color = "black")
  for (i in 1:nmods){
    # print(i)

      p1 <- p1 + geom_line(aes_string(  listq[[i]]$TrvalQ, listq[[i]]$qmean),  linetype="solid", color=colors[[i]])
      p1 <- p1 + geom_line(aes_string(  listq[[i]]$TrvalQ, listq[[i]]$qupper), linetype="dashed",color=colors[[i]])
      p1 <- p1 + geom_line(aes_string(  listq[[i]]$TrvalQ, listq[[i]]$qlower), linetype="dashed",color=colors[[i]])
  }
  p1 = p1 + coord_trans(x="log10")
  p1 = p1 + labs(y = "Quantile [mm/day]", x= "Return Time [years]")
  p1 = p1 + theme_bw()
p2 <- ggplot() +
  geom_point(aes( listq[[1]]$Xival, listq[[1]]$Fival),  size = 1.5, color = "black")
    for (i in 1:nmods){
          modelcolors[[listq[[i]]$model]] = colors[i]
      p2 = p2 + annotate("text", x =150, y= 0.01 + 0.04*(i-1), label = listq[[i]]$model, color = colors[[i]])
      p2 <- p2 + geom_line(aes_string(  listq[[i]]$Xival, listq[[i]]$fmean),  linetype="solid", color=colors[[i]])
      p2 <- p2 + geom_line(aes_string(  listq[[i]]$Xival, listq[[i]]$fupper), linetype="dashed",color=colors[[i]])
      p2 <- p2 + geom_line(aes_string(  listq[[i]]$Xival, listq[[i]]$flower), linetype="dashed",color=colors[[i]])
    }

  p2 = p2 + coord_trans(x="log10")
  p2 = p2 + labs(x = "Quantile [mm/day]", y = "Non exceedance probability")
  p2 = p2 + theme_bw()
fig <- grid.arrange(p1, p2, ncol=2)
return(fig)
}



#' @export
plot_ppd_pdfs <- function(model_fit, datamat, ndraws2plot = 100){
 " plot posterior predictive distributions for N, maxima and daily data xij
   ndraws2plot -> number of draws to limit at for plotting only
  # Mgen must be equal in size to the observed sample
  datamat -> matrix with obs data (nyears*Nt) "

# test posterior predictive distribution for maxima and ordinary values
model_sim = rstan::extract(model_fit$model_fit)
simdata = hbev_xij_ppd(model_fit)

maxima = apply(datamat, 1, max)
Mgen = model_fit$Mgen
exceedances = pmax(datamat - model_fit$thresh_hbev, 0)
excesses = datamat[datamat > model_fit$thresh_hbev] - model_fit$thresh_hbev
nyears = dim(exceedances)[1]
N = rep(0, nyears)
for (i in 1:nyears){
  samplei = exceedances[i,]
  N[i] = length(samplei[samplei > 0])
}
fig1 <- ppc_dens_overlay(log(maxima[1:Mgen]),log(simdata$max[1:ndraws2plot,])) + labs(x="log(max xij)")
fig2 <- ppc_dens_overlay(N[1:(Mgen)],model_sim$Ngen[1:ndraws2plot,]) + labs(x="nj")

fig3 <- ggplot()
for (s in 1:ndraws2plot){
  sample_s = simdata$xij[s, ]
  wets_s = sample_s[sample_s > 0]
  # print(sum(is.na(log(wets_s))))
  fig3 = fig3 + geom_line( aes_string(log(wets_s)), stat = 'density', adjust=1,
                           col='lightblue', lwd = 0.2)

}
# print(max(log(excesses)))
# print(min(log(excesses)))
fig3 = fig3 + geom_line( aes(log(excesses)), stat = 'density', adjust=1, col = 'black', lwd=1)
# fig3 = fig3 +  coord_trans(x="log10")
fig3 = fig3 + scale_x_continuous(limits = c(-4, 6))
fig3 = fig3 + labs( y = "",  x= "log(xij)")
fig3 = fig3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.y=element_blank(), axis.ticks=element_blank())
fig2 = fig2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.y=element_blank(), axis.ticks=element_blank())
fig1 = fig1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.y=element_blank(), axis.ticks=element_blank())

fig <- grid.arrange(fig1, fig2, fig3,ncol=3)
return(fig)
}
