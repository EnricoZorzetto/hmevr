    hbev_priors <- function(model, empirical=FALSE,  data = FALSE, thresh_hbev=0,
                        draw_rng=FALSE, ndraws = 1000, priorpar = NULL){
  "---------------------------------------------------------------------
  Prior Elicitation for the hbev -family models.

  ARGUMENTS:
  model -> name of the model (e.g., 'wei_dyn_sim') see fit_ev_model
  empirical -> if TRUE, use information (sample mean only) to elicit priors
  data -> matrix of dim nyears*Nt with the data
  thresh_hbev -> threshold used in the hbev model (defualt 0)
  draw_rng -> if true sample from the prior distributions
  ndraws -> number of draws from the prior (default 1000)

  RETURN: list with following elements:
  prior_par -> list of arrays with the parameters of the prior distributions.
  prior_name -> list of string with the name of the prior distr for each par.
  prior_rng -> list with the arrays of draws for each parameter in the model
  ----------------------------------------------------------------------"

  print(sprintf('defining priors for the the %s model', model))
  mod_names = unlist(strsplit(model, "_"))
  n_mod = mod_names[[3]]
  x_mod = paste(mod_names[[1]], mod_names[[2]], sep="_")

  if (empirical == TRUE){
      exceedances = pmax(data-thresh_hbev, 0.0)
      meanwet = mean(exceedances[exceedances > 0])
      nyears = dim(exceedances)[1]
      N = rep(0, nyears)
      for (i in 1:nyears){
        samplei = exceedances[i,]
        N[i] = length(samplei[samplei > 0])
      }
      meanN = mean(N)
  } else {
    meanwet <- 10 # prior belief on charactristic event magnitude
    meanN   <- 100 # prior belief on average number events / year
  }

  prior_w <- 0.7 # prior belief on shape parameter (Wilson and Tuomi, 2005, GRL)
  prior_c = meanwet/gamma(1+1/prior_w)
  prior_a <- 1.1 # gamma shape -> slighter heavier tail than exponential

  # priorpar = list(inf_w0 = 100)
  # prior_c = 1
  # prior_w <- 0.66 # prior belief on shape parameter (Wilson and Tuomi, 2005, GRL)

  # set DEFAULT values for prior parameters:
  defp = list(
    prior_w = prior_w,
    prior_c = prior_c,

    # weibull models
      inf_mc0 = 10, # the larger the more informative the prior
      inf_sc0 = 10,
      inf_mw0 = 10,
      inf_sw0 = 10,
      exp_sc0 = 0.2, # as a fraction of expected mc
      exp_sw0 = 0.1, # as a fraction of expected mw
      inf_c0 = 1,
      inf_w0 = 10,
      exp_c0 =  1,
      exp_w0 =  1,
      exp_mc0 = 1, # as a fraction of expected mc
      exp_mw0 = 1, # as a fraction of expected mc
      # dynamic Gumbel model:
      # the other are like weibull dynamic
      gu_exp_sc0 = 0.3,
      gu_exp_sw0 = 0.1,

      # mixed weibull models:
      inf_mc10 = 100, # the larger the more informative the prior
      inf_sc10 = 100,
      inf_mw10 = 100,
      inf_sw10 = 100,
      exp_sc10 = 0.1, # as a fraction of expected mc
      exp_sw10 = 0.05, # as a fraction of expected mw
      exp_mc10 = 0.8,
      exp_mw10 = 1.0,

      inf_mc20 = 100, # the larger the more informative the prior
      inf_sc20 = 100,
      inf_mw20 = 100,
      inf_sw20 = 100,
      exp_sc20 = 0.1, # as a fraction of expected mc
      exp_sw20 = 0.05, # as a fraction of expected mw
      exp_mc20 = 1.2,
      exp_mw20 = 1.2,

      # rho0a = 60,
      # rho0b = 60,
      # mrho0a = 60,
      # mrho0b = 60,
      # srho0a = 10,
      # srho0b = 60,

      rho0a = 100,
      rho0b = 100,
      mrho0a = 200,
      mrho0b = 100,
      srho0a = 20,
      srho0b = 100,

      inf_c10 = 100*1,
      inf_w10 = 100*10,
      exp_c10 = 0.8,
      exp_w10 = 1,
      inf_c20 = 100*1,
      inf_w20 = 100*10,
      exp_c20 = 1.2,
      exp_w20 = 1.2,

      # n model prior parameters::
      pn0_an = 4,
      pn0_bn = 4,
      mu0_an = 10,
      mu0_bn = 0.1,
      omega0_an = 0.2,
      omega0_bn = 0.2,

    # hbev ar1 model
  infrho = 1000,
  m_nun = -0.75,
  m_nuc = 3,
  m_nuw = -0.1,
  s_nun = 0.1,
  s_nuc = 1,
  s_nuw = 0.1,
  a_sign = 30,
  b_sign = 0.05*30,
  a_sigc = 5,
  a_sigw = 5,
  b_sigc = 0.2,
  b_sigw = 0.2,

  # wei_dex_bin model
  inf_exploc = 10,
  inf_expsc = 200,
  exp_cloc = 0.8, # as a fraction of expected values - not needed for now
  exp_wloc = 0.9, # as a fraction of expected values - not needed for now
  exp_charvar_c = 0.2,
  exp_charvar_w = 0.1,

  # wei_dsn_bin model
  inf_snloc = 10,
  inf_snsc = 100,
  inf_snsh = 100,
  sn_cloc = 0.7, # as a fraction of expected values - not needed for now
  sn_wloc = 0.9, # as a fraction of expected values - not needed for now
  sn_charvar_c = 0.3,
  sn_charvar_w = 0.1,
  sn_shape_c = 4,
  sn_shape_w = 4


  )

  # if provided, substitute this values
  if (!is.null(priorpar)){
    namesp = names(priorpar)
    nnamesp = length(namesp)
    for (i in 1:nnamesp){
      mypname = namesp[i]
      if (mypname %in% names(defp)){
      # print(mypname)
      defp[[mypname]] <- priorpar[[mypname]]
      }
    }
  }

  # wei:
  # wei dynamic
  sc0prior = c(defp[['inf_sc0']], defp[['inf_sc0']]*defp[['exp_mc0']]*defp[['exp_sc0']]*defp[['prior_c']])       # a, b of an inverse gamma
  sw0prior = c(defp[['inf_sw0']], defp[['inf_sw0']]*defp[['exp_mw0']]*defp[['exp_sw0']]*defp[['prior_w']])       # a, b of an inverse gamma
  mc0prior = c(defp[['inf_mc0']], defp[['inf_mc0']]*defp[['exp_mc0']]*defp[['prior_c']])               # a, b of an inverse gamma
  mw0prior = c(defp[['inf_mw0']], defp[['inf_mw0']]*defp[['exp_mw0']]*defp[['prior_w']])               # a, b of an inverse gamma

  # wei static
  C0prior = c(defp[['inf_c0']]*defp[['exp_c0']]*defp[['prior_c']], defp[['inf_c0']])    # gamma a, b
  w0prior = c(defp[['inf_w0']]*defp[['exp_w0']]*defp[['prior_w']], defp[['inf_w0']])    # gamma a, b

  # # mix dyn / dyp / sta:
  sc10prior = c(defp[['inf_sc10']], defp[['inf_sc20']]*defp[['exp_mc20']]*defp[['exp_sc20']]*defp[['prior_c']])       # a, b of an inverse gamma
  sw10prior = c(defp[['inf_sw10']], defp[['inf_sw20']]*defp[['exp_mw20']]*defp[['exp_sw20']]*defp[['prior_w']])       # a, b of an inverse gamma
  mc10prior = c(defp[['inf_mc10']], defp[['inf_mc20']]*defp[['exp_mc20']]*defp[['prior_c']])               # a, b of an inverse gamma
  mw10prior = c(defp[['inf_mw10']], defp[['inf_mw20']]*defp[['exp_mw20']]*defp[['prior_w']])               # a, b of an inverse gamma
  sc20prior = c(defp[['inf_sc20']], defp[['inf_sc20']]*defp[['exp_mc20']]*defp[['exp_sc20']]*defp[['prior_c']])       # a, b of an inverse gamma
  sw20prior = c(defp[['inf_sw20']], defp[['inf_sw20']]*defp[['exp_mw20']]*defp[['exp_sw20']]*defp[['prior_w']])       # a, b of an inverse gamma
  mc20prior = c(defp[['inf_mc20']], defp[['inf_mc20']]*defp[['exp_mc20']]*defp[['prior_c']])               # a, b of an inverse gamma
  mw20prior = c(defp[['inf_mw20']], defp[['inf_mw20']]*defp[['exp_mw20']]*defp[['prior_w']])               # a, b of an inverse gamma

  rhoj0prior = c(defp[['rho0a']], defp[['rho0b']])     # a, b of a beta
  mrho0prior = c(defp[['mrho0a']], defp[['mrho0b']])  # a, b of a beta
  srho0prior = c(defp[['srho0a']], defp[['srho0b']])  # a, b of a beta

  C10prior = c(defp[['inf_c10']]*defp[['exp_c10']]*defp[['prior_c']], defp[['inf_c10']])    # gamma a, b
  w10prior = c(defp[['inf_w10']]*defp[['exp_w10']]*defp[['prior_w']], defp[['inf_w10']])    # gamma a, b
  C20prior = c(defp[['inf_c20']]*defp[['exp_c20']]*defp[['prior_c']], defp[['inf_c20']])    # gamma a, b
  w20prior = c(defp[['inf_w20']]*defp[['exp_w20']]*defp[['prior_w']], defp[['inf_w20']])    # gamma a, b

  # wei_ar1_bin model:


rhon0prior  = c(defp[['infrho']], defp[['infrho']]) # a, b of a beta
rhoc0prior  = c(defp[['infrho']], defp[['infrho']]) # a, b of a beta
rhow0prior  = c(defp[['infrho']], defp[['infrho']]) # a, b of a beta

nun0prior = c(defp[['m_nun']], defp[['s_nun']]) # mu, sigma of a normal
nuc0prior = c(defp[['m_nuc']], defp[['s_nuc']]) # mu, sigma of a normal
nuw0prior = c(defp[['m_nuw']], defp[['s_nuw']]) # mu, sigma of a normal

sign0prior = c(defp[['a_sign']], defp[['b_sign']]) # a, b of a gamma
sigc0prior = c(defp[['a_sigc']], defp[['b_sigc']]) # a, b of a gamma
sigw0prior = c(defp[['a_sigw']], defp[['b_sigw']]) # a, b of a gamma

# wei_dex_bin model::
cexploc0prior = c(defp[['inf_exploc']]*defp[['prior_c']]*defp[['exp_cloc']], defp[['inf_exploc']])    # gamma a, b
wexploc0prior = c(defp[['inf_exploc']]*defp[['prior_w']]*defp[['exp_wloc']], defp[['inf_exploc']])    # gamma a, b
cexpsc0prior  = c(defp[['inf_expsc']]*defp[['exp_charvar_c']]*defp[['prior_c']], defp[['inf_expsc']])    # gamma a, b
wexpsc0prior  = c(defp[['inf_expsc']]*defp[['exp_charvar_w']]*defp[['prior_w']], defp[['inf_expsc']])    # gamma a, b

# wei_dsn_bin model::
csnloc0prior = c(defp[['inf_snloc']]*defp[['prior_c']]*defp[['sn_cloc']], defp[['inf_snloc']])    # gamma a, b
wsnloc0prior = c(defp[['inf_snloc']]*defp[['prior_w']]*defp[['sn_wloc']], defp[['inf_snloc']])    # gamma a, b
csnsc0prior  = c(defp[['inf_snsc']]*defp[['sn_charvar_c']]*defp[['prior_c']], defp[['inf_snsc']])    # gamma a, b
wsnsc0prior  = c(defp[['inf_snsc']]*defp[['sn_charvar_w']]*defp[['prior_w']], defp[['inf_snsc']])    # gamma a, b
csnsh0prior  = c(defp[['inf_snsh']]*defp[['sn_shape_c']], defp[['inf_snsh']])    # gamma a, b
wsnsh0prior  = c(defp[['inf_snsh']]*defp[['sn_shape_w']], defp[['inf_snsh']])    # gamma a, b

  # wei-dgu:: Weibull dynamic with latent gumbels
  gusc0prior = c(defp[['inf_sc0']], defp[['inf_sc0']]*defp[['exp_mc0']]*defp[['gu_exp_sc0']]*defp[['prior_c']])       # a, b of an inverse gamma
  gusw0prior = c(defp[['inf_sw0']], defp[['inf_sw0']]*defp[['exp_mw0']]*defp[['gu_exp_sw0']]*defp[['prior_w']])       # a, b of an inverse gamma
  gumc0prior = c(defp[['inf_mc0']], defp[['inf_mc0']]*defp[['exp_mc0']]*defp[['prior_c']])               # a, b of an inverse gamma
  gumw0prior = c(defp[['inf_mw0']], defp[['inf_mw0']]*defp[['exp_mw0']]*defp[['prior_w']])               # a, b of an inverse gamma


  # inf_exploc = 10,
  # inf_expsc = 200,
  # exp_cloc = 0.9, # as a fraction of expected values - not needed for now
  # exp_wloc = 0.9, # as a fraction of expected values - not needed for now
  # exp_charvar_c = 0.2,
  # exp_charvar_w = 0.1,
  #
  # # wei_dsn_bin model
  # inf_snloc = 10,
  # inf_snsc = 200,
  # inf_snsh = 200,
  # sn_cloc = 0.9, # as a fraction of expected values - not needed for now
  # sn_wloc = 0.9, # as a fraction of expected values - not needed for now
  # sn_charvar_c = 0.2,
  # sn_charvar_w = 0.1,
  # sn_shape_c = 10,
  # sn_shape_w = 10

  # wei lornormal
  prior_logc = log(meanwet/gamma(1+1/prior_w))
  prior_logw = log(prior_w)
  muc0prior = c(prior_logc, 1)   # normal mean, sigma2
  muw0prior = c(prior_logw, 0.3) # normal mean, sigma2
  sigmac0prior = c(5, 0.2) # parameters a, b of an inverse gamma
  sigmaw0prior = c(5, 0.2) # parameters a, b of an inverse gamma

  # gamma
  prior_b = prior_a/meanwet
  ma0prior = c(200, 200*prior_a) # (a, b) of inverse gamma distributions:
  sa0prior = c(0.5, 0.5*prior_a/50)
  mb0prior = c(20, 20*prior_b)
  # sb0prior = c(10, 0.1*20)
  sb0prior = c(2, 2/10*prior_b)
  a0prior = c(10*prior_a, 10 ) # gamma a, b
  b0prior = c(30*prior_b, 30)    # gamma a, b

  # priors for the N model::
  # pn0prior = c(2,2) # a, b of a beta
  # mu0prior = c(10, 0.1) # a, b of a gammma
  # omega0prior = c(0.2, 0.2) # a, b of a gammma

  pn0prior = c(defp[['pn0_an']],defp[['pn0_bn']]) # a, b of a beta
  mu0prior = c(defp[['mu0_an']],defp[['mu0_bn']]) # a, b of a beta
  omega0prior = c(defp[['omega0_an']],defp[['omega0_bn']]) # a, b of a beta

  prior_par = list()
  prior_name = list()
  # prior for the Nj model
  if (n_mod == 'bbn'){
    prior_par[['mun0prior']] = mu0prior
    prior_par[['omegan0prior']] = omega0prior
    prior_name[['mun0prior']] = 'gamma'
    prior_name[['omegan0prior']] = 'gamma'
  } else if (n_mod == 'bin'){
     prior_par[['pn0prior']]  = pn0prior
    prior_name[['pn0prior']] = 'beta'
  }
  # prior for the xij model
  if (x_mod == 'wei_dyn' || x_mod == 'gga_dyn' || x_mod == 'gga_sta' || x_mod == 'cop_dyn' || x_mod == 'wei_dyl' || x_mod == 'wei_dig'){
    prior_par[['sc0prior']] = sc0prior
    prior_par[['mc0prior']] = mc0prior
    prior_par[['sw0prior']] = sw0prior
    prior_par[['mw0prior']] = mw0prior
    prior_name[['sc0prior']] = 'invgamma'
    prior_name[['mc0prior']] = 'invgamma'
    prior_name[['sw0prior']] = 'invgamma'
    prior_name[['mw0prior']] = 'invgamma'
  }
   else if (x_mod == 'wei_dgu'){ # new asymmetric model
    prior_par[['sc0prior']] = gusc0prior
    prior_par[['mc0prior']] = gumc0prior
    prior_par[['sw0prior']] = gusw0prior
    prior_par[['mw0prior']] = gumw0prior
    prior_name[['sc0prior']] = 'invgamma'
    prior_name[['mc0prior']] = 'invgamma'
    prior_name[['sw0prior']] = 'invgamma'
    prior_name[['mw0prior']] = 'invgamma'

 } else if (x_mod == 'wei_dex'){ # new asymmetric model
    prior_par[['cexploc0prior']] = cexploc0prior
    prior_par[['wexploc0prior']] = wexploc0prior
    prior_par[['cexpsc0prior']] =  cexpsc0prior
    prior_par[['wexpsc0prior']] =  wexpsc0prior
    prior_name[['cexploc0prior']] = 'gamma'
    prior_name[['wexploc0prior']] = 'gamma'
    prior_name[['cexpsc0prior']]  = 'gamma'
    prior_name[['wexpsc0prior']]  = 'gamma'
  } else if (x_mod == 'wei_dsn'){ # new asymmetric model
    prior_par[['csnloc0prior']] = csnloc0prior
    prior_par[['wsnloc0prior']] = wsnloc0prior
    prior_par[['csnsc0prior']] =  csnsc0prior
    prior_par[['wsnsc0prior']] =  wsnsc0prior
    prior_par[['csnsh0prior']] =  csnsh0prior
    prior_par[['wsnsh0prior']] =  wsnsh0prior
    prior_name[['csnloc0prior']] = 'gamma'
    prior_name[['wsnloc0prior']] = 'gamma'
    prior_name[['csnsc0prior']]  = 'gamma'
    prior_name[['wsnsc0prior']]  = 'gamma'
    prior_name[['csnsc0prior']]  = 'gamma'
    prior_name[['wsnsc0prior']]  = 'gamma'

  } else if (x_mod == 'wei_sta'){
    prior_par[['C0prior']] = C0prior
    prior_par[['w0prior']] = w0prior
    prior_name[['C0prior']] = 'gamma'
    prior_name[['w0prior']] = 'gamma'
  } else if (x_mod == 'wei_dsc'){
    prior_par[['sc0prior']] = sc0prior
    prior_par[['mc0prior']] = mc0prior
    prior_par[['w0prior']] = w0prior
    prior_name[['sc0prior']] = 'invgamma'
    prior_name[['mc0prior']] = 'invgamma'
    prior_name[['w0prior']] = 'gamma'
  # } else if (x_mod == 'wei_dyl'){
  #   prior_par[['muc0prior']] = muc0prior
  #   prior_par[['muw0prior']] = muw0prior
  #   prior_par[['sigmac0prior']] = sigmac0prior
  #   prior_par[['sigmaw0prior']] = sigmaw0prior
  #   prior_name[['muc0prior']] = 'normal'
  #   prior_name[['muw0prior']] = 'normal'
  #   prior_name[['sigmac0prior']] = 'invgamma'
  #   prior_name[['sigmaw0prior']] = 'invgamma'
  } else if (x_mod == 'wei_dsl'){
    prior_par[['muc0prior']] = muc0prior
    prior_par[['sigmac0prior']] = sigmac0prior
    prior_par[['w0prior']] = w0prior
    prior_name[['muc0prior']] = 'normal'
    prior_name[['sigmac0prior']] = 'invgamma'
    prior_name[['w0prior']] = 'gamma'
  } else if (x_mod == 'gam_dyn'){
    prior_par[['ma0prior']] = ma0prior
    prior_par[['mb0prior']] = mb0prior
    prior_par[['sa0prior']] = sa0prior
    prior_par[['sb0prior']] = sb0prior
    prior_name[['ma0prior']] = 'invgamma'
    prior_name[['mb0prior']] = 'invgamma'
    prior_name[['sa0prior']] = 'invgamma'
    prior_name[['sb0prior']] = 'invgamma'
  } else if (x_mod == 'gam_dsc'){
    prior_par[['mb0prior']] = mb0prior
    prior_par[['sb0prior']] = sb0prior
    prior_par[['a0prior']] = a0prior
    prior_name[['mb0prior']] = 'invgamma'
    prior_name[['sb0prior']] = 'invgamma'
    prior_name[['a0prior']] = 'gamma'
  } else if (x_mod == 'gam_sta'){
    prior_par[['a0prior']] = a0prior
    prior_par[['b0prior']] = b0prior
    prior_name[['a0prior']] = 'gamma'
    prior_name[['b0prior']] = 'gamma'

    } else if (x_mod == 'mix_dyn' || x_mod == 'mix_fix'){

    prior_par[['sc10prior']] = sc10prior
    prior_par[['mc10prior']] = mc10prior
    prior_par[['sw10prior']] = sw10prior
    prior_par[['mw10prior']] = mw10prior
    prior_name[['sc10prior']] = 'invgamma'
    prior_name[['mc10prior']] = 'invgamma'
    prior_name[['sw10prior']] = 'invgamma'
    prior_name[['mw10prior']] = 'invgamma'

    prior_par[['sc20prior']] = sc20prior
    prior_par[['mc20prior']] = mc20prior
    prior_par[['sw20prior']] = sw20prior
    prior_par[['mw20prior']] = mw20prior
    prior_name[['sc20prior']] = 'invgamma'
    prior_name[['mc20prior']] = 'invgamma'
    prior_name[['sw20prior']] = 'invgamma'
    prior_name[['mw20prior']] = 'invgamma'

    prior_par[['rhoj0prior']] = rhoj0prior
    prior_name[['rhoj0prior']] = 'beta'

   }  else if (x_mod == 'mgw_dyn'){ # mixture gamma / Weibull

    prior_par[['sc0prior']] = sc0prior
    prior_par[['mc0prior']] = mc0prior
    prior_par[['sw0prior']] = sw0prior
    prior_par[['mw0prior']] = mw0prior
    prior_name[['sc0prior']] = 'invgamma'
    prior_name[['mc0prior']] = 'invgamma'
    prior_name[['sw0prior']] = 'invgamma'
    prior_name[['mw0prior']] = 'invgamma'

    prior_par[['ma0prior']] = ma0prior
    prior_par[['mb0prior']] = mb0prior
    prior_par[['sa0prior']] = sa0prior
    prior_par[['sb0prior']] = sb0prior
    prior_name[['ma0prior']] = 'invgamma'
    prior_name[['mb0prior']] = 'invgamma'
    prior_name[['sa0prior']] = 'invgamma'
    prior_name[['sb0prior']] = 'invgamma'

    prior_par[['rhoj0prior']] = rhoj0prior
    prior_name[['rhoj0prior']] = 'beta'

  } else if (x_mod == 'mix_dyp'){

    prior_par[['C10prior']] = C10prior
    prior_par[['w10prior']] = w10prior
    prior_name[['C10prior']] = 'gamma'
    prior_name[['w10prior']] = 'gamma'
    prior_par[['C20prior']] = C20prior
    prior_par[['w20prior']] = w20prior
    prior_name[['C20prior']] = 'gamma'
    prior_name[['w20prior']] = 'gamma'
    prior_par[['mrho0prior']] = mrho0prior
    prior_name[['mrho0prior']] = 'beta'
    prior_par[['srho0prior']] = srho0prior
    prior_name[['srho0prior']] = 'beta'

  } else if (x_mod == 'mix_sta'){

    prior_par[['C10prior']] = C10prior
    prior_par[['w10prior']] = w10prior
    prior_name[['C10prior']] = 'gamma'
    prior_name[['w10prior']] = 'gamma'
    prior_par[['C20prior']] = C20prior
    prior_par[['w20prior']] = w20prior
    prior_name[['C20prior']] = 'gamma'
    prior_name[['w20prior']] = 'gamma'
    prior_par[['rhoj0prior']] = rhoj0prior
    prior_name[['rhoj0prior']] = 'beta'
  }  else if (x_mod == 'wei_ar1'){

    prior_par[['rhon0prior']]  =rhon0prior
    prior_par[['rhoc0prior']]  =rhoc0prior
    prior_par[['rhow0prior']]  =rhow0prior
    prior_par[['nun0prior']]  =nun0prior
    prior_par[['nuc0prior']]  =nuc0prior
    prior_par[['nuw0prior']]  =nuw0prior
    prior_par[['sign0prior']]  =sign0prior
    prior_par[['sigc0prior']]  =sigc0prior
    prior_par[['sigw0prior']]  =sigw0prior

    prior_name[['rhon0prior']]  = 'beta'
    prior_name[['rhoc0prior']]  = 'beta'
    prior_name[['rhow0prior']]  = 'beta'
    prior_name[['nun0prior']]  = 'normal'
    prior_name[['nuc0prior']]  = 'normal'
    prior_name[['nuw0prior']]  = 'normal'
    prior_name[['sign0prior']]  = 'gamma'
    prior_name[['sigc0prior']]  = 'gamma'
    prior_name[['sigw0prior']]  = 'gamma'

  }
  prior_rng = list()
  if (draw_rng==TRUE){
      # print(prior_par)
    for (namei in names(prior_par)){
      # print(namei)
      namep = unlist(strsplit(namei, "0"))
      # print(namep)
      name = prior_name[[namei]]
      # print(name)
      # print(prior_name)
      if (name == 'gamma'){
        prior_rng[[namei]] = rgamma(ndraws, prior_par[[namei]][1], prior_par[[namei]][2])
      } else if (name == 'invgamma'){
        prior_rng[[namei]] = extraDistr::rinvgamma(ndraws, prior_par[[namei]][1], prior_par[[namei]][2])
      } else if (name == 'normal'){
        prior_rng[[namei]] = rnorm(ndraws, prior_par[[namei]][1], prior_par[[namei]][2])
      } else if (name == 'beta'){
        prior_rng[[namei]] = rbeta(ndraws, prior_par[[namei]][1], prior_par[[namei]][2])
      }
    }
  }
  prior = list(prior_par=prior_par,
               prior_name = prior_name,
               prior_rng=prior_rng)
  return(prior)
}


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


#' @export
hbev_wei_cdf <- function(x, C = C, W = W, N = N){
  # only for scalar x!
  cdf = mean(   pweibull(x, W, C )^N  )
return(cdf)
}


#' @export
hbev_wei_pdf <- function(x, C = C, W = W, N = N){
  # only for scalar x!
  pdf = mean( N*pweibull(x, W, C)^(N-1)*dweibull(x, W, C))
return(pdf)
}


#' @export
hbev_wei_quant <- function(p, C = C, W = W, N = N){
  # only for scalar x!
  myfunq <-function(x) hbev_wei_cdf(x, C=C, W=W, N=N) - p
  F0 <- p
  x0 <- mean(C)*(log(mean(N)/(1-F0)))^(1/mean(W))
  optim <- nleqslv(x0, myfunq)
  quant <- optim$x
  return(quant)
}


#' @export
dhbev <- function(x, ptrue = ptrue, ntrue = ntrue,
                  pdistr = 'wei_dgu' , ndistr = 'bin',
                  Nt=366, nsamples = 50){
  # compute theoretical hbev cdf for the values in x (or Tr)
  # if nsamples = 1, draw from the static model with C=mc, W=mw

  if (ndistr == 'bin'){
      Ngen = rbinom(nsamples, Nt, ntrue$pn)
  } else if (ndistr == 'bbn') {
      Ngen = extraDistr::rbbinom(nsamples, Nt, alpha = ntrue$an, beta = ntrue$bn)
  }

  if (pdistr == 'wei_dgu'){
      Cgen = extraDistr::rgumbel(nsamples, mu = mctrue, sigma = sctrue)
      Wgen = extraDistr::rgumbel(nsamples, mu = mwtrue, sigma = swtrue)
  } else if (pdistr == 'wei_dyn'){
      Cgen = rgamma(nsamples, ptrue$ac, ptrue$bc)
      Wgen = rgamma(nsamples, ptrue$aw, ptrue$bw)
  } else if (pdistr == 'wei'){
      Cgen = rep(nsamples, ptrue$C)
      Wgen = rep(nsamples, ptrue$w)
  }
    Nx = length(x)
    res = rep(0, Nx)
    for (i in 1:Nx){
      res[i] <-  hbev_wei_pdf(x[i], C=Cgen, W=Wgen, N=Ngen)
    }
    return(res)
}

#' @export
rhbev <- function(n, ptrue = ptrue, ntrue = ntrue,
                  pdistr = 'wei_dgu' , ndistr = 'bin',
                  Nt=366, nsamples = 50){
  # draw n annual maxima from hbev distribution
  u <- runif(n)
  res <- qhbev(u, ptrue = ptrue, ntrue = ntrue,
               pdistr = pdistr, ndistr = ndistr,
               Nt = Nt, nsamples = nsamples)

  return(res)
}

#' @export
phbev <- function(x, ptrue = ptrue, ntrue = ntrue,
                  pdistr = 'wei_dgu' , ndistr = 'bin',
                  Nt=366, nsamples = 50, Tr = FALSE){
  # compute theoretical hbev cdf for the values in x (or Tr)
  # if nsamples = 1, draw from the static model with C=mc, W=mw

  if (ndistr == 'bin'){
      Ngen = rbinom(nsamples, Nt, ntrue$pn)
  } else if (ndistr == 'bbn') {
      Ngen = extraDistr::rbbinom(nsamples, Nt, alpha = ntrue$an, beta = ntrue$bn)
  }
  if (pdistr == 'wei_dgu'){
      Cgen = extraDistr::rgumbel(nsamples, mu = mctrue, sigma = sctrue)
      Wgen = extraDistr::rgumbel(nsamples, mu = mwtrue, sigma = swtrue)
  } else if (pdistr == 'wei_dyn'){
      Cgen = rgamma(nsamples, ptrue$ac, ptrue$bc)
      Wgen = rgamma(nsamples, ptrue$aw, ptrue$bw)
  } else if (pdistr == 'wei'){
      Cgen = rep(nsamples, ptrue$C)
      Wgen = rep(nsamples, ptrue$w)
  }
    Nx = length(x)
    res = rep(0, Nx)
    for (i in 1:Nx){
      res[i] <-  hbev_wei_cdf(x[i], C=Cgen, W=Wgen, N=Ngen)
    }
    if (Tr == TRUE){
      res = 1/(1-res)
    }
    return(res)
}

#' @export
qhbev <- function(fi, ptrue = ptrue, ntrue = ntrue,
                  pdistr = 'wei_dgu' , ndistr = 'bin',
                  Nt=366, nsamples = 50, Tr = FALSE){
  # compute theoretical hbev quantiles for the values in fi.
  # values in fi must be non exceedance probabilities or return times
  # (second option only if fromTr= FALSE)
  # if nsamples = 1, draw from the static model with C=mc, W=mw
  if (Tr == FALSE){
    Fi = fi
    if (!(max(Fi)<1)){
      print('phbev ERROR: check if Fi or Tr')
    }
  } else {
    Tr = fi
    Fi = 1-1/Tr
    if (!(min(Tr)>1)){
      print('phbev ERROR: check if Fi or Tr')
    }
  }
  if (ndistr == 'bin'){
      Ngen = rbinom(nsamples, Nt, ntrue$pn)
  } else if (ndistr == 'bbn') {
      Ngen = extraDistr::rbbinom(nsamples, Nt, alpha = ntrue$an, beta = ntrue$bn)
  }
  if (pdistr == 'wei_dgu'){
      Cgen = extraDistr::rgumbel(nsamples, mu = mctrue, sigma = sctrue)
      Wgen = extraDistr::rgumbel(nsamples, mu = mwtrue, sigma = swtrue)
  } else if (pdistr == 'wei_dyn'){
      Cgen = rgamma(nsamples, ptrue$ac, ptrue$bc)
      Wgen = rgamma(nsamples, ptrue$aw, ptrue$bw)
  } else if (pdistr == 'wei'){
      Cgen = rep(nsamples, ptrue$C)
      Wgen = rep(nsamples, ptrue$w)
  }
    Nx = length(fi)
    res = rep(0, Nx)
    for (i in 1:Nx){
      res[i] <- hbev_wei_quant(Fi[i], C=Cgen, W=Wgen, N=Ngen)
    }
    return(res)
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



#' @export
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
