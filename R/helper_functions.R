####
### Functions that can be used by any pomp models
####

pull_outbreak_data <- function(outbreak=c("Yambuku","Kikwit","Mweka2007","Mweka2008","Isiro","Boende", "Equator"),
                               data = NULL){
  ## Sets up the data based on the given outbreak name
  outbrk <- match.arg(outbreak)
  data %>%
    ungroup() %>%
    filter(outbreak==outbrk) %>%
    select(times,cases) -> data
  
  data <- as.data.frame(data)  
}

randomize_parms <- function(pomp_obj, model_used){
  ## Randomizes parameters for mif2 starting iteration
  if(model_used == "combo"){
    parm_box <- rbind(
      beta0 = c(1e-4, 5),
      p0 = c(1e-8, 1),
      k = c(1e-5, 1),
      tau1 = c(1, max(pomp_obj@times))
    )
  } else if (model_used == "ss"){
    parm_box <- rbind(
      beta0 = c(1e-4, 10),
      p0 = c(1e-8, 1)
    )
  } else if (model_used == "int"){
    parm_box <- rbind(
      beta0 = c(1e-4, 10),
      k = c(1e-5, 1),
      tau1 = c(1, max(pomp_obj@times))
    )
  } else{
    stop("Do not recognize the model_used parameter. Either a typo, or need to add new parameter randomization for new model")
  }
  rand_parms <- apply(parm_box, 1, function(x) exp(runif(1, log(x[1]), log(x[2]))))
  sub_parms(rand_parms, coef(pomp_obj))
}

mif2_single_run <- function(parms,
                            pomp_obj, 
                            settings) {
  ## Function to run single mif
  ## With the specified parameters
  
  num_particles <- settings$mif_nparticles
  num_mif_iterations <- settings$mif_niter
  
  mif2(pomp_obj,
       start=parms,
       Np=num_particles,
       Nmif=num_mif_iterations,
       cooling.type="geometric",
       cooling.fraction.50=0.6,
       transform=TRUE,
       rw.sd= settings$parms_sd)
}

mif2_multirun <- function(pomp_obj, 
                          settings,
                          refresh = FALSE,
                          n = 10){
  ## Function to run multiple mifs
  
  ## First setup the location for the .rda file
  dest <- paste0("data_produced/outbreak_rda/", settings$model_used,"_", settings$outbreak, "_mif.Rda")
  
  ## Generate list of parms for running models
  parm_list <- n %>% rerun(randomize_parms(pomp_obj, settings$model_used))
  # parm_list <- n %>% rerun(coef(pomp_obj))

  if(file.exists(dest) & refresh == FALSE){
    load(dest)
  } else{
    mifs_global <- parm_list %>% map(mif2_single_run, 
                                     pomp_obj=pomp_obj, 
                                     settings=settings)
    mifs_global <- as(mifs_global, "mif2List")
    save(mifs_global, file = dest)
  }
  
  return(mifs_global)
}
  
find_max_ll_mif <- function(mif2_list){
  ## Returns the MLE mif run
  mif2_list[[which.max(map(mif2_list, logLik) %>% flatten_dbl())]]
}


get_parm_bounds <- function(parm, mif2_obj, settings){
  if (settings$intensive==FALSE) {
    if(parm == "p0"){
      lower <- mif2_obj@params[parm] / 5
      upper <- 0.999
    } else{
        lower <- mif2_obj@params[parm] / 5
        upper <- mif2_obj@params[parm] * 5  
        }
  } else {
    bounds <- unname(unlist(settings$bounds))
    if (parm == "beta0") {
      lower_ind = 1
      upper_ind = 2
    } else if (parm == "p0" | parm == "k") {
      lower_ind = 3
      upper_ind = 4
    } else if (parm == "tau1") {
      lower_ind = 5
      upper_ind = 6
    }
    lower = bounds[lower_ind]
    upper = bounds[upper_ind]
  }

  c(lower, upper)
}

get_list_parm_bounds <- function(parms, mif2_obj, settings){
  parms %>% map(get_parm_bounds, mif2_obj=mif2_obj, settings) %>% setNames(parms)
}

get_parm_slice_input <- function(parm, parm_bounds, settings){
  rep(seq(from=parm_bounds[[parm]][1], to=parm_bounds[[parm]][[2]], length=settings$slice_length), each=settings$slice_reps)
}

get_lik_slice <- function(mif2_obj, settings){
  ## Gives the slice design used for the likelihood profile
  est_parms <- settings$est_parms
  model_used <-settings$model_used
  
  parm_bounds <- get_list_parm_bounds(est_parms, mif2_obj, settings)

  if(model_used == "combo"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"],coef(mif2_obj)["p0"], 
               coef(mif2_obj)["tau1"], coef(mif2_obj)["k"],
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      p0 = get_parm_slice_input("p0", parm_bounds, settings),
      k = get_parm_slice_input("k", parm_bounds, settings),
      tau1 = get_parm_slice_input("tau1", parm_bounds, settings)
    )
  } else if (model_used == "ss"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"],coef(mif2_obj)["p0"], 
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      p0 = get_parm_slice_input("p0", parm_bounds, settings)
    )
  } else if (model_used == "int"){
    sliceDesign(
      center=c(coef(mif2_obj)["beta0"], 
               coef(mif2_obj)["tau1"], coef(mif2_obj)["k"],
               sigma = 1/9.312799, gamma = 1/7.411374, ff = 49/69),
      beta0 = get_parm_slice_input("beta0", parm_bounds, settings),
      k = get_parm_slice_input("k", parm_bounds, settings),
      tau1 = get_parm_slice_input("tau1", parm_bounds, settings)
    )
  } else{
    stop("Do not recognize the model_used parameter. Either a typo, or need to add new parameter randomization for new model")
  }
}

get_single_lik <- function(pomp_obj, num_particles, ...){
  ## Completes the particle filtering for a certain set of parameters
  ## The parameters are specified by the ... and called through the prof_like_run function
  
  pfilter(pomp_obj, params=unlist(list(...)), Np = num_particles)
}

prof_lik_run <- function(mif2_obj, est_parms, settings, outbreak, model_used, refresh=FALSE){
  ## This function calculates the likelihood profile for the MLE parameters, and
  ## Returns a data_frame with the results (and saves it for later)
  
  dest <- paste0("data_produced/outbreak_rda/", settings$model_used, "_", settings$outbreak, "_prof.rda")
  if(file.exists(dest) & refresh == FALSE){
    load(dest)
  } else{
    parms <- get_lik_slice(mif2_obj, settings) %>% as_data_frame
    

    
    slice_runs <- parms %>% 
      pmap(.f = get_single_lik, mif2_obj, settings$prof_lik_nparticles) 
    
    prof_lik <- parms %>% mutate(ll = map(slice_runs, logLik) %>% unlist)
    save(prof_lik, file = dest)
  }
  prof_lik
}

plot_prof_lik <- function(df, max_mif, settings){
  max_df <- as_data_frame(matrix(coef(max_mif)[settings$est_parms], ncol = length(settings$est_parms)))
  names(max_df) <- settings$est_parms
  max_df <- max_df %>% gather(key, val)
  
  df %>% gather(key,val, settings$est_parms) %>% 
    filter(key == slice) %>% 
    group_by(key, val) %>% 
    summarize(avg_ll = mean(ll)) %>%
    mutate(avg_ll = avg_ll - max(avg_ll)) %>%
    # mutate(avg_ll = ll-max(ll)) %>%
    ggplot(aes(val, avg_ll)) + 
    facet_wrap(~key, scales="free_x") + 
    geom_point(alpha=1,size=.5) +
    coord_cartesian(ylim = c(-30,0)) +
    geom_vline(data=max_df, aes(xintercept=val), color = "darkred", lty=2) -> pl_plot
  paste0("data_produced/figs/", settings$model_used, "_", settings$outbreak, "_prof_plot.pdf") -> dest
  save_plot(dest,pl_plot,base_height = 4, base_aspect_ratio = 1.6)
}



calc_k <- function(parms) {
  p <- unname(parms["p0"])
  R <- calc_rnot(parms)
  p <- 1-p
  f <- function(x)  ((1+(R/x))^(-x)-p)
  uniroot(f, lower=0.00000001, upper= 2)$root
}

conf_interval <- function(prof_lik, max_mif, settings) {
  prof_lik %>%
    gather(key,val, settings$est_parms) %>%
    filter(key == slice) -> prof_lik
  if (settings$model_used == "ss") {
    p_bounds <- conf_int_calc('p0', coef(max_mif)["p0"], prof_lik, -1.96, settings)
    b_bounds <- conf_int_calc('beta0', coef(max_mif)["beta0"], prof_lik, -1.96, settings)
    bounds <- c(p_bounds,b_bounds)
  }
  else {
    b_bounds <- conf_int_calc('beta0', coef(max_mif)["beta0"], prof_lik, -1.96, settings)
    k_bounds <- conf_int_calc('k', coef(max_mif)["k"], prof_lik, -1.96, settings)
    t_bounds <- conf_int_calc('tau1', coef(max_mif)["tau1"], prof_lik, -1.96, settings)
    bounds <- c(b_bounds,k_bounds,t_bounds)
  }
  return(bounds)
}

conf_int_calc <- function(param, mle, df, likelihood_cutoff, settings) {
  require(scam)

  df <- df %>% filter(slice == param)
  vals_to_consider <- df %>% group_by(val) %>% 
    summarize(avg_ll = mean(ll)) %>% 
    mutate(avg_ll = avg_ll - max(avg_ll) ) %>% 
    filter(avg_ll > -30) %>% 
    pull(val)
  
  lower_df <- df %>% 
    filter(val %in% vals_to_consider, 
           val <= mle)
  
  upper_df <- df %>% 
    filter(val %in% vals_to_consider, 
           val >= mle)
  
  lower_mod <- scam(ll ~ s(val, bs = "mpi"), data = lower_df)
  upper_mod <- scam(ll ~ s(val, bs = "mpd"), data = upper_df)
  
  if(sum(vals_to_consider < mle) < 2){
    lower_bound = vals_to_consider[1]
  } else{
    lower_bound <- approx(lower_mod$fitted.values - max(lower_mod$fitted.values), 
                              y = lower_df$val, 
                              xout = likelihood_cutoff)$y
  }
  
  if(sum(vals_to_consider > mle) < 2){
    upper_bound = vals_to_consider[length(vals_to_consider)]
  } else{
    x_vals <- upper_mod$fitted.values - max(upper_mod$fitted.values)
    if (contains_unique(x_vals)) {
    upper_bound <- approx(upper_mod$fitted.values - max(upper_mod$fitted.values), 
                          y = upper_df$val, 
                          xout = likelihood_cutoff)$y  
    } else {
      upper_bound <- 100
    }
    }
    
  
  plot_scam_fit(lower_df, lower_mod, upper_df, upper_mod, param, settings)
  
  if(is.na(upper_bound)){
    upper_bound <- max(upper_df$val)
  }
  if(is.na(lower_bound)){
    lower_bound <- max(lower_df$val)
  }
  
  lower_name <- paste0(param, '_lower')
  upper_name <- paste0(param, '_upper')
  
  bounds <- c(lower_bound,upper_bound)
  names(bounds) <- c(lower_name, upper_name)
  return(bounds)
  

}



plot_scam_fit <- function(df_lower, mod_lower, df_upper, mod_upper, param, settings) {
  if(!dir.exists("figs/scam_fits")){
    dir.create("figs/scam_fits")
  }
  
  pdf(paste0("figs/scam_fits/", settings$model_used, "_", settings$outbreak,"_", param, "_scam_plot.pdf"), width = 8, height=4)
  par(mfrow=c(1,2))
  plot(df_lower$val, df_lower$ll, xlab = "Value", ylab = "Likelihood")
  lines(df_lower$val, predict(mod_lower, newdata = data.frame(val=df_lower$val)), col = "red")
  title(main = paste0(param, "_lower"))
  plot(df_upper$val, df_upper$ll, xlab = "Value", ylab = "Likelihood")
  lines(df_upper$val, predict(mod_upper, newdata = data.frame(val=df_upper$val)), col = "red")
  title(main = paste0(param, "_upper"))
  dev.off()
}


sub_parms <- function(sub_parms=NULL,
                      ref_parms) {
  for(nm in names(sub_parms)) {
    ref_parms[nm] <- sub_parms[nm]
  }
  ref_parms 
}
calc_rnot <- function(fit_parms){
  unname(fit_parms["beta0"] * fit_parms["p0"] / fit_parms["gamma"])
}

calc_cv <- function(fit_parms){
  p <- unname(fit_parms["p0"])
  rnot <- calc_rnot(fit_parms)
  sqrt( (rnot * (2 / p - 1) + 1) / rnot)
}

ss_results <- function(outbreak, max_mif, conf_int) {
  pars <- max_mif@params
  outbreak <- c(outbreak, outbreak)
  parameter <- c('beta', 'p')
  estimate <- c(unname(pars['beta0']), unname(pars['p0'])) 
  lower <- c(unname(conf_int['beta0_lower']), unname(conf_int['p0_lower']))
  upper <- c(unname(conf_int['beta0_upper']), unname(conf_int['p0_upper']))
  results <- data_frame(outbreak, parameter, estimate, lower, upper)
  return(results)
}

int_results <- function(outbreak, max_mif, conf_int) {
  pars <- max_mif@params
  outbreak <- c(outbreak, outbreak, outbreak)
  parameter <- c('beta', 'k', 'tau')
  estimate <- c(unname(pars['beta0']), unname(pars['k']), unname(pars['tau1']))
  lower <- c(unname(conf_int['beta0_lower']), unname(conf_int['k_lower']), unname(conf_int['tau1_lower']))
  upper <- c(unname(conf_int['beta0_upper']), unname(conf_int['k_upper']), unname(conf_int['tau1_upper']))
  results <- data_frame(outbreak, parameter, estimate, lower, upper)
  return(results)
}

contains_unique <- function(x) {
    if (length(x) == 1) return(TRUE)
    x <- range(x) / mean(x)
    !isTRUE(all.equal(x[1], x[2], 0))
  }