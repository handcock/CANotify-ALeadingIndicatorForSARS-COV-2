#=======================================================================
# Loading Packages
#=======================================================================
library(tidyverse)
library(rjags)
library(coda)
library(bayesplot)
library(MCMCvis)
library(runjags, warn.conflicts = F)
library(DescTools)
library(magrittr, warn.conflicts = F)
library(grid)
library(gridExtra, warn.conflicts = F)
library(corrplot)
library(matrixStats)

source("util.R")
#=======================================================================
# Set variables specific to this run
#=======================================================================

skip_run <- TRUE
skip_run <- FALSE
use_raw_data <- FALSE
run_with_en <- FALSE
run_with_en <- TRUE
n_test <- 1
lag_vec <- c(1:20)
lag_poly_order <- 8
lag_vec <- c(1:14)
lag_poly_order <- 5
lag_vec <- c(1:10)
lag_poly_order <- 4
#lag_vec <- c(1:5)
#lag_poly_order <- 3
#=======================================================================
# finding files
# - assume run directory is parallel to mcmc and data directories
#=======================================================================
get.mcmc.file <- function(fname) {
    return(sprintf("mcmc/%s", fname))
}
get.data.file <- function(fname) {
    return(sprintf("data/%s", fname))
}
out.file.sfx <- function(jags_data) {
    if(all(is.null(jags_data$en_lags))) {
        return("no_lag")
    }
    return(sprintf("en%dx%d", ncol(jags_data$en_lags), jags_data$N_lag_basis))
}

#=======================================================================
# Load dataset from file
#   Here we only care about a few of the many columns:
# 1  - date
# 2  - codes issued
# 3  - codes used
# 4  - en webpage visits
# 5  - reported_cases
# 6  - hospitalizations
# 7  - (7day) codes issued
# 8  - (7day) codes used
# 9  - (7day) en webpage visits
# 10 - (7day) reported_cases
# 11 - (7day) hospitalizations  --->>> basic data stops here
# 12 - iEN
# 13 - icases
# Here we'll just read the data and the 3 counts we're interested in.
# We set all the ones we don't care about to "NULL" to leave out
#
# Prep basic data frame
# - convert date to a numeric day offset
#
#=======================================================================

pull_cases <- function(en_col, cases_col, fname="canotify.tsv", extended=FALSE) {
    all_cols = c(
        "date",
        "codes issued",
        "codes used",
        "en webpage visits",
        "reported_cases",
        "hospitalizations",
        "(7day) codes issued",
        "(7day) codes used",
        "(7day) en webpage visits",
        "(7day) reported_cases",
        "(7day) hospitalizations"
        )
    if(extended) {
        all_cols = c(all_cols,
                     "iEN",
                     "icases")
    }
    setClass("num.with.commas")
    setAs("character", "num.with.commas", 
          function(from) as.numeric(gsub(",", "", from) ) )
    colClasses = replace(
        replace(rep("NULL", length(all_cols)), 1, "Date"),
        which(all_cols %in% c(en_col, cases_col)), "num.with.commas"
    )

    raw_data <- read.csv(
        file=get.data.file(fname),
        sep='\t',
        colClasses=colClasses) %>%
        rename(time=date,
               en:=!!str_replace_all(en_col, " ", "."),
               cases:=!!str_replace_all(cases_col, " ", "."))
    idx_col = !is.na(raw_data$time)
    if (en_col == 'en webpage visits') {
        idx_col = idx_col & (raw_data$time > as.Date("2020-12-11"))
    }
    raw_data = raw_data[idx_col,]

    return(raw_data)
}

#=======================================================================
# Polynomial Basis Coefficients
#=======================================================================

get.basis <- function(n_lags, poly_order) {
    df_poly = read.csv(
        file=get.data.file(
            sprintf("basis_%dx%d.csv", n_lags, poly_order)
        )
    )[,-1]
    return(as.matrix(df_poly))
}

#=======================================================================
# Effect vectors
#
# - dow (day of week)
#   - N x 6 matrix with 1 in row i, col j if $time is on weekday j
#     - NOTE: skip one day of the week as that is the base element
#   - Add an extra row for each day of forecast
# - en_lag
#   - N x Np matrix for EN's, with p being lag
#   - these can't be NA, so we just use the last value for now?
#   - technically we should be estimating those guys or just not
#     estimating beyond the first lag
#=======================================================================
get.dow <- function(df_in, n_ahead){
    weekday <- factor(
        weekdays(c(df_in$time, tail(df_in$time,1) + c(1:n_ahead))),
        levels=c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday",
                 "Friday", "Saturday")
    )
    return(as.matrix(model.matrix(~weekday))[,-1] %>%
           set_colnames(sub("weekday","", colnames(.))))
}


## Return matrix whose columns are lagged versions of the given column
## from df_in (repeating first value).  lags arg can be either:
## - vector of integer lags to use
## - max lag (in which case, 0:max will be used)
get.lags <- function(df_in, var, lags, n_ahead) {
    if(is.atomic(lags) && length(lags) == 1L){
        lags = c(0:lags)
    }
    names = str_glue("{var}.Lag{lags}")
    vec = c(df_in[,var], rep(tail(df_in[,var],1), n_ahead))
    Nv = length(vec)
    return(sapply(lags, function(x)head(c(rep(vec[1], x), vec), Nv)) %>%
           set_colnames(names))
}

#=======================================================================
# Gather parameters and vectors for jags run
#=======================================================================

get.jags.data <- function(df_in, n_ahead, lags, lag_poly_order, beta_poly_order=2) {
    df_in$Time <- as.numeric(df_in$time) - min(as.numeric(df_in$time)) + 1
    dow <- get.dow(df_in, n_ahead)

    n_days = nrow(df_in)

    beta_poly <- poly(1:(n_days+n_ahead), beta_poly_order)
    mu_beta <- rep(0, beta_poly_order+1)
    tau_beta <- diag(rep(0.001, beta_poly_order+1))

    dat = list(
            "Y" = c(df_in$cases, rep(NA, n_ahead)),
            "EN" = c(log(df_in$en), rep(NA, n_ahead)),
            "N_data" = n_days,
            "mu_beta" = mu_beta,
            "tau_beta" = tau_beta,
            "beta_basis" = cbind(rep(1, nrow(beta_poly)), beta_poly),
            "dow" = dow,
            "N_ahead" = n_ahead
        )

    if(all(is.na(lags))){
        return(dat)
    }
    dat$en_lags <- get.lags(df_in, "en", lags, n_ahead)
    f <- (get.lags(df_in, "en", lags, n_ahead))
   #f <- cbind(f[,3], f[,-1] - f[,-ncol(f)])
    dat$en_lags <- log(f)
    dat$lag_basis <- get.basis(length(lags), lag_poly_order)
    dat$N_lag_basis <- ncol(dat$lag_basis)
    return(dat)
}

#=======================================================================
# initialization
#=======================================================================
initfunction <- function(chain) {
  return(switch(chain,
                "1" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=1),
                "2" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=2),
                "3" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=3),
                ))
}


#=======================================================================
# Conditional setup
#=======================================================================

nChains <- 3
nAdaptSteps <- 2000
if (run_with_en) {
    nBurninSteps <- 200000
    nThinSteps <- 50
    nUseSteps <- 50000

    nBurninSteps <- 200000
    nThinSteps <- 50
    nUseSteps <- 10000

    nBurninSteps <- 50000
    nThinSteps <- 25
    nUseSteps <- 15000
} else {
    nBurninSteps <- 50000
    nThinSteps <- 100
    nUseSteps <- 10000
}

if (use_raw_data) {
    work_data <- pull_cases('en webpage visits',
                            'reported_cases',
                            fname='canotify_ts_2021-12-20.tsv')
} else {
    work_data <- pull_cases('iEN', 'icases', extended=TRUE, fname='canotify_augmented_2022-1-11.tsv') # changes the variable here
    work_data[,"cases"] <- log(work_data[,"cases"])
  # work_data <- work_data[51:nrow(work_data),]
}

n_train = nrow(work_data) - n_test
train_data = work_data[1:n_train,]
test_data = work_data[n_train+1:length(work_data),]

jags_params <- c("beta",
                 "nu",
                 "dw",
                 "dw_en",
                 "sd_ar",
                 "sd_ar_en",
                 "pi",
                 "sd_omega",
                 "sd_ln",
                 "sd_ln_en",
                 "phi",
                 "ar",
                 "ar_en",
                 "Y_pred",
                 "lambda")
if (run_with_en) {
    jags_data <- get.jags.data(train_data, n_ahead=n_test, lags=lag_vec,
                               lag_poly_order=lag_poly_order)
    jags_params <- append(jags_params, "sin")
    model_file <- get.mcmc.file("model_ln_with_en4.jag")
} else {
    jags_data <- get.jags.data(train_data, n_ahead=n_test, lags=NA)
    model_file <- get.mcmc.file("model_ln_no_en.jag")
}

if(skip_run){
    load(file=(ifelse(run_with_en,sprintf("runJagsOut.en4_ln_%dx%d.RData", length(lag_vec), lag_poly_order),
                                  sprintf("runJagsOut.noen_ln_%dx%d.RData", length(lag_vec), lag_poly_order))))
} else {
    runJagsOut <- run.jags(method = "parallel",
                       model = model_file,
                       monitor = jags_params,
                       data = jags_data,
                       n.chains = nChains,
                       adapt = nAdaptSteps,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunction)
    save(runJagsOut, file=(ifelse(run_with_en,sprintf("SavedResults/runJagsOut.en2_ln_%dx%d.RData", length(lag_vec), lag_poly_order),
                                              sprintf("SavedResults/runJagsOut.noen_ln_%dx%d.RData", length(lag_vec), lag_poly_order))))
}

#=======================================================================
# coda samples - MCMC
#=======================================================================
coda_samples <- as.mcmc.list(runJagsOut)
str(coda_samples[[1]])
colnames(coda_samples[[1]])

#=======================================================================
# Y_pred posterior as matrix
# - rows correspond to mcmc iterations
# - columns correspond to dates
#=======================================================================
get.Ypred.posterior <- function(coda_samples, jags_data, chain=1) {
    pred_mat = NULL
    for(i in 1:(jags_data$N_data + jags_data$N_ahead)) {
        pred_mat = cbind(pred_mat, coda_samples[[chain]][,sprintf("Y_pred[%d]", i)])
    }
    return(pred_mat)
}

#=======================================================================
# mcmc draws for cin based on actual sin walk
#=======================================================================
get.cin.posterior <- function(coda_samples, jags_data, chain=1:3) {
    asin = NULL
    sin = NULL
    for(ch in chain){
      for(i in 1:jags_data$N_lag_basis) {
        sin = cbind(sin, coda_samples[[ch]][,sprintf("sin[%d]", i)])
      }
      asin = rbind(asin, sin %*% t(jags_data$lag_basis))
    }
    return(asin)
}
get.nu.posterior <- function(coda_samples, jags_data, chain=1:3) {
    asin = NULL
    for(ch in chain){
      sin = NULL
      for(i in 1:7) {
        sin = cbind(sin, coda_samples[[ch]][,sprintf("nu[%d]", i)])
      }
      asin = rbind(asin, sin)
    }
    return(asin)
}

#=======================================================================
# forecast on a given day using samples up to the day before
#=======================================================================
get.one.day.fcast.vec <- function(coda_samples, jags_data, day, chain=1:3, run_with_en=TRUE) {
    if(is.na(day)){
        day <- jags_data$N_data + 1
    }
    T <- day - 1
    acases_fore <- NULL
    for(ch in chain){
    s <- coda_samples[[ch]]
    nums <- nrow(s)

    omega <- apply(matrix(s[,"sd_omega"], ncol=1), 1,
                   function(s){rnorm(n=1, mean=0, sd=s)})
    lambda <- apply(matrix(s[,"pi"], ncol=1), 1,
                   function(s){rbernoulli(n=1, p=s)})
    Omega <- lambda*omega

    dow_sum <- 0
     for(i in 1:6) {
        dow_sum <- dow_sum + s[, sprintf("dw[%d]", i)] * jags_data$dow[day, i]
     }

    en_sum <- rep(0, nums)
    if(run_with_en){
     if(!all(is.null(jags_data$en_lags))) {
        cin <- get.cin.posterior(coda_samples, jags_data, chain=ch)
        for(i in 1:ncol(jags_data$en_lags)) {
            en_sum <- en_sum + cin[, i] * jags_data$en_lags[day, i]
        }
     }
    }

    phi_day <- rep(0, nrow=nums)
    for(i in 1:ncol(jags_data$beta_basis)){
      phi_day <- phi_day + s[,sprintf("beta[%d]", i)]*jags_data$beta_basis[day,i]
    }
    phi_day <- phi_day + en_sum/1000

    gamma_day <- phi_day*s[,sprintf("ar[%d]", day-1)] + rnorm(n=nums, mean=0, sd=s[, "sd_ar"])

    mu <- exp(gamma_day + Omega + dow_sum)

    nsim <- 10
    cases_fore <- matrix(NA, ncol=nsim, nrow=length(mu))
    for(i in seq(along=mu)){
     cases_fore[i,] = rnorm(n=nsim, sd=s[i,"sd_ln"], mean=mu[i])
    }

    acases_fore <- rbind(acases_fore, cases_fore)
    }
    return(acases_fore)
}

#=======================================================================
# forecast on a given day using samples up to the day before
#=======================================================================
get.multi.day.fcast.vec <- function(coda_samples, jags_data, day, chain=1:3, run_with_en=TRUE, ahead=7) {
    if(is.na(day)){
        day <- jags_data$N_data + 1
    }
    acases_fore <- NULL
    T <- day - ahead
    for(ch in chain){
    s <- coda_samples[[ch]]
    sel <- 1:nrow(s)
  # sel <- sample.int(nrow(s), size=10)
    s <- s[sel,]
    nums <- nrow(s)
    cases_fore <- matrix(NA, ncol=ahead, nrow=nums)
    gamma_en <- matrix(NA, ncol=ahead, nrow=nums)

    for(ah in (1:ahead)){ 

    omega <- apply(matrix(s[,"sd_omega"], ncol=1), 1,
                   function(s){rnorm(n=1, mean=0, sd=s)})
    lambda <- apply(matrix(s[,"pi"], ncol=1), 1,
                   function(s){rbernoulli(n=1, p=s)})
    Omega <- lambda*omega

    dow_sum <- 0
    for(i in 1:6) {
       dow_sum <- dow_sum + s[, sprintf("dw[%d]", i)] * jags_data$dow[T+ah, i]
    }

    en_sum <- rep(0, nums)
    if(run_with_en){
     if(!all(is.null(jags_data$en_lags))) {
        cin <- get.cin.posterior(coda_samples, jags_data, chain=ch)[sel,]
        for(i in 1:ncol(jags_data$en_lags)) {
            en_sum <- en_sum + cin[, i] * jags_data$en_lags[T+ah, i]
        }
     }
    }

    phi_day <- rep(0, nrow=nums)
    for(i in 1:ncol(jags_data$beta_basis)){
      phi_day <- phi_day + s[,sprintf("beta[%d]", i)]*jags_data$beta_basis[T+ah,i]
    }
   #phi_day <- phi_day + en_sum/1000

    if(run_with_en){
     if(ah == 1){
      gamma_day <- s[,sprintf("ar[%d]", T)]
      for(i in 2:7){
        gamma_en[,i] <- s[,sprintf("ar_en[%d]", T-i+6)]
      }
      gamma_en[,1]  <- s[,sprintf("ar_en[%d]", T+5)]
     }else{
      gamma_day <- phi_day*gamma_day + rnorm(n=nums, mean=0, sd=s[, "sd_ar"])
      for(i in ah:2){
       gamma_en[,i] <-  gamma_en[,i-1]
      }
      gamma_en[,1] <-  phi_day*gamma_en[,2] + rnorm(n=nums, mean=0, sd=s[, "sd_ar_en"])
     }
    }else{
     if(ah == 1){
      gamma_day <- s[,sprintf("ar[%d]", T)]
     }else{
      gamma_day <- phi_day*gamma_day + rnorm(n=nums, mean=0, sd=s[, "sd_ar"])
     }
    }

    mu <- gamma_day + Omega + dow_sum + en_sum/1000
    if(run_with_en){
     for(i in 1:7){
      mu <- mu + s[,sprintf("nu[%d]", i)]*gamma_en[,i]
     }
    }

   #for(i in seq(along=mu)){
   # cases_fore[i,ah] = rnorm(n=1, sd=s[i,"sd_ln"], mean=mu[i])
   ##cases_fore[i,ah] = mu[i]
   #}
    cases_fore[,ah] = mu
    }

    acases_fore <- rbind(acases_fore, cases_fore)
    }
    return(acases_fore)
}

get.one.day.fcast <- function(coda_samples, jags_data, day, chain=1, run_with_en=TRUE) {
    cases_fore = get.one.day.fcast.vec(coda_samples, jags_data, day, chain, run_with_en=run_with_en)
    return(data.frame("cases" = median(cases_fore),
                      "Lower" = quantile(cases_fore,probs=0.025, na.rm=T),
                      "Upper" = quantile(cases_fore,probs=0.975, na.rm=T)))
}

# MSE for the last numPred days
get.pred.rmse <- function(coda_samples, jags_data, numPred, chain=1:3, run_with_en=TRUE, ahead=7){
    df_rmse = data.frame(
        PredictDay=c(-(numPred+1):0 + length(work_data$cases))
    )
    calc.rmse <- function(idx) {
        cases_fore = get.multi.day.fcast.vec(coda_samples, jags_data, idx, chain, run_with_en=run_with_en)
       #t.cases_fore <- cases_fore[cases_fore < quantile(cases_fore,0.95) & cases_fore > quantile(cases_fore,0.05)]
       #t.cases_fore <- apply(cases_fore,2,median)
        t.cases_fore <- mode.density(cases_fore[,7-4])
       #return(sqrt(mean(t.cases_fore-work_data$cases[idx-(6:0)])^2))
        return(sqrt(mean(t.cases_fore-work_data$cases[idx-4])^2))
    }
    df_rmse$MSE <- sapply(df_rmse$PredictDay, calc.rmse)
    return(df_rmse)
}

print("One Day Forecast Reconstructed:")
get.one.day.fcast(coda_samples, jags_data, day=NA, chain=1)


#=======================================================================
# parameter estimates
#=======================================================================
export.estimates <- function(work_data, coda_samples, jags_data, chain=1) {
    # Table of scalar and fixed-length vector parameters
    fixed_pnames <- c(str_glue("beta[{c(1:3)}]"),
                      str_glue("dw[{c(1:6)}]"),
                      str_glue("dw_en[{c(1:6)}]"),
                      str_glue("nu[{c(1:7)}]"),
                      "ar_en[150]","ar_en[151]",
                      "sd_ar",
                      "pi",
                      "sd_omega",
                      "sd_ln")
    df_fixed <- MCMCsummary(coda_samples, params=fixed_pnames, ISB=FALSE, round = 4) %>%
        rename(Median="50%")

    if(!all(is.null(jags_data$en_lags))) {
        sin_pnames <- c(str_glue("sin[{c(1:jags_data$N_lag_basis)}]"))
        df_sin <- MCMCsummary(coda_samples, params=sin_pnames, ISB=FALSE, round = 4) %>%
            rename(Median="50%")
    }
            
    # Table of forecasts for last data day plus test data
    pred_days <- jags_data$N_data + c(0:jags_data$N_ahead)
    pred_names <- str_glue("Y_pred[{pred_days}]")
    df_fcast <- MCMCsummary(coda_samples, params=pred_names, ISB=FALSE, round = 4) %>%
        rename(Median="50%")

    # mean squared prediction error
    df_rmse <- get.pred.rmse(coda_samples, jags_data, 100, chain=chain)
    save(df_rmse, file=ifelse(run_with_en,"SavedResults/df_rmse.en4_ln.RData","SavedResults/df_rmse.noen_ln.RData"))

    # Export to pdf
    pdf(sprintf("Output/mcmc_estimates_ln4_%s.pdf", out.file.sfx(jags_data)))
    grid.table(df_fixed)
    # Couldn't make variable number of grobs work here
    if(!all(is.null(jags_data$en_lags))) {
        grid.arrange(tableGrob(df_sin),
                     tableGrob(df_fcast))
    } else {
        grid.newpage()
        grid.table(df_fcast)
    }

    # NOTE: only around 20 at a time fit on a page.
   #grid.newpage()
   #grid.table(df_rmse[1:15,], row=NULL)
   #grid.newpage()
   #grid.table(df_rmse[16:30,], row=NULL)

    dev.off()
}

export.estimates(work_data, coda_samples, jags_data, chain=1)


#=======================================================================
# Plots
#=======================================================================

plot.trace.param.vec <- function(coda_samples, params, labels) {
    npar <- length(params)
    max_par <- 6
    if(npar > max_par){
        if(mod(npar, max_par) == 1) {
            max_par <- max_par - 1
        }
        par_chunks <- split(params, ceiling(seq_along(params)/max_par))
        lab_chunks <- split(labels, ceiling(seq_along(labels)/max_par))
        for(i in 1:length(par_chunks)){
            plot.trace.param.vec(coda_samples, par_chunks[[i]], lab_chunks[[i]])
        }
        return(NA)
    }
    posterior <- as.matrix(coda_samples[,params])
    plots <- list(NA)
    cur_plot=1
    for(i in 1:npar){
        pars <- c(params[i], params[i])
        plots[[cur_plot]] <- mcmc_trace(posterior[, c(i, i)], pars=pars)
            labs(y=labels[i])
        cur_plot<-cur_plot+1
        plots[[cur_plot]] <- mcmc_areas(posterior[, c(i,i)],
                                        pars=pars,
                                        prob = 0.8) +
            labs(y=labels[i]) + yaxis_text(on=FALSE)
        cur_plot<-cur_plot+1
    }
    if(npar < 6) {
        for(i in (npar+1):6) {
            for(j in 1:2){
                plots[[cur_plot]] <- grid.rect(gp=gpar(col="white"))
                cur_plot<-cur_plot+1
            }
        }
    }
    do.call("grid.arrange", c(plots, ncol=2))
}

plot.intervals.param.vec <- function(coda_samples, params) {
    posterior <- as.matrix(coda_samples[,params]) 
    mcmc_intervals(posterior)
}

plot.cin.errbar <- function(coda_samples, jags_data, chain=1, minmax=FALSE) {
    c_post <- get.cin.posterior(coda_samples, jags_data, chain)
    if(minmax) {
        probs <- c(0, .5, 1)
        title <- "Median with Min and Max values"
        df_q <- as.data.frame(colQuantiles(c_post, probs=probs)) %>%
            setNames(c("Lower", "Median", "Upper"))
    } else {
        sds <- colSds(c_post)
        meds <- colQuantiles(c_post, probs=c(0.5))
        df_q <- as.data.frame(cbind(meds-1.96*sds, meds, meds+1.96*sds)) %>%
            setNames(c("Lower", "Median", "Upper"))
        probs <- c(.05, .5, .95)
        title <- "Median +/- 1.96*std deviation"
    }
    ggplot(df_q, aes(x=row.names(df_q), y=Median)) +
        geom_point(size=2) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
        scale_x_discrete(limits=c(1:dim(df_q)[1])) + 
        labs(x="Lag (Days)", y="Coefficient",
             title=title)
}

plot.nu.errbar <- function(coda_samples, jags_data, chain=1, minmax=FALSE) {
    c_post <- get.nu.posterior(coda_samples, jags_data, chain)
    if(minmax) {
        probs <- c(0, .5, 1)
        title <- "Median with Min and Max values"
        df_q <- as.data.frame(colQuantiles(c_post, probs=probs)) %>%
            setNames(c("Lower", "Median", "Upper"))
    } else {
        sds <- colSds(c_post)
        meds <- colQuantiles(c_post, probs=c(0.5))
        df_q <- as.data.frame(cbind(meds-1.96*sds, meds, meds+1.96*sds)) %>%
            setNames(c("Lower", "Median", "Upper"))
        probs <- c(.05, .5, .95)
        title <- "Median +/- 1.96*std deviation"
        title <- "Coefficients of EN for each lag (+/- 1.96*SD)"
    }
    ggplot(df_q, aes(x=0:6, y=Median)) +
        geom_point(size=2) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
        scale_x_discrete(limits=c(0:dim(df_q)[1])) + 
        labs(x="Lag (Days)", y="Coefficient",
             title=title)
}

 
#-----------------------------------------------------------------------
## MCMC diagnostics
#-----------------------------------------------------------------------

pdf(sprintf("Output/mcmc_plots_ln4_%s.pdf", out.file.sfx(jags_data)))

# Plot beta posterior
params <- c(sprintf("beta[%d]", c(1:3)), sprintf("nu[%d]", c(1:7)), "sd_ln", "pi", "sd_ar", "sd_omega", "sd_ln_en","sd_ar_en")
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

# Plot dw posterior
params <- sprintf("dw[%d]", c(1:6))
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

# Plot nu posterior
params <- sprintf("nu[%d]", c(1:7))
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

# Plot dw posterior
params <- sprintf("dw_en[%d]", c(1:6))
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

# Plot dw posterior
params <- sprintf("ar[%d]", c(40, 150))
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

# Plot sin posterior
if (run_with_en) {
    params <- sprintf("sin[%d]", c(1:jags_data$N_lag_basis))
    labels <- params
    plot.intervals.param.vec(coda_samples, params)
    plot.trace.param.vec(coda_samples, params, labels)
}

# Plot posterior predictive points at and beyond t
times <- c(0:jags_data$N_ahead)
params <- sprintf("Y_pred[%d]", jags_data$N_data+times)
labels <- sprintf("pred t+%d", times)
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

dev.off()


#-----------------------------------------------------------------------
# Reconstruction plots
#
# For some reason if run in the script, a warning in plot.cin.errbar
# ends up delaying plot until after dev.off call.
#-----------------------------------------------------------------------
pdf(sprintf("Output/coefficient_plot_ln4_%s.pdf", out.file.sfx(jags_data)))
if (run_with_en) {
    plot.nu.errbar(coda_samples, jags_data, chain=1:3)
}
dev.off()
