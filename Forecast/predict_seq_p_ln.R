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

library(compiler)
enableJIT(3)
cmpfile("predict_cases_fn.R", "predict_cases_fn.cmpR")

loadcmp("predict_cases_fn.cmpR")
#=======================================================================
# Set variables specific to this run
#=======================================================================

skip_run <- TRUE
skip_run <- FALSE
use_raw_data <- FALSE
run_with_en <- FALSE
run_with_en <- TRUE
n_test <- 7
lag_vec <- c(1:20)
lag_poly_order <- 8
lag_vec <- c(1:14)
lag_poly_order <- 5
lag_vec <- c(1:10)
lag_poly_order <- 4
#lag_vec <- c(1:5)
#lag_poly_order <- 3

#=======================================================================
# Conditional setup
#=======================================================================

nChains <- 1
nAdaptSteps <- 2000
if (run_with_en) {
    nBurninSteps <- 200000
    nThinSteps <- 50
    nUseSteps <- 50000

    nBurninSteps <- 100000
    nThinSteps <- 100
    nUseSteps <- 10000

  # nBurninSteps <- 20000
  # nThinSteps <- 10
  # nUseSteps <- 1000

  # nBurninSteps <- 10000
  # nThinSteps <- 3
  # nUseSteps <- 1000

    nBurninSteps <- 120000
    nThinSteps <- 10
    nUseSteps <- 5000

} else {
    nBurninSteps <- 50000
    nThinSteps <- 100
    nUseSteps <- 10000
 
    nBurninSteps <- 20000
    nThinSteps <- 10
    nUseSteps <- 10000
}

if (use_raw_data) {
    work_data <- pull_cases('en webpage visits',
                            'reported_cases',
                            fname='canotify_ts_2021-12-20.tsv')
} else {
    work_data <- pull_cases('iEN', 'icases', extended=TRUE) # changes the variable here
    work_data[,"cases"] <- log(work_data[,"cases"])
   #work_data <- work_data[51:nrow(work_data),]
}

jags_params <- c("beta",
                 "dw",
                 "sd_ar",
                 "pi",
                 "sd_omega",
                 "sd_ln",
                 "phi",
                 "ar",
                 "Y_pred",
                 "lambda")

library(future)
library(doRNG)  # parallel random number generation streams
ncores <- 30 # number of physical cores
### parallel commands to set up cluster
 doFuture::registerDoFuture()
#cl <- parallel::makeCluster(3)
#future::plan(cluster, workers = parallel.ncores)
#plan(multisession)  ## on MS Windows
 future::plan(multicore, workers = ncores)     ## on Linux, Solaris, and macOS
 ### start each virtual machine with libraries loaded
 # The machines start up empty otherwise (like a new R session)
 packagenames <- c("rjags")
 message(sprintf("Starting parallel using %d cores.",ncores))
# ok, good to go now

fore.at = c(150, 151)
fore.at = c(150:360)
fore <- matrix(0,ncol=8,nrow=length(fore.at))
#for(i in seq_along(fore.at)){
fore <-
     foreach::foreach (i=seq_along(fore.at), .packages=packagenames, .combine=rbind
     ) %dorng% {

#train_data = work_data[51:(fore.at[i]-n_test),]
train_data = work_data[1:(fore.at[i]-n_test),]

if (run_with_en) {
    jags_data <- get.jags.data(train_data, n_ahead=n_test, lags=lag_vec,
                               lag_poly_order=lag_poly_order)
    jags_params <- append(jags_params, "sin")
    model_file <- get.mcmc.file("model_ln_with_en.jag")
} else {
    jags_data <- get.jags.data(train_data, n_ahead=n_test, lags=NA)
    model_file <- get.mcmc.file("model_ln_no_en.jag")
}

if(skip_run){
    load(file=(ifelse(run_with_en,sprintf("runJagsOut.en_ln_%dx%d.RData", length(lag_vec), lag_poly_order),
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
                       silent.jags=TRUE,
                       inits = initfunction)
  # save(runJagsOut, file=(ifelse(run_with_en,sprintf("runJagsOut.en_ln_%dx%d.RData", length(lag_vec), lag_poly_order),
  #                                           sprintf("runJagsOut.noen_ln_%dx%d.RData", length(lag_vec), lag_poly_order))))
}

#=======================================================================
# coda samples - MCMC
#=======================================================================
coda_samples <- as.mcmc.list(runJagsOut)

print("Multi Day Forecast Reconstructed:")
#a <- get.multi.day.fcast.vec(coda_samples, jags_data, day=fore.at[i]-50, chain=1, run_with_en=run_with_en)
a <- get.multi.day.fcast.vec(coda_samples, jags_data, day=fore.at[i], chain=1, run_with_en=run_with_en)
#fore[i,] <- c(fore.at[i], apply(a$cases,2,mode.density))
#a$cases
# c(fore.at[i], apply(a$cases,2,median))
  c(fore.at[i], apply(a[[1]],2,mode.density))
}

print(fore)
save(fore, file= ifelse(run_with_en,
              sprintf("fore.en_ln_%dx%d.RData", length(lag_vec), lag_poly_order),
              sprintf("fore.noen_ln_%dx%d.RData", length(lag_vec), lag_poly_order)))

#export.estimates(work_data, coda_samples, jags_data, chain=1)
