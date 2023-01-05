[<img src="../Misc/images/UCLADepartmentofStatisticsSmall.png" align="left" width=500 alt="UCLA STAT Logo"/>](http://statistics.ucla.edu/)<br/>  
   
<p/>      <br/> 
                            
    
Analysis of the Californian Exposure Notification System (CA Notify)
==========

These are code and output related to the analysis of data from the [Californian Exposure Notification System](https://canotify.ca.gov/)

The materials relate to the paper: 
*Exposure Notification System activity as a leading indicator for SARS-COV-2 caseload forecasting* by Sepideh Mazrouee1, Eliah Aronoff-Spencer, Rishi Graham, Mark S. Handcock, Kevin Nguyen, Camille Nebeker, Mohsen Malekinejad, and Christopher A. Longhurst, to appear in PLOS ONE in 2023.


This directory (`Forecast`) contain the primary code used to build a model for forecasting CA SARS-COV-2 cases based on the CA exposure notification (EN) system activity. It uses a log-normal model fitted in `JAGS`. There are two variants of the model. The *baseline* model is similar to the state-of-the-art model in 
[Oliveira and Moral (2021)](https://doi.org/10.1038/s41598-021-87230-x) except that we use a log-normal variant of their negative binomial model. We do this as the negative binomial model is over-dispersed for our situation.
The *EN* model includes the history of exposure notification (EN) to help forecast cases up to seven days into the future.

* This primary code is in the top-level of the directory. The code is mostly in `R`. 
  * `predict_cases_ln.R` is the canonical code. It fits the model using the log-normal model.  The variable `run_with_en=TRUE` toggles fitting the
  model with *EN* verses the *baseline* model. By default it uses 10 lags of history of `EN` and a 4th degree spline for the 10 coefficients of
  those lags. The lagged regressive model is on `sqrt(EN)` rather than the natural `log(EN)` or raw `EN`.
  * `predict_cases_ln4.R` is a variant of `predict_cases_ln.R` that models `EN` with a similar AR(1) latent model to that of cases. It then regresses cases on 7 lags the latent `EN`. This has the advantage of allowing cases to depend on contemporaneous and lagged version of the latent `EN` improving accuracy. Which is does.
  * `predict_seq_p_ln.R` is a complication of `predict_cases_ln.R` the sequentially fits the model. It first fits the model to days 1 through 143 and forecasts days 144 through 150 based on that data. That is, it does 1:7 day ahead forecasts. The code then repeats this for data from days 144, 145, ..., 360. The results are saved in `SavedResults/fore.en_ln_10x4.RData`.
  * `makeRMSE_ln4_seq.R` reads in the results from `predict_seq_p_ln4.R` and produces a simple analysis. These are in `Output/makeRMSE_ln4_seq.R.Rout`. In brief, for this example model the EN improves the RMSE of forecasting by up to 9% for 7-day ahead forecast compared to the baseline model.
  
* Other directories are:
  * `Output`: Saved output from the code, typically `PDF` files.
  * `SavedResults`: Saved output for re-analysis, typically `.RData` files.
  * `data`: The data plus basis files.
  
**Brief description of the results**

* The lagged variables in the *EN*, as a set, are probabilistically greater than zero. Specifically, the posterior probability that at least one of the lagged coefficients is positive is 99.8%. The pattern can be seen in. This can be seen in [Output/coefficient_plot_ln4_en10x4.pdf](Output/coefficient_plot_ln4_en10x4.pdf), which plots 95% ranges for each coefficient. While many of the coefficients are possibly zero, it is unlikely that they all are.
* A secondary question is if the EN variables can significantly improve the forecasting of cases. To assess this we consider forecasting days based on history. Specifically, we fit the model to data from day 1 to day `T` and then forecast days `T+1, T+2, T+7`. Each of these forecasts can then be compared to the recorded cases for those days. We did this for days `T=143, 144, ..., 360`. Two measures are considered: The mean absolute deviations of the forecast from the recorded and the mean-squared error of the forecast from the recorded. The results are given in [Output/compare_ln4_seq.pdf](Output/compare_ln4_seq.pdf). As can be seen, the reductions are modest (5% to 17%), but real.

 **References**

Oliveira, T.P., Moral, R.A. Global short-term forecasting of COVID-19 cases. Sci Rep 11, 7555 (2021).
[https://doi.org/10.1038/s41598-021-87230-x](https://doi.org/10.1038/s41598-021-87230-x)
