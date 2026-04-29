# STAT6990-Final-Barnabas-Rory-Tripp
This is the code and data used for our Causal Inference final paper and presentation. It consists of a complete R file, two Excel datasets, and supporting documentation.

Evaluating the Effect of New Zealand's Earthquake Commission Act
Authors: Barnabas Amlalo, Rory McDermott, Tripp Modzelewski
Last Updated: April 29, 2026
Course: [STAT 6990 - CAUSAL INFERENCE]

Research Questions
RQ1: To what extent did the introduction of the EQC Act in New Zealand affect the claims severity of natural disaster events?
RQ2: Is the EQC effective in reducing the financial impact of these events on households?

Data Sources
New Zealand catastrophe history (CASdatasets R package, 1968-2014)
Australian catastrophe history (CASdatasets R package)
Population data: Australian Bureau of Statistics (2026), Grimes & Tarrant (2013)

Methods Used                          Method	Purpose
Segmented regression	                Estimate immediate and gradual effects assuming linearity
Bayesian structural time series      	Estimate causal effect without linearity assumption
Difference-in-Differences (unmatched)	Compare NZ vs Australia before/after 1994
Propensity score matching	            Create comparable events between countries
Difference-in-Differences (matched)	  DiD on matched sample (caliper = 0.3)
Durbin-Watson test	                  Test for autocorrelation

Required R Packages
```r install.packages(c( "CASdatasets", "dplyr", "readxl", "ggplot2", "MatchIt", "tidyr", "stringr", "lubridate", "car", "CausalImpact" ))
