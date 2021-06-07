# causal-eval-rct repository
Code for evaluating causal modeling algorithms using data from RCTs, and assessing the correctness of the OSRCT algorithm

**RCT-analysis-pipline.R** is the main analysis file. It contains a single function (OSRCT_analysis) that does the core biasing and analysis pipeline
* OSRCT_analysis takes in a data set and a bunch of parameters, performs a biased sub-sampling, applies various causal modeling methods, and evaluates and returns those results
* PC, propensity score matching, regression, and BART are implemented as inference algorithms, and ATE and outcome estimation are implemented as evaluation measures
* if an id variable is specified, OSRCT_analysis can also be used for OSAPO analysis

**RCT-analysis-helper-functions.R** provides various utility functions used by RCT-analysis-pipeline, to help keep things clean.  Some of this code was taken from Dan's software systems eval code.  Includes functions that:
* apply covariate biasing and generate biased sample
* create a blacklist for structure learning
* compute ground truth ATE from RCT data
* get Markov equivalence class from PDAG
* discretize data
* convert between an odds ratio and relative risk

**software-rct-analysis.R** sets up and runs the experiments to test the equivalence of APO and RCT sampling, using the software data.
* postgres, jdk, and networking are all implemented
* This file sub-samples the software data to create experimental RCT or APO data, passes the data to OSRCT_analysis, aggregates the results, and creates graphics to compare APO and RCT sampling

**RCT-analysis-multirun.R** performs the OSRCT analysis on a set of data sets in a provided directory
* This file reads in all data sets and sets up the data frames, passes the data to OSRCT_analysis, aggregates the results, and creates graphics
* the data sets for analysis must be provided in a folder, containing data files (of the form "name_data.csv") and config files (of the form "name_config.txt").  Config files must be of this form, where X = N if that value is numeric, and X = F if that value is a discrete factor:
	* treatment X
	* outcome X
	* bias_covar X
	* covariate X
	* covariate X
	* ...
* Each data set provided is checked to confirm it fits the requirements of OSRCT_analysis (treatment is binary, outcome is either continuous or binary, and sub-samples data if it's too big)

**synthetic-eval-test.R** was an initial proof-of-concept test, to check the pipeline idea with synthetic data and confirm the intuiting that bias strength doesn't affect biaesd sample size

