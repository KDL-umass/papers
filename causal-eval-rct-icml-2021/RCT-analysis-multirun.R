library(egg)
library(stringr)

setwd("~/Documents/workspace/causal-eval-rct")
# setwd("~/research/causal_evaluation_codebase/causal-eval-rct/")
source("RCT-analysis-pipeline.R")

#assumes treatment is binary and 0,1 (0 = control, 1 = treatment)
#if outcome is binary, assumes outcome is 0,1 (1 = 'outcome of interest')

all_files <- list.files("data")
dataset_names <- unlist(unique(sapply(all_files, function(file){
  if (any(!is.na(str_locate(file, "_config.txt"))))
    return(substr(file, 1, (str_locate(file, "_config.txt")[1,"start"]-1)))
  if (any(!is.na(str_locate(file, "_data.csv"))))
    return(substr(file, 1,(str_locate(file, "_data.csv")[1,"start"]-1)))
})))
#dataset_names <- dataset_names[c(19, 23:26, 31)]
dataset_names <- dataset_names[c(23:25, 27:29, 35:37)]
#dataset_names <- "UK_Data_Service_852874"
#dataset_names <- "Whynot_opioid"
#dataset_names <- "Nemo-breeding"
#dataset_names <- "acic_4"
#dataset_names <- c("Whynot_opioid", "Whynot_world2", "Whynot_zika")
#dataset_names <- c("Nemo-delet-model", "Nemo-breeding")
#dataset_names <- "KNB_f1qf8r51t"
# algorithms_to_run <- c("naive", "dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle")
algorithms_to_run <- c("neural_network")
#algorithms_to_run <- c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "neural_network")
#algorithms_to_run <- "naive"

measures=c("ATE", "outcome")
num_runs <- 30 #will estimate ATE and outcome this many times per algorithm
#algorithms.to.run <- "tree-based"
bias_strength <- 1
max_data_size <- 6000
two_bias_covars <- FALSE
skip_pre_biased <- FALSE
fix_bias_seed <- FALSE

#indexed by dataset name, true if binary outcome, false if continuous outcome
binary_datasets <- list()
results <- list()
for (dataset in dataset_names){
  skip <- FALSE
  print(dataset)

  df <- read.csv(paste0("data/", dataset, "_data.csv"))
  config <- read.csv(paste0("data/", dataset, "_config.txt"), header=FALSE, sep=" ")
  
  apo <- FALSE
  if (substr(dataset, 1, 4) == "Nemo" || dataset %in% c("postgres", "jdk", "networking")){#nemo data sets are APO
    apo <- TRUE
  }
  
  if (substr(dataset, 1, 3) == "ibm" || substr(dataset,1,4) == "acic"){#post-biased data set, need to read in counterfactual outcomes
    ibm_acic <- TRUE
    skip <- skip_pre_biased
    true_ATE = mean(df$counterfactual_outcome_1) - mean(df$counterfactual_outcome_0)

    treatment <- config[1,1]
    outcome <- config[2,1]
    bias_covar <- "" #already biased
    covars <- c()
    if (nrow(config) > 2)#other covariates to read in
      covars <- config[3:nrow(config),1]
    
    df <- df[, c(treatment, outcome, covars)]
  } else{#data set with bias covar specified, will need to be biased later
    ibm_acic <- FALSE

    if (apo){#if apo, then the first line of the config file is actually the index variable
      id_var <- config[1,1]
      config <- config[2:nrow(config),]
    }
    treatment <- config[1,1]
    outcome <- config[2,1]
    bias_covar <- config[3,1]
    covar_start <- 4
    
    if (two_bias_covars){
      if (nrow(config) < 4){
        print(paste0("unable to use two bias covars for ", dataset))
        skip <- TRUE
      } else{
        bias_covar <- c(bias_covar, config[4,1])
        covar_start <- 5
      }
    }
    
    covars <- c()
    if (nrow(config) >= covar_start)#other covariates to read in
      covars <- config[covar_start:nrow(config),1]
    
    if (apo)#if apo, need to keep the id_var around
      df <- df[, c(treatment, outcome, bias_covar, covars, id_var)]
    else
      df <- df[, c(treatment, outcome, bias_covar, covars)]
  }
  
  if (!skip){
    #cast each column to the correct file type
    df[,treatment] <- as.factor(as.character(df[,treatment]))
    for (i in 2:nrow(config)){
      if (!(config[i,1] %in% colnames(df)))
        print(paste0("ERROR: invalid column ", config[i,1], " in config file"))
      if (config[i,2] == 'f'){
        df[,config[i,1]] <- as.factor(as.character(df[,config[i,1]]))
      } else if (config[i,2] == 'n'){
        df[,config[i,1]] <- as.numeric(df[,config[i,1]])
      } else{
        print(paste0("ERROR: invalid data type ", config[i,2], " for column ", config[i,1]))
      }
    }
    
    #a few sanity checks before we get a more subtle error in the OSRCT analysis
    if (length(unique(df[,outcome])) < 3){
      binary_datasets[[dataset]] <- TRUE
      if (!(0 %in% df[,outcome]) || !(1 %in% df[,outcome])){# need to make outcome {0,1}
        print("Warning: binary outcome but levels aren't (0,1)!")
        levels(df[,outcome]) <- c(0,1)
        }
    } else{
      binary_datasets[[dataset]] <- FALSE
    }
    if (length(unique(df[,treatment])) > 3)
      print("Warning: treatment is not binary!")
    if (!(0 %in% df[,treatment]) || !(1 %in% df[,treatment])){# need to make outcome {0,1}
      print("Warning: treatment levels aren't (0,1)!")
      levels(df[,outcome]) <- c(0,1)
    }
    
    if (dataset %in% c("postgres", "jdk", "networking")){#need to subsample one value for each potential outcome for each individual
      df.apo <- c()
      treatment_vals <- unique(df[,treatment])
      for (id in unique(df[,id_var])){
        for (treatment_val in treatment_vals){
          treatment_rows <- df[df[,id_var] == id & df[,treatment] == treatment_val,]
          df.apo <- rbind(df.apo, treatment_rows[sample(1:nrow(treatment_rows), 1),])
        }
      }
      df <- df.apo
    }
    
    #if df is too big to run efficiently, subsample
    if (nrow(df) > max_data_size && !apo)
      df <- df[sample(1:nrow(df), max_data_size),]
    if (nrow(df) > max_data_size*2 && apo){
      sample_ids <- sample(unique(df[,id_var]), max_data_size)
      df <- df[df[,id_var] %in% sample_ids,]
    }
    
    
    
    
    #if we want the same biased samples for all 30 runs, we choose a random seed here
    random_seed <- -1
    if (fix_bias_seed)
      random_seed <- sample(2000000000, 1)
    results[[dataset]] <- list()
    ate_count <- 1
    while (ate_count <= num_runs){
      print(paste0("Dataset ", dataset, ", run ", ate_count))
      
      tryCatch({
      if (ibm_acic){#pre-biased
        curr_results <- OSRCT_analysis(df, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures=measures, already_biased = TRUE, two_bias_covars=two_bias_covars)
      } else if (apo){#apo
        curr_results <- OSRCT_analysis(df, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures=measures, hold_out_test_set = TRUE, id_var = id_var, two_bias_covars=two_bias_covars, random_seed = random_seed)
      } else{#rct
        curr_results <- OSRCT_analysis(df, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures=measures, two_bias_covars=two_bias_covars, random_seed = random_seed)
      }
        
        #save ATE for the current run
        for (alg in algorithms_to_run){
          if (ate_count == 1){
            #only get outcome once, so save outcome for that one ron now
            results[[dataset]][[alg]][["ATE"]] <- c()
            if (ibm_acic){
              results[[dataset]][["true_ATE"]] <- true_ATE
            } else{
              results[[dataset]][["true_ATE"]] <- curr_results$true_ATE
            }
            results[[dataset]][[alg]][["outcome"]] <- c()
          }
          results[[dataset]][[alg]][["ATE"]] <- c(results[[dataset]][[alg]][["ATE"]], curr_results$ATE[[alg]])
          results[[dataset]][[alg]][["outcome"]] <- c(results[[dataset]][[alg]][["outcome"]], curr_results$outcome[[alg]])
        }
      }, 
      error=function(cond){
        message(cond)
        print(cond)
        assign("ate_count", ate_count-1, env=globalenv())
      })
      
      ate_count <- ate_count + 1
    }
    print(results[[dataset]])
  }
}

save(results, file=paste0("all-data-results/final-results/simulator-sample-sizes/nn/run-30-sample-", max_data_size, "-bias-", bias_strength, "-nn-results.RData"))
save(binary_datasets, file=paste0("all-data-results/final-results/simulator-sample-sizes/nn/run-30-sample-", max_data_size, "-bias-", bias_strength, "-nn-binary_datasets.Rdata"))

#save(results, file=paste0("all-data-results/final-results/2000-sample-size-with-acic/50-runs/two-bias/bias-", bias_strength, "/run-50-bias-", bias_strength, "-nn-results.RData"))
#save(binary_datasets, file=paste0("all-data-results/final-results/2000-sample-size-with-acic/50-runs/two-bias/bias-", bias_strength, "/run-50-bias-", bias_strength, "-nn-binary_datasets.Rdata"))


#save(results, file=paste0("all-data-results/neuropathic-pain-1-100000-nn-results.RData"))
#save(binary_datasets, file=paste0("all-data-results/neuropathic-pain-1-100000-nn-binary_datasets.Rdata"))
#true_ATEs <- c()
#for(dataset in dataset_names){
#  df <- read.csv(paste0("data/", dataset, "_data.csv"))
#  config <- read.csv(paste0("data/", dataset, "_config.txt"), header=FALSE, sep=" ")
#  
#  if (substr(dataset, 1, 4) == "Nemo" || dataset %in% c("postgres", "jdk", "networking")){#nemo data sets are APO
#    id_var <- config[1,1]
#    config <- config[2:nrow(config),]
#  }
#  
#  treatment <- config[1,1]
#  outcome <- config[2,1]
#  
#  true_ATEs <- rbind(true_ATEs, c(dataset, get_ATE_from_RCT(df, treatment, outcome), max(df[,outcome])-min(df[,outcome])))
#}
#write.csv(true_ATEs, file="all-ATE.csv")

#outcome_ranges <- c()
#for (dataset in dataset_names){
#  df <- read.csv(paste0("data/", dataset, "_data.csv"))
#  config <- read.csv(paste0("data/", dataset, "_config.txt"), header=FALSE, sep=" ")
#  
#  if (substr(dataset, 1, 4) == "Nemo" || dataset %in% c("postgres", "jdk", "networking")){#nemo data sets are APO
#    id_var <- config[1,1]
#    config <- config[2:nrow(config),]
#  }
#  
#  treatment <- config[1,1]
#  outcome <- config[2,1]
#  
#  outcome_range <- rbind(outcome_ranges, c(dataset, abs(diff(range(df[,outcome])))))
#}