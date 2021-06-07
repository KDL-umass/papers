#library(pcalg)
library(bnlearn)
#library(MXM)
#library(gRain)
library(MatchIt)
library(dplyr)
library(ggplot2)
options(java.parameters = "-Xmx3g")
library(bartMachine)
library(bartCause)
library(stringr)
library(grf)
#library(drtmle)
#library(SuperLearner)
library(ipw)
library(survey)
#library(drgee)
library(reticulate) # R-python integration
library(gbm3)
library(fastDR)

# set up steps for python r integration
# To run neural network approaches, we create conda environment to install the required packages
# and their compatible versions. We use reticulate package for R-python integration. 

# create conda environment and install tensorflow, keras, numpy, scikit-learn and pandas. 
# ONLY TO BE RUN ONCE FOR CREATION and INSTALLATION OF PACKAGES. After that only run lines 34-38.
conda_lists = conda_list()
if ("causal_eval_nn" %in% conda_lists[["name"]] == FALSE){
  conda_create("causal_eval_nn", packages = "python")
  conda_install("causal_eval_nn", "tensorflow==1.14")
  conda_install("causal_eval_nn", "keras==2.2.4")
  conda_install("causal_eval_nn", "numpy ==1.16.4")
  conda_install("causal_eval_nn", "scikit-learn")
  conda_install("causal_eval_nn", "pandas")
}

  # use created conda environment and make sure python path is sourced from there. 
use_condaenv(condaenv = "causal_eval_nn", required = TRUE)
#  use_python(" /Users/agentzel/opt/anaconda3/envs/causal_eval_nn/bin/python3.7")
  # print current path to python to cross-check if environment is properly set 
print(py_config())
nn_model <- import("dragonnet.src.experiment.dragonnet")
#  nn_algs <- c("dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle")
#}

#calculate ATE for structure learning, propensity score, regression, and BART
#calculate outcome prediction for structure learning, regression, and BART
#ATE calculation for binary outcome is actually relative risk

setwd("~/Documents/workspace/causal-eval-rct")
# setwd("~/research/causal_evaluation_codebase/causal-eval-rct/")
source("RCT-analysis-helper-functions.R")

#assumes treatment is binary and 0,1 (0 = control, 1 = treatment)
#if outcome is binary, assumes outcome is 0,1 (1 = 'outcome of interest')


#read in dataset
#df <- read.csv("~/Documents/workspace/causal-eval-rct/Datasets/ISPS d113/wi_ballot_secrecy.csv")
#treatment <- "treat"
#outcome <- "vs_12_update"
#bias_covar <- "milwaukee"
#covars <- c("female", "black", "latino", "mideastern", "nativeamer", "asian", "yearssincereg", "age_corrected")
#covar_discrete <- TRUE

#df <- df[, c(treatment, outcome, bias_covar, covars)]
#df <- df[, c(treatment, outcome, bias_covar)]

#df$treat <- as.factor(df$treat)
#df$vs_12_update <- as.factor(df$vs_12_update)
#df$black <- as.factor(df$black)
#df$yearssincereg <- as.numeric(df$yearssincereg)
#df$age_corrected <- as.numeric(df$yearssincereg)
#df$female <- as.factor(df$female)
#df$milwaukee <- as.factor(df$milwaukee)
#df$latino <- as.factor(df$latino)
#df$mideastern <- as.factor(df$mideastern)
#df$nativeamer <- as.factor(df$nativeamer)
#df$asian <- as.factor(df$asian)


#df <- read.csv("~/Documents/workspace/causal-eval-rct/Datasets/KNB 1596312/BLWarmingExperiment_welldata.csv")
#treatment <- "Treatment"
#outcome <- "CO2_ppm"
#bias_covar <- "Depth"
#covars <- c("Plot")
#df <- df[, c(treatment, outcome, bias_covar, covars)]
#df$Treatment <- as.factor(df$Treatment)
#df$Depth <- as.factor(df$Depth)
#df$Plot <- as.factor(df$Plot)
#df$CO2_ppm <- as.numeric(df$CO2_ppm)
#df[,treatment] <- as.factor(as.integer(df[,treatment] == "Heated"))
#df <- df[!is.na(df[,outcome]),]


alpha = 0.05
cp_queries_to_average_ATE <- 100
cp_queries_to_average_outcome <- 10



#df: experimental data
#treatment: column in df corresponding to treatment
#outcome: column in df corresponding to outcome
#bias_covar: column in df corresponding to biasing covariate
#covars: column(s) in df corresponding to additional pre-treatment covariates
#bias_strength: value denoting strength of bias (typically, -3 through 3)
#algorithms_to_run: vector of names of algorithsm to run (currently supports "structure_learning", "propensity_score", "regression", and "bart")
#measures: a vector of measures to calculate (currently supports "ATE" and "outcome") - if "", estimates both
#hold_out_test_set: whether or not to set aside part of the biased data as a test set
#id_var: column in df corresponding to the id.  If provided, APO-style biasing is used
#already_biased: if true, the data set already has confounding bias, and we should skip the bias step
#df.rct, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures = measures, hold_out_test_set=FALSE, id_var = ""

#results <- OSRCT_analysis(df, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures=measures)

OSRCT_analysis <- function(df, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures=NULL, hold_out_test_set = FALSE, id_var = "", already_biased = FALSE, two_bias_covars = FALSE, random_seed = -1){
  if (is.null(measures))
    measures <- c("ATE", "outcome")

  if (!already_biased){#do biased subsampling
    covar_discrete <- is.factor(df[,bias_covar])
    print("constructing observational data")
    if (random_seed != -1)#if a random seed was passed in, we want to ensure that the biasing protocol produces the same biased samples every time
      set.seed(random_seed)
    samples <- get_biased_samples(df, treatment, bias_covar, covar_discrete, bias_strength, hold_out_test_set, id_var, two_bias_covars=two_bias_covars)
    if (random_seed != -1)#we don't want the same random seeds for the algorithms, so randomize it here
      sample(2000000000, 1)
    df.biased <- samples[["primary"]]
    df.complement <- samples[["complement"]]
    
    covars <- c(covars, bias_covar)
  } else{#if already biased, just need to hold out a test set
    df <- df[sample(1:nrow(df)),]
    midpoint <- as.integer(nrow(df)/2)
    df.complement <- df[1:midpoint, ]
    df.biased <- df[(midpoint+1):nrow(df), ]
    df.complement$p <- .5 #no weighting necessary
  }
  print(paste0("size of df.biased: ", nrow(df.biased)))
  print(paste0("size of df.complement: ", nrow(df.complement)))
 # samples <- get_biased_samples(df.rct, treatment, bias_covar, covar_discrete, bias_strength, FALSE, "")
#  samples.rct <- samples
#  samples <- get_biased_samples(df.apo, treatment, bias_covar, covar_discrete, bias_strength, TRUE, id_var)
#  samples.apo <- samples
#  df.biased.rct <- samples.rct[["primary"]]
#  df.biased.apo <- samples.apo[["primary"]]
#  df.complement.rct <- samples.rct[["complement"]]
#  df.complement.apo <- samples.apo[["complement"]]
#  df.complement.apo$repo_name <- NULL
#  df.biased.apo$repo_name <- NULL
#  regress.rct <- eval(parse(text=paste0("lm(", outcome, " ~ ", paste(c(treatment, bias_covar, covars), collapse=" + "), ", data=df.biased.rct)")))
#  regress.apo <- eval(parse(text=paste0("lm(", outcome, " ~ ", paste(c(treatment, bias_covar, covars), collapse=" + "), ", data=df.biased.apo)")))
#  outcome.apo <- predict(regress.apo, df.complement.apo)
#  outcome.rct <- predict(regress.rct, df.complement.rct)
  
  
  if (!is.null(id_var) && id_var != ""){#remove id
    df.biased <- df.biased[,!(colnames(df.biased) %in% id_var)]
    df.complement <- df.complement[,!(colnames(df.complement) %in% id_var)]
  }

  
  ATEs <- list()
  outcome_estimates <- list()
  binary_outcome <- length(unique(df.biased[,outcome])) < 3
  if (!binary_outcome && is.factor(df.biased[,outcome]))
    print("ERROR: non-binary factor outcomes not supported.  Convert to binary or numeric")
  #run causal modeling algorithms
  if ("structure_learning" %in% algorithms_to_run){
    print("structure learning")
    
    #need to discretize the data to perform cp queries
    df.disc <- discretize_data(df.biased, exclude=outcome)
    
    if (!all(df.disc == df.biased)){
      df.complement_p <- df.complement$p
      df.complement.disc <- apply_discretization(df.complement[,!(colnames(df.complement) == "p")], df.disc)
      df.complement.disc$p <- df.complement_p
    }else {
      df.complement.disc <- df.complement#if nothing changed, we don't have to discretize df.complement}
    }
  #  pc_result <- pc(suffStat = list(dataset=df), indepTest = ci.mm2, u2pd="rand", skel.method="stable.fast",
  #            alpha = alpha, labels=colnames(df), verbose=FALSE)
    
    blacklist <- make_blacklist(c(covars), treatment, outcome)
    
    if (binary_outcome){
      pc_result <- pc.stable(df.disc, blacklist = blacklist)
    }else{
      pc_result <- pc.stable(df.disc, test="mi-cg", blacklist = blacklist)
    }  
    #choose a Markov equivalence class member (fix this)
    mec <- amat(pc_result)
    class.members <- expand.class(mec)
    class.member <- class.members[[sample(1:length(class.members), 1)]]
    #remove incoming edges to treatment
    intervened <- empty.graph(rownames(class.member))
    class.member[, treatment] <- 0
    amat(intervened) <- class.member
    
    print(intervened$arcs)
    
    #fit parameters
    fit <- bn.fit(intervened, df.disc, method="mle")
 #   fit <- remove_zero_probs(fit, outcome)
    
    if ("ATE" %in% measures){
      #to calculate ATE, intervene to set T = 1, and then to set T = 0.  For each, calculate outcome (P(outcome) for binary, average outcome for continuous)
      if (binary_outcome){
        treatment1_query <- paste0("cpquery(fit, event=(", outcome, " == '1'), evidence = (", treatment, " == '1'))")
        treatment1 <- mean(sapply(1:cp_queries_to_average_ATE, function(i){eval(parse(text=treatment1_query))}))
        treatment0_query <- paste0("cpquery(fit, event=(", outcome, " == '1'), evidence = (", treatment, " == '0'))")
        treatment0 <- mean(sapply(1:cp_queries_to_average_ATE, function(i){eval(parse(text=treatment0_query))}))
        #structure_learning_ATE <- treatment1-treatment0
        structure_learning_ATE <- treatment1/treatment0
      } else{#continuous outcome
        #get the P(Y=1|T=1) and P(Y=1|T=0) using cpdist.  Because response varies with repeat runs, take the mean of many runs for each 
        treatment1_query <- paste0("cpdist(fit, nodes=c(outcome), evidence = (", treatment, " == '1'))")
        treatment1 <- mean(sapply(1:cp_queries_to_average_ATE, function(i){mean(eval(parse(text=treatment1_query))[[1]])}))
        treatment0_query <- paste0("cpdist(fit, nodes=c(outcome), evidence = (", treatment, " == '0'))")
        treatment0 <- mean(sapply(1:cp_queries_to_average_ATE, function(i){mean(eval(parse(text=treatment0_query))[[1]])}))
        structure_learning_ATE <- treatment1-treatment0
      }
      print(paste0("ATE calculated across ", nrow(df.disc), " observations"))
      ATEs[["structure_learning"]] <- structure_learning_ATE
    }
    
    if ("outcome" %in% measures){
      #outcome estimation
      structure_learning_outcome <- c()
      for (i in 1:nrow(df.complement.disc)){
#        if (i %% 100 == 0)
#          print(i)
        #create the string the denotes conditioning on treatment+biasing covariates
        evidence_str <- paste(sapply(c(covars, treatment), function(var){paste0("(", var, " == '", df.complement.disc[i,var], "')")}), collapse = " & ")
    
        if (binary_outcome){
          outcome_query <- paste0("cpquery(fit, event=(", outcome, " == '1'), evidence = ", evidence_str, ")")
          structure_learning_outcome <- c(structure_learning_outcome, mean(sapply(1:cp_queries_to_average_outcome, function(i){eval(parse(text=outcome_query))})))
        } else{#continuous outcome
          outcome_query <- paste0("cpdist(fit, nodes=c(outcome), evidence = ", evidence_str, ")")
      #    evidence_list <- as.list(df.complement.disc[i,-8])
      #    outcome_query <- paste0("cpdist(fit, nodes=c('", outcome, "'), evidence=evidence_list, method='lw')")
          vals_to_average <- sapply(1:cp_queries_to_average_outcome, function(i){mean(eval(parse(text=outcome_query))[[1]])})
          vals_to_average <- vals_to_average[!is.na(vals_to_average)]
  #        if (any(is.na(vals_to_average)))
  #          print(all(is.na(vals_to_average)))
          structure_learning_outcome <- c(structure_learning_outcome, mean(vals_to_average))
        }
      }
      
      outcome_estimates[["structure_learning"]] <- structure_learning_outcome
    }
  }
  
  if ("naive" %in% algorithms_to_run){
    ATEs[["naive"]] <- get_ATE_from_RCT(df.biased, treatment, outcome, risk_ratio=FALSE)
  }
  
  if ("propensity_score" %in% algorithms_to_run){
    print("propensity score")
  #  ps_model <- eval(parse(text=paste0("glm(", treatment, " ~ ", paste(c(covars, bias_covar), collapse=" + "), ", family = binomial(), data = df.biased)")))
  #  df.ps <- df.biased
  #  df.ps$ps <- predict(ps_model, type="response")
    
    ps_matched <- eval(parse(text=paste0("matchit(", treatment, " ~ ", paste(covars, collapse=" + "), ", method = 'nearest', replace=TRUE, data = df.biased)")))
    df.matched <- match.data(ps_matched)
    matched.treatment <- df.matched[df.matched[,treatment] == 1, c(outcome, "weights")]
    matched.control <- df.matched[df.matched[,treatment] == 0, c(outcome, "weights")]
    #currently doing matching with weighting, so need to take that into account when estimating ATE (ATT??)
    
    if ("ATE" %in% measures){
      #Currently, binary and continuous look exactly the same.  That may not be true in the future, though, so keeping them separate
      if (binary_outcome){
        matched.treatment[,outcome] <- as.numeric(as.character(matched.treatment[,outcome]))
        matched.control[,outcome] <- as.numeric(as.character(matched.control[,outcome]))
        #P(O=1 | T = 1) - P(O==1 | T = 0), probabilities are weighted counts
        treatment1 <- sum(apply(matched.treatment, 1, function(row){row[1]*row[2]}))/sum(matched.treatment$weights)
        treatment0 <- sum(apply(matched.control, 1, function(row){row[1]*row[2]}))/sum(matched.control$weights)
        #ps_ATE <- treatment1 - treatment0
        ps_ATE <- treatment1-treatment0
      } else{#continuous outcome
        #ATT
        #t-test approach doesn't work if matching returns a weighted sample
        #t_test_result <- eval(parse(text=paste0("with(df.matched, t.test(", outcome, " ~ ", treatment, "))")))
        #ps_ATE <- diff(t_test_result$estimate)[[1]]
        treatment1 <- sum(apply(matched.treatment, 1, function(row){row[1]*row[2]}))/sum(matched.treatment$weights)
        treatment0 <- sum(apply(matched.control, 1, function(row){row[1]*row[2]}))/sum(matched.control$weights)
        ps_ATE <- treatment1 - treatment0
    #    mean(df.matched[df.matched$Treatment == 1, "CO2_ppm"])
    #    mean(df.matched[df.matched$Treatment == 0, "CO2_ppm"])
      }
      ATEs[["propensity_score"]] <- ps_ATE
    }
  }
  
  #the commented-out version if ipw should be equivalent to the version that's currently running
  #but ipwpoint seems to break for data sets with ~50+ variables, so I'm not using it here.
  #However, ipwpoint seems to handle non-numeric variables better, so not a perfect choice
  #(that's why I'm now casting to numeric)
  if ("ipw" %in% algorithms_to_run && "ATE" %in% measures){
    print("ipw")
    #need to cast to numeric
    df.relevant <- df.biased[,c(treatment, outcome, covars)]
    df.cont <- make_numeric(df.relevant)
    #training_command <- paste0("ipwpoint(exposure = ", treatment, ", family = \"binomial\", link = \"logit\", numerator = ~ 1, denominator = ~ ", paste(colnames(df.cont)[!(colnames(df.cont) %in% c(treatment, outcome))], collapse = " + "), ", data = df.cont)")
    #df.weights <- df.cont
    #df.weights$sw <- eval(parse(text=training_command))$ipw.weights
    
    #if (binary_outcome)
    #  df.weights[,outcome] <- as.numeric(df.weights[,outcome])
    #inference_command = paste0("svyglm(", outcome, " ~ ", treatment, ", design = svydesign(~ 1, weights = ~ sw, data = df.weights))")
    #msm <- eval(parse(text=inference_command))
    
    #if (binary_outcome){
    #  ATEs[["ipw"]] <- (msm$coefficients[1]+msm$coefficients[2])-msm$coefficients[1]
    #} else{
    #  ATEs[["ipw"]] <- msm$coefficients[2]
    #}
    
    #turns out the ipwpoint command breaks if you have ~50 variables, so just do it more manually
    ps_command <- paste0("glm(", treatment, " ~ ", paste(colnames(df.cont)[!(colnames(df.cont) %in% c(treatment, outcome))], collapse = " + "), ", family = binomial(link=\"logit\"), data=df.cont)")
    ps_preds <- predict(eval(parse(text=ps_command)), type="response")
    weights <- ifelse(df.biased[,treatment] == 1, 1/ps_preds, 1/(1-ps_preds))

    if (binary_outcome){
      glm_command <- paste0("glm(", outcome, " ~ ", treatment, ", weights=weights, data=df.cont, family = binomial(link=\"identity\"))")
      glm <- eval(parse(text=glm_command))
      ATEs[["ipw"]] <- glm$coefficients[2]
    } else{
      glm_command <- paste0("glm(", outcome, " ~ ", treatment, ", weights=weights, data=df.biased, family = gaussian)")
      glm <- eval(parse(text=glm_command))
      ATEs[["ipw"]] <- glm$coefficients[2]
    }
    
  #  glm <- eval(parse(text=glm_command))
    
  }
  
  if ("regression" %in% algorithms_to_run){
    print("regression")
    
    if (binary_outcome){
      regress <- eval(parse(text=paste0("glm(", outcome, " ~ ", paste(c(treatment, covars), collapse=" + "), ", data=df.biased, family='binomial')")))
      #regress_ATE <- regress$coefficients[2]
      df.treatment1.interv <- df.biased
      df.treatment1.interv[,treatment] <- 1
      df.treatment1.interv[,treatment] <- factor(df.treatment1.interv[,treatment], levels=levels(df.biased[,treatment]))
      df.treatment0.interv <- df.biased
      df.treatment0.interv[,treatment] <- 0
      df.treatment0.interv[,treatment] <- factor(df.treatment0.interv[,treatment], levels=levels(df.biased[,treatment]))
      #relative_risk <- convert_odds_ratio_to_relative_risk(regress, df, treatment, outcome)
      risk_difference <- mean(predict(regress, df.treatment1.interv, type="response")) -  mean(predict(regress, df.treatment0.interv, type="response"))
      regress_ATE <- risk_difference
      regress_outcome <- predict(regress, df.complement, type="response")
    } else{#continuous outcome
      regress <- eval(parse(text=paste0("lm(", outcome, " ~ ", paste(c(treatment, covars), collapse=" + "), ", data=df.biased)")))
      regress_ATE <- regress$coefficients[2]
      regress_outcome <- predict(regress, df.complement)
    }
    
    if ("ATE" %in% measures)
      ATEs[["regression"]] <- regress_ATE
    if ("outcome" %in% measures)
      outcome_estimates[["regression"]] <- regress_outcome
  }
  
  if ("bart" %in% algorithms_to_run){
    print("bart")
    bart_model <- bartMachine(df.biased[,c(treatment, covars)], df.biased[,outcome], verbose=FALSE)
    df.treatment1.interv <- df.biased
    df.treatment1.interv[,treatment] <- 1
    df.treatment0.interv <- df.biased
    df.treatment0.interv[,treatment] <- 0
    
    if ("ATE" %in% measures){
      if (binary_outcome){
  #      treatment1_preds <- as.integer(as.character(predict(bart_model, df.treatment1.interv[,c(treatment, bias_covar, covars)], type="class")))
  #      treatment0_preds <- as.integer(as.character(predict(bart_model, df.treatment0.interv[,c(treatment, bias_covar, covars)], type="class")))
        #bart_ATE <- sum(treatment1_preds)/length(treatment1_preds) - sum(treatment0_preds)/length(treatment0_preds)
        #the '1-' is because, for some reason, predict treats '0' as the target level, so need to adjust for that
        treatment1_probs <- 1-predict(bart_model, df.treatment1.interv[,c(treatment, covars)], type="prob")
        treatment0_probs <- 1-predict(bart_model, df.treatment0.interv[,c(treatment, covars)], type="prob")
        bart_ATE <- mean(treatment1_probs)-mean(treatment0_probs)
      } else{#continuous outcome
        treatment1_preds <- predict(bart_model, df.treatment1.interv[,c(treatment, covars)])
        treatment0_preds <- predict(bart_model, df.treatment0.interv[,c(treatment, covars)])
        bart_ATE <- mean(treatment1_preds) - mean(treatment0_preds)
      }
      
      ATEs[["bart"]] <- bart_ATE
    }
    
    if ("outcome" %in% measures){
      if (binary_outcome){
        #the predictions from bartMachine appear to predict the probability of being in the first class in bart_model$y_levels, which it gets from levels(df.biased[,outcome])
        bart_outcome <- predict(bart_model, df.complement[,c(treatment, covars)])
        if (!is.null(bart_model$y_levels) && bart_model$y_levels[1] == "0")#if it's actually giving us the probability that Y=0, need to switch it
          bart_outcome <- 1 - bart_outcome
      }else{#continuous outcome
        bart_outcome <- predict(bart_model, df.complement[,c(treatment, covars)])
      }
      outcome_estimates[["bart"]] <- bart_outcome
    }
  }
  
  if ("bartCause" %in% algorithms_to_run && "ATE" %in% measures){
    print("bartCause")
    df.temp <- df.biased
    df.temp[,treatment] <- as.numeric(as.character(df.biased[,treatment]))
    df.temp[,outcome] <- as.numeric(as.character(df.biased[,outcome]))
    bart_model <- bartc(df.temp[,outcome], df.temp[,treatment], df.temp[,covars], method.rsp="bart", method.trt="bart", keepTrees=TRUE, keepCall=TRUE, estimand="ate")
    ATEs[["bartCause"]] <- fitted(bart_model)
  }
  
  if ("causal-forest" %in% algorithms_to_run){
    print("causal forests")
    #need to make data continuous for use with causal forests
    if ("ATE" %in% measures){
      df.relevant <- df.biased[,c(treatment, outcome, covars)]
      df.cont <- make_numeric(df.relevant)
      c.forest <- causal_forest(df.cont[,!(colnames(df.cont) %in% c(treatment, outcome))], df.cont[,outcome], df.cont[,treatment])

      ATEs[["causal-forest"]] <- average_treatment_effect(c.forest)["estimate"]
      
      #tau.hat <- predict(c.forest)$predictions
      #if (binary_outcome){
        #ATEs[["causal-forest"]] <- mean(c.forest$Y.hat + (1-c.forest$W.hat)*tau.hat)-mean(c.forest$Y.hat + c.forest$W.hat*tau.hat)
      #}else{
          #ATEs[["causal-forest"]] <- mean(tau.hat)
      #  ATEs[["causal-forest"]] <- average_treatment_effect(c.forest)["estimate"]
      #}
    }
  }
  
  if ("doubly-robust" %in% algorithms_to_run && "ATE" %in% measures){
    print("doubly-robust")
    df.temp <- df.biased
    df.temp[,treatment] <- as.numeric(as.character(df.temp[,treatment]))
    df.temp[,outcome] <- as.numeric(as.character(df.temp[,outcome]))
    df.temp$id_var <- 1:nrow(df.temp)
    df.temp$weights <- rep(1, nrow(df.temp))
    
    if (binary_outcome){
      dr.method <- "quasibinomial"
    } else
      dr.method <- "gaussian"
    fastDR.command <- paste0("fastDR(list(y.form=~", outcome, ", t.form=~", treatment, ", x.form=~", paste(c(bias_covar,covars), collapse="+"), ", weights.form=~weights, key.form=~id_var), data=df.temp, y.dist=\"", dr.method, "\", n.trees=3000, interaction.dept=3, shrinkage=0.005, verbose=FALSE, smooth.lm=0)")
    fastDR.model <- eval(parse(text=fastDR.command))
    
    ATEs[["doubly-robust"]] <- fastDR.model$effects[[outcome]]["dr","TE"]
#    drgee.command <- paste0("drgee(oformula=formula(", outcome, "~", paste(covars, collapse="+"), "), eformula=formula(", treatment, "~", paste(covars, collapse="+"), "), data=df.biased)")
#    drgee.model <- eval(parse(text=drgee.command))
#    ATEs[["doubly-robust"]] <- drgee.model$coefficients

#    df.tmle <- df.biased
#    df.tmle[,treatment] <- as.numeric(as.character(df.tmle[,treatment]))
#    w <- df.tmle[,covars]
#    if (length(covars) < 2)
#      w <- as.data.frame(w)
    
    #drtmle freaks out if there is a relatively rare class, so we're gonna remove those here....
#    imbalanced_cols <- c()
#    for (col in colnames(w))
#      if (is.factor(w[,col]) && min(plyr::count(w[,col])[,2]) < 10)
#        imbalanced_cols <- c(imbalanced_cols, col)
#    if (length(imbalanced_cols) > 0)
#      print(paste0("NOTE: Removing columns ", paste(imbalanced_cols, collapse=", "), " for drtmle due to low-frequency factors"))
   # w <- w[,!(colnames(w) %in% imbalanced_cols)]
#    if (binary_outcome){
#      df.tmle[,outcome] <- as.numeric(as.character(df.tmle[,outcome]))
      
#      dr.fit <- drtmle(Y = df.tmle[,outcome], A = df.tmle[,treatment], W = w, family = binomial(), SL_Q = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#                     SL_g = c("SL.glm", "SL.mean", "SL.glm.interaction"), SL_Qr = "SL.glm", SL_gr = "SL.glm", maxIter=1)
#    } else {
#      dr.fit <- drtmle(Y = df.tmle[,outcome], A = df.tmle[,treatment], W = w, family = gaussian(), SL_Q = c("SL.glm", "SL.mean", "SL.glm.interaction"),
#                      SL_g = c("SL.glm", "SL.mean", "SL.glm.interaction"), SL_Qr = "SL.glm", SL_gr = "SL.glm", maxIter=1)
#    }
#    if (df.tmle[1,treatment] == 1)#drtmle figures out what to call 'Treatment 0' based on which treatment value it sees first, so have to adjust for that
#      ATEs[["doubly-robust"]] <- dr.fit$drtmle$est[1]-dr.fit$drtmle$est[2]
#    else
#      ATEs[["doubly-robust"]] <- dr.fit$drtmle$est[2]-dr.fit$drtmle$est[1]
  }
#  nn_check = all(nn_algs %in% algorithms_to_run)
  if("neural_network" %in% algorithms_to_run){
    print("neural network based approaches")
    
    df.brelevant <- df.biased[,c(treatment, outcome, covars)]
    df.biased_cont <- make_numeric(df.brelevant)
    
    df.crelevant <- df.complement[,c(treatment, outcome, covars)]
    df.complement_cont <- make_numeric(df.crelevant)
    
    n_train = nrow(df.biased_cont)
    n_test = nrow(df.complement_cont)
    y_tr = array(df.biased_cont[, outcome], dim = c(n_train, 1))
    y_te = array(df.complement_cont[, outcome], dim = c(n_test, 1))

    t_tr = array(df.biased_cont[,treatment], dim = c(n_train, 1))
    t_te = array(df.complement_cont[,treatment], dim = c(n_test, 1))
    x_tr = df.biased_cont[,!(colnames(df.biased_cont) %in% c(treatment, outcome))]
    x_te = df.complement_cont[,!(colnames(df.biased_cont) %in% c(treatment, outcome))]
   
    #tarnet_out = nn_model$run(y_tr, t_tr, x_tr, y_te, t_te, x_te, knob = 'tarnet')
    dgnet_out = nn_model$run(y_tr, t_tr, x_tr, y_te, t_te, x_te, knob = 'dragonnet', binary_flag = binary_outcome)
    #tarnet_out = nn_model$run(y_tr, t_tr, x_tr, y_te, t_te, x_te, knob = 'tarnet', binary_flag = binary_outcome)

    #ATEs[["dragonnet"]] <- dgnet_out$t_reg_false$train_ate
    #ATEs[["dragonnet_tmle"]] <- dgnet_out$t_reg_true$train_ate
    #ATEs[["tarnet"]] <- tarnet_out$t_reg_false$train_ate
    #ATEs[["tarnet_tmle"]] <- tarnet_out$t_reg_true$train_ate
    ATEs[["neural_network"]] <- dgnet_out$t_reg_true$train_ate
    
    #outcome_estimates[["dragonnet"]] <- dgnet_out$t_reg_false$test_outcome
    #outcome_estimates[["dragonnet_tmle"]] <- dgnet_out$t_reg_true$test_outcome
    #outcome_estimates[["tarnet"]] <- tarnet_out$t_reg_false$test_outcome
    #outcome_estimates[["tarnet_tmle"]] <- tarnet_out$t_reg_true$test_outcome
    outcome_estimates[["neural_network"]] <- dgnet_out$t_reg_true$test_outcome
   }
  
  
  #calculate ground truth effects
  true_ATE <- get_ATE_from_RCT(df, treatment, outcome, risk_ratio=FALSE)
  
  #for outcome prediction for continuous outcome, do we calculate error for outcome estimate, and weight those errors by inverse probability of inclusion?
  #and what do we do for binary outcome?  count binary errors (correct vs incorrect) and weight those by inverse probability of inclusion?
  
  weighted_errors <- list()
  if (binary_outcome){
    for (method in names(outcome_estimates)){
      print(paste0("Evaluating ", method))
      estimates <- outcome_estimates[[method]]
      actual_outcomes <- as.integer(as.character(df.complement[,outcome]))
      #calculate weighted probability of incorrect estimate
#      weighted_errors[[method]] <- sum(sapply(1:length(estimates), function(i){abs(actual_outcomes[i]-estimates[i])*df.complement[i,"p"]}))/sum(df.complement[,"p"])
      weighted_errs <- sapply(1:length(estimates), function(i){abs(actual_outcomes[i]-estimates[i])*df.complement[i,"p"]/(1-df.complement[i,"p"])})
      weighted_errors[[method]] <- sum(weighted_errs)/sum(sapply(1:length(estimates), function(i){df.complement[i,"p"]/(1-df.complement[i,"p"])}))
    }
  }else{
    for (method in names(outcome_estimates)){
      print(paste0("Evaluating ", method))
      estimates <- outcome_estimates[[method]]
      actual_outcomes <- df.complement[,outcome]
      #calculate mean weighted absolute error
      weighted_errs_alt <- sapply(1:length(estimates), function(i){abs(estimates[i]-actual_outcomes[i])*df.complement[i,"p"]})
      weighted_errs <- sapply(1:length(estimates), function(i){abs(estimates[i]-actual_outcomes[i])*df.complement[i,"p"]/(1-df.complement[i,"p"])})
      valid.indices <- !(1:length(weighted_errs) %in% which(is.na(weighted_errs)))
      weighted_errors[[method]] <- sum(weighted_errs[valid.indices])/sum(sapply(1:length(estimates), function(i){df.complement[i,"p"]/(1-df.complement[i,"p"])})[valid.indices])
 #     print(paste0(weighted_errors[[method]], "=new, ", weighted_errors[[method]] <- sum(weighted_errs_alt[valid.indices])/sum(df.complement[valid.indices,"p"]), "=original"))
      #valid.indices <- !(1:length(weighted_errs) %in% which(is.na(weighted_errs)))
      #weighted_errors[[method]] <- sum(weighted_errs[valid.indices])/sum(df.complement[valid.indices,"p"])
      if (any(is.na(weighted_errs)))
        print(paste0("WARNING: Some estimates were NA for method = ", method))
      
    }
  }
  
  return(list(ATE=ATEs, true_ATE=true_ATE, outcome=weighted_errors))
}





