library(plyr)
library(arules)
library(caret)

#takes in a vector for a single covariate (either in one-hot encoding or scaled to be probabilities) and biasing coefficients
#applies the biasing coefficients to the covariate, applies a logistic function, and samples from a binomial
logistic_biasing <- function(bias_vals, coefficients, include_prob = FALSE){#if include_prob, return a dataframe including probability of the selected value
  swept <- sweep(bias_vals, MARGIN=2, coefficients, FUN="*")
  p <- 1 / (1 + exp(-rowSums(swept)))
  p[p > .999] <- .999
  p[p < .001] <- .001
  desired_t <- rbinom(n=nrow(bias_vals), size=1, p)
  if (include_prob){
    return(data.frame(t=desired_t, p=p))  
  }else
    return(desired_t)
}
#samples <- get_biased_samples(df, treatment, bias_covar, covar_discrete, bias_strength, hold_out_test_set, id_var, two_bias_covars=two_bias_covars)
#get_biased_samples(df, treatment, bias_covar, covar_discrete, bias_strength, hold_out_test_set, id_var)
get_biased_samples <- function(df, treatment, bias_covar, covar_discrete, bias_strength, hold_out_test_set = FALSE, id_var = "", fixed.direction = TRUE, two_bias_covars = FALSE){
  #testing the hypothesis that the increased data size for APO biasing is leading to lower variance, gonna try removing that size difference before biasing
  if (id_var != ""){
    print("biasing, APO-style")
    #each instance is in the data twice by default for APO, so get every other row
    df.for.biasing <- df[seq(1, nrow(df), by=2),]
  } else{
    print("biasing, RCT-style")
    df.for.biasing <- df
  }
  if (!two_bias_covars){#normal version, 1 bias covar
    if (covar_discrete){#single discrete covar
      df.for.biasing[,bias_covar] <- as.factor(df.for.biasing[,bias_covar])
      bias_val_matrix <- model.matrix(~ df.for.biasing[,bias_covar] + 0)
      signs <- ifelse(1:ncol(bias_val_matrix) %% 2 == 0, -1, 1)
      bias_val_matrix <- sweep(bias_val_matrix, MARGIN=2, signs, FUN="*")
    }else{#single numeric covar
      bias_cdf <- ecdf(df.for.biasing[,bias_covar[1]])(df.for.biasing[,bias_covar[1]])
      bias_val_matrix <- scale(bias_cdf, center=TRUE, scale=FALSE)
    }
  } else{#funky new version, two bias covars
    if (is.factor(df.for.biasing[,bias_covar[1]]) && is.factor(df.for.biasing[,bias_covar[2]])){
      #if both are factors, combine their values
      #start by getting a matrix with a column for each level of the first bias covar
      bias_val_matrix <- model.matrix(~ df.for.biasing[,bias_covar[1]] + 0)
      signs <- ifelse(1:ncol(bias_val_matrix) %% 2 == 0, -1, .5)
      bias_val_matrix <- sweep(bias_val_matrix, MARGIN=2, signs, FUN="*")
      
      #get an augment based on the value of the second bias covar (.5 or 0) to add to the -1 or .5 from the first bias covar
      bias_val_matrix_2 <- model.matrix(~ df.for.biasing[,bias_covar[2]] + 0)
      augment <- ifelse(1:ncol(bias_val_matrix_2) %% 2 == 0, .5, 0)
      augment <- rowSums(sweep(bias_val_matrix_2, MARGIN=2, augment, FUN="*"))
      
      #combine the two components
      #this is absurdly inefficient - we should just rowSum bias_val_matrix and add augment, but that doesn't work with the fact that coefficient isn't multiplied until later
      #would require more of a rework than I have time for right now
      for (row in 1:nrow(bias_val_matrix)){
        df.for.biasing[row, bias_covar[2]]
        for (col in 1:ncol(bias_val_matrix)){
          if (bias_val_matrix[row, col] != 0)
            bias_val_matrix[row,col] <- bias_val_matrix[row,col] + augment[row]
        }
      }
            
    } else if (is.numeric(df.for.biasing[,bias_covar[1]]) && is.numeric(df.for.biasing[,bias_covar[2]])){
      #if both are numeric, combine their values
      temp <- scale(df.for.biasing[,bias_covar[1]]) + scale(df.for.biasing[,bias_covar[2]]) + scale(scale(df.for.biasing[,bias_covar[1]])*scale(df.for.biasing[,bias_covar[2]]))
      bias_val_matrix <- scale(ecdf(temp)(temp), scale=FALSE)
#      bias_val_matrix <- scale(scale(df.for.biasing[,bias_covar[1]]) + scale(df.for.biasing[,bias_covar[2]]) + scale(scale(df.for.biasing[,bias_covar[1]])*scale(df.for.biasing[,bias_covar[2]])), scale=FALSE)
    } else {
      #one is factor and one is numeric - need to determine which is which first
      if (is.factor(df.for.biasing[,bias_covar[1]])){
        factor_covar <- bias_covar[1]
        numeric_covar <- bias_covar[2]
      } else {
        factor_covar <- bias_covar[2]
        numeric_covar <- bias_covar[1]
      }
      #setup numeric covar as in single variable case, then modify based on value of factor covar
      bias_cdf <- scale(ecdf(df.for.biasing[,numeric_covar])(df.for.biasing[,numeric_covar]), center=TRUE, scale=FALSE)
      factor_mult <- c(.5,1,-.5, -1)
      factor_mult <- rep(factor_mult, ceiling(length(unique(df.for.biasing[,factor_covar]))/4))
      factor_mult <- factor_mult[1:length(unique(df.for.biasing[,factor_covar]))]
      bias_val_matrix <- bias_cdf*rowSums(sweep(model.matrix(~ df.for.biasing[,factor_covar] + 0), MARGIN=2, factor_mult, FUN="*"))
    }
  }
#  #randomly assign either positive or negative bias to each value of a discrete variable (or overall bias for a continuous variable)
#  if (sample(2,1) == 1 | fixed.direction){#fixed direction lets us run multiple trials and compare directly, without extra unnecessary variability
#    signs <- ifelse(1:ncol(bias_val_matrix) %% 2 == 0, -1, 1)
#  }else{
#    signs <- ifelse(1:ncol(bias_val_matrix) %% 2 == 0, 1, -1)
#    print("NEGATIVE BIASING, NOOOOOO")
#  }  
  
  #multiply the bias strength into the desired biasing direction
  coefficients <- rep(bias_strength, ncol(bias_val_matrix))
  
  #get the desired treatment given the biasing covariate
  desired_t <- logistic_biasing(bias_val_matrix, coefficients, include_prob = TRUE)
#  print(paste0("length of desired_t = ", nrow(desired_t)))
  if (id_var != "")
    desired_t$id_var <- df.for.biasing[,id_var]
  
  #loop through every instance, adding it to the primary sample if t == desired_t (and the complementary sample otherwise)
  biased_df <- c()
  complement_df <- c()
  
  if (id_var != ""){#sub-sample by ID, not by row - only works in the APO framework
    print("performing ID sub-sampling")
    for (id in unique(df[,id_var])){
      selected_t <- desired_t[desired_t$id_var == id, "t"]
      biased_row <- df[df[,id_var] == id & df[,treatment] == selected_t,]
      biased_df <- rbind(biased_df, biased_row)
    }
  } else{
    for (i in 1:nrow(df)){
      if (df[i,treatment] == desired_t[i,"t"]){
        biased_df <- rbind(biased_df, df[i,])
      }else{
       row <- df[i,]
       row["p"] <- desired_t[i,"p"]
       complement_df <- rbind(complement_df, row)
      }
    }
  }
  dim(biased_df)
  dim(complement_df)
  
#  if (id_var == ""){
#    write(mean(biased_df[biased_df$index_level == 1, "runtime"]) - mean(biased_df[biased_df$index_level == 0, "runtime"]), file="author response data sets/biased_subsample_ATE_RCT.txt", append=TRUE)
#  }else{
#    write(mean(biased_df[biased_df$index_level == 1, "runtime"]) - mean(biased_df[biased_df$index_level == 0, "runtime"]), file="author response data sets/biased_subsample_ATE_APO.txt", append=TRUE)
#  }
  
  
  #if hold_out_test_set, we want to replace complement_df with a separate sample from biased_df
  #this mostly makes sense in the APO setting, where complement_df is empty
  if (hold_out_test_set){
    print("creating held-out test set")
    #randomize and create a held-out test set
    biased_df <- biased_df[sample(1:nrow(biased_df)),]
    midpoint <- as.integer(nrow(biased_df)/2)
    complement_df <- biased_df[1:midpoint, ]
    biased_df <- biased_df[(midpoint+1):nrow(biased_df), ]
    complement_df$p <- .5
  }
  
#  if (id_var == ""){
#    write(mean(biased_df[biased_df$index_level == 1, "runtime"]) - mean(biased_df[biased_df$index_level == 0, "runtime"]), file="author response data sets/biased_equal_size_ATE_RCT.txt", append=TRUE)
#  }else{
#    write(mean(biased_df[biased_df$index_level == 1, "runtime"]) - mean(biased_df[biased_df$index_level == 0, "runtime"]), file="author response data sets/biased_equal_size_ATE_APO.txt", append=TRUE)
#  }
  
  return(list(primary = biased_df, complement = complement_df))
}

#for compatibility with bnlearn, create a dataframe with a 'from' and a 'to' column of blacklisted edges
#blacklist edges from treatment to covars, from outcome to covars, and from outcome to treatment
make_blacklist <- function(covars, treatments, outcomes){
  blacklist <- c()
  
  blacklist <- rbind(blacklist, adply(treatments, 1, function(treatment){
    t(sapply(covars, function(covar){c(treatment, covar)}))
  }))
  
  blacklist <- rbind(blacklist, adply(outcomes, 1, function(outcome){
    t(sapply(covars, function(covar){c(outcome, covar)}))
  }))
  
  blacklist <- rbind(blacklist, adply(outcomes, 1, function(outcome){
    t(sapply(treatments, function(treatment){c(outcome, treatment)}))
  }))
  
  blacklist$X1 <- NULL
  colnames(blacklist) <- c("from", "to")
  
  return(blacklist)
}


#outcome is passed in because it's the only non-factor field
remove_zero_probs <- function(fit, outcome){
  zero_placeholder <- .1
  for (var in names(fit)){
    if (var != outcome){
      cpt <- fit[[var]]$prob
      if (any(cpt == 0)){#need to remove 0 probs
        print(paste0(var, " has 0 probabilities in cpt: setting those to ", zero_placeholder))
        new_cpt <- ifelse(cpt == 0, zero_placeholder, cpt)
        for (col in 1:ncol(new_cpt)){
          new_cpt[,col] <- new_cpt[,col]/sum(new_cpt[,col])
        }
        fit[[var]] <- new_cpt
      }
    }
  }
  return(fit)
}


#this function only works with a single treatment and single outcome value - if there are multiple, only pass one of each into this function
#also assumes treatment is binary (0,1)
#checks for outcome dimensionality - if binary, assumes 1 = outcome of interest
get_ATE_from_RCT <- function(df, treatment, outcome, risk_ratio = FALSE){
  if (length(unique(df[,outcome])) < 3){#binary outcome
    #not sure this is right
    #P(O=1|T=1) - P(O=1|T=0) = P(O=1 & T=1)/P(T=1) - P(O=1 & T=0)/P(T=0)
    treatment1 <- nrow(df[df[,outcome] == 1 & df[,treatment] == 1,])/nrow(df[df[,treatment] == 1,])
    treatment0 <- nrow(df[df[,outcome] == 1 & df[,treatment] == 0,])/nrow(df[df[,treatment] == 0,])
  } else{#numeric outcome
    treatment1 <- 1/nrow(df[df[,treatment] == 1,])*sum(df[df[,treatment] == 1,outcome])
    treatment0 <- 1/nrow(df[df[,treatment] == 0,])*sum(df[df[,treatment] == 0,outcome])
  }
  if (risk_ratio)
    return(treatment1/treatment0)
  return(treatment1-treatment0)
  
}

#discretize df, without discretizing columns in 'exclude'
discretize_data <- function(df, exclude=c()){
  num_bins <- 5
  discretization_scheme <- list()
  #for each variable, check if it needs to be discretized, and if so, discretize it and add to df
  for (var in colnames(df)){
    if (!(var %in% exclude)){
      if (is.numeric(df[,var])){
        if (length(unique(df[,var])) > num_bins){
          discretization_scheme[[var]] <- arules::discretize(df[,var], method="frequency", breaks=num_bins, infinity=TRUE)
          df[,var] <- discretization_scheme[[var]]
        } else{
          df[,var] <- as.factor(df[,var])
          print(paste0(var, " has only ", length(unique(df[,var])), " unique values, converting into a factor"))
        }
      }
    }
  }
  return(df)
}

#apply the discretization scheme in df.disc to the data in df
apply_discretization <- function(df, discretization_scheme){
  return(discretizeDF(df, methods=discretization_scheme))
}

#make all columns in df numeric, except for any columns listed in 'exclude'
#this is done through one-hot encoding
make_numeric <- function(df, exclude=c()){
  #don't want to make dummy vars for any factor that already just has two levels
  for (col in colnames(df)[!(colnames(df) %in% exclude)]){
    if (is.factor(df[,col]) && length(unique(df[,col])) < 3){
      exclude <- c(exclude, col)
      df[,col] <- as.character(df[,col])
    }
  }
  
  factor_colnames <- colnames(df)[which(sapply(colnames(df), function(col){is.factor(df[,col])}))]
  factor_colnames <- factor_colnames[!(factor_colnames %in% exclude)]
  if (length(factor_colnames) == 0){
    for (col in exclude)
      df[,col] <- as.numeric(df[,col])
    return(df)
  }
  
  factor_cols <- df[,factor_colnames]
  if (length(factor_colnames) == 1){#R does stupid things with single-column dataframes, so handle that here
    factor_cols <- as.data.frame(factor_cols)
    colnames(factor_cols) <- factor_colnames
  }
  dummies <- dummyVars("~ .", factor_cols)
  dummies <- data.frame(predict(dummies, newdata = factor_cols))
  
  #now paste together the exclude columns, the columns that were already numeric, and the new dummy columns
  df.numeric <- df[, !(colnames(df) %in% factor_colnames)]
  df.numeric[,colnames(dummies)] <- dummies
  
  #everything in 'exclude' is just cast to numeric, rather than converted into dummy variables
  for (col in exclude){
    if (is.logical(df.numeric[,col]))
      df.numeric[,col] <- as.numeric(df.numeric[,col])
    else if (is.factor(df.numeric[,col]))
      df.numeric[,col] <- as.numeric(df.numeric[,col])-1
    else if (is.character(df.numeric[,col]) && any(is.na(as.numeric(df.numeric[,col])))){#if it's a character that isn't just a number stored as a character
      df.numeric[,col] <- as.numeric(as.factor(df.numeric[,col]))-1
    } else
      df.numeric[,col] <- as.numeric(as.character(df.numeric[,col]))
  }
  return(df.numeric)
}

#converts coefficient coef in the regression from an odds ratio to a relative risk
#calculates p0 (baseline risk) by counting in the dataframe used to train the regression
convert_odds_ratio_to_relative_risk <- function(regress, df, treatment, outcome){
  coef <- paste0(treatment, "1")
#  control <- regress$coefficients["(Intercept)"]
#  p0 <- exp(control)/(1+exp(control))
  p0 = nrow(df[df[,treatment] == 0 & df[,outcome] == 1,])/nrow(df[df[,treatment] == 0,])
  return(exp(regress$coefficients[coef])/(1-p0+(p0*exp(regress$coefficients[coef]))))
}

convert_odds_ratio_to_risk_difference <- function(regress, df, treatment, outcome){
  coef <- paste0(treatment, "1")
  p0 = nrow(df[df[,treatment] == 0 & df[,outcome] == 1,])/nrow(df[df[,treatment] == 0,])
  
}

expand.class <- function(adj.mat, return.pair.names=F, max_members=100){  
  dags <- list()
  pair.names <- list()
  
  #find all undirected edges
  poss.edges <- combn(colnames(adj.mat), 2, simplify=F)
  undirected.edges <- llply(poss.edges, function(e){
    if(adj.mat[e[[1]], e[[2]]] == 1 & adj.mat[e[[2]], e[[1]]] == 1){
      return(e)
    }
  })
  undirected.edges <- compact(undirected.edges)
  
  #get all permutations of orientations
  poss.orient <- expand.grid(undirected.edges)
  if(ncol(poss.orient) > 0) {
    colnames(poss.orient) <- 1:ncol(poss.orient)
    
    # if we need to stop early, make sure we selected a random sample of members
    if(nrow(poss.orient) > max_members) {
      poss.orient <- poss.orient[sample(nrow(poss.orient)), ]
    }
    #get pdag of all possible orientations of each undirected edge
    l_ply(1:nrow(poss.orient),function(orient.num){
      if(length(dags) <= max_members) {
        adj.mat.temp <- adj.mat
        temp <- c()
        l_ply(colnames(poss.orient), function(edge){
          cause <- which(undirected.edges[[as.integer(edge)]] == poss.orient[[orient.num,edge]])[[1]]
          effect <- which(undirected.edges[[as.integer(edge)]] != poss.orient[[orient.num,edge]])[[1]]
          adj.mat.temp[[ undirected.edges[[as.integer(edge)]][[effect]], undirected.edges[[as.integer(edge)]][[cause]] ]] <<- 0
        })
        
        if(is.acyclic.amat(adj.mat.temp)) {    
          dags[[length(dags)+1]] <<- adj.mat.temp #add our new dag orientation to the list
        }
      }
    })
  } else {
    dags[[1]] <- adj.mat
  }
  
  return(dags)
}


# rebuild the network structure using a new adjacency matrix.
is.acyclic.amat <- function(value) {
  bn <- empty.graph(rownames(value))
  
  if (missing(value))
    stop("no adjacency matrix specified.")
  
  # check the adjacency matrix.
  value = bnlearn:::check.amat(amat = value, nodes = rownames(value))
  
  # update the arcs of the network.
  bn$arcs = bnlearn:::amat2arcs(value, nodes = rownames(value))
  
  return(bnlearn:::is.acyclic(nodes = names(bn$nodes), arcs = bn$arcs))
}