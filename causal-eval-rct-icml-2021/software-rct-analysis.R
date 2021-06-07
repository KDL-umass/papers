library(gridExtra)
library(egg)

setwd("~/Documents/workspace/causal-eval-rct")
source("RCT-analysis-pipeline.R")

styles <- c("RCT", "APO") #either APO or RCT
dataset <- "postgres"
algorithms_to_run <- c("structure_learning", "propensity_score", "regression", "tree-based")
#algorithms_to_run <- "regression"
measures <- c("ATE", "outcome")
bias_strength <- 1
num_reps <- 100
log_outcome <- FALSE

if (dataset == "postgres"){
  df <- read.csv("../causal-eval/postgres/results_abbrv_all.csv")
  id_var <- "url"
  treatment = "index_level"
  outcome = "runtime"
  if (log_outcome){
    df[,paste0("log_", outcome)] <- log(df[,outcome])
    outcome <- paste0("log_", outcome)
  }
  bias_covar <- "rows"
  covars <- c("creation_year", "num_ref_tables", "num_joins", "num_group_by", "queries_by_user", "length_chars", "total_ref_rows")
  df[,treatment] <- as.factor(df[,treatment])
  df[,outcome] <- as.numeric(df[,outcome])
  df[,bias_covar] <- as.numeric(df[,bias_covar])
  df$creation_year <- as.numeric(df$creation_year)
  df$num_ref_tables <- as.numeric(df$num_ref_tables)
  df$num_joins <- as.numeric(df$num_joins)
  df$num_group_by <- as.numeric(df$num_group_by)
  df$queries_by_user <- as.numeric(df$queries_by_user)
  df$length_chars <- as.numeric(df$length_chars)
  df$total_ref_rows <- as.numeric(df$total_ref_rows)
  
  #make sure each url has both potential outcomes - if one doesn't, get rid of it
  cursed_ids <- c()
  for (id in unique(df[,id_var])){
    if (length(unique(df[df[,id_var] == id,treatment])) < 2)
      cursed_ids <- c(cursed_ids, id)
  }
  df <- df[!(df[,id_var] %in% cursed_ids), ]
}

if (dataset == "jdk"){
  df <- read.csv("../causal-eval/jdk/results.csv")
  id_var <- "repo_name"
  treatment <- "obfuscate"
  outcome <- "num_bytecode_ops"
  bias_covar <- "test_javadocs"
  covars <- c("source_ncss", "test_classes", "test_functions", "test_ncss")
  df[,treatment] <- as.factor(df[,treatment])
  df[,outcome] <- as.numeric(df[,outcome])
  df[,bias_covar] <- as.numeric(df[,bias_covar])
  df$source_ncss <- as.numeric(df$source_ncss)
  df$test_classes <- as.numeric(df$test_classes)
  df$test_functions <- as.numeric(df$test_functions)
  df$test_ncss <- as.numeric(df$test_ncss)
}

if (dataset == "networking"){
  df <- read.csv("../causal-eval/networking/results2.dat")
  id_var <- "site"
  treatment <- "proxy"
  outcome <- "elapsed"
  bias_covar <- "server.class"
  covars <- c()
  df[,treatment] <- as.factor(df[,treatment])
  df[,outcome] <- as.numeric(df[,outcome])
  df[,bias_covar] <- as.factor(df[,bias_covar])
}


df <- df[,c(id_var, treatment, outcome, bias_covar, covars)]

#df <- df[df$url %in% sample(unique(df$url), 5000),]

results <- c()

random_sample_ATE_RCT <- c()
random_sample_ATE_APO <- c()

#after initial biasing, when APO has ~double the instances of RCT
#biased_subsample_ATE_RCT <- c(read.csv("author response data sets/biased_subsample_ATE_RCT.txt"))[[1]]
#biased_subsample_ATE_APO <- c(read.csv("author response data sets/biased_subsample_ATE_APO.txt"))[[1]]

#after random subsampling APO
#biased_equal_size_ATE_RCT <- c(read.csv("author response data sets/biased_equal_size_ATE_RCT.txt"))[[1]]
#biased_equal_size_ATE_APO <- c(read.csv("author response data sets/biased_equal_size_ATE_APO.txt"))[[1]]

true_ATE <- mean(df[df[,treatment] == 1, outcome]) - mean(df[df[,treatment] == 0, outcome])

for (rep in 1:num_reps){
  print(paste0("STARTING REP ", rep))

  if ("APO" %in% styles){
    #sub-sample to have both potential outcomes once each for every individual (APO style)
    df.apo <- c()
    treatment_vals <- unique(df[,treatment])
    for (id in unique(df[,id_var])){
      for (treatment_val in treatment_vals){
        treatment_rows <- df[df[,id_var] == id & df[,treatment] == treatment_val,]
        df.apo <- rbind(df.apo, treatment_rows[sample(1:nrow(treatment_rows), 1),])
      }
    }
    
    print(paste0("APO ATE = ", (mean(df.apo[df.apo$index_level == 1, "runtime"]) - mean(df.apo[df.apo$index_level == 0, "runtime"]))))
    random_sample_ATE_APO <- c(random_sample_ATE_APO, (mean(df.apo[df.apo$index_level == 1, "runtime"]) - mean(df.apo[df.apo$index_level == 0, "runtime"])))
    
    APO_result <- OSRCT_analysis(df.apo, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures = measures, hold_out_test_set=TRUE, id_var = id_var)
    names(APO_result$ATE$regression) <- NULL
  }
  if ("RCT" %in% styles){
    #sub-sample to only have a single treatment for every individual (RCT style)

    df.rct <- c()
    treatment_vals <- unique(df.apo[,treatment])
    if ("APO" %in% styles){ #if we're also doing APO sampling, sub-sample that data to create RCT data
      print("sampling RCT data from APO data")
      for (id in unique(df.apo[,id_var])){
        treatment_val <- sample(treatment_vals,1)
        df.rct <- rbind(df.rct, df.apo[df.apo[,id_var] == id & df.apo[,treatment] == treatment_val,])
      }
    } else{#if we're not also doing APO sampling, sub-sample the full data as normal
      print("sampling RCT data from original data")
      for (id in unique(df[,id_var])){
        treatment_val <- sample(treatment_vals,1)
        id_rows <- df[df[,id_var] == id & df[,treatment] == treatment_val,]
        df.rct <- rbind(df.rct, id_rows[sample(1:nrow(id_rows), 1),])
      }
    }
    df.rct[,id_var] <- NULL
    
    print(paste0("RCT ATE = ", (mean(df.rct[df.rct$index_level == 1, "runtime"]) - mean(df.rct[df.rct$index_level == 0, "runtime"]))))
    random_sample_ATE_RCT <- c(random_sample_ATE_RCT, (mean(df.rct[df.rct$index_level == 1, "runtime"]) - mean(df.rct[df.rct$index_level == 0, "runtime"])))
    
    RCT_result <- OSRCT_analysis(df.rct, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, measures = measures, hold_out_test_set=FALSE, id_var = "")
    names(RCT_result$ATE$regression) <- NULL
  }
  
  #get results out of the mess returned by OSRCT_analysis and append to the master 'results' dataframe
  for (style in styles){
    curr_result <- get(paste0(style, "_result"))
    #   true_ATE <- curr_result[["true_ATE"]]
    print(paste0("true_ATE = ", true_ATE))
    for (measure in names(curr_result)[!(names(curr_result) == "true_ATE")]){
      algs <- names(curr_result[[measure]])
      for (alg in algs){
        if (measure == "ATE"){
          results <- rbind(results, c(style, measure, alg, curr_result[[measure]][[alg]]-true_ATE))
          print(c(style, measure, alg, curr_result[[measure]][[alg]]-true_ATE))
        }else{
          results <- rbind(results, c(style, measure, alg, curr_result[[measure]][[alg]]))
          print(c(style, measure, alg, curr_result[[measure]][[alg]]))
        }
      }
    }
  }
  
  
  #print(curr_result)
  

    #OSRCT_analysis(df.rct, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, FALSE, "")
    #OSRCT_analysis(df.apo, treatment, outcome, bias_covar, covars, bias_strength, algorithms_to_run, TRUE, id_var)
}
save(results, file="author response data sets/postgres_apo_vs_rct_200_runs_reducing_variance_final.RData")

#results.postgres
#load("author response data sets/postgres_apo_vs_rct_100_runs_reducing_variance_final.RData")
#load("results/postgres_apo_vs_rct_19_runs_one_direction.RData")

results <- data.frame(results)
colnames(results) <- c("style", "measure", "method", "error")
results$style <- as.factor(results$style)
results$measure <- as.factor(results$measure)
results$method <- as.factor(results$method)
results$error <- as.numeric(results$error)
results$method <- factor(results$method, levels=c("structure_learning", "regression", "tree-based", "propensity_score"))
levels(results$method) <- c("Structure\nlearning", "Regression", "Tree-based", "Propensity\nscore")

#write.csv(results, "results/postgres_apo_vs_rct_40_runs_log_outcome.csv")
#results <- read.csv("results/postgres_apo_vs_rct_40_runs_log_outcome.csv")
p1 <- ggplot(results[results$measure == "ATE",]) + geom_violin(aes(x=method, y=error, fill=style), trim=FALSE) + ggtitle("MAE of ATE estimates") + 
  ylab("Mean error") + xlab("Inference method") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2 <- ggplot(results[results$measure == "outcome",]) + geom_violin(aes(x=method, y=error, fill=style), trim=FALSE) + ggtitle("Mean error of outcome estimates") +
  ylab("Mean absolute error") + xlab("Inference method") + theme_bw() + theme(legend.position = c(0.85, 0.8), legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
          plot.title = element_text(hjust = 0.5),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- egg::ggarrange(p1, p2, nrow = 1)
p
ggsave("results/postgres-100-runs-one-direction-violin-fixed-variance-final.png", p, width = 18, height = 7, units = "cm", dpi = 600)

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- cbind(g1, g2, size="first")
g$widths <- unit.pmax(g1$widths, g2$widths)
grid.newpage()
grid.draw(g)

png("results/postgres-40-runs-any-direction.png")
grid.arrange(p1, p2, nrow=1, widths=c(1,1))
dev.off()

png("results/postgres-ATE-40-runs-any-direction.png")
p1
dev.off()

png("results/posgres-outcome-40-runs-any-direction.png")
p2
dev.off()