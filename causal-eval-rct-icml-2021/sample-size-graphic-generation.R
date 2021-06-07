outcome_ranges <- read.csv("~/Documents/workspace/causal-eval-rct/outcome_ranges.csv")

#plotting and analysis of results
#dataset_names <- c("Dryad_6d4p06g", "Dryad_B8KG77", "ICPSR_23980", "ISPS_d037", "ISPS_d084", "ISPS_d113", "KNB_1596312", "KNB_f1qf8r51t", "NIDA_P1S1", "UK_Data_Service_853369", "UK_Data_Service_854092")
#dataset_name_mapping <- c("D-1", "D-2", "C-1", "I-1", "I-2", "I-3", "K-1", "K-2", "N-1", "U-1", "U-2")

dataset_names <- c("acic_27", "acic_4", "acic_47", "acic_65", "acic_71", "Dryad_4f4qrfj95", "Dryad_B8KG77", "HDV WT4I9N", "ibm_19e667b985624159bae940919078d55f", "ibm_1b50aae9f0e34b03bdf03ac195a5e7e9", "ibm_2b6d1d419de94f049d98c755beea4ae2", "ibm_7510d73712fe40588acdb129ea58339b", "ibm_c55cbee849534815ba80980975c4340b", "ICPSR_20160213", "ICPSR_23980", "ISPS_d037", "ISPS_d084", "ISPS_d113", "jdk", "KNB_1596312", "KNB_f1qf8r51t", "musiclab", "Nemo-breeding", "Nemo-delet-model", "Nemo-dispersal-rate", "networking", "neuropathic-pain-10000-sample1", "neuropathic-pain-10000-sample2", "neuropathic-pain-10000-sample3", "NIDA_P1S1", "postgres", "UK_Data_Service_852874", "UK_Data_Service_853369", "UK_Data_Service_854092", "Whynot_opioid", "Whynot_world2", "Whynot_zika")
dataset_name_mapping <- c("SR-1", "SR-2", "SR-3", "SR-4", "SR-5", "RCT-1", "RCT-2", "RCT-3", "SR-6", "SR-7", "SR-8", "SR-9", "SR-10", "RCT-4", "RCT-5", "RCT-6", "RCT-7", "RCT-8", "APO-1", "RCT-9", "RCT-10", "RCT-11", "Sim-1", "Sim-2", "Sim-3", "APO-2", "Sim-4", "Sim-5", "Sim-6", "RCT-12", "APO-3", "RCT-13", "RCT-14", "RCT-15", "Sim-7", "Sim-8", "Sim-9")


#columns are dataset, algorithm, measure
df.risk_ratio <- c()
df.ATE <- c()
df.outcome.binary <- c()
df.outcome.continuous <- c()
for (sample_size in seq(1000,10000,1000)){
  load(paste0("all-data-results/final-results/simulator-sample-sizes/run-30-sample-", sample_size, "-bias-1-results.RData"))
  load(paste0("all-data-results/final-results/simulator-sample-sizes/run-30-sample-", sample_size, "-bias-1-binary_datasets.RData"))
  
  for (dataset in names(results)){
    algs <- names(results[[dataset]])
    true_ATE <- results[[dataset]][["true_ATE"]]
    algs <- algs[!(algs == "true_ATE")]
    
    new_dataset_name <- dataset_name_mapping[which(dataset_names == dataset)]
    print(paste0("Mapping ", dataset, " to ", new_dataset_name))
    for (alg in algs){
      ATEs <- results[[dataset]][[alg]][["ATE"]]
      outcome_estimates <- results[[dataset]][[alg]][["outcome"]]
      
      for (ATE in ATEs){
        if (binary_datasets[[dataset]]){
          df.risk_ratio <- rbind(df.risk_ratio, c(new_dataset_name, alg, ATE-true_ATE, sample_size))
        }else{
          df.ATE <- rbind(df.ATE, c(new_dataset_name, alg, ATE-true_ATE, sample_size))
        }
      }
      for (outcome in outcome_estimates){
        if(binary_datasets[[dataset]]){
          df.outcome.binary <- rbind(df.outcome.binary, c(new_dataset_name, alg, outcome, sample_size))
        }else{
          df.outcome.continuous <- rbind(df.outcome.continuous, c(new_dataset_name, alg, outcome, sample_size))
        }
      }
    }
  }
  colnames(df.risk_ratio) <- c("dataset", "algorithm", "error", "sample_size")
  colnames(df.ATE) <- c("dataset", "algorithm", "error", "sample_size")
  colnames(df.outcome.binary) <- c("dataset", "algorithm", "error", "sample_size")
  colnames(df.outcome.continuous) <- c("dataset", "algorithm", "error", "sample_size")
  df.ATE <- data.frame(df.ATE)
  df.risk_ratio <- data.frame(df.risk_ratio)
  df.outcome.binary <- data.frame(df.outcome.binary)
  df.outcome.continuous <- data.frame(df.outcome.continuous)
  df.risk_ratio$error <- as.numeric(df.risk_ratio$error)
  df.ATE$error <- as.numeric(df.ATE$error)
  df.outcome.binary$error <- as.numeric(df.outcome.binary$error)
  df.outcome.continuous$error <- as.numeric(df.outcome.continuous$error)
}



df.ATE <- df.ATE[df.ATE$algorithm != "bartCause",]
df.risk_ratio <- df.risk_ratio[df.risk_ratio$algorithm != "bartCause",]
df.outcome.binary <- df.outcome.binary[df.outcome.binary$algorithm != "bartCause",]
df.outcome.continuous <- df.outcome.continuous[df.outcome.continuous$algorithm != "bartCause",]

df.risk_ratio$dataset <- factor(df.risk_ratio$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6"))
df.ATE$dataset <- factor(df.ATE$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9"))
df.outcome.binary$dataset <- factor(df.outcome.binary$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6"))
df.outcome.continuous$dataset <- factor(df.outcome.continuous$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9"))

df.risk_ratio$algorithm <- as.factor(df.risk_ratio$algorithm)
df.ATE$algorithm <- as.factor(df.ATE$algorithm)
df.outcome.binary$algorithm <- as.factor(df.outcome.binary$algorithm)
df.outcome.continuous$algorithm <- as.factor(df.outcome.continuous$algorithm)


new_names <- c("Naive", "PSM", "IPTW", "OR", "BART", "BART-cause", "CF", "DRE")
df.risk_ratio$algorithm <- factor(df.risk_ratio$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
levels(df.risk_ratio$algorithm) <- new_names
df.ATE$algorithm <- factor(df.ATE$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
levels(df.ATE$algorithm) <- new_names
df.outcome.binary$algorithm <- factor(df.outcome.binary$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
levels(df.outcome.binary$algorithm) <- new_names
df.outcome.continuous$algorithm <- factor(df.outcome.continuous$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
levels(df.outcome.continuous$algorithm) <- new_names

df.ATE <- df.ATE[!is.na(df.ATE$error),]

for (dataset in unique(df.ATE$dataset)){
  unmapped_dataset <- dataset_names[which(dataset_name_mapping == dataset)]
  outcome_range <- outcome_ranges[outcome_ranges$V1 == unmapped_dataset, "V2"]
  df.ATE[df.ATE$dataset == dataset, "error"] <- df.ATE[df.ATE$dataset == dataset, "error"]/outcome_range
}

for (dataset in unique(df.outcome.continuous$dataset)){
  unmapped_dataset <- dataset_names[which(dataset_name_mapping == dataset)]
  outcome_range <- outcome_ranges[outcome_ranges$V1 == unmapped_dataset, "V2"]
  df.outcome.continuous[df.outcome.continuous$dataset == dataset, "error"] <- df.outcome.continuous[df.outcome.continuous$dataset == dataset, "error"]/outcome_range
}

#need to get the lower and upper bounds within each set of dataset/algorithm/sample size
df.ATE$error_lower <- 0
df.ATE$error_higher <- 0
for (curr_dataset in unique(df.ATE$dataset)){
  for (alg in unique(df.ATE$algorithm)){
    for (size in unique(df.ATE$sample_size)){
      df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error_lower"] <- min(df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error"])
      df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error_higher"] <- max(df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error"])
      df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error_mean"] <- mean(df.ATE[df.ATE$dataset == curr_dataset & df.ATE$algorithm == alg & df.ATE$sample_size == size, "error"])
    }
  }
}
df.ATE <- df.ATE[,colnames(df.ATE) != "error"]
df.ATE <- unique(df.ATE)
df.ATE$sample_size <- as.numeric(df.ATE$sample_size)

df.risk_ratio$error_lower <- 0
df.risk_ratio$error_higher <- 0
for (curr_dataset in unique(df.risk_ratio$dataset)){
  for (alg in unique(df.risk_ratio$algorithm)){
    for (size in unique(df.risk_ratio$sample_size)){
      df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error_lower"] <- min(df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error"])
      df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error_higher"] <- max(df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error"])
      df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error_mean"] <- mean(df.risk_ratio[df.risk_ratio$dataset == curr_dataset & df.risk_ratio$algorithm == alg & df.risk_ratio$sample_size == size, "error"])
    }
  }
}
df.risk_ratio <- df.risk_ratio[,colnames(df.ATE) != "error"]
df.risk_ratio <- unique(df.risk_ratio)
df.risk_ratio$sample_size <- as.numeric(df.risk_ratio$sample_size)

for (alg in unique(df.ATE$algorithm[df.ATE$algorithm != "Naive"])){
#  alg <- "CF"
  p <- ggplot(df.ATE[df.ATE$algorithm == alg,]) + geom_ribbon(aes(x=sample_size, ymin=error_lower, ymax=error_higher, fill=dataset), alpha=0.4) + 
    geom_line(aes(x=sample_size, y=error_mean, color=dataset)) + geom_hline(yintercept=0, linetype="dashed") + xlab("sample size") + ylab("Error") + guides(fill=guide_legend(title="Data set"))
#  p
  ggsave(paste0("all-data-results/final-results/simulator-sample-sizes/ATE_", alg, ".png"), p, width=22, height=18, units="cm", dpi=600)
  
  p <- ggplot(df.risk_ratio[df.risk_ratio$algorithm == alg,]) + geom_ribbon(aes(x=sample_size, ymin=error_lower, ymax=error_higher, fill=dataset), alpha=0.4) + 
    geom_line(aes(x=sample_size, y=error_mean, color=dataset)) + geom_hline(yintercept=0, linetype="dashed") + xlab("sample size") + ylab("Error") + guides(fill=guide_legend(title="Data set"))
#  p
  ggsave(paste0("all-data-results/final-results/simulator-sample-sizes/rd_", alg, ".png"), p, width=22, height=18, units="cm", dpi=600)
}



#for (dataset in names(results.apo)){
#  results[[dataset]] <- results.apo[[dataset]]
#}





#for (sample_size in seq(1000,10000,1000)){
#  load(paste0("all-data-results/final-results/simulator-sample-sizes/replacement-apo-subsampling/run-30-sample-", sample_size, "-bias-1-Nemo-results.RData"))
#  results.nemo <- results
#  
#  load(paste0("all-data-results/final-results/simulator-sample-sizes/invalid-apo-subsampling/run-30-sample-", sample_size, "-bias-1-results.RData"))
#  
#  for (dataset in names(results.nemo)){
#    results[[dataset]] <- results.nemo[[dataset]]
#  }
#  save(results, file=paste0("all-data-results/final-results/simulator-sample-sizes/run-30-sample-", sample_size, "-bias-1-results.RData"))
#}
