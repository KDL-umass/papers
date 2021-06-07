library(RColorBrewer)
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/bias-1/initial-30-run-bias-1-with-nn-results.RData"))
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/bias-1/initial-30-run-bias-1-with-nn-binary_datasets.Rdata"))

#bias <- 1
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/two-bias/bias-", bias, "/run-30-bias-", bias, "-two-bias-with-nn-results.RData"))
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/two-bias/bias-", bias, "/run-30-bias-", bias, "-two-bias-with-nn-binary_datasets.Rdata"))

#bias <- 1
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/bias-1/run-30-bias-", bias, "-with-nn-results.RData"))
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/bias-1/run-30-bias-", bias, "-with-nn-binary_datasets.Rdata"))

#bias <- 1
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/1000-sample-size/bias-", bias, "/run-30-bias-", bias, "-1000-samples-results.RData"))
#load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/1000-sample-size/bias-", bias, "/run-30-bias-", bias, "-1000-samples-binary_datasets.RData"))
load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/final-results/run-30-bias-5-apo-two-bias-results.RData"))
results.apo <- results


simplified_naming <- TRUE #if true, abbreviate all algorithms far shorter
final_paper_graphics <- TRUE #whether to do crazy formatting for the paper
remove_bartCause <- TRUE
bias <- 5
nn_included <- FALSE
remove_nn <- TRUE #use when the data you read in contains the nn method, but we don't want to plot it
includes_acic_ibm <- FALSE
folder <- paste0("final-results/5000-sample-size/two-bias/bias-", bias)
#folder <- paste0("final-results/2000-sample-size-with-acic/bias-", bias)
suffix <- ""
load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/", folder, "/run-30-bias-", bias, suffix, "-results.RData"))
load(file=paste0("~/Documents/workspace/causal-eval-rct/all-data-results/", folder, "/run-30-bias-", bias, suffix, "-binary_datasets.RData"))

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
        df.risk_ratio <- rbind(df.risk_ratio, c(new_dataset_name, alg, ATE-true_ATE))
      }else{
        df.ATE <- rbind(df.ATE, c(new_dataset_name, alg, ATE-true_ATE))
      }
    }
    for (outcome in outcome_estimates){
      if(binary_datasets[[dataset]]){
        df.outcome.binary <- rbind(df.outcome.binary, c(new_dataset_name, alg, outcome))
      }else{
        df.outcome.continuous <- rbind(df.outcome.continuous, c(new_dataset_name, alg, outcome))
      }
    }
  }
}
colnames(df.risk_ratio) <- c("dataset", "algorithm", "error")
colnames(df.ATE) <- c("dataset", "algorithm", "error")
colnames(df.outcome.binary) <- c("dataset", "algorithm", "error")
colnames(df.outcome.continuous) <- c("dataset", "algorithm", "error")
df.ATE <- data.frame(df.ATE)
df.risk_ratio <- data.frame(df.risk_ratio)
df.outcome.binary <- data.frame(df.outcome.binary)
df.outcome.continuous <- data.frame(df.outcome.continuous)
df.risk_ratio$error <- as.numeric(df.risk_ratio$error)
df.ATE$error <- as.numeric(df.ATE$error)
df.outcome.binary$error <- as.numeric(df.outcome.binary$error)
df.outcome.continuous$error <- as.numeric(df.outcome.continuous$error)


if (remove_bartCause){
  df.ATE <- df.ATE[df.ATE$algorithm != "bartCause",]
  df.risk_ratio <- df.risk_ratio[df.risk_ratio$algorithm != "bartCause",]
  df.outcome.binary <- df.outcome.binary[df.outcome.binary$algorithm != "bartCause",]
  df.outcome.continuous <- df.outcome.continuous[df.outcome.continuous$algorithm != "bartCause",]
}

#use these when all data sets are included
if (includes_acic_ibm){
  df.risk_ratio$dataset <- factor(df.risk_ratio$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6", "RCT-1", "RCT-3", "RCT-4", "RCT-5", "RCT-6", "RCT-7", "RCT-8", "RCT-10", "RCT-13", "RCT-14"))
  df.ATE$dataset <- factor(df.ATE$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9", "SR-1", "SR-2", "SR-3", "SR-4", "SR-5", "SR-6", "SR-7", "SR-8", "SR-9", "SR-10", "RCT-2", "RCT-9", "RCT-11", "RCT-12", "RCT-15", "APO-1", "APO-2", "APO-3"))
  df.outcome.binary$dataset <- factor(df.outcome.binary$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6", "RCT-1", "RCT-3", "RCT-4", "RCT-5", "RCT-6", "RCT-7", "RCT-8", "RCT-10", "RCT-13", "RCT-14"))
  df.outcome.continuous$dataset <- factor(df.outcome.continuous$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9", "SR-1", "SR-2", "SR-3", "SR-4", "SR-5", "SR-6", "SR-7", "SR-8", "SR-9", "SR-10", "RCT-2", "RCT-9", "RCT-11", "RCT-12", "RCT-15", "APO-1", "APO-2", "APO-3"))
} else{
  #use these when IBM and ACIC aren't included
  df.risk_ratio$dataset <- factor(df.risk_ratio$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6", "RCT-1", "RCT-3", "RCT-4", "RCT-5", "RCT-6", "RCT-7", "RCT-8", "RCT-10", "RCT-13", "RCT-14"))
  df.ATE$dataset <- factor(df.ATE$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9", "RCT-2", "RCT-9", "RCT-11", "RCT-12", "RCT-15", "APO-1", "APO-2", "APO-3"))
  df.outcome.binary$dataset <- factor(df.outcome.binary$dataset, levels=c("Sim-3", "Sim-4", "Sim-5", "Sim-6", "RCT-1", "RCT-3", "RCT-4", "RCT-5", "RCT-6", "RCT-7", "RCT-8", "RCT-10", "RCT-13", "RCT-14"))
  df.outcome.continuous$dataset <- factor(df.outcome.continuous$dataset, levels=c("Sim-1", "Sim-2", "Sim-7", "Sim-8", "Sim-9", "RCT-2", "RCT-9", "RCT-11", "RCT-12", "RCT-15", "APO-1", "APO-2", "APO-3"))
}

df.risk_ratio$algorithm <- as.factor(df.risk_ratio$algorithm)
df.ATE$algorithm <- as.factor(df.ATE$algorithm)
df.outcome.binary$algorithm <- as.factor(df.outcome.binary$algorithm)
df.outcome.continuous$algorithm <- as.factor(df.outcome.continuous$algorithm)

if (!nn_included){
  if (simplified_naming){
    new_names <- c("Naive", "PSM", "IPTW", "OR", "BART", "BART-cause", "CF", "DRE")
  }else{
    new_names <- c("Naive", "PS\nmatching", "IPTW", "Outcome\nregression", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust")
  }
  df.risk_ratio$algorithm <- factor(df.risk_ratio$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
  levels(df.risk_ratio$algorithm) <- new_names
  df.ATE$algorithm <- factor(df.ATE$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
  levels(df.ATE$algorithm) <- new_names
  df.outcome.binary$algorithm <- factor(df.outcome.binary$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
  levels(df.outcome.binary$algorithm) <- new_names
  df.outcome.continuous$algorithm <- factor(df.outcome.continuous$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust"))
  levels(df.outcome.continuous$algorithm) <- new_names
} else{
  if (simplified_naming){
    new_names <- c("Naive", "PSM", "IPTW", "OR", "BART", "BART-cause", "CF", "DRE", "NN")
  }else{
    new_names <- c("Naive", "PS\nmatching", "IPTW", "Outcome\nregression", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust", "Neural\nnetwork")
  }
  df.risk_ratio$algorithm <- factor(df.risk_ratio$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "neural_network"))
  levels(df.risk_ratio$algorithm) <- new_names
  df.ATE$algorithm <- factor(df.ATE$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "neural_network"))
  levels(df.ATE$algorithm) <- new_names
  df.outcome.binary$algorithm <- factor(df.outcome.binary$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "neural_network"))
  levels(df.outcome.binary$algorithm) <- new_names
  df.outcome.continuous$algorithm <- factor(df.outcome.continuous$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "neural_network"))
  levels(df.outcome.continuous$algorithm) <- new_names
}



#df.risk_ratio$algorithm <- factor(df.risk_ratio$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle"))
#levels(df.risk_ratio$algorithm) <- c("Naive", "PS matching", "IPTW", "Outcome\nmodeling", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust", "Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle")
#df.ATE$algorithm <- factor(df.ATE$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle"))
#levels(df.ATE$algorithm) <- c("Naive", "PS matching", "IPTW", "Outcome\nmodeling", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust", "Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle")
#df.outcome.binary$algorithm <- factor(df.outcome.binary$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle"))
#levels(df.outcome.binary$algorithm) <- c("Naive", "PS matching", "IPTW", "Outcome\nmodeling", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust", "Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle")
#df.outcome.continuous$algorithm <- factor(df.outcome.continuous$algorithm, levels=c("naive", "propensity_score", "ipw", "regression", "bart", "bartCause", "causal-forest", "doubly-robust", "dragonnet", "dragonnet_tmle", "tarnet", "tarnet_tmle"))
#levels(df.outcome.continuous$algorithm) <- c("Naive", "PS matching", "IPTW", "Outcome\nmodeling", "BART", "BART-cause", "Causal\nforest", "Doubly\nrobust", "Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle")
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


#ggplot(df.risk_ratio) + geom_bar(aes(x=dataset,y=estimate,fill=algorithm), stat="identity", position=position_dodge())
#png("results/Justin-results/6-datasets-30-risk-ratio.png")
#ggplot(df.risk_ratio) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm))+ ylab("mean error") + xlab("data set")
#dev.off()
#png("results/Justin-results/6-datasets-30-ATE.png")
#ggplot(df.ATE) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm)) + ylab("normalized mean error") + xlab("data set")
#dev.off()
#png("results/Justin-results/6-datasets-30-binary-outcome.png")
#ggplot(df.outcome.binary) + geom_boxplot(aes(x=dataset, y=error, fill=algorithm)) + ylab("mean absolute error") + xlab("data sets")
#dev.off()
#png("results/Justin-results/6-datasets-30-continuous-outcome.png")
#ggplot(df.outcome.continuous) + geom_boxplot(aes(x=dataset, y=error, fill=algorithm)) + xlab("data sets") + ylab("normalized mean absolute error")
#dev.off()

#p1 <- ggplot(df.risk_ratio) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Mean error") + xlab("Data set") + ggtitle("ATE for binary outcome") + theme_bw() + theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
#                                                                                                                                                                                                                        plot.title = element_text(hjust = 0.5),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p2 <- ggplot(df.ATE) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Normalized mean error") + xlab("Data set") + ggtitle("ATE for continuous outcome") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p3 <- ggplot(df.outcome.binary) + geom_violin(aes(x=dataset, y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Mean absolute error") + xlab("Data sets") + ggtitle("Outcome estimation for binary outcome") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                                      panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())# + scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF"))
#p4 <- ggplot(df.outcome.continuous) + geom_violin(aes(x=dataset, y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + xlab("Data sets") + ylab("Normalized mean absolute error") + ggtitle("Outcome estimation for continuous outcome") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                                                         panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())# + scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF"))

#p <- egg::ggarrange(p1, p2, p3, p4, nrow = 2)
#ggsave("results/datasets-all-11-violin-one-direction-correct-normalizing.png", p, width = 22, height = 18, units = "cm", dpi = 600)

if (remove_nn){
  df.risk_ratio <- df.risk_ratio[df.risk_ratio$algorithm != "Neural\nnetwork" & df.risk_ratio$algorithm != "NN",]
  df.ATE <- df.ATE[df.ATE$algorithm != "Neural\nnetwork" & df.ATE$algorithm != "NN",]
  df.outcome.binary <- df.outcome.binary[df.outcome.binary$algorithm != "Neural\nnetwork" & df.outcome.binary$algorithm != "NN",]
  df.outcome.continuous <- df.outcome.continuous[df.outcome.continuous$algorithm != "Neural\nnetwork" & df.outcome.continuous$algorithm != "NN",]
}

clear_gg_formatting <- function(p, algorithms=c(), zero_line=FALSE, legend=FALSE, simplified_naming=TRUE){
  if (simplified_naming){
    palette_map <- c("Naive", "PSM", "IPTW", "OR", "BART", "CF", "DRE", "NN")
  } else{
    palette_map <- c("Naive", "PS\nmatching", "IPTW", "Outcome\nregression", "BART", "Causal\nforest", "Doubly\nrobust", "Neural\nnetwork")
  }
  palette <- c("#45C19C", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#990000")
  
  if (length(algorithms) > 1){#make sure the color-algorithm mapping is consistent, no matter which algorithms are actually in this graphic
    palette <- palette[which(palette_map %in% algorithms)]
  }
  if (legend){
    p <- p + theme_bw() + theme(axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=14), axis.text.x = element_text(angle = 90), legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }else{
    p <- p + theme_bw() + theme(axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=14), axis.text.x = element_text(angle = 90), legend.position = "none", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  if (zero_line)
    p <- p + geom_hline(yintercept=0, linetype="dashed")
#  p <- p + scale_fill_brewer(palette="Dark2") + guides(fill=guide_legend(nrow=1,byrow=TRUE))
  p <- p + scale_fill_manual(values=palette) + guides(fill=guide_legend(nrow=1,byrow=TRUE))
  return(p)
}

if (final_paper_graphics){
  if (includes_acic_ibm){
      p <- ggplot(df.risk_ratio) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Error") + xlab("Data set")
      p <- clear_gg_formatting(p, algorithms=unique(df.risk_ratio$algorithm), zero_line=TRUE, legend=TRUE, simplified_naming=simplified_naming)
      p <- p + geom_vline(xintercept=4.5, linetype="dashed")
      p <- p + coord_cartesian(ylim=c(-0.4, .27))#one bias covar
      p <- p + coord_cartesian(ylim=c(-0.3, .27))#with nn
      p <- p + coord_cartesian(ylim=c(-0.18, .23))#2000 sample size
      ggsave(paste0("all-data-results/", folder, "/run-30-risk-diff-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
      
      p <- ggplot(df.ATE) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Normalized error") + xlab("Data set")
      p <- clear_gg_formatting(p, algorithms=unique(df.ATE$algorithm), zero_line=TRUE, legend=TRUE, simplified_naming=simplified_naming)
      p <- p + geom_vline(xintercept=5.5, linetype="dashed") + geom_vline(xintercept=15.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed")
      if (!("SR-6" %in% df.ATE$dataset))#need different lines if it's just acic, not ibm
        p <- p + geom_vline(xintercept=10.5, linetype="dashed")
      p <- p + coord_cartesian(ylim=c(-.15, .15))#one bias covar
      p <- p + coord_cartesian(ylim=c(-0.05, 0.06))#2000 sample size
      p <- p + coord_cartesian(ylim=c(-.8,.8))#with nn
      ggsave(paste0("all-data-results/", folder, "/run-30-ATE-bias-", bias, ".png"), p, width=44, height=18, units="cm", dpi=600)
      
      p <- ggplot(df.outcome.binary) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Absolute error") + xlab("Data set")
      p <- clear_gg_formatting(p, algorithms=unique(df.outcome.binary$algorithm), zero_line=FALSE, legend=TRUE, simplified_naming=simplified_naming)
      p <- p + geom_vline(xintercept=4.5, linetype="dashed")
      ggsave(paste0("all-data-results/", folder, "/run-30-outcome-binary-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
      
      p <- ggplot(df.outcome.continuous) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Normalized absolute error") + xlab("Data set")
      p <- clear_gg_formatting(p, algorithms=unique(df.outcome.continuous$algorithm), zero_line=FALSE, legend=TRUE, simplified_naming=simplified_naming)
      p <- p + geom_vline(xintercept=5.5, linetype="dashed") + geom_vline(xintercept=15.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed")
      p <- p + coord_cartesian(ylim=c(0, .6))#with nn
      ggsave(paste0("all-data-results/", folder, "/run-30-outcome-continuous-bias-", bias, ".png"), p, width=44, height=18, units="cm", dpi=600)
  } else{
    p <- ggplot(df.risk_ratio) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Error") + xlab("Data set")
    p <- clear_gg_formatting(p, algorithms=unique(df.risk_ratio$algorithm), zero_line=TRUE, legend=TRUE, simplified_naming=simplified_naming)
    p <- p + geom_vline(xintercept=4.5, linetype="dashed")
    p <- p + coord_cartesian(ylim=c(-0.35, 0.3))#2 bias covarm, bias=1
    p <- p + coord_cartesian(ylim=c(-0.5, 0.55))#2 bias covars, bias=5
    ggsave(paste0("all-data-results/", folder, "/run-30-risk-diff-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
    
    p <- ggplot(df.ATE) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Normalized error") + xlab("Data set")
    p <- clear_gg_formatting(p, algorithms=unique(df.ATE$algorithm), zero_line=TRUE, legend=TRUE, simplified_naming=simplified_naming)
    p <- p + geom_vline(xintercept=5.5, linetype="dashed") + geom_vline(xintercept=10.5, linetype="dashed")
    p <- p + coord_cartesian(ylim=c(-.05, .1))#2 bias covar, bias=1
    p <- p + coord_cartesian(ylim=c(-.15, .19))#2 bias covar, bias=5
    ggsave(paste0("all-data-results/", folder, "/run-30-ATE-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
    
    p <- ggplot(df.outcome.binary) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Absolute error") + xlab("Data set")
    p <- clear_gg_formatting(p, algorithms=unique(df.outcome.binary$algorithm), zero_line=FALSE, legend=TRUE, simplified_naming=simplified_naming)
    p <- p + geom_vline(xintercept=4.5, linetype="dashed")
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-binary-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
    
    p <- ggplot(df.outcome.continuous) + geom_boxplot(aes(x=dataset,y=error, fill=algorithm), width=.7) + ylab("Normalized absolute error") + xlab("Data set")
    p <- clear_gg_formatting(p, algorithms=unique(df.outcome.continuous$algorithm), zero_line=FALSE, legend=TRUE, simplified_naming=simplified_naming)
    p <- p + geom_vline(xintercept=5.5, linetype="dashed") + geom_vline(xintercept=10.5, linetype="dashed")
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-continuous-bias-", bias, ".png"), p, width=22, height=18, units="cm", dpi=600)
    
  }
} else{#original versions, without ggplot formatting removed
  if (includes_acic_ibm){
    p1 <- ggplot(df.risk_ratio) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-risk-diff-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.ATE[df.ATE$dataset %in% levels(df.ATE$dataset)[c(1:12)],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    p2 <- ggplot(df.ATE[df.ATE$dataset %in% levels(df.ATE$dataset)[12:23],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-ATE-1-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    ggsave(paste0("all-data-results/", folder, "/run-30-ATE-2-bias-", bias, ".png"), p2, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.outcome.binary) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-binary-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.outcome.continuous[df.outcome.continuous$dataset %in% levels(df.outcome.continuous$dataset)[1:11],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    p2 <- ggplot(df.outcome.continuous[df.outcome.continuous$dataset %in% levels(df.outcome.continuous$dataset)[12:23],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-continuous-1-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-continuous-2-bias-", bias, ".png"), p2, width=22, height=18, units="cm", dpi=600)
  } else{
    p1 <- ggplot(df.risk_ratio) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-risk-diff-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.ATE) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-ATE-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.outcome.binary) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-binary-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
    
    p1 <- ggplot(df.outcome.continuous) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
    ggsave(paste0("all-data-results/", folder, "/run-30-outcome-continuous-bias-", bias, ".png"), p1, width=22, height=18, units="cm", dpi=600)
  }
}




#    p1 <- ggplot(df.risk_ratio[df.risk_ratio$algorithm %in% c("Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle"),]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    ggsave("all-data-results/bias-1/initial-30-risk-diff-bias-1-nn.png", p1, width=22, height=18, units="cm", dpi=600)
#    
#    df.ATE.nn <- df.ATE[df.ATE$algorithm %in% c("Dragonnet", "Dragonnet\ntmle", "Tarnet", "Tarnet\ntmle"),]
#    p1 <- ggplot(df.ATE.nn[df.ATE.nn$dataset %in% levels(df.ATE.nn$dataset)[1:11],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    p2 <- ggplot(df.ATE.nn[df.ATE.nn$dataset %in% levels(df.ATE.nn$dataset)[12:23],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    ggsave("all-data-results/bias-1/initial-30-ATE-1-bias-1-nn.png", p1, width=22, height=18, units="cm", dpi=600)
#    ggsave("all-data-results/bias-1/initial-30-ATE-2-bias-1-nn.png", p2, width=22, height=18, units="cm", dpi=600)
#    
#    p1 <- ggplot(df.outcome.binary) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    ggsave("all-data-results/bias-1/initial-30-outcome-binary-bias-1-nn.png", p1, width=22, height=18, units="cm", dpi=600)
#    
#    p1 <- ggplot(df.outcome.continuous[df.outcome.continuous$dataset %in% levels(df.outcome.continuous$dataset)[1:11],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    p2 <- ggplot(df.outcome.continuous[df.outcome.continuous$dataset %in% levels(df.outcome.continuous$dataset)[12:23],]) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE)+theme(axis.text.x = element_text(angle = 90))
#    ggsave("all-data-results/bias-1/initial-30-outcome-continuous-1-bias-1-nn.png", p1, width=22, height=18, units="cm", dpi=600)
#    ggsave("all-data-results/bias-1/initial-30-outcome-continuous-2.bias-1-nn.png", p2, width=22, height=18, units="cm", dpi=600)




#p1 <- ggplot(df.risk_ratio) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Mean error") + xlab("Data set") + ggtitle("ATE for binary outcome") + theme_bw() + theme(legend.position = "none",
#                                                                                                                                                                                                                        plot.title = element_text(hjust = 0.5),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p2 <- ggplot(df.ATE) + geom_violin(aes(x=dataset,y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Normalized mean error") + xlab("Data set") + ggtitle("ATE for continuous outcome") + theme_bw() + theme(legend.position = c(0.2, 0.8), legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p3 <- ggplot(df.outcome.binary) + geom_violin(aes(x=dataset, y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + ylab("Mean absolute error") + xlab("Data sets") + ggtitle("Outcome estimation for binary outcome") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                                      panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())# + scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF"))
#p4 <- ggplot(df.outcome.continuous) + geom_violin(aes(x=dataset, y=error, fill=algorithm), width=.7, scale="width", trim=FALSE) + xlab("Data sets") + ylab("Normalized mean absolute error") + ggtitle("Outcome estimation for continuous outcome") + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
#                                                                                                                                                                                                                                                                         panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())# + scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF"))

#p <- egg::ggarrange(p1, p2, p3, p4, nrow = 2)
#ggsave("results/datasets-all-11-violin-one-direction-correct-normalizing.png", p, width = 22, height = 18, units = "cm", dpi = 600)