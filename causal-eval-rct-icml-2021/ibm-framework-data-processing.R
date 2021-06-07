params <- read.csv("IBM causal inference benchmarking framework/scaling/params.csv")

params <- params[params$size == 10000,]

#we'll try to get diversity on X._effect_size, effect_size, link_type, tx_importance, deg.y, deg.z, and dgp

#X_effect size: [-401.90, 611.4]
#effect size: [-12.8, 26]
#link type: (exp, log, poly)
#tx importance: (0.2, 1.4) or NA
#deg.y: (1,101)
#deg.z: (1.96)
#dgp: (0,63)

params[sample(1:nrow(params), 10, ), c("X._effect_size", "effect_size", "link_type", "tx_importance", "deg.y.", "deg.z.", "dgp")]

#param nums:
#1461
#1389
#1720
#1695
#1517
paramIDs <- params[c("1461", "1389", "1720", "1695", "1517"), "ufid"]

for (paramID in paramIDs){
  curr_params <- params[params$ufid == paramID, ]
  factuals <- read.csv(paste0("IBM causal inference benchmarking framework/scaling/factuals/", paramID, ".csv"))
  counterfactuals <- read.csv(paste0("IBM causal inference benchmarking framework/scaling/counterfactuals/", paramID, "_cf.csv"))
  covars <- read.csv("IBM causal inference benchmarking framework/x.csv")
  df <- data.frame(sample_id = factuals$sample_id, treatment = factuals$z, outcome = factuals$y)
  df <- merge(df, counterfactuals, by="sample_id")
  df <- rename(df, c("y0"="counterfactual_outcome_0", "y1"="counterfactual_outcome_1"))
  df <- merge(df, covars, by="sample_id")
  df <- df[,(colnames(df) != "sample_id")]
  
  write.csv(df, paste0("data/ibm/", paramID, ".csv"))
}