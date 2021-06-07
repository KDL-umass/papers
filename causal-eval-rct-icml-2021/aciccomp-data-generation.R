if (require("remotes", quietly = TRUE) == FALSE) {
  install.packages("remotes")
  require("remotes")
}
remotes::install_github("vdorie/aciccomp/2016")

library(aciccomp2016)

parameterNum <- 1 #1-77
simulationNum <- 1 #1-100
data <- dgp_2016(input_2016, parameterNum, simulationNum, extraInfo=TRUE)

#we'll use parameter sets:
#4 (polynomial treatment, 35% treated, exponential outcome, 75% alignment, high te hetero)
#27 (polynomial treatment, 35% treated, step outcome, 25% alignment, medium te hetero)
#47 (polynomial treatment, 65% treated, expoential outcome, 75% alignment, high te hetero)
#65 (step treatment, 65% treated, step outcome, 75% alignment, medium te hetero)
#71 (step treatment, 65% treated, step outcome, 25% alignment, high te hetero)
params <- c(4, 27, 47, 65, 71)

for (paramNum in params){
  data <- dgp_2016(input_2016, paramNum, 1, extraInfo=TRUE)
  df <- data.frame(treatment=data$z, outcome=data$y)
  df[,names(data$x)] <- data$x
  df$treatment <- ifelse(df$treatment == "trt", 1, 0)
  df$counterfactual_outcome_1 <- data$y.1
  df$counterfactual_outcome_0 <- data$y.0
  write.csv(df, paste0("data/acic/acic_", paramNum, ".csv"))
}
