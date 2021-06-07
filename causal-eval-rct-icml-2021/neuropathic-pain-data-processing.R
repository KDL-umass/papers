setwd("~/Documents/workspace/causal-eval-rct")
df <- read.csv("neuropathic-pain-10000-sample2.csv")
network <- read.csv("neuropathic-pain-full-network.csv")

translation <- read.csv("neuropathic-pain-translation.csv")

tier_1_vars <- colnames(df)[1:26]
tier_2_vars <- colnames(df)[27:79]
tier_3_vars <- colnames(df)[80:222]

#translate the column names from Swedish to English (translations from Google Translate)
#english_cols <- c()
#for (col in colnames(df)){
#  english_cols <- c(english_cols, translation[translation$Swedish == col, "English"])
#}
#colnames(df) <- english_cols

num_top_level <- c()
for (row in 1:nrow(df)){
  num_top_level <- c(num_top_level, sum(df[row, 1:26]))
}

for (i in 1:26){
  print(colnames(df)[i])
  print(count(df[,i]))
}
#reasonably even treatment groups:
  #DLS.C5.C6
  #DLS.L4.L5 (t1)
  #DLS.L5.S1 (b1)
#at least 1000 in each group:
  #DLS.C2.C3
#at least 2000 in each group:
  #DLS.C3.C4
  #DLS.C4.C5
  #DLS.C6.C7
for (col in colnames(df)[1:26])
  print(paste0(col, ": ", sum(df[,col])))

#if we want top-level variables that result in approximately equal treatment group sizes, could use "DLS.L4.L5 or DLS.L5.S1

#unroll the network from DLS.L4.L5 and DLS.L5.S1
get_effects <- function(root, network){
  descendents <- network[network$From == root, "To"]
  return(c(descendents, network[network$From %in% descendents, "To"]))
}

DLS.L4.L5_effects <- get_effects("DLS.L4.L5", network)
DLS.L5.S1_effects <- get_effects("DLS.L5.S1", network)

DLS.C5.C6_effects <- get_effects("DLS.C5.C6", network)
DLS.C3.C4_effects <- get_effects("DLS.C3.C4", network)
DLS.C4.C5_effects <- get_effects("DLS.C4.C5", network)
DLS.C6.C7_effects <- get_effects("DLS.C6.C7", network)
DLS.C2.C3_effects <- get_effects("DLS.C2.C3", network)


common_effects <- union(DLS.L4.L5_effects, DLS.L5.S1_effects)
tier_3_common_effects <- common_effects[common_effects %in% tier_3_vars]

common_effects <- union(DLS.C5.C6_effects, DLS.C3.C4_effects)
tier_3_common_effects <- common_effects[common_effects %in% tier_3_vars]

common_effects <- union(DLS.C6.C7_effects, DLS.C4.C5_effects)
tier_3_common_effects <- common_effects[common_effects %in% tier_3_vars]

#test the correlation between each of the tier 1 vars and their common effects
cor1s <- c()
cor2s <- c()
for (curr_var in tier_3_common_effects){
  if (var(df[,curr_var]) != 0){
    cor1 <- cor(df$DLS.C6.C7, df[,curr_var])
    cor2 <- cor(df$DLS.C4.C5, df[,curr_var])
    cor1s <- c(cor1s, cor1)
    cor2s <- c(cor2s, cor2)
    if (cor1 > .1 && cor2 > .1)
      print(paste(curr_var, cor1, cor2))
  }
}

#treatment = DLS.L4.L5
#outcome = Lombago
#bias covar = DLS.L5.S1
#other covars = all other tier 1 vars


#treatment = DLS.C5.C6
#outcome = R.Skulderbesvär
#bias_covar = DLS.C3.C4
#other covars = all other tier 1 vars

#treatment = DLS.C4.C5
#outcome = R.Axelbesvär
#bias_covar = DLS.C6.C7
#other covars = all other tier 1 vars