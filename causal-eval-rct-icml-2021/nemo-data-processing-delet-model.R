library(ggplot2)
setwd("~/Documents/workspace/causal-eval-rct/Nemo\ data/")

mono <- read.csv("Nemo-2.3.52/example/test/data/MONO_dcost01_ISM.txt", sep="\t")
poly <- read.csv("Nemo-2.3.52/example/test/data/POLY_dcost01_ISM.txt", sep="\t")


mono <- mono[mono$generation == 100,]
poly <- poly[poly$generation == 100,]



ggplot() + geom_density(aes(mono$off.delfreq), color="blue") + geom_density(aes(poly$off.delfreq), color="red")

#possible outcomes:
#off.nbr = number of offspring **
#off.density = average offspring density ***
#extrate = proportion of extinct patches in the popluation *
#off.prop.fsibs = mean viability of inbred individuals between full-sib parents ***
#off.delfreq = mean deleterious mutation frequency ****


all_files <- list.files("results-breeding-on-adlt.delhmz/results")
poly_missing <- c()
mono_missing <- c()
nums <- 1:1000
for (i in nums){
  if (!paste0("POLY_iter", i, ".txt") %in% all_files)
    poly_missing <- c(poly_missing, i)
  if (!paste0("MONO_iter", i, ".txt") %in% all_files)
    mono_missing <- c(mono_missing, i)
}



results_dir <- "results-delet-model"
#results_dir <- "results-breeding-on-adlt.delhmz"

all_files <- list.files(paste0(results_dir, "/results"))
parameters <- read.csv(paste0(results_dir, "/parameters.csv"))
#read in the second row of each file (the final generation) and collect in a dataframe (with an index for matching POs)
df <- c()
one_filenames <- all_files[substring(all_files, 13, 13) == "1"]
params_dataset <- c()
for (i in 1:length(one_filenames)){
  one_filename <- one_filenames[i]
  two_filename <- paste0("delet_model_2_", substring(one_filename, 15, nchar(one_filename)))
  iter_num <- substring(one_filename, 19, nchar(one_filename)-4)
  params <- parameters[parameters$iter == iter_num,]
  one_file <- read.csv(paste0(results_dir, "/results/", one_filename), sep="\t")
  one_file <- one_file[one_file$generation == 100, ]
  two_file <- read.csv(paste0(results_dir, "/results/", two_filename), sep="\t")
  two_file <- two_file[two_file$generation == 100, ]
  if (nrow(one_file) != 0 & nrow(two_file) != 0){
    df <- rbind(df, c(unlist(one_file), unlist(params), "1", i))
    df <- rbind(df, c(unlist(two_file), unlist(params), "2", i))
    params_dataset <- rbind(params_dataset, c(unlist(params), "fit"))
  } else{
    params_dataset <- rbind(params_dataset, c(unlist(params), "unfit"))
  }
  if(length(c(unlist(one_file), unlist(params), "mono", i)) < 98)
    print(i)
}
colnames(df)[(length(colnames(df))-1):length(colnames(df))] <- c("delet_model", "index")
df <- as.data.frame(df)

colnames(params_dataset)[length(colnames(params_dataset))] <- "survive"
params_dataset <- as.data.frame(params_dataset)
for (col in 1:13)
  params_dataset[,col] <- as.numeric(params_dataset[,col])
ggplot(params_dataset) + geom_density(aes(total_capacity, color=survive))

possible_outcomes <- c()
for (outcome in colnames(df)[3:83]){
  for (covar in colnames(df)[85:96]){
    correlation <- abs(cor(as.numeric(df[,outcome]), as.numeric(df[,covar])))
    if (!is.na(correlation) && correlation > 0.1){
      print(paste0(covar, " -> ", outcome, ": ", correlation))
      possible_outcomes <- rbind(possible_outcomes, c(covar, outcome, correlation))
    }
  }
}
colnames(possible_outcomes) <- c("covar", "outcome", "correlation")

print(unique(possible_outcomes[,"outcome"]))
for (val in possible_outcomes[,"outcome"])
  df[,val] <- as.numeric(df[,val])
ggplot(df) + geom_density(aes(adlt.delsegr, color=delet_model))
#adlt.delfreq, adlt.delhmz



#treatment = delet_model
#outcome = adlt.delfreq
#covar = delet_mutation_rate
ggplot(df) + geom_point(aes(delet_mutation_rate, adlt.delfreq))
ggplot(df) + geom_density(aes(adlt.delfreq, color=delet_model))
