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



results_dir <- "results-breeding-10000"

all_files <- list.files(paste0(results_dir, "/results"))
parameters <- read.csv(paste0(results_dir, "/parameters.csv"))
#read in the second row of each file (the final generation) and collect in a dataframe (with an index for matching POs)
df <- c()
mono_filenames <- all_files[substring(all_files, 1, 4) == "MONO"]
params_dataset <- c()
for (i in 1:length(mono_filenames)){
 mono_filename <- mono_filenames[i]
 poly_filename <- paste0("POLY", substring(mono_filename, 5, nchar(mono_filename)))
 iter_num <- substring(substring(mono_filename, 10, nchar(mono_filename)), 1, nchar(mono_filename)-13)
 params <- parameters[parameters$iter == iter_num,]
 mono_file <- read.csv(paste0(results_dir, "/results/", mono_filename), sep="\t")
 mono_file <- mono_file[mono_file$generation == 100, ]
 poly_file <- read.csv(paste0(results_dir, "/results/", poly_filename), sep="\t")
 poly_file <- poly_file[poly_file$generation == 100, ]
 if (nrow(mono_file) != 0 & nrow(poly_file) != 0){
   df <- rbind(df, c(unlist(mono_file), unlist(params), "mono", i))
   df <- rbind(df, c(unlist(poly_file), unlist(params), "poly", i))
   params_dataset <- rbind(params_dataset, c(unlist(params), "fit"))
 } else{
   params_dataset <- rbind(params_dataset, c(unlist(params), "unfit"))
 }
 if(length(c(unlist(mono_file), unlist(params), "mono", i)) < 98)
   print(i)
}
colnames(df)[(length(colnames(df))-1):length(colnames(df))] <- c("breeding", "index")
df <- as.data.frame(df)
colnames(params_dataset)[length(colnames(params_dataset))] <- "survive"
params_dataset <- as.data.frame(params_dataset)
for (col in 1:13)
  params_dataset[,col] <- as.numeric(params_dataset[,col])
ggplot(params_dataset) + geom_density(aes(dispersal_rate, color=survive))


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

print(unique(possible_outcomes))
for (val in possible_outcomes)
  df[,val] <- as.numeric(df[,val])
ggplot(df) + geom_density(aes(off.dvar, color=breeding))

df$off.nbr <- as.numeric(df$off.nbr)
df$off.density <- as.numeric(df$off.density)
df$extrate <- as.numeric(df$extrate)
df$off.prop.fsibs <- as.numeric(df$off.prop.fsibs)
df$off.prop.hsibs <- as.numeric(df$off.prop.hsibs)
df$off.delfreq <- as.numeric(df$off.delfreq)
ggplot(df) + geom_density(aes(off.nbr, color=breeding))
ggplot(df) + geom_density(aes(off.density, color=breeding))
ggplot(df) + geom_density(aes(extrate, color=breeding))
ggplot(df) + geom_density(aes(off.prop.fsibs, color=breeding))
ggplot(df) + geom_density(aes(off.prop.hsibs, color=breeding))
ggplot(df) + geom_density(aes(off.delfreq, color=breeding))
#possible outcomes:
#7, 11, 12, 15**, 16*, 19, 20*, 21*, 22, 26, 27, 30, 31*, 32, 33, 42*, 43*, 45, 53*, 54*, 58*, 59*, 60*, 61*, 62*, 63*, 64**, 66*, 72, 73, 74**, 76, 77*, 81 
possible_outcomes <- colnames(df)[c(7, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 26, 27, 30, 31, 32, 33, 42, 43, 45, 53, 54, 58, 59, 60, 61, 62, 63, 64, 66, 72, 73, 74, 76, 77, 81)]
strong_possible_outcomes <- colnames(df)[c(15, 16, 20, 21, 31, 42, 43, 53, 54, 58, 59, 60, 61, 62, 63, 64, 66, 74, 77)]
#favorites: off.delhtz, off.delsegrp, adlt.delsegrp, off.prop.outb, off,prop.outw, adlt.prop.outb, adlt.prop.outw, off.disp, off.fdisp, off.mdisp, adlt.disp, adlt.fdisp, adlt.mdisp, adlt.allnbp, off.allnbp, off.fixlocp, adlt.fixlocp
df$off.fit <- as.numeric(df$off.fit)
df$off.fst <- as.numeric(df$off.fst)
df$off.fis <- as.numeric(df$off.fis)
ggplot(df) + geom_density(aes(off.fst, color=breeding))

possible_covars <- colnames(parameters)[c(2:6, 8:11)]
for (outcome in possible_outcomes){
  for (covar in possible_covars){
    correlation <- cor(df[,outcome], df[,covar])
    if (!is.na(correlation) &&correlation > 0.05)
      print(paste0(covar, " -> ", outcome, ": ", correlation))
  }
}

#patch_capacity
#mating_proportion
#dispersal_rate
#dispersal_model
#dispersal_cost | dispersal_model
#delet_loci
#(delet_mutation_rate, delet_effect_mean) | delet_loci
for (col in colnames(parameters)){
  df[,col] <- as.numeric(df[,col])
}
df$dispersal_model <- as.factor(df$dispersal_model)

ggplot(df) + geom_point(aes(delet_effects_mean, off.fst))
ggplot(df) + geom_point(aes(delet_loci, off.fst))
ggplot(df) + geom_density(aes(off.fis, color=dispersal_model))
ggplot



#treatment = breeding
#outcome = adlt.delhmz - mean deletirious mutation homozygosity
#covar = dispersal_rate, dispersal_model, delet_mutation_rate, or delet_mutation_model

#treatment = breeding
#outcome = adlt.delhtz - mean deletirious mutation heterozygosity
#covar = dispersal_rate, dispersal_model, delet_mutation_rate, delet_effects_mean, delet_mutation_model

#treatment = breeding
#outcome = off.ho - heterozygosity
#covar = dispersal_rate

#treatment = breeding
#outcome = adlt.ho - heterozygosity
#covar = dispersal_rate

#treatment = breeding
#outcome = adlt.allnbp - mean number of alleles within demes
#covar = disperssal_rate

###treatment = breeding
###outcome = adlt.delhmz
###covar = delet_mutation_rate

df$adlt.delhmz <- as.numeric(df$adlt.delhmz)
df$delet_mutation_rate <- as.numeric(df$delet_mutation_rate)
ggplot(df) + geom_point(aes(adlt.delhmz, delet_mutation_rate))
ggplot(df) + geom_density(aes(adlt.delhmz, color=breeding))

write.csv(df, paste0(results_dir, "/Nemo-breeding_data.csv"))
