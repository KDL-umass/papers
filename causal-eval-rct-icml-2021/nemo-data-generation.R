setwd("~/Documents/workspace/causal-eval-rct/Nemo\ data")
library(MASS)

generate_sample_breeding <- function(iter, results_dir){
  #we'll always do 10 patches, but randomize patch capacity within that
  patch_capacities <- as.integer(rnorm(16, mean=120, sd=60))
  for (i in 1:length(patch_capacities)){
    while (patch_capacities[i] < 10)
      patch_capacities[i] <- as.integer(rnorm(1, mean=120, sd=60))
  }
  
  mating_proportion <- rnorm(1, mean=0.8, sd=0.05)
  if (mating_proportion >= 1)
    mating_proportion <- 0.99
  
  dispersal_rate <- rlnorm(1, mean=-2, sd=1)
  while (dispersal_rate > 0.6)
    dispersal_rate <- rlnorm(1, mean=-2, sd=)
  
  dispersal_model <- sample(1:4, 1)
  
  dispersal_cost <- -1
  while (dispersal_cost < 0){
    if (dispersal_model == 1){
      dispersal_cost <- rnorm(1, mean= 0.1, sd=0.05)
    } else if (dispersal_model == 2){
      dispersal_cost <- rnorm(1, mean= 0.05, sd=0.01)
    } else if (dispersal_model == 3){
      dispersal_cost <- rnorm(1, mean= 0.03, sd=0.01)
    } else{
      dispersal_cost <- rnorm(1, mean= 0.005, sd=0.0006)
    }
  }
  
  delet_loci <- as.integer(rnorm(1, mean=80, sd=15))
  while (delet_loci <= 0 || delet_loci >= 128)
    delet_loci <- as.integer(rnorm(1, mean=80, sd=15))
  
  
  #mat <- matrix(c(0.005, 0, 0, 0.05), 2, 2)
  #mat.second <- matrix(c(1, .3, .3, 1), 2, 2)
  #mat.final <- (mat.second %*% mat) %*% solve(mat.second)
  
  delet_vals <- mvrnorm(1000, mu=c(0.005, 0.05), Sigma=matrix(c(1, .3, .3, 1),2,2))
  delet_mutation_rate <- ((delet_vals[1]-.005)/2000)+0.005
  delet_effects_mean <- ((delet_vals[2]-.05)/50)+0.05
  
  while (delet_mutation_rate <= 0 || delet_effects_mean <= 0){
    delet_vals <- mvrnorm(1, mu=c(0.005, 0.05), Sigma=matrix(c(1, .3, .3, 1),2,2))
    delet_mutation_rate <- ((delet_vals[1]-.005)/2000)+0.005
    delet_effects_mean <- ((delet_vals[2]-.05)/50)+0.05
  }
  
  delet_mutation_model <- sample(1:2, 1)
  disp_mutation_rate <- -1
  while (disp_mutation_rate < 0)
    disp_mutation_rate <- rnorm(1, mean=.01, sd=0.05)
  
  #now write the .ini file
  
  ini_text <- paste0("logfile                 example.log
  
  run_mode                overwrite
  
  random_seed             486532
  
  root_dir                ", results_dir, "   ## the directory in which all results and log-files will be saved
  
  replicates              1
  generations             100
                     
  filename                %'4[POLYMONO]'1_iter", iter, "
  
  ###         POPULATION          ###
  patch_capacity          {{", paste(patch_capacities, collapse=", "), "}}
  
  ###      LIFE CYCLE EVENTS      ###
  
  # the number after the LCE name is its rank in the life cycle
  
  breed_selection         1
  save_stats              2   # placed here, we will have access to the adults and the juveniles
  disperse                3   # it is imperative that dispersal is before aging!
  viability_selection     4
  aging                   5
  save_files              6
  extinction              7
  
  ## SELECTION ##
  selection_trait         delet    # fitness (offspring survival) determined by the deleterious mutations
  selection_model         direct
  
  ## EXTINCTION ##
  extinction_rate         0.05
  
  ## MATING SYSTEM (BREED) ##/
    1 = promiscuity 
    2 = polygyny 
    3 = monogyny 
    4 = selfing
    5 = cloning    
  /#
  mating_system            2 3        ### <--- sequential parameter no. 3
  mean_fecundity           4
  mating_proportion        ", mating_proportion, "
  
  ## DISPERSAL LCE ##
  ## dispersal models:
  ## 1 = Island Model Migrant Pool
  ## 2 = Island Model Propagule Pool
  ## 3 = Stepping Stone Model 1D
  ## 4 = Lattice Model
  dispersal_model          ", dispersal_model, "
  dispersal_rate           ", dispersal_rate, "
  dispersal_cost           ", dispersal_cost, "
  ", ifelse(dispersal_model == 4, "dispersal_border_model   3", ""), "
  ", ifelse(dispersal_model == 4, "dispersal_lattice_range   1", ""), "
  ", ifelse(dispersal_model == 2, "dispersal_propagule_prob  0.3", ""), "
  
  ###          TRAITS          ###
  
  ## NEUTRAL MARKERS ##
  ntrl_loci                20
  ntrl_all                 12
  ntrl_mutation_rate       0.0001
  ntrl_recombination_rate  0.5
  ntrl_mutation_model      1
  # output #
  ntrl_save_genotype       fstat
  ntrl_save_freq           allele
  ntrl_output_dir          ntrl
  ntrl_output_logtime      100
  
  ## DELETERIOUS MUTATIONS ##
  delet_loci               ", delet_loci, "
  delet_init_freq          0
  delet_mutation_rate      ", delet_mutation_rate, "
  delet_mutation_model     ", delet_mutation_model, "
  delet_recombination_rate 0.5
  delet_effects_mean       ", delet_effects_mean, "
  delet_dominance_mean     0.36
  delet_effects_distribution lognormal
  delet_effects_dist_param1 -1.4 #was -6.4
  delet_effects_dist_param2 5.3
  # output #
  delet_save_genotype
  delet_genot_dir          delet
  delet_genot_logtime      100
  
  ## DISPERSAL TRAIT ##
  disp_mutation_rate       ", disp_mutation_rate, "
  disp_mutation_mean       0.2
  
  
  ###           OUTPUT           ###
  
  ## STATS ##
  # the parameters for the save_stats LCE #
  stat fstat delet viability disp demography extrate
  stat_log_time 100
  stat_dir results
  ")
  
  write(ini_text, "./Nemo2.ini")
  
  return(c(-1, sum(patch_capacities), max(patch_capacities), mean(patch_capacities), mating_proportion, dispersal_rate, dispersal_model, dispersal_cost, delet_loci, delet_mutation_rate, delet_effects_mean, delet_mutation_model, disp_mutation_rate))
}

generate_sample_delet_model <- function(iter, results_dir){
  #we'll always do 10 patches, but randomize patch capacity within that
  patch_capacities <- as.integer(rnorm(16, mean=120, sd=60))
  for (i in 1:length(patch_capacities)){
    while (patch_capacities[i] < 10)
      patch_capacities[i] <- as.integer(rnorm(1, mean=120, sd=60))
  }
  
  mating_proportion <- rnorm(1, mean=0.8, sd=0.05)
  if (mating_proportion >= 1)
    mating_proportion <- 0.99
  
  dispersal_rate <- rlnorm(1, mean=0.2, sd=0.1)
  while (dispersal_rate > 0.6 || dispersal_rate < 0.01)
    dispersal_rate <- rnorm(1, mean=0.2, sd=0.1)
  
  dispersal_model <- sample(1:4, 1)
  
  dispersal_cost <- -1
  while (dispersal_cost < 0){
    if (dispersal_model == 1){
      dispersal_cost <- rnorm(1, mean= 0.1, sd=0.05)
    } else if (dispersal_model == 2){
      dispersal_cost <- rnorm(1, mean= 0.05, sd=0.01)
    } else if (dispersal_model == 3){
      dispersal_cost <- rnorm(1, mean= 0.03, sd=0.01)
    } else{
      dispersal_cost <- rnorm(1, mean= 0.005, sd=0.0006)
    }
  }
  
  delet_loci <- as.integer(rnorm(1, mean=80, sd=13))
  while (delet_loci <= 0 || delet_loci >= 128)
    delet_loci <- as.integer(rnorm(1, mean=80, sd=13))
  
  
  #mat <- matrix(c(0.005, 0, 0, 0.05), 2, 2)
  #mat.second <- matrix(c(1, .3, .3, 1), 2, 2)
  #mat.final <- (mat.second %*% mat) %*% solve(mat.second)
  
  delet_vals <- mvrnorm(1000, mu=c(0.005, 0.05), Sigma=matrix(c(1, .3, .3, 1),2,2))
  delet_mutation_rate <- ((delet_vals[1]-.005)/2000)+0.005
  delet_effects_mean <- ((delet_vals[2]-.05)/50)+0.05
  
  while (delet_mutation_rate <= 0 || delet_effects_mean <= 0){
    delet_vals <- mvrnorm(1, mu=c(0.005, 0.05), Sigma=matrix(c(1, .3, .3, 1),2,2))
    delet_mutation_rate <- ((delet_vals[1]-.005)/2000)+0.005
    delet_effects_mean <- ((delet_vals[2]-.05)/50)+0.05
  }
  
  delet_mutation_model <- -1
  disp_mutation_rate <- -1
  while (disp_mutation_rate < 0)
    disp_mutation_rate <- rnorm(1, mean=.01, sd=0.05)
  
  #now write the .ini file
  
  ini_text <- paste0("logfile                 example.log
  
  run_mode                overwrite
  
  random_seed             486532
  
  root_dir                ", results_dir, "   ## the directory in which all results and log-files will be saved
  
  replicates              1
  generations             100
                     
  filename                delet_model_%'1[12]'1_iter", iter, "
  
  ###         POPULATION          ###
  patch_capacity          {{", paste(patch_capacities, collapse=", "), "}}
  
  ###      LIFE CYCLE EVENTS      ###
  
  # the number after the LCE name is its rank in the life cycle
  
  breed_selection         1
  save_stats              2   # placed here, we will have access to the adults and the juveniles
  disperse                3   # it is imperative that dispersal is before aging!
  viability_selection     4
  aging                   5
  save_files              6
  extinction              7
  
  ## SELECTION ##
  selection_trait         delet    # fitness (offspring survival) determined by the deleterious mutations
  selection_model         direct
  
  ## EXTINCTION ##
  extinction_rate         0.05
  
  ## MATING SYSTEM (BREED) ##/
    1 = promiscuity 
    2 = polygyny 
    3 = monogyny 
    4 = selfing
    5 = cloning    
  /#
  mating_system            2        ### <--- sequential parameter no. 3
  mean_fecundity           4
  mating_proportion        ", mating_proportion, "
  
  ## DISPERSAL LCE ##
  ## dispersal models:
  ## 1 = Island Model Migrant Pool
  ## 2 = Island Model Propagule Pool
  ## 3 = Stepping Stone Model 1D
  ## 4 = Lattice Model
  dispersal_model          ", dispersal_model, "
  dispersal_rate           ", dispersal_rate, "
  dispersal_cost           ", dispersal_cost, "
  ", ifelse(dispersal_model == 4, "dispersal_border_model   3", ""), "
  ", ifelse(dispersal_model == 4, "dispersal_lattice_range   1", ""), "
  ", ifelse(dispersal_model == 2, "dispersal_propagule_prob  0.3", ""), "
  
  ###          TRAITS          ###
  
  ## NEUTRAL MARKERS ##
  ntrl_loci                20
  ntrl_all                 12
  ntrl_mutation_rate       0.0001
  ntrl_recombination_rate  0.5
  ntrl_mutation_model      1
  # output #
  ntrl_save_genotype       fstat
  ntrl_save_freq           allele
  ntrl_output_dir          ntrl
  ntrl_output_logtime      100
  
  ## DELETERIOUS MUTATIONS ##
  delet_loci               ", delet_loci, "
  delet_init_freq          0
  delet_mutation_rate      ", delet_mutation_rate, "
  delet_mutation_model     1 2
  delet_recombination_rate 0.5
  delet_effects_mean       ", delet_effects_mean, "
  delet_dominance_mean     0.36
  delet_effects_distribution lognormal
  delet_effects_dist_param1 -1.4 #was -6.4
  delet_effects_dist_param2 5.3
  # output #
  delet_save_genotype
  delet_genot_dir          delet
  delet_genot_logtime      100
  
  ## DISPERSAL TRAIT ##
  disp_mutation_rate       ", disp_mutation_rate, "
  disp_mutation_mean       0.2
  
  
  ###           OUTPUT           ###
  
  ## STATS ##
  # the parameters for the save_stats LCE #
  stat fstat delet viability disp demography extrate
  stat_log_time 100
  stat_dir results
  ")
  
  write(ini_text, "./Nemo2.ini")
  
  return(c(2, sum(patch_capacities), max(patch_capacities), mean(patch_capacities), mating_proportion, dispersal_rate, dispersal_model, dispersal_cost, delet_loci, delet_mutation_rate, delet_effects_mean, delet_mutation_model, disp_mutation_rate))
}


generate_sample_dispersal_rate <- function(iter, results_dir){
  #we'll always do 10 patches, but randomize patch capacity within that
  
  breeding_system <- sample(1:3, 1, prob=c(0.15, 0.55, 0.30))
  patch_capacities <- as.integer(rnorm(16, mean=120, sd=60))
  for (i in 1:length(patch_capacities)){
    while (patch_capacities[i] < 10)
      patch_capacities[i] <- as.integer(rnorm(1, mean=120, sd=60))
  }
  
  mating_proportion <- rnorm(1, mean=0.8, sd=0.05)
  if (mating_proportion >= 1)
    mating_proportion <- 0.99
  
#  dispersal_rate <- rlnorm(1, mean=-2.5, sd=1)
#  while (dispersal_rate > 0.6)
#    dispersal_rate <- rlnorm(1, mean=-2.5, sd=1)
  dispersal_rate <- -1
  
  dispersal_model <- sample(1:4, 1)
  
  dispersal_cost <- -1
  while (dispersal_cost < 0){
    if (dispersal_model == 1){
      dispersal_cost <- rnorm(1, mean= 0.1, sd=0.05)
    } else if (dispersal_model == 2){
      dispersal_cost <- rnorm(1, mean= 0.07, sd=0.03)
    } else if (dispersal_model == 3){
      dispersal_cost <- rnorm(1, mean= 0.05, sd=0.02)
    } else{
      dispersal_cost <- rnorm(1, mean= 0.005, sd=0.0006)
    }
  }
  
  delet_loci <- as.integer(rnorm(1, mean=90, sd=25))
  while (delet_loci <= 0 || delet_loci >= 140)
    delet_loci <- as.integer(rnorm(1, mean=90, sd=25))
  
  
  #mat <- matrix(c(0.005, 0, 0, 0.05), 2, 2)
  #mat.second <- matrix(c(1, .3, .3, 1), 2, 2)
  #mat.final <- (mat.second %*% mat) %*% solve(mat.second)
  
  delet_vals <- mvrnorm(1000, mu=c(0.005, 0.08), Sigma=matrix(c(1, .3, .3, 1),2,2))
  delet_mutation_rate <- ((delet_vals[1]-.005)/500)+0.005
  delet_effects_mean <- ((delet_vals[2]-.08)/30)+0.08
  
  while (delet_mutation_rate <= 0 || delet_effects_mean <= 0){
    delet_vals <- mvrnorm(1, mu=c(0.005, 0.05), Sigma=matrix(c(1, .3, .3, 1),2,2))
    delet_mutation_rate <- ((delet_vals[1]-.005)/500)+0.005
    delet_effects_mean <- ((delet_vals[2]-.08)/30)+0.08
  }
  
  delet_mutation_model <- sample(1:2, 1)
  disp_mutation_rate <- -1
  while (disp_mutation_rate < 0)
    disp_mutation_rate <- rnorm(1, mean=.02, sd=0.08)
  
  #now write the .ini file
  
  ini_text <- paste0("logfile                 example.log
  
  run_mode                overwrite
  
  random_seed             486532
  
  root_dir                ", results_dir, "   ## the directory in which all results and log-files will be saved
  
  replicates              1
  generations             100
                     
  filename                dispersal_%'4[0.050.15]'1_iter", iter, "
  
  ###         POPULATION          ###
  patch_capacity          {{", paste(patch_capacities, collapse=", "), "}}
  
  ###      LIFE CYCLE EVENTS      ###
  
  # the number after the LCE name is its rank in the life cycle
  
  breed_selection         1
  save_stats              2   # placed here, we will have access to the adults and the juveniles
  disperse                3   # it is imperative that dispersal is before aging!
  viability_selection     4
  aging                   5
  save_files              6
  extinction              7
  
  ## SELECTION ##
  selection_trait         delet    # fitness (offspring survival) determined by the deleterious mutations
  selection_model         direct
  
  ## EXTINCTION ##
  extinction_rate         0.05
  
  ## MATING SYSTEM (BREED) ##/
    1 = promiscuity 
    2 = polygyny 
    3 = monogyny 
    4 = selfing
    5 = cloning    
  /#
  mating_system            ", breeding_system, "
  mean_fecundity           4
  mating_proportion        ", mating_proportion, "
  
  ## DISPERSAL LCE ##
  ## dispersal models:
  ## 1 = Island Model Migrant Pool
  ## 2 = Island Model Propagule Pool
  ## 3 = Stepping Stone Model 1D
  ## 4 = Lattice Model
  dispersal_model          ", dispersal_model, "
  dispersal_rate           0.005 0.15
  dispersal_cost           ", dispersal_cost, "
  ", ifelse(dispersal_model == 4, "dispersal_border_model   3", ""), "
  ", ifelse(dispersal_model == 4, "dispersal_lattice_range   1", ""), "
  ", ifelse(dispersal_model == 2, "dispersal_propagule_prob  0.3", ""), "
  
  ###          TRAITS          ###
  
  ## NEUTRAL MARKERS ##
  ntrl_loci                20
  ntrl_all                 12
  ntrl_mutation_rate       0.0001
  ntrl_recombination_rate  0.5
  ntrl_mutation_model      1
  # output #
  ntrl_save_genotype       fstat
  ntrl_save_freq           allele
  ntrl_output_dir          ntrl
  ntrl_output_logtime      100
  
  ## DELETERIOUS MUTATIONS ##
  delet_loci               ", delet_loci, "
  delet_init_freq          0
  delet_mutation_rate      ", delet_mutation_rate, "
  delet_mutation_model     ", delet_mutation_model, "
  delet_recombination_rate 0.5
  delet_effects_mean       ", delet_effects_mean, "
  delet_dominance_mean     0.36
  delet_effects_distribution lognormal
  delet_effects_dist_param1 -1.4 #was -6.4
  delet_effects_dist_param2 5.3
  # output #
  delet_save_genotype
  delet_genot_dir          delet
  delet_genot_logtime      100
  
  ## DISPERSAL TRAIT ##
  disp_mutation_rate       ", disp_mutation_rate, "
  disp_mutation_mean       0.2
  
  
  ###           OUTPUT           ###
  
  ## STATS ##
  # the parameters for the save_stats LCE #
  stat fstat delet viability disp demography extrate
  stat_log_time 100
  stat_dir results
  ")
  
  write(ini_text, "./Nemo2.ini")
  
  return(c(breeding_system, sum(patch_capacities), max(patch_capacities), mean(patch_capacities), mating_proportion, dispersal_rate, dispersal_model, dispersal_cost, delet_loci, delet_mutation_rate, delet_effects_mean, delet_mutation_model, disp_mutation_rate))
}

results_dir <- "results-dispersal-rate"

num_files <- length(list.files(paste0(results_dir, "/results")))
all_params <- c()
for (iter in 9001:10000){
  print(iter)
  all_params <- rbind(all_params, c(iter, generate_sample_dispersal_rate(iter, results_dir)))
  system("nemo")
  if (length(list.files(paste0(results_dir, "/results"))) == num_files){
    print("ABORT!")
    print(iter)
    system("cp ./Nemo2.ini ./Nemo2-failed.ini")
  }
  num_files <- length(list.files(paste0(results_dir, "/results")))
}
colnames(all_params) <- c("iter", "breeding", "total_capacity", "max_patch_capacity", "mean_patch_capacity", "mating_proportion", "dispersal_rate", "dispersal_model", "dispersal_cost", "delet_loci", "delet_mutation_rate", "delet_effects_mean", "delet_mutation_model", "disp_mutation_rate")
write.csv(all_params, paste0(results_dir, "/parameters9.csv"), row.names = FALSE)

#let's set:
#patch_capacity
#mating_proportion
#dispersal_rate
#dispersal_model
#dispersal_cost | dispersal_model
#delet_loci
#(delet_mutation_rate, delet_effect_mean) | delet_loci
#disp_mutation_rate
