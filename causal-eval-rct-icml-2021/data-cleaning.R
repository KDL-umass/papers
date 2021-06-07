data <- read.csv("~/Documents/workspace/causal-eval-rct/data/Dryad_B8KG77_data.csv")


for (col in colnames(data)){
  if (is.character(data[,col]) && length(grep('\'', data[,col])) > 0 || length(grep('\"', data[,col])) > 0){#at least one ' or " in this column
    print(col)
    fixed <-sapply(data[,col], function(val){gsub('\'', "", val)})
    fixed <- sapply(fixed, function(val){gsub('\"', "", val)})
    data[,col] <- fixed
  }
}

write.csv(data, "~/Documents/workspace/causal-eval-rct/data/Dryad_B8KG77_data.csv")


data <- read.csv("~/Documents/workspace/causal-eval-rct/data/UK_Data_service_852874_data.csv")

data[data$EmploymentStatus > 2, "EmploymentStatus"] <- 2
write.csv(data, "~/Documents/workspace/causal-eval-rct/data/UK_Data_service_852874_data.csv")




data <- read.csv("~/Documents/workspace/causal-eval-rct/data/Dryad_B8KG77_data.csv")
data[data$primary_role == "Visiting scholar","primary_role"] <- "Faculty"
write.csv(data, "~/Documents/workspace/causal-eval-rct/data/Dryad_B8KG77_data.csv")
