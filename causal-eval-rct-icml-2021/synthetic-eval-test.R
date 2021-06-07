library(bnlearn)
library(pcalg)
library(gtools)
library(igraph)

#generate data from specified model
num.levels <- 3
dependence.strength <- 2
num.subjects <- 10000
dag <- as.graphNEL(model2network("[T][C1][C2][C3][O|T:C1:C2:C3]"))
graph <- graph.adjacency(as(dag, "matrix"))
nodes <- V(graph)$name
parent.map <- list()
for(node in nodes){
  parent.map[[node]] <- neighbors(graph, node, mode="in")$name
}


chickering_generate <- function(node, parents, setup) {
  node.cpts <- setup$indexed.cpts[[node]]
  if(ncol(parents) > 0) {
    parent.key <- do.call(paste, c(parents, sep="."))
    split.parents <- split(parents, parent.key)
    
    allvals <- data.frame()
    for(key in names(split.parents)) {
      cur.entries <- node.cpts[[key]]
      vals <- sample(cur.entries[, node], nrow(split.parents[[key]]), prob=cur.entries$prob, replace=TRUE)
      cur.frame <- data.frame(val=vals)
      rownames(cur.frame) <- rownames(split.parents[[key]])
      allvals <- rbind(allvals, cur.frame)
    }
    return(allvals[as.character(1:nrow(parents)), ])
  } else {
    return(sample(node.cpts[, node], nrow(parents), prob = node.cpts$prob, replace=TRUE))
  }##NOTE!!!! When I grabbed this code from Dan's codebase, the prob = node.cpts$prob was missing, so every node w/no parents ignored their cpt and sampled uniformlyinstal from their domain
}

nodes <- V(graph)$name
var.levels <- list()
for(node in nodes) {
  if(grepl("^T", node) || grepl("^O", node)) {
    var.levels[[node]] <- c(0, 1)
  } else {
    var.levels[[node]] <- 1:num.levels
  }
}

cpts <- list()
for(node in nodes) {
  parents <- neighbors(graph, node, mode="in")$name
  node.settings <- var.levels[c(node, parents)]
  cpt.frame <- expand.grid(node.settings)
  setting.number <- 0
  cpt.frame <- ddply(cpt.frame, parents, function(df) {
    setting.number <<- setting.number + 1
    nlevels <- length(var.levels[[node]])
    base.weights <- (1 / (1:nlevels))
    base.weights <- base.weights / sum(base.weights)
    shifted.weights <- base.weights[((0:(nlevels-1) - setting.number) %% nlevels) + 1]
    df$prob <- rdirichlet(1, (1/dependence.strength) * 10 * shifted.weights)[1, ]
    return(df)
  })
  cpts[[node]] <- cpt.frame
}

indexed.cpts <- list()
for(node in nodes) {
  parents <- neighbors(graph, node, mode="in")$name
  if(length(parents) > 0) {
    indexed.distrs <- dlply(cpts[[node]], parents, function(df) {
      return(df)
    })
    indexed.cpts[[node]] <- indexed.distrs
  } else {
    indexed.cpts[[node]] <- cpts[[node]]
  }
}

topo.vertices <- topological.sort(graph)$name
generated <- data.frame(subject.id = 1:num.subjects)
for(node in topo.vertices) {
  parents <- parent.map[[node]]
  generated[, node] <- chickering_generate(node, generated[, parents, drop=FALSE], list(indexed.cpts=indexed.cpts))
}
generated <- generated[,!(colnames(generated) == "subject.id")]

generated$C1 <- as.factor(generated$C1)
generated$C2 <- as.factor(generated$C2)
generated$C3 <- as.factor(generated$C3)
generated$T <- as.factor(generated$T)
generated$O <- as.factor(generated$O)


train = generated[1:5000,]
test = generated[5001:10000,]

#bias data
bias.strength <- 2

#Bias using C1
bias_vals <- with(train, model.matrix(~ C1 + 0))
coefficients <- bias.strength * c(-1,1,1)
swept <- sweep(bias_vals, MARGIN=2, coefficients, FUN="*")
p <- 1 / (1 + exp(-rowSums(swept)))
desired_t <- rbinom(n=nrow(bias_vals), size=1, p)

biased_train <- c()
complement_train <- c()
for (i in 1:nrow(train)){
  if (train[i,"T"] == desired_t[i])
    biased_train <- rbind(biased_train, train[i,])
  else
    complement_train <- rbind(complement_train, train[i,])
}
dim(biased_train)
dim(complement_train)

#learn model using data
mmhc.result <- mmhc(biased_train)
mmhc.mec <- amat(mmhc.result)
mmhc.fit <- bn.fit(mmhc.result, biased_train)
grain.model <- as.grain(mmhc.fit)

num_correct_train <- 0
for (i in 1:nrow(train)){
  setEvidence(grain.model, "T", train[i,"T"])
  outcome.probs <- querygrain(grain.model, nodes = "O")$O
  pred.outcome <- names(outcome.probs)[which.max(outcome.probs)]
  print(paste0(pred.outcome + " vs " train[i,"O"]))
  if (pred.outcome == train[i,"O"])
    num_correct_train = num_correct_train + 1
  retractEvidence(grain.model)
}
num_correct_test <- 0
for (i in 1:nrow(test)){
  setEvidence(grain.model, "T", test[i,"T"])
  outcome.probs <- querygrain(grain.model, nodes = "O")$O
  pred.outcome <- names(outcome.probs)[which.max(outcome.probs)]
  if (pred.outcome == test[i,"O"])
    num_correct_test = num_correct_test + 1
  retractEvidence(grain.model)
}

num_correct_train
num_correct_test

#iamb.result <- iamb(biased_train)
#iamb.mec <- amat(iamb.result)
#iamb.fit <- bn.fit(iamb.result, biased_train)

#use model to predict outcome on training data
#intervene (remove edges into T), and get predicted value of treatment for pre-biased data
#mmhc.mec[,'T'] <- 0
#amat(mmhc) <- mmhc.mec



#use model to predict outcome on test data
