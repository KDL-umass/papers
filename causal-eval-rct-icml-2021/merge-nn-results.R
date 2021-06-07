results <- ""
binary_datasets <- ""

#load("~/Downloads/nn_initial-30-run-bias-1-results_whynot_world_normalized.RData")
#results.world2 <- results

#load("~/Downloads/nn_initial-30-run-bias-1-binary_datsets_whynot.Rdata")
#binary_datasets.world2 <- binary_datasets

#load("all-data-results/1000-sample-size/bias-1/run-30-bias-1-nn-results.RData")
load("all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-nn-results.RData")
load("all-data-results/bias-1/run-30-run-bias-1-nn-results.RData")
results.nn <- results

load("all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-nn-binary_datasets.RData")
#load("all-data-results/bias-1/run-30-run-bias-1-nn-binary_datasets.Rdata")
binary_datasets.nn <- binary_datasets



#load("all-data-results/1000-sample-size/bias-1/run-30-bias-1-1000-samples-results.Rdata")
#load("all-data-results/1000-sample-size/bias-1/run-30-bias-1-1000-samples-binary_datasets.RData")
load("all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-two-bias-results.RData")
load("all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-two-bias-binary_datasets.RData")


for (dataset in names(results)){
  new_methods <- names(results.nn[[dataset]])[!(names(results.nn[[dataset]]) %in% c("true_ATE", "naive"))]
  for (new_method in new_methods){
    results[[dataset]][[new_method]] <- results.nn[[dataset]][[new_method]]
  }
}

#save(results, file="all-data-results/1000-sample-size/bias-1/run-30-bias-1-1000-samples-with-nn-results.RData")
#save(binary_datasets, file="all-data-results/1000-sample-size/bias-1/run-30-bias-1-1000-samples-with-nn-binary_datasets.RData")
save(results, file="all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-two-bias-with-nn-results.RData")
save(binary_datasets, file="all-data-results/two-bias-no-NAs/bias-1/run-30-bias-1-two-bias-with-nn-binary_datasets.RData")

