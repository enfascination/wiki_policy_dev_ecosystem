# initialize
#results <- list()
# load
results <- readRDS("generated_data/wikiprod.Rds")

# dummy data
results$test <- 3.1415926
results$mean1 <- 3.5
results$mean2 <- 35

# section Empirial Setting
results$nWikisiFoundation = 342 # number of language editions in the foundation
results$nWikisIncubating = 15
results$nWikisAnalysis = 245 #number we analyze
results$yearData = 2023
results$nAdmins = 3400
results$nEditors = 287000
results$nUsersMillions = 120
results$nSharedPolicies = 60
results$nSharedPolicyAdoptions = 3611

# section Study One: Institutional diversity or homogeneity
results$nFullyConnectedCoreCommunities = 30
results$nFullyConnectedCorePolicies = 37
results$nKCore = 34
results$nPeripheralPoliciesUnderDegree10 = 4
results$pctSequencesLength1 = 25

# section Study 2: One path or many?
# sequence statistics Will i be reporting Chi^2 % p-value or just bayes factor estimates (k=1 to k=2, k=2 to k=3, k=3 to k=4)
results$pctSequences = 25

# section Appendix


saveRDS(results, file="generated_data/wikiprod.Rds")

## test
#results2 <- readRDS("generated_data/wikiprod.Rds")
#results2 
