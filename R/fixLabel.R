# Fixing label switching problem

library(label.switching)
library(MASS)


set.seed(123)
n <- 100  # Number of data points
k <- 3    # Number of clusters
true_labels <- sample(1:k, n, replace = TRUE)  # True cluster assignments

data <- NULL
for (i in 1:k) {
  cluster_mean <- runif(2, min = -5, max = 5)  # Generate random cluster mean
  cluster_sd <- runif(1, min = 0.5, max = 2)   # Generate random cluster standard deviation
  data <- rbind(data, mvrnorm(n = sum(true_labels == i), mu = cluster_mean, Sigma = diag(2)))
}


m = 2000 # of draws
K = 2 # of classes
J = 5 # of component-wise parameters

# initialize mcmc arrays
mcmc <- array(data = NA, dim = c(m = m, K = K, J = J))

# assign posterior draws to the array
mcmc[, , 1] <- post_par$alpha
for (i in 1:(J - 1)) {
  mcmc[, , i + 1] <- post_par$p[, , i]
}

# set of selected relabeling algorithm
set <-
  c("PRA",
    "ECR",
    "ECR-ITERATIVE-1",
    "AIC",
    "ECR-ITERATIVE-2",
    "STEPHENS",
    "DATA-BASED")

# find the MAP draw as a pivot

ls_lcm <-
  label.switching(
    method = set,
    zpivot = post_class[mapindex,],
    z = post_class,
    K = K,
    prapivot = mcmc[mapindex, ,],
    constraint = 1,
    mcmc = mcmc,
    p = post_class_p,
    data = lca_data$y
  )

mcmc_permuted <- permute.mcmc(mcmc, ls_lcm$permutations$ECR)

# change dimension for each parameter defined as in the Stan code
mcmc_permuted <-
  array(
    data = mcmc_permuted$output,
    dim = c(2000, 1, 10),
    dimnames = list(NULL, NULL, pars)
  )


fit_permuted <-
  monitor(mcmc_permuted, warmup = 0,  digits_summary = 3)


df <- data.frame(c('e[1]', 'e[2]', 'e[3]', 'env[1]', 'env[2]','env[4]'),
                 Q5 = c(-0.323,0.351, -1.293, 1.000, 1.000, 3.000),
                 Q50 = c(0.421, 1.129, -0.549, 1.000,1.000,3.000),
                 Q95 = c( 0.674, 1.395,-0.271, 1.000,1.000,3.000),
                 Mean = c(0.427, 1.135, -0.533, 1.000,1.000,3.000),
                 SD = c(0.396, 0.407, 0.401, 0.000, 0.000, 0.000),
                 Rhat = c(1.002, 1.003, 1.002, 1.000, 1.000,1.000))

colnames(df)[1] <- ' '
df


# Estimated cluster assignments
estimated_labels <- results$estimated.labels
cat("Estimated Cluster Assignments:\n")
print(estimated_labels)

# Mode of the estimated permutation matrix
mode_permutation <- results$mode.permutation
cat("Mode of Permutation Matrix:\n")
print(mode_permutation)

# Posterior probabilities of true labels
posterior_probs <- results$posterior.probabilities
cat("Posterior Probabilities of True Labels:\n")
print(posterior_probs)
