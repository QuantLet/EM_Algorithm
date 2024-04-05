
rm(list = ls(all = TRUE))

# change to your owner folder
wdir = "/Users/ruting/Documents/Github/FRM/FRM_Quantlet/FRM_All"

# Parameter setting
# Initialize parameters and storage for history
mu1_init <- 0
mu2_init <- 2
sigma1_init <- 1
sigma2_init <- 1
pi_init <- 0.5
iterations <- 250

n <- 250
mu1_true <- 0
mu2_true <- 3
sigma1_true <- 1
sigma2_true <- 2
pi1_true <- 0.6


# EM algorithm for a mixture of two normals with convergence plot
# Expectation step
e_step <- function(data, pi, mu1, mu2, sigma1, sigma2) {
  tau1 <- pi * dnorm(data, mean = mu1, sd = sigma1)
  tau2 <- (1 - pi) * dnorm(data, mean = mu2, sd = sigma2)
  gamma <- tau1 / (tau1 + tau2)
  return(gamma)
}

# Maximization step
m_step <- function(data, gamma) {
  mu1 <- sum(gamma * data) / sum(gamma)
  mu2 <- sum((1 - gamma) * data) / sum(1 - gamma)
  sigma1 <- sqrt(sum(gamma * (data - mu1)^2) / sum(gamma))
  sigma2 <- sqrt(sum((1 - gamma) * (data - mu2)^2) / sum(1 - gamma))
  pi <- mean(gamma)
  return(list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, pi = pi))
}

# Simulate data from a mixture of two normals
set.seed(0)
x1 <- rnorm(pi1_true * n, mu1_true, sigma1_true)
x2 <- rnorm((1 - pi1_true) * n, mu2_true, sigma2_true)
data <- sample(c(x1, x2)) # Shuffling the data

history <- matrix(NA, nrow = iterations, ncol = 5)  # To store the parameter history

# Run EM algorithm and store the parameter history
for (i in 1:iterations) {
  gamma <- e_step(data, pi_init, mu1_init, mu2_init, sigma1_init, sigma2_init)
  params <- m_step(data, gamma)
  history[i, ] <- c(params$mu1, params$mu2, params$sigma1, params$sigma2, params$pi)
  mu1_init <- params$mu1
  mu2_init <- params$mu2
  sigma1_init <- params$sigma1
  sigma2_init <- params$sigma2
  pi_init <- params$pi
}

# Assuming the 'history' matrix is already populated with parameter values from the EM algorithm

# Open a PDF device to save the plots
pdf("convergence_plots.pdf", width = 12, height = 4)

# Setup the layout for the plots
par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 1), las = 1)

# Plot the convergence of means
plot(1:iterations, history[, 1], type = 'l', col = 'blue', ylim = range(c(history[, 1:2], mu1_true, mu2_true)), xlab = 'Iteration', ylab = 'Mean value', main = 'Convergence of Means')
lines(1:iterations, history[, 2], col = 'red')
abline(h = mu1_true, col = 'blue', lty = 2)
abline(h = mu2_true, col = 'red', lty = 2)


# Plot the convergence of proportion
plot(1:iterations, history[, 5], type = 'l', col = 'orange', ylim = range(c(history[, 5], pi_true)), xlab = 'Iteration', ylab = 'Proportion', main = 'Convergence of Proportion')
abline(h = pi1_true, col = 'orange', lty = 2)

# Close the PDF device
dev.off()

