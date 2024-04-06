
rm(list = ls(all = TRUE))
library(ggplot2)

# change to your owner folder
wdir = "/Users/ruting/Documents/Github/EM_Algorithm/"

# Parameter setting
# Initialize parameters and storage for history
mu1_init_0 <- 1
mu2_init_0 <- 2
pi_init_0 <- 0.5
sigma <- 1
iterations <- 100

mu1_true <- 0
mu2_true <- 3
pi1_true <- 0.6


# EM algorithm for a mixture of two normals with convergence plot
# Expectation step
e_step <- function(data, pi, mu1, mu2, sigma) {
  tau1 <- pi * dnorm(data, mean = mu1, sd = sigma)
  tau2 <- (1 - pi) * dnorm(data, mean = mu2, sd = sigma)
  gamma <- tau1 / (tau1 + tau2)
  return(gamma)
}

# Maximization step
m_step <- function(data, gamma) {
  mu1 <- sum(gamma * data) / sum(gamma)
  mu2 <- sum((1 - gamma) * data) / sum(1 - gamma)
  pi <- mean(gamma)
  return(list(mu1 = mu1, mu2 = mu2, pi = pi))
}

# Simulate data from a mixture of two normals
set.seed(0)
n = 250
x1 <- rnorm(pi1_true * n, mu1_true, sigma)
x2 <- rnorm((1 - pi1_true) * n, mu2_true, sigma)
data <- sample(c(x1, x2)) # Shuffling the data

history <- matrix(NA, nrow = iterations, ncol = 3)  # To store the parameter history

# Run EM algorithm and store the parameter history
mu1_init = mu1_init_0
mu2_init = mu2_init_0
pi_init = pi_init_0
for (i in 1:iterations) {
  gamma <- e_step(data, pi_init, mu1_init, mu2_init, sigma)
  params <- m_step(data, gamma)
  history[i, ] <- c(params$mu1, params$mu2, params$pi)
  mu1_init <- params$mu1
  mu2_init <- params$mu2
  pi_init <- params$pi
}

# Assuming the 'history' matrix is already populated with parameter values from the EM algorithm

# Open a PDF device to save the plots

png(paste0(wdir,paste0("Convergence_N_",iterations,"_Mu1_",mu1_init_0,"Mu2_",mu2_init_0,"Pro_",pi_init_0,".png")), 
    width = 1500, height = 600, bg = "transparent")
# Setup the layout for the plots
par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 1), las = 1)

# Plot the convergence of means
plot(1:iterations, history[, 1], type = 'l', col = 'blue', ylim = range(c(history[, 1:2], mu1_true, mu2_true)), 
     xlab = 'Iteration', ylab = '', main = 'Convergence of Means',
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.4, lwd = 2)
lines(1:iterations, history[, 2], col = 'red',lwd = 2)
abline(h = mu1_true, col = 'blue', lty = 2)
abline(h = mu2_true, col = 'red', lty = 2)


# Plot the convergence of proportion
plot(1:iterations, history[, 3], type = 'l', col = 'orange', ylim = range(c(history[, 3], pi1_true)),
     xlab = 'Iteration', ylab = '', main = 'Convergence of Proportion',
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.4, lwd = 2)
abline(h = pi1_true, col = 'orange', lty = 2)

# Close the device
dev.off()

# Plot density
# Simulate data with the final estimated parameters
set.seed(1)  # For reproducibility
final_mu1 <- history[iterations, 1]
final_mu2 <- history[iterations, 2]
final_pi <- history[iterations, 3]

x1_final <- rnorm(final_pi * n, final_mu1, sigma)
x2_final <- rnorm((1 - final_pi) * n, final_mu2, sigma)
data_final <- sample(c(x1_final, x2_final))

# Combine the actual and simulated data for plotting
data_to_plot <- data.frame(value = c(data, data_final),
                           type = factor(c(rep("Actual", length(data)), rep("Simulated", length(data_final)))))


# Plot the densities with transparent background and no gridlines, retaining axes
ggplot(data_to_plot, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), # makes panel background transparent
        plot.background = element_rect(fill = "transparent", colour = NA), # makes plot background transparent
        panel.grid.major = element_blank(), # removes major gridlines
        panel.grid.minor = element_blank(), # removes minor gridlines
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18), # Adjust global text size
        axis.title = element_text(size = 18), # Adjust axis titles size
        plot.title = element_text(size = 20, hjust = 1.5)) # Adjust plot title size  
        + scale_fill_manual(values = c("Actual" = "blue", "Simulated" = "red"))

# When saving the plot with transparent background
ggsave(paste0(wdir,"density_plot_EM.png"), bg = "transparent",width = 15, height = 8)

# Plot original data density
ggplot(data.frame(value = data, type = 'Actural'), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), # makes panel background transparent
        plot.background = element_rect(fill = "transparent", colour = NA), # makes plot background transparent
        panel.grid.major = element_blank(), # removes major gridlines
        panel.grid.minor = element_blank(), # removes minor gridlines
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18), # Adjust global text size
        axis.title = element_text(size = 18), # Adjust axis titles size
        plot.title = element_text(size = 20, hjust = 1.5) # Adjust plot title size
        )

# When saving the plot with transparent background
ggsave(paste0(wdir,"density_plot_Original.png"), bg = "transparent",width = 15, height = 8)

