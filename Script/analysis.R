# Set working directory -------------------

setwd("C:\\Documents\\Bi-cluster\\")


# Load functions -------------------

source("function.R") 


# Read data and impirical values -------------------

X <- matrix(unlist(read.table("SimulationData.txt")), nrow = 100, byrow = F)   
n <- nrow(X)
p <- ncol(X)
X_expand <- cbind(c(1:n), X)
mu0 <- apply(X, 2, mean)
sigma0sq <- apply(X, 2, var)
sigma0 <- sqrt(sigma0sq)


# Initial Values -------------------

iter <- 150000   ## set the number of iterations in the MCMC
k_max <- 5   ## set the value of K
alpha_init <- 1   ## set the initial value of alpha
gamma_a <- 2   ## set the shape parameter in the prior distribution of alpha
gamma_b <- 2   ## set the scale parameter in the prior distribution of alpha
sigmasq_a_init <- 0.5  ## set the shape parameter in the prior distribution of sigmasq
sigmasq_b_init <- 0.5 * var(as.vector(X))   ## set the scale parameter in the prior distribution of sigmasq
m_ini <- n   ## set the initial value of number of subject-level clusters
C_ini <- seq(1, n)   ## set the initial value of C's
sigmasq_ini <- 1   ## set the initial value of sigmasq
rho_ini <- 0   ## set the initial value of rho
k_ini <- rep(0, n)   ## set the initial value of k's
s_ini <- matrix(rep(c(1, p+1, rep(NA, k_max)), n), nrow = n, byrow = T)   ## extended S's
mu1 <- unlist(apply(X, 1, mean))   ## set the initial value of mu's
mu_star_ini <- matrix(c(mu1, rep(NA, k_max * n)), nrow = n, byrow = F)   ## extended mu's
current <- list(M = m_ini,
                Cluster = C_ini,
                Heights = mu_star_ini,
                K = k_ini,
                S = s_ini,
                Var_X = sigmasq_ini,
                Corr = rho_ini,
                Alpha = alpha_init,
                Info = "Initial")


# Run MCMC iterations and output results -------------------

for ( rept in 1:iter ) {
  current <- Update.C(current)
  current <- Update.RJMCMC(current)
  current <- Update.Sigmasq(current)
  current <- Update.Rho(current,method = "uniform")
  current <- Update.Alpha(current)
  # output
  C.out <- paste(current$Cluster, collapse = " ")
  cat("\n", C.out, "\n", file = "C.txt", append = TRUE)
  
  m.out <- paste(current$M, collapse = " ")
  cat("\n", m.out, "\n", file = "m.txt", append = TRUE)		
  
  K.out <- paste(current$K, collapse = " ")
  cat("\n", K.out, "\n", file = "K.txt", append = TRUE)
  
  S.out <- unlist(apply(current$S, 1, paste, collapse = " "))
  S.out <- paste(S.out, collapse = "\n", sep = "")
  cat("\n", S.out, "\n", file = "S.txt", append = TRUE)
  
  mu_star.out <- unlist(apply(current$Heights, 1, paste, collapse = " "))
  mu_star.out <- paste(mu_star.out, collapse = "\n", sep = "")
  cat("\n", mu_star.out, "\n", file = "mu_star.txt", append = TRUE)
  
  sigmasq.out <- paste(current$Var_X, collapse = " ")
  cat("\n", sigmasq.out, "\n", file = "sigmasq.txt", append = TRUE)
  
  rho.out <- paste(current$Corr, collapse = " ")
  cat("\n", rho.out, "\n", file = "rho.txt", append = TRUE)
  
  alpha.out <- paste(current$Alpha, collapse = " ")
  cat("\n", alpha.out, "\n", file = "alpha.txt", append = TRUE)
}