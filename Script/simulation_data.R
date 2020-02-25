# Set working directory -------------------

setwd("C:\\Documents\\Bi-cluster\\")


# Set parameters -------------------
 
p <- 24   ## set the length of the longitudinal profile p
rho <- 0.5   ## set the correlation coefficient rho
sigmasq <- 30   ## set the variance sigma^2
mu <- matrix(c(rep(60, 12), rep(80, 12),
               rep(55, 12), rep(72, 6), rep(88, 6)),
             nrow = 2, byrow = T)    ## set the subject-level cluster means mu*
n1 <- 50   ## set the number of subjects in the first cluster
n2 <- 50   ## set the number of subjects in the second cluster
### If more subject-level clusters are desired to simulate, there should be n3, n4, etc. And mu should be modified accordingly.  ###
N <- n1 + n2
Create.RhoMat <- function(r) {         ## function to create a AR(1) matrix
  Rho <- matrix(NA, nrow = p, ncol = p)
  for ( i in 1:p ) {
    for ( j in 1:p ) {
      Rho[i,j] <- r ^ abs(i - j)
    }
  }
  return(Rho)
}
Sigma_x <- Create.RhoMat(rho) * sigmasq


# Simulate data -------------------

library(MASS)

X <- matrix(0, nrow = N, ncol = p)
X[(1:n1), ] <- round(mvrnorm(n1, mu = mu[1, ], Sigma = Sigma_x), 2) 
X[((n1+1):N), ] <- round(mvrnorm(n2, mu = mu[2, ], Sigma = Sigma_x), 2)
### If more subject-level clusters are desired to simulate, X should be modified accordingly.  ###


# Output simulated data -------------------

X.out <- unlist(apply(X, 1, paste, collapse = " "))
X.out <- paste(X.out, collapse = "\n", sep = "")
cat(X.out, file = "SimulationData.txt")
