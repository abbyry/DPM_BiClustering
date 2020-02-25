
library(MASS)
library(combinat)



# preliminary functions -------------------

quad1 <- function(v1, v2) crossprod(as.vector(v1) - as.vector(v2))
quad2 <- function(v, M1, M2) matrix(v[2:length(v)] - M1[current$Cluster[v[1]], ], nrow = 1) %*%
                             solve(M2) %*%
                             matrix(v[2:length(v)] - M1[current$Cluster[v[1]], ], ncol = 1)
quad3 <- function(v1,v2,M) matrix(v1 - v2, nrow = 1) %*% solve(M) %*% matrix(v1 - v2, ncol = 1)

LRatio <- function(para_curr, para_prop, M, Sig) {
  n_j <- nrow(M)
  crossp <- rep(NA, n_j)
  for (i in 1:n_j) {
    crossp[i] <- quad3(M[i,], para_prop, Sig) - quad3(M[i,], para_curr, Sig)
  }
  out <- exp(-sum(crossp)/2) 
  return(out)
}

Create.MuForIntegral <- function() {
  K <- sample(k_max + 1, 1) - 1
  if ( K > 0 ) {
    S_all <- t(combn(c(2:p), K))
    pick <- sample(nrow(S_all), 1)
    S <- c(1, S_all[pick, ], p + 1)
  } else {
    S <- c(1, p + 1)
  }
  mu_s <- rep(NA, K + 1)
  for (pointer in 1: (K + 1)) {
    mu0_s <- mean(mu0[S[pointer]:(S[pointer + 1] - 1)])
    sig0_s <- mean(sigma0sq[S[pointer]:(S[pointer + 1] - 1)])
    mu_s[pointer] <- rnorm(1, mu0_s, sd = sqrt(sig0_s))
  }
  mu_vec <- Create.Mu(mu_s, S, K)
  return(mu_vec)
}

Create.Mu <- function(in_v, in_s, in_k) {
  mu_vec <- NULL
  for (pointer in 1:(in_k + 1)) {
    mu_vec <- c(mu_vec, rep(in_v[pointer], in_s[pointer + 1] - in_s[pointer]))
  }
  return(mu_vec)
}

Create.Prob <- function(k) {
  if (k==0) {
    eta <- 0.1
    pi <- 0
    b <- 0.9
    d <- 0
  } else {
    if (k==k_max) {
      eta <- 0.05
      pi <- 0.05
      b <- 0
      d <- 0.9
    } else {
      eta <- 0.05
      pi <- 0.05
      b <- 0.45
      d <- 0.45
    }
  }
  prob <- c(eta, pi, b, d)
  return(prob)
}

Create.RhoMat <- function(r) {
  mat <- diag(p)
  mat <- r ^ abs(row(mat) - col(mat))
  return(mat)
}

probf <- function(x, pr1, pr2) x[pr1] * exp(x[pr2])


# RJMCMC transit functions - H move
H.move <- function(curr, m, Sig) {	
  k <- curr$K
  s <- curr$S
  mu_star_cur <- curr$Heights
  j <- sample(k + 1,1)
  sig <- mean(sigma0sq[s[j]:(s[j+1] - 1)])
  mu <- mean(mu0[s[j]:(s[j+1] - 1)])
  mu_prop <- rnorm(1, mu_star_cur[j], sd = sqrt(sig))
  mu_star_pri <- replace(mu_star_cur, j, mu_prop)
  mu_vec_pri <- Create.Mu(mu_star_pri, s, k)
  mu_vec_cur <- Create.Mu(mu_star_cur, s, k)
  H_prob <- LRatio(mu_vec_cur, mu_vec_pri, m, Sig) * 
            exp(-((mu_prop + mu_star_cur[j] - 2 * mu) * (mu_prop - mu_star_cur[j])) / (2 * sig))
  u_star <- runif(1)
  if (u_star < H_prob) {
    mu_star_cur <- mu_star_pri
  }
  out_state <- list(K = k, S = s, Heights = mu_star_cur)
  return(out_state)
}

# RJMCMC move functions - P move
P.move <- function(curr, m, Sig) {
  k <- curr$K
  s_cur <- curr$S
  mu_star <- curr$Heights
  lag2 <- (s_cur - rep(c(0, 0, s_cur), length = k+2))[-(1:2)]
  Num <- sum(lag2 > 2)
  Select <- sample(Num, 1)
  j <- (1:k)[((lag2 > 2) == 1)][Select] + 1
  s_prop <- s_cur[j-1] + sample(s_cur[j+1] - s_cur[j-1] - 1, 1)
  s_pri <- replace(s_cur, j, s_prop)
  mu_vec_pri <- Create.Mu(mu_star, s_pri, k)
  mu_vec_cur <- Create.Mu(mu_star, s_cur, k)
  P_prob <- LRatio(mu_vec_cur, mu_vec_pri, m, Sig)
  u_star <- runif(1)
  if (u_star < P_prob) {
    s_cur <- s_pri
  } 
  out_state <- list(K = k, S = s_cur, Heights = mu_star)
  return(out_state)
}

# RJMCMC transit functions - B move
B.move <- function(curr, m, Sig) {
  k_cur <- curr$K
  s_cur <- curr$S
  mu_star_cur <- curr$Heights
  k_pri <- k_cur + 1
  s_star <- (1:(p+1))[-s_cur][sample(p - k_cur - 1, 1)]
  j <- (1:(k_cur+2))[(s_star < s_cur)==1][1] - 1
  s_pri <- c(s_cur[1:j], s_star, s_cur[-(1:j)])
  sig0_j_pri <- mean(sigma0sq[s_cur[j]:(s_star - 1)])
  U <- rnorm(1,0,sd = sqrt(sig0_j_pri))
  mu_j_pri <- mu_star_cur[j] + U
  mu_j1_pri <- mu_star_cur[j] - (s_star - s_cur[j]) / (s_cur[j+1] - s_star) * U
  if (j==1) {
    mu_star_pri <- c(mu_j_pri, mu_j1_pri, mu_star_cur[-1])
  } else {
    mu_star_pri <- c(mu_star_cur[1:(j-1)], mu_j_pri, mu_j1_pri, mu_star_cur[-(1:j)])	
  }
  mu_vec_pri <- Create.Mu(mu_star_pri, s_pri, k_pri)
  mu_vec_cur <- Create.Mu(mu_star_cur, s_cur, k_cur)
  LR <- LRatio(mu_vec_cur, mu_vec_pri, m, Sig)
  sig0_j1_pri <- mean(sigma0sq[s_star:(s_cur[j+1] - 1)])
  sig0_j_cur <- mean(sigma0sq[s_cur[j]:(s_cur[j+1] - 1)])
  mu0_j_pri <- mean(mu0[s_cur[j]:(s_star - 1)])
  mu0_j1_pri <- mean(mu0[s_star:(s_cur[j+1] - 1)])
  mu0_j_cur <- mean(mu0[s_cur[j]:(s_cur[j+1] - 1)])
  A <- (mu_j_pri - mu0_j_pri) ^ 2 / sig0_j_pri + (mu_j1_pri - mu0_j1_pri) ^ 2 / sig0_j1_pri - (mu_star_cur[j] - mu0_j_cur) ^ 2 / sig0_j_cur
  pri_ratio <- k_pri * exp(-A / 2) / (p - k_pri) / sqrt(2 * pi) * sqrt(sig0_j_cur / sig0_j_pri / sig0_j1_pri)
  b_k <- Create.Prob(k_cur)[3]
  d_k1 <- Create.Prob(k_pri)[4]
  pro_ratio <- d_k1 * (p - k_pri) * sqrt(pi * 2 * sig0_j_pri) / b_k / k_pri * exp(U ^ 2 / 2 / sig0_j_pri)
  Jacobian <- (s_cur[j+1] - s_cur[j]) / (s_cur[j+1] - s_star)
  B_prob <- LR * pri_ratio * pro_ratio * Jacobian
  u_star <- runif(1)
  if (u_star < B_prob) {
    s_cur <- s_pri
    mu_star_cur <- mu_star_pri
    k_cur <- k_pri
  } 
  out_state <- list(K = k_cur, S = s_cur, Heights = mu_star_cur)
  return(out_state)
}

# RJMCMC transit functions - D move
D.move <- function(curr, m, Sig) {
  k_cur <- curr$K
  s_cur <- curr$S
  mu_star_cur <- curr$Heights
  k_pri <- k_cur - 1
  SJ1 <- s_cur[2:(k_cur + 1)][sample(k_cur, 1)]
  j1 <- (1:(k_cur + 2))[(s_cur == SJ1)]
  s_pri <- s_cur[-j1]
  mu_j_pri <- (mu_star_cur[j1-1] * (SJ1 - s_cur[j1-1]) + mu_star_cur[j1] * (s_cur[j1+1] - SJ1)) /
            (s_cur[j1+1] - s_cur[j1-1])
  if (j1==2) {
    mu_star_pri <- c(mu_j_pri,mu_star_cur[-(1:2)])
  } else {
    mu_star_pri <- c(mu_star_cur[1:(j1 - 2)], mu_j_pri, mu_star_cur[-(1:j1)])
  }
  mu_vec_pri <- Create.Mu(mu_star_pri, s_pri, k_pri)
  mu_vec_cur <- Create.Mu(mu_star_cur, s_cur, k_cur)
  LR <- LRatio(mu_vec_cur, mu_vec_pri, m, Sig)
  sig0_j_pri <- mean(sigma0sq[s_cur[j1-1]:(s_cur[j1+1] - 1)])
  sig0_j_cur <- mean(sigma0sq[s_cur[j1-1]:(SJ1 - 1)]) 
  sig0_j1_cur <- mean(sigma0sq[SJ1:(s_cur[j1+1] - 1)])
  mu0_j_pri <- mean(mu0[s_cur[j1-1]:(s_cur[j1+1] - 1)])
  mu0_j_cur <- mean(mu0[s_cur[j1-1]:(SJ1 - 1)]) 
  mu0_j1_cur <- mean(mu0[SJ1:(s_cur[j1+1] - 1)])
  B <- (mu_star_cur[j1-1] - mu0_j_cur) ^ 2 / sig0_j_cur + (mu_star_cur[j1] - mu0_j1_cur) ^ 2 / sig0_j1_cur - 
       (mu_j_pri - mu0_j_pri) ^ 2 / sig0_j_pri
  pri_ratio <- (p - k_cur) * sqrt(2 * pi) * sqrt(sig0_j_cur * sig0_j1_cur / sig0_j_pri) / k_cur * exp(B / 2)
  b_k <- Create.Prob(k_pri)[3]
  d_k1 <- Create.Prob(k_cur)[4]
  U <- mu_star_cur[j1-1] - mu_j_pri
  pro_ratio <- sqrt(1 / pi / 2) * b_k * k_cur / d_k1 / (p - k_cur) / sqrt(sig0_j_cur) * exp(-U ^ 2 / 2 / sig0_j_cur)
  Jacobian <- (s_cur[j1+1] - SJ1) / (s_cur[j1+1] - s_cur[j1-1])
  D_prob <- LR * pri_ratio * pro_ratio * Jacobian
  u_star <- runif(1)
  if (u_star < D_prob) {
    s_cur <- s_pri
    mu_star_cur <- mu_star_pri
    k_cur <- k_pri
  } 
  out_state <- list(K = k_cur, S = s_cur, Heights = mu_star_cur)
  return(out_state)
}



# Gibbs sampling update functions
Update.C <- function(curr) {
  C <- curr$Cluster
  m <- curr$M
  mu_star <- curr$Heights
  k <- curr$K
  s <- curr$S
  sigmasq <- curr$Var_X
  rho <- curr$Corr
  alpha <- curr$Alpha
  Sigma <- sigmasq * Create.RhoMat(rho)
  p0 <- c(1, 1)
  p0[1] <- alpha
  for (i in 1:n){
    mu <- matrix(NA, nrow = m, ncol = p)
    for (i_mu in 1:m) {
      mu[i_mu, ] <- Create.Mu(mu_star[i_mu, (1:(k[i_mu] + 1))], s[i_mu, (1:(k[i_mu] + 2))], k[i_mu])
    }
    C0 <- C[-i]
    p0[2] <- -quad3(X[i, ], Create.MuForIntegral(), Sigma) / 2
    if (sum(C0==C[i]) == 0) {
      prob <- matrix(nrow = m, ncol = 2, 0)
      prob[C[i], ] <- p0
    } else {
      prob <- matrix(nrow = m + 1, ncol = 2, 0)
      prob[m+1,] <- p0
      prob[C[i],1] <- sum(C0 == C[i])
      prob[C[i],2] <- -quad3(X[i, ],matrix(mu, ncol = p, byrow = F)[C[i], ], Sigma) / 2
    }
    for (j in c(1:m)[-C[i]]) {
      prob[j,1] <- sum(C0==j)
      prob[j,2] <- -quad3(X[i, ], matrix(mu, ncol = p, byrow = F)[j, ], Sigma) / 2	
    }
    maxi <- max(prob[ ,2])
    prob[,2] <- prob[ ,2] - maxi
    dimnames(prob)[[2]] <- letters[1:2]
    p1 <- apply(prob, 1, probf, pr1 = "a", pr2 = "b")
    u_c <- runif(1, 0, 1)
    p2 <- p1 / sum(p1)
    j <- (1:length(p2))[((u_c <= cumsum(p2)) == 1)][1]
    if (j<=m & sum(C0==j)==0) {
      mu_star[j, ] <- c(mean(X[i, ]), rep(NA, k_max))
      k[j] <- 0
      s[j,] <- c(1, p+1, rep(NA, k_max))
    } else if (j>m) {
      mu_star <- rbind(mu_star, c(mean(X[i, ]), rep(NA, k_max)))
      k <- c(k, 0)
      s <- rbind(s, c(1, p + 1, rep(NA, k_max)))			
    }
    C[i] <- j
    mu_star <- matrix(mu_star[sort(unique(C)), ], ncol = k_max + 1)
    k <- k[sort(unique(C))]
    s <- matrix(s[sort(unique(C)), ], ncol = k_max + 2)
    C <- as.integer(as.factor(C))
    m <- max(C)
  }
  out <- list(M = m, Cluster = C, Heights = mu_star, K = k, S = s, Var_X = sigmasq, Corr = rho, Alpha = alpha, Info = "Update C")
  return(out)
}

Update.RJMCMC <- function(curr) {
  C <- curr$Cluster
  m <- curr$M
  mu_star <- curr$Heights
  k <- curr$K
  s <- curr$S
  sigmasq <- curr$Var_X
  rho <- curr$Corr
  alpha <- curr$Alpha
  Sigma <- sigmasq * Create.RhoMat(rho)
  for (i in 1:m) {
    curr_temp <- list(K = k[i], S = s[i,(1:(k[i]+2))], Heights = mu_star[i,(1:(k[i]+1))])
    data <- matrix(X[(C==i), ], ncol = p)
    # type of move
    move_prob <- Create.Prob(curr_temp$K)
    u_move <- runif(1)
    move_type <- (1:4)[((u_move < cumsum(move_prob))==1)][1]   # 1=H, 2=P, 3=B, 4=D
    # update theta(k)
    if (move_type==1) curr_temp <- H.move(curr_temp, data, Sigma)
    if (move_type==2) curr_temp <- P.move(curr_temp, data, Sigma)
    if (move_type==3) curr_temp <- B.move(curr_temp, data, Sigma)
    if (move_type==4) curr_temp <- D.move(curr_temp, data, Sigma)
    k[i] <- curr_temp$K
    s[i,] <- c(curr_temp$S, rep(NA, k_max - k[i]))
    mu_star[i,] <- c(curr_temp$Heights, rep(NA, k_max - k[i]))		
  }
  s <- matrix(s, ncol = 2 + k_max)
  mu_star <- matrix(mu_star, ncol = 1 + k_max)
  out <- list(M = m, Cluster = C, Heights = mu_star, K = k, S = s, Var_X = sigmasq, Corr = rho, Alpha = alpha, Info = "Update RJMCMC")
  return(out)
}

Update.Sigmasq <- function(curr) {
  C <- curr$Cluster
  m <- curr$M
  mu_star <- curr$Heights
  k <- curr$K
  s <- curr$S
  sigmasq <- curr$Var_X
  rho <- curr$Corr
  alpha <- curr$Alpha
  Rho <- Create.RhoMat(rho)
  mu <- matrix(NA, nrow = m, ncol = p)
  for (i in 1:m) {
    mu[i,] <- Create.Mu(mu_star[i,(1:(k[i]+1))], s[i,(1:(k[i]+2))], k[i])
  }
  sigmasq_a <- n * p / 2 + sigmasq_a_init
  sigmasq_b <- sigmasq_b_init + sum(apply(X_expand, 1, quad2, mu, Rho)) / 2
  sigmasq <- 1 / rgamma(n = 1, shape = sigmasq_a, scale = 1 / sigmasq_b)
  out <- list(M = m, Cluster = C, Heights = mu_star, K = k, S = s, Var_X = sigmasq, Corr = rho, Alpha = alpha, Info = "Update sigmasq")
  return(out)
}

Update.Rho <- function(curr,method = c("uniform", "normal")) {
  C <- curr$Cluster
  m <- curr$M
  mu_star <- curr$Heights
  k <- curr$K
  s <- curr$S
  sigmasq <- curr$Var_X
  rho <- curr$Corr
  alpha <- curr$Alpha
  Sigma <- sigmasq * Create.RhoMat(rho)
  mu <- matrix(NA, nrow = m, ncol = p)
  for (i in 1:m) {
    mu[i,] <- Create.Mu(mu_star[i,(1:(k[i]+1))], s[i,(1:(k[i]+2))], k[i])
  }
  if (method == "uniform") {
   rho_pri <- runif(1, -1, 1)	
   Sigma_pri <- sigmasq * Create.RhoMat(rho_pri)	
   alpha_0 <- exp((1 - n * (p - 1) / 2) * log((1 - rho_pri ^ 2) / (1 - rho ^ 2)) - 
                    sum(apply(X_expand, 1, quad2, mu, Sigma_pri)) / 2 + sum(apply(X_expand, 1, quad2, mu, Sigma)) / 2)
  } else {
   rho_pri <- 2
   mu_pro <- 0.26
   sigma_pro <- 0.05
   while (rho_pri > 1 || rho_pri < -1 ) rho_pri <- rnorm(1, mu_pro, sd = sigma_pro)
   Sigma_pri <- sigmasq * Create.RhoMat(rho_pri)
   alpha_0 <- exp((1 - n * (p - 1) / 2) * log((1 - rho_pri ^ 2) / (1 - rho ^ 2)) - 
                    sum(apply(X_expand, 1, quad2, mu, Sigma_pri)) / 2 + sum(apply(X_expand, 1, quad2, mu, Sigma)) / 2 + 
                    (rho_pri - rho) * (rho_pri + rho - 2 * mu_pro) / 2 / sigma_pro ^ 2)
  }
  u <- runif(1)
  if ( u < alpha_0 ) {
    rho <- rho_pri
  } 
  out <- list(M = m, Cluster = C, Heights = mu_star, K = k, S = s, Var_X = sigmasq, Corr = rho, Alpha = alpha, Info = "Update rho")
  return(out)
}

Update.Alpha <- function(curr) {
  C <- curr$Cluster
  m <- curr$M
  mu_star <- curr$Heights
  k <- curr$K
  s <- curr$S
  sigmasq <- curr$Var_X
  rho <- curr$Corr
  alpha <- curr$Alpha
  x_pi <- rbeta(1, alpha + 1, n)
  pi_x <- (gamma_a - 1 + m) / (gamma_a - 1 + m + n * (1 / gamma_b + log(x_pi)))
  u_pi <- runif(1)
  if (u_pi < pi_x) {
    alpha <- rgamma(1, shape = gamma_a + m, rate = 1 / gamma_b - log(x_pi))
  } else {
    alpha <- rgamma(1, shape = gamma_a + m - 1, rate = 1 / gamma_b - log(x_pi))
  }
  out <- list(M = m, Cluster = C, Heights = mu_star, K = k, S = s, Var_X = sigmasq, Corr = rho, Alpha = alpha, Info = "Update alpha")
  return(out)
}