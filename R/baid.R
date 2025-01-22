## ========================================================================== ##
##                BAYESIAN ANALYSIS OF INDIVIDUAL DIFFERENCES
## ========================================================================== ##

## ==============================
## LOAD/INSTALL REQUIRED PACKAGES

# Install/load pacman
if (!require(pacman, quietly = TRUE)) {
  
  install.packages("pacman")
  require(pacman)
  
}

# Required Packages
pkgs <- c("dplyr", "purrr", "rstan", "MCMCpack", "msm")
p_load(pkgs, character.only = T)


## ==============================
## (INVISIBLE) HELP FUNCTIONS

# Restructure list of vectors into matrix  
.make_matrix <- function (x) {
  
  max_n <- purrr::map_dbl(x, length) |> max()
  out <- purrr::map(x, \(x) {
    length(x) = max_n
    return(x)
  }) |> 
    reduce(rbind)
  dimnames(out) <- NULL
  return(out)
  
}

# Restructure data frame for efficient sampling (requires columns sub, cond, y)
.transform_data <- function (df) {
  
  if (!is.data.frame(df))
    stop("Please provide data frame with columns 'sub', 'cond', 'y'.")
  
  if (!all(c("sub", "cond", "y") %in% names(df))) 
    stop("Please provide data frame with columns 'sub', 'cond', 'y'.")
  
  list(
    N = length(unique(df$sub)),
    k = tapply(df$sub, df$sub, length) |> as.vector(),
    y = .make_matrix(split(df$y, df$sub)),
    x = .make_matrix(split(df$cond, df$sub)),
    Ntot = nrow(df)
  )
  
}

# Calculate conditional posterior rho in mixture-model with 3 latent classes
.calc_rho_3 <- function (y, x, t1, t2, alpha, J, s2, p) {
  
  out <- numeric(3)
  
  tmp1 <- exp(- (2 * t1 * sum((y - alpha) / x, na.rm = T) - J * t1^2) / (8 * s2))
  tmp2 <- exp(- (2 * (t1 - t2) * sum((y - alpha) / x, na.rm = T) + J * (t2^2 - t1^2)) / (8 * s2))
  out[1] <- 1 / (1 + p[3]/p[1] * tmp1 + p[2]/p[1] * tmp2)
  
  tmp1 <- exp(- (J * t1^2 - 2 * t1 * sum((y - alpha) / x, na.rm = T)) / (8 * s2))
  tmp2 <- exp(- (J * t2^2 - 2 * t2 * sum((y - alpha) / x, na.rm = T)) / (8 * s2))
  out[3] <- 1 / (1 + p[1]/p[3] * tmp1 + p[2]/p[3] * tmp2)
  
  out[2] <- 1 - sum(out[c(1, 3)])
  out <- ifelse(out < 0, 0, out)
  
  return (out)
}

# Calculate convergence statistics (Rhat and effective sample sizes)
.check_convergence <- function (samps, chains) {
  
  tmp <- as.matrix(samps)
  
  # Note: To calculate statistics, samples are retransformed 
  # into n_samples x n_chains matrix
  
  rhat <- apply(tmp, 2, \(x) {
    m = matrix(x, ncol = chains)
    rstan::Rhat(m)
  }, simplify = TRUE)
  
  bulk_ess <- apply(tmp, 2, \(x) {
    m = matrix(x, ncol = chains)
    rstan::ess_bulk(m)
  }, simplify = TRUE)
  
  tail_ess <- apply(tmp, 2, \(x) {
    m = matrix(x, ncol = chains)
    rstan::ess_tail(m)
  }, simplify = TRUE)
  
  # Return df with convergence statistics (rows = pars, cols = stats)
  data.frame(
    rhat = rhat,
    bulk_ess = bulk_ess,
    tail_ess = tail_ess
  )
  
}

# log of multivariate beta function
.lmbeta <- function (alpha) {
  
  sum(lgamma(alpha)) - lgamma(sum(alpha))
  
}

# Calculate Bayes factors
.calc_bfs <- function (data, priors, mcmc_samples, h1) {
  
  # Prior probability of constraint coded as numeric vector h1
  par_a <- sum(priors$q * h1)
  par_b <- priors$q[!h1]
  prior_prob <- exp(.lmbeta(c(par_a + data$N, par_b)) - .lmbeta(c(par_a, par_b)))
  
  # Posterior probability of constraint in population
  post_prob_population <- 
    exp(data$N * log(rowSums(mcmc_samples$rho[, h1 * 1:3, drop = F]))) |> mean()
  
  # Posterior probability of constraint in sample
  p_sum <- apply(mcmc_samples$p, c(1, 3), \(x) sum(x[h1 * c(1, 2, 3)]))
  post_prob_sample <- apply(p_sum, 1, \(x) exp(sum(log(x)))) |> mean()
  
  # Bayes Factor (sample)
  bf_samp <- post_prob_sample / prior_prob
  
  # Bayes Factor (population)
  bf_pop <- post_prob_population / prior_prob
  
  # return
  bf <- c(
    sample = bf_samp,
    population = bf_pop
  )
  
  return(bf)
}

# MCMC for mixture model with 3 classes (1 chain)
.mcmc_mixture_3_single <- function (data, priors, n_samples, burnin) {
  
  ### Get data
  N    <- data$N
  k    <- data$k
  y    <- data$y
  x    <- data$x
  Ntot <- data$Ntot
  
  ### Get priors
  q         <- priors$q
  am        <- priors$mu[1]
  bm        <- priors$mu[2]
  a1        <- priors$nu1[1]
  b1        <- priors$nu1[2]
  a2        <- priors$nu2[1]
  b2        <- priors$nu2[2]
  s1_eta    <- .5
  s2_eta    <- .5 * priors$eta2^2
  s1_delta1 <- .5
  s2_delta1 <- .5 * priors$delta2_1^2
  s1_delta2 <- .5
  s2_delta2 <- .5 * priors$delta2_2^2
  s1_sigma  <- .5
  s2_sigma  <- .5 * priors$sigma2^2
  
  ### Initialize output objects
  alpha  <- matrix(0, n_samples, N)
  theta1 <- matrix(0, n_samples, N)
  theta2 <- matrix(0, n_samples, N)
  rho    <- matrix(0, n_samples, 3)
  z      <- array(0, dim = c(n_samples, 3, N))
  p      <- array(0, dim = c(n_samples, 3, N))
  eta    <- numeric(n_samples)
  mu     <- numeric(n_samples)
  # pp_p   <- numeric(n_samples)
  # pp_rho <- numeric(n_samples)
  nu1    <- numeric(n_samples)
  nu2    <- numeric(n_samples)
  delta1 <- numeric(n_samples)
  delta2 <- numeric(n_samples)
  sigma  <- numeric(n_samples)
  
  ### Create random starting values
  alpha[1,]  <- runif(N, -.1, .1)
  mu[1]      <- runif(1, -.1, .1)
  eta[1]     <- runif(1, 0, 1)
  theta1[1,] <- runif(N, .1, .5)
  theta2[1,] <- runif(N, -.5, -.1)
  rho[1,]    <- rdirichlet(1, q)
  z[1,,]     <- rmultinom(N, 1, rho[1,])
  p[1,,]     <- rho[1,]
  # pp_p[1]    <- runif(1, 0, 1)
  # pp_rho[1]  <- runif(1, 0, 1)
  nu1[1]     <- runif(1, .1, .5)
  nu2[1]     <- runif(1, -.5, -.1)
  delta1[1]  <- runif(1, 0, 1)
  delta2[1]  <- runif(1, 0, 1)
  sigma[1]   <- runif(1, 0, 1)
  
  ### Start Gibbs sampler
  for (m in 2:n_samples) {
    
    # alpha
    v <- 1 / ((k / sigma[m-1]) + 1 / eta[m-1])
    c <- rowSums(y - x * (theta1[m-1,]^z[(m-1),1,])*(theta2[m-1,]^z[(m-1),2,])*(0^z[(m-1),3,]), na.rm = T) / 
      sigma[m-1] +  mu[m-1] / eta[m-1]
    alpha[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # theta1
    v <- 1 / (z[(m-1),1,] * (k / (4 * sigma[m-1])) + 1 / delta1[m-1])
    c <- z[(m-1),1,] * rowSums((y - alpha[m,]) / x, na.rm = T) / 
      (4 * sigma[m-1]) +  nu1[m-1] / delta1[m-1]
    theta1[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), lower = 0)
    
    # theta2
    v <- 1 / (z[(m-1),2,] * (k / (4 * sigma[m-1])) + 1 / delta2[m-1])
    c <- z[(m-1),2,] * rowSums((y - alpha[m,]) / x, na.rm = T) / 
      (4 * sigma[m-1]) +  nu2[m-1] / delta2[m-1]
    theta2[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), upper = 0)
    
    # z
    for (i in 1:N) {
      p[m,,i] <- .calc_rho_3(y[i,], x[i,], theta1[m,i], theta2[m,i], alpha[m,i],
                             k[i], sigma[m-1], rho[m-1,])
      z[m,,i] <- rmultinom(1, 1, p[m,,i])
    }
    
    # posterior prob does everybody (sample)
    # pp_p[m] <- exp(sum(log(p[m,1,])))
    
    # rho
    rho[m,] <- rdirichlet(1, q + rowSums(z[m,,]))
    
    # posterior prob does everybody (population)
    # pp_rho[m] <- exp(N * log(rho[m,1]))
    
    # nu1
    v <- 1 / ((N / delta1[m-1]) + 1 / b1)
    c <- sum(theta1[m,]) / delta1[m-1] +  a1 / b1
    nu1[m] <- rtnorm(1, mean = v*c, sd = sqrt(v), lower = 0)
    
    # nu2
    v <- 1 / ((N / delta2[m-1]) + 1 / b2)
    c <- sum(theta2[m,]) / delta2[m-1] +  a2 / b2
    nu2[m] <- rtnorm(1, mean = v*c, sd = sqrt(v), upper = 0)
    
    # mu
    v <- 1 / ((N / eta[m-1]) + 1 / bm)
    c <- sum(alpha[m,]) / eta[m-1] +  am / bm
    mu[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # eta
    s1 <- N / 2 + s1_eta
    tmp <- sum((alpha[m,] - mu[m])^2)
    s2 <- tmp / 2 + s2_eta
    eta[m] <- rinvgamma(1, s1, s2)
    
    # delta1
    s1 <- N / 2 + s1_delta1
    tmp <- sum((theta1[m,] - nu1[m])^2)
    s2 <- tmp / 2 + s2_delta2
    delta1[m] <- rinvgamma(1, s1, s2)
    
    # delta2
    s1 <- N / 2 + s1_delta2
    tmp <- sum((theta2[m,] - nu2[m])^2)
    s2 <- tmp / 2 + s2_delta2
    delta2[m] <- rinvgamma(1, s1, s2)
    
    # sigma2
    s1 <- Ntot / 2 + s1_sigma
    tmp <- y - alpha[m,] - x * theta1[m,]^z[m,1,] * theta2[m,]^z[m,2,] * 0^z[m,3,]
    s2 <- sum(tmp^2, na.rm = T) / 2 + s2_sigma
    sigma[m] <- rinvgamma(1, s1, s2)
    
  }
  
  ### Discard burn-in samples
  alpha <- alpha[(burnin+1):n_samples,]
  theta1 <- theta1[(burnin+1):n_samples,]
  theta2 <- theta2[(burnin+1):n_samples,]
  rho <- rho[(burnin+1):n_samples,]
  z <- z[(burnin+1):n_samples,,]
  p <- p[(burnin+1):n_samples,,]
  # pp_p <- pp_p[(burnin+1):n_samples]
  # pp_rho <- pp_rho[(burnin+1):n_samples]
  mu <- mu[(burnin+1):n_samples]
  nu1 <- nu1[(burnin+1):n_samples]
  nu2 <- nu2[(burnin+1):n_samples]
  eta <- eta[(burnin+1):n_samples]
  delta1 <- delta1[(burnin+1):n_samples]
  delta2 <- delta2[(burnin+1):n_samples]
  sigma <- sigma[(burnin+1):n_samples]
  
  ### Samples for combined theta
  theta <- z[,1,] * theta1 + z[,2,] * theta2 + z[,3,] * 0
  
  ### Return chains
  list(
    # arrays/matrices
    alpha = alpha,
    p = p,
    # pp_p = pp_p,
    # pp_rho = pp_rho,
    rho = rho,
    theta = theta,
    theta1 = theta1,
    theta2 = theta2,
    z = z,
    # vectors
    delta1 = delta1,
    delta2 = delta2,
    eta = eta,
    mu = mu,
    nu1 = nu1,
    nu2 = nu2,
    sigma = sigma
  )
  
}

# MCMC for mixture model with 3 classes (multiple chains)
.mcmc_mixture_3 <- function (data, priors, n_samples, burnin, chains, progress) {
  
  mcmc_samples_list <- purrr::map(
    .x = 1:chains,
    .f = \(x) {
      .mcmc_mixture_3_single(data, priors, n_samples, burnin)
    },
    .progress = progress
  )
  
  mcmc_samples <- do.call(dplyr::bind_rows, args = list(mcmc_samples_list))
  
  return(mcmc_samples)
  
}


## ==============================
## MAIN FUNCTIONS


# Main function for Bayesian analysis of individual differences
baid <- function (
    data,       # data frame; must contain columns sub, cond, y
    priors,     # list of priors
    n_samples,  # No of posterior samples per chain
    burnin,     # No of burnin samples to be discarded from the chains
    chains,     # No of chains 
    check_convergence = FALSE, # logical; calculate convergence statistics?
    bayes_factors = FALSE,     # logical; calculate Bayes factors?
    h1 = c(1, 0, 0),           # single vector or list of vectors representing class structures to be tested against general model
    progress = FALSE           # logical; display progress bar?
    ) {
  
  data <- .transform_data(data)
  
  if (progress)
    message("Running MCMC for latent-mixture model with ", chains, " chains...")
  
  mcmc_samples <- .mcmc_mixture_3(data, priors, n_samples, burnin, chains, progress)
  
  if (check_convergence) {
    
    if (progress)
      message("\nCalculating convergence criteria...")
    
    if (chains > 1) {
      
      checks <- mcmc_samples |> 
        dplyr::select(-c(p, z)) |> 
        .check_convergence(chains = chains)
      
    } else {
      
      warning("Chains must be > 1 to run convergence tests.")
      checks <- NULL
      
    }
    
  } else {
    
    checks <- NULL
    
  }
  
  if (bayes_factors) {
    
    if (progress)
      message("\nCalculating Bayes factors...")
    
    if (is.list(h1)) {
      
      bf <- list()
      
      for (i in h1) {
        
        bf[[toString(i)]] <- .calc_bfs(data, priors, mcmc_samples, h1 = i)
        
      }
      
      
    } else {
      
      bf <- .calc_bfs(data, priors, mcmc_samples, h1)
      
    }
    
    
  } else {
    
    bf <- NULL
    
  }
  
  return(
    list(
      mcmc_samples = mcmc_samples,
      checks = checks,
      bf = bf
    )
  )
  
}

# naive approach
naid <- function (data, ci = .90) {
  
  data <- .transform_data(data)
  
  ### Get data
  N <- data$N
  y <- data$y
  x <- data$x
  
  y_obs <- numeric(N)
  ll <- numeric(N)
  ul <- numeric(N)
  for (i in 1:N) {
    # mean
    tmp <- tapply(y[i,], x[i,], mean)
    y_obs[i] <- as.vector(tmp[2] - tmp[1])
    # vars
    tmp <- tapply(y[i,], x[i,], var)
    k <- tapply(y[i,], x[i,], length)
    se <- sqrt(tmp[1] / k[1] + tmp[2] / k[2])
    ll[i] <- y_obs[i] - qt(1 - (1 - ci) / 2, (sum(k) - 2)) * se
    ul[i] <- y_obs[i] + qt(1 - (1 - ci) / 2, (sum(k) - 2)) * se
  }
  
  class <- case_when(
    ll > 0 ~ 1,
    ul < 0 ~ -1,
    TRUE ~ 0
  )
  
  list(
    theta = y_obs,
    ll = ll,
    ul = ul,
    class = class
  )
  
}

