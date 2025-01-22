## ========================================================================= ##
##  HELP FUNCTIONS
## ========================================================================= ##

# MCMC sampling Mixture with 3 components (1 chain)
.mcmc_mixture_3_single = function (data, priors, n_samples, burnin) {
  
  ### Get data
  N <- data$N
  k <- data$k
  y <- data$y
  x <- data$x
  Ntot <- data$Ntot
  
  ### Get priors
  q     <- priors$q
  am    <- priors$mu[1]
  bm    <- priors$mu[2]
  a1    <- priors$nu1[1]
  b1    <- priors$nu1[2]
  a2    <- priors$nu2[1]
  b2    <- priors$nu2[2]
  s1_eta <- .5
  s2_eta <- .5 * priors$eta2^2
  s1_delta1 <- .5
  s2_delta1 <- .5 * priors$delta2_1^2
  s1_delta2 <- .5
  s2_delta2 <- .5 * priors$delta2_2^2
  s1_sigma <- .5
  s2_sigma <- .5 * priors$sigma2^2
  
  ### Initialize output objects
  alpha <- matrix(0, n_samples, N)
  theta1 <- matrix(0, n_samples, N)
  theta2 <- matrix(0, n_samples, N)
  rho <- matrix(0, n_samples, 3)
  z <- array(0, dim = c(n_samples, 3, N))
  p <- array(0, dim = c(n_samples, 3, N))
  
  eta <- numeric(n_samples)
  mu <- numeric(n_samples)
  pp_p <- numeric(n_samples)
  pp_rho <- numeric(n_samples)
  nu1 <- numeric(n_samples)
  nu2 <- numeric(n_samples)
  delta1 <- numeric(n_samples)
  delta2 <- numeric(n_samples)
  sigma <- numeric(n_samples)
  
  ### Starting values
  alpha[1,] <- runif(N, -.1, .1)
  mu[1] <- runif(1, -.1, .1)
  eta[1] <- runif(1, 0, 1)
  theta1[1,] <- runif(N, .1, .5)
  theta2[1,] <- runif(N, -.5, -.1)
  rho[1,] <- rdirichlet(1, q)
  z[1,,] <- rmultinom(N, 1, rho[1,])
  p[1,,] <- rho[1,]
  pp_p[1] <- runif(1, 0, 1)
  pp_rho[1] <- runif(1, 0, 1)
  nu1[1] <- runif(1, .1, .5)
  nu2[1] <- runif(1, -.5, -.1)
  delta1[1] <- runif(1, 0, 1)
  delta2[1] <- runif(1, 0, 1)
  sigma[1] <- runif(1, 0, 1)
  
  ### Sampling
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
      p[m,,i] <- .calc_rho_3(y[i,], x[i,], theta1[m,i], theta2[m,i], alpha[m,i], k[i], sigma[m-1], rho[m-1,])
      z[m,,i] <- rmultinom(1, 1, p[m,,i])
    }
    
    # posterior prob does everybody (sample)
    pp_p[m] <- exp(sum(log(p[m,1,])))
    
    # rho
    rho[m,] <- rdirichlet(1, q + rowSums(z[m,,]))
    
    # posterior prob does everybody (population)
    pp_rho[m] <- exp(N * log(rho[m,1]))
    
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
  pp_p <- pp_p[(burnin+1):n_samples]
  pp_rho <- pp_rho[(burnin+1):n_samples]
  mu <- mu[(burnin+1):n_samples]
  nu1 <- nu1[(burnin+1):n_samples]
  nu2 <- nu2[(burnin+1):n_samples]
  eta <- eta[(burnin+1):n_samples]
  delta1 <- delta1[(burnin+1):n_samples]
  delta2 <- delta2[(burnin+1):n_samples]
  sigma <- sigma[(burnin+1):n_samples]
  
  #Samples for combined theta
  theta <- z[,1,] * theta1 + z[,2,] * theta2 + z[,3,] * 0
  
  ### Return chains
  list(
    # arrays/matrices
    alpha = alpha,
    p = p,
    pp_p = pp_p,
    pp_rho = pp_rho,
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

# Mixture with multiple chains
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

# MCMC sampling Mixture with 2 components (1 chain)
.mcmc_mixture_2_single <- function (data, priors, n_samples, burnin) {
  
  ### Get data
  N <- data$N
  k <- data$k
  y <- data$y
  x <- data$x
  Ntot <- data$Ntot
  
  ### Get priors
  l <- priors$rho[1]
  u <- priors$rho[2]
  am <- priors$mu[1]
  bm <- priors$mu[2]
  a1 <- priors$nu1[1]
  b1 <- priors$nu1[2]
  a2 <- priors$nu2[1]
  b2 <- priors$nu2[2]
  
  s1_eta <- .5
  s2_eta <- .5 * priors$eta2^2
  s1_delta1 <- .5
  s2_delta1 <- .5 * priors$delta2_1^2
  s1_delta2 <- .5
  s2_delta2 <- .5 * priors$delta2_2^2
  s1_sigma <- .5
  s2_sigma <- .5 * priors$sigma2^2
  
  ### Initialize output objects
  alpha <- matrix(0, n_samples, N)
  theta1 <- matrix(0, n_samples, N)
  theta2 <- matrix(0, n_samples, N)
  z <- matrix(0, n_samples, N)
  rho <- numeric(n_samples)
  p <- matrix(0, n_samples, N)
  
  mu <- numeric(n_samples)
  pp_rho <- numeric(n_samples)
  pp_p <- numeric(n_samples)
  nu1 <- numeric(n_samples)
  nu2 <- numeric(n_samples)
  eta <- numeric(n_samples)
  delta1 <- numeric(n_samples)
  delta2 <- numeric(n_samples)
  sigma <- numeric(n_samples)
  
  ### Starting values
  alpha[1,] <- runif(N, -.5, .5)
  theta1[1,] <- runif(N, .1, .5)
  theta2[1,] <- runif(N, -.5, -.1)
  z[1,] <- rbinom(N, 1, .5)
  p[1,] <- rep(.5, N)
  rho[1] <- runif(1, .25, .75)
  nu1[1] <- runif(1, .1, .5)
  nu2[1] <- runif(1, -.5, -.1)
  eta[1] <- runif(1, .1, .5)
  delta1[1] <- runif(1, 0, 1)
  delta2[1] <- runif(1, 0, 1)
  sigma[1] <- runif(1, 0, 1)
  
  ### Sampling
  for (m in 2:n_samples) {
    
    # alpha
    v <- 1 / ((k / sigma[m-1]) + 1 / eta[m-1])
    c <- rowSums(y - x*(theta1[m-1,]^z[m-1,])*(theta2[m-1,]^(1-z[m-1,])), na.rm = T) / 
      sigma[m-1] +  mu[m-1] / eta[m-1]
    alpha[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # theta1
    v <- 1 / (z[(m-1),] * (k / (4 * sigma[m-1])) + 1 / delta1[m-1])
    c <- z[(m-1),] * rowSums((y - alpha[m,]) / x, na.rm = T) / 
      (4 * sigma[m-1]) +  nu1[m-1] / delta1[m-1]
    theta1[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), lower = 0)
    
    # theta2
    v <- 1 / ((1 - z[(m-1),]) * (k / (4 * sigma[m-1])) + 1 / delta2[m-1])
    c <- (1 - z[(m-1),]) * rowSums((y - alpha[m,]) / x, na.rm = T) / 
      (4 * sigma[m-1]) + nu2[m-1] / delta2[m-1]
    theta2[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), upper = 0)
    
    # z
    for (i in 1:N) {
      p[m,i] <- .calc_rho(y[i,], theta1[m,i], theta2[m,i], alpha[m,i], 
                         x[i,], sigma[m-1], rho[m-1])
    }
    z[m,] <- rbinom(N, 1, p[m,])
    
    # Does everybody -- sample level
    pp_p[m] <- prod(p[m,])
    
    # rho
    rho[m] <- rbeta(1, l + sum(z[m,]), u + sum(1 - z[m,]))
    
    # Does everybody -- population level
    pp_rho[m] <- rho[m]^N
    
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
    s2 <- tmp / 2 + s2_delta1
    delta1[m] <- rinvgamma(1, s1, s2)
    
    # delta2
    s1 <- N / 2 + s1_delta2
    tmp <- sum((theta2[m,] - nu2[m])^2)
    s2 <- tmp / 2 + s2_delta2
    delta2[m] <- rinvgamma(1, s1, s2)
    
    # sigma2
    s1 <- Ntot / 2 + s1_sigma
    tmp <- y - alpha[m,] - x * (theta1[m,]^z[m,]) * (theta2[m,]^(1-z[m,]))
    s2 <- sum(tmp^2, na.rm = T) / 2 + s2_sigma
    sigma[m] <- rinvgamma(1, s1, s2)
    
  }
  
  ### Discard burn-in samples
  alpha <- alpha[(burnin+1):n_samples,]
  theta1 <- theta1[(burnin+1):n_samples,]
  theta2 <- theta2[(burnin+1):n_samples,]
  z <- z[(burnin+1):n_samples,]
  rho <- rho[(burnin+1):n_samples]
  p <- p[(burnin+1):n_samples,]
  mu <- mu[(burnin+1):n_samples]
  pp_rho <- pp_rho[(burnin+1):n_samples]
  pp_p <- pp_p[(burnin+1):n_samples]
  nu1 <- nu1[(burnin+1):n_samples]
  nu2 <- nu2[(burnin+1):n_samples]
  eta <- eta[(burnin+1):n_samples]
  delta1 <- delta1[(burnin+1):n_samples]
  delta2 <- delta2[(burnin+1):n_samples]
  sigma <- sigma[(burnin+1):n_samples]
  
  #Samples for combined theta
  theta <- z * theta1 + (1 - z) * theta2
  
  ### Return chains
  list(
    # matrices
    alpha = alpha,
    p = p,
    theta = theta,
    theta2 = theta2,
    theta1 = theta1,
    z = z,
    # vectors
    delta2 = delta2,
    delta1 = delta1,
    eta = eta,
    mu = mu,
    nu2 = nu2,
    nu1 = nu1,
    pp_p = pp_p,
    pp_rho = pp_rho,
    rho = rho,
    sigma = sigma
  )
  
}

# MCMC sampling Normal (1 chain)
.mcmc_normal_single <- function (data, priors, n_samples, burnin) {
  
  ### Get data
  N <- data$N
  k <- data$k
  y <- data$y
  x <- data$x
  Ntot <- data$Ntot
  
  ### Get priors
  am <- priors$mu[1]
  bm <- priors$mu[2]
  an <- priors$nu[1]
  bn <- priors$nu[2]
  qa <- .5
  sa <- .5 * priors$eta2^2
  qt <- .5
  st <- .5 * priors$delta2^2
  qs <- .5
  ss <- .5 * priors$sigma2^2
  
  ### Initialize output objects
  alpha <- matrix(0, n_samples, N)
  theta <- matrix(0, n_samples, N)
  p <- matrix(0, n_samples, N)
  
  mu <- numeric(n_samples)
  nu <- numeric(n_samples)
  eta <- numeric(n_samples)
  delta <- numeric(n_samples)
  sigma <- numeric(n_samples)
  rho <- numeric(n_samples)
  
  ### Starting values
  alpha[1,] <- runif(N, -.5, .5)
  theta[1,] <- runif(N, -.5, .5)
  p[1,] <- rep(.5, N)
  mu[1] <- runif(1, -.5, -.5)
  nu[1] <- runif(1, -.5, -.5)
  eta[1] <- runif(1, 0, 1)
  delta[1] <- runif(1, 0, 1)
  sigma[1] <- runif(1, 0, 1)
  
  ### Sampling
  for (m in 2:n_samples) {
    
    # alpha
    v <- 1 / ((k / sigma[m-1]) + 1 / eta[m-1])
    c <- rowSums(y - x * theta[m-1,], na.rm = T) / sigma[m-1] +  mu[m-1] / eta[m-1]
    alpha[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # theta
    v <- 1 / ((k / (4 * sigma[m-1])) + 1 / delta[m-1])
    c <- rowSums((y - alpha[m,]) / x, na.rm = T) / (4 * sigma[m-1]) + nu[m-1] / delta[m-1]
    theta[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # p
    p[m,] <- pnorm(0, mean = v*c, sd = sqrt(v), lower.tail = FALSE)
    
    # nu
    v <- 1 / ((N / delta[m-1]) + 1 / bn)
    c <- sum(theta[m,]) / delta[m-1] +  an / bn
    nu[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # mu
    v <- 1 / ((N / eta[m-1]) + 1 / bm)
    c <- sum(alpha[m,]) / eta[m-1] +  am / bm
    mu[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # eta
    q_ <- N / 2 + qa
    tmp <- sum((alpha[m,] - mu[m])^2)
    s_ <- tmp / 2 + sa
    eta[m] <- rinvgamma(1, q_, s_)
    
    # delta
    q_ <- N / 2 + qt
    tmp <- sum((theta[m,] - nu[m])^2)
    s_ <- tmp / 2 + st
    delta[m] <- rinvgamma(1, q_, s_)
    
    # sigma2
    q_ <- Ntot / 2 + qs
    tmp <- sum((y - alpha[m,] - x * theta[m,])^2, na.rm = T)
    s_ <- tmp / 2 + ss
    sigma[m] <- rinvgamma(1, q_, s_)
    
  }
  
  ### Discard burn-in samples
  alpha <- alpha[(burnin+1):n_samples,]
  theta <- theta[(burnin+1):n_samples,]
  p <- p[(burnin+1):n_samples,]
  mu <- mu[(burnin+1):n_samples]
  nu <- nu[(burnin+1):n_samples]
  eta <- eta[(burnin+1):n_samples]
  delta <- delta[(burnin+1):n_samples]
  sigma <- sigma[(burnin+1):n_samples]
  
  ### Return chains
  list(
    # matrices
    alpha = alpha,
    p = p,
    theta = theta,
    # vectors
    delta = delta,
    eta = eta,
    mu = mu,
    nu = nu,
    sigma = sigma
  )
  
}

# Mixture with multiple chains
.mcmc_mixture_2 <- function (data, priors, n_samples, burnin, chains, progress) {
  
  mcmc_samples_list <- purrr::map(
    .x = 1:chains,
    .f = \(x) {
      .mcmc_mixture_2_single(data, priors, n_samples, burnin)
    },
    .progress = progress
  )
  
  mcmc_samples <- do.call(dplyr::bind_rows, args = list(mcmc_samples_list))
  
  return(mcmc_samples)
  
}

# Normal with multiple chains
.mcmc_normal <- function (data, priors, n_samples, burnin, chains, progress) {
  
  mcmc_samples_list <- purrr::map(
    .x = 1:chains,
    .f = \(x) {
      .mcmc_normal_single(data, priors, n_samples, burnin)
    },
    .progress = progress
  )
  
  mcmc_samples <- do.call(dplyr::bind_rows, args = list(mcmc_samples_list))
  
  return(mcmc_samples)
  
}

# General function
baid <- function (data, method, priors, n_samples, burnin, chains,
                  check_convergence = FALSE, bayes_factors = FALSE,
                  h1 = c(1, 0, 0), progress = FALSE) {
  
  if (is.data.frame(data)) {
    
    data <- .transform_data(data)
    
  }
  
  if (progress)
    message("Running MCMC for ", method, " model with ", chains, " chains...")
  
  if (method == "mixture_2") {
    
    mcmc_samples <- .mcmc_mixture_2(data, priors, n_samples, burnin, chains, progress)
    
  } else if (method == "normal") {
    
    mcmc_samples <- .mcmc_normal(data, priors, n_samples, burnin, chains, progress)
    
  } else if (method == "mixture_3") {
    
    mcmc_samples <- .mcmc_mixture_3(data, priors, n_samples, burnin, chains, progress)
  }
  
  if (check_convergence) {
    
    if (progress)
      message("\nCalculating convergence criteria...")
    
    if (chains > 1) {
      
      if (method == "mixture_3") {
        checks <- mcmc_samples |> 
          dplyr::select(-c(p, z)) |> 
          .check_convergence(chains = chains)
      } else {
        
        checks <- .check_convergence(mcmc_samples, chains = chains)
        
      }
      
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
        
        bf[[toString(i)]] <- .calc_bfs(data, priors, mcmc_samples, method, h1 = i)
        
      }
      
      
    } else {
      
      bf <- .calc_bfs(data, priors, mcmc_samples, method, h1)
      
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
  
  p <- if_else(ul < 0, 0, 1)
  
  list(
    theta = y_obs,
    ll = ll,
    ul = ul,
    p = p
  )
  
}

