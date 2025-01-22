# Load Data (Copied from Thiele et al., 2017)

data_1a = function ()
{
  inDat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data1/master/sysfactorial/SysFac02/SysFac02.all'))
  
  names=c("sub","trial","block","trialBlock","sizeChange","angleChange","changeType","bothDiff",
          "size1","angle1","size2","angle2","resp","acc","rt","bad","badCount",
          "tS","tF","tT","tB","tP","tR","percMem","feedback")
  colnames(inDat)=names
  
  bad1 = inDat$sub %in% c(4,21,33)
  bad2 = inDat$rt<.3 | inDat$rt>5
  bad3 = inDat$block==0
  bad4 = inDat$acc==0
  bad=(bad1 | bad2 | bad3 |bad4)
  return(inDat[!bad,])	
}



# Restructure Data
.transform_data_thiele <- function (df) {
  
  dat <- df[df$bothDiff==1,]
  y <- dat$rt
  sub <- as.integer(as.factor(dat$sub))
  sj <- 2 * dat$sizeChange - 5
  sk <- 2 * dat$angleChange - 5
  sjk <- sj * sk
  
  data <- list(
    N = max(sub),
    Ji = tapply(sub, sub, length) |> as.vector(),
    y = .make_matrix(split(y, sub)),
    sj = .make_matrix(split(sj, sub)),
    sk = .make_matrix(split(sk, sub)),
    sjk = .make_matrix(split(sjk, sub)),
    Ntot = length(sub)
  )
  
  return(data)
  
}

# Calculate Rho
.calc_rho_thiele <- function (y, sj, sk, sjk, mu, a, b, g, t1, t2, J, s2, p) {
  
  out <- numeric(3)
  
  y_tmp <- (y - mu - a - sj*b - sk*g) / sjk
  
  tmp1 <- exp(- (2 * t1 * sum(y_tmp, na.rm = T) - J * t1^2) / (2 * s2))
  tmp2 <- exp(- (2 * (t1 - t2) * sum(y_tmp, na.rm = T) + J * (t2^2 - t1^2)) / (2 * s2))
  out[1] <- 1 / (1 + p[3]/p[1] * tmp1 + p[2]/p[1] * tmp2)
  
  tmp1 <- exp(- (J * t1^2 - 2 * t1 * sum(y_tmp, na.rm = T)) / (2 * s2))
  tmp2 <- exp(- (J * t2^2 - 2 * t2 * sum(y_tmp, na.rm = T)) / (2 * s2))
  out[3] <- 1 / (1 + p[1]/p[3] * tmp1 + p[2]/p[3] * tmp2)
  
  out[2] <- 1 - sum(out[c(1, 3)])
  out <- ifelse(out < 0, 0, out)
  
  return (out)
}

# MCMC sampling Mixture with 3 components (1 chain)
.mcmc_mixture_thiele <- function (data, priors, n_samples, burnin) {
  
  ### Get data
  N    <- data$N
  Ji   <- data$Ji
  y    <- data$y
  sj   <- data$sj
  sk   <- data$sk
  sjk  <- data$sjk
  Ntot <- data$Ntot
  
  ### Get priors
  q     <- priors$q
  am    <- priors$mu[1]
  bm    <- priors$mu[2]
  ab    <- priors$phiBeta[1]
  bb    <- priors$phiBeta[2]
  ag    <- priors$phiGamma[1]
  bg    <- priors$phiGamma[2]
  a1    <- priors$nu1[1]
  b1    <- priors$nu1[2]
  a2    <- priors$nu2[1]
  b2    <- priors$nu2[2]
  qa    <- .5
  sa    <- .5 * priors$eta2alpha^2
  qb    <- .5
  sb    <- .5 * priors$eta2beta^2
  qg    <- .5
  sg    <- .5 * priors$eta2gamma^2
  q1    <- .5
  s1    <- .5 * priors$delta2_1^2
  q2    <- .5
  s2    <- .5 * priors$delta2_2^2
  qs    <- .5
  ss    <- .5 * priors$sigma2^2
  
  # initiate objects
  alpha <- matrix(0, n_samples, N)
  mu <- numeric(n_samples)
  etaAlpha <- numeric(n_samples)
  beta <- matrix(0, n_samples, N)
  phiBeta <- numeric(n_samples)
  etaBeta <- numeric(n_samples)
  gamma <- matrix(0, n_samples, N)
  phiGamma <- numeric(n_samples)
  etaGamma <- numeric(n_samples)
  
  rho <- matrix(0, n_samples, 3)
  z <- array(0, dim = c(n_samples, 3, N))
  p <- array(0, dim = c(n_samples, 3, N))
  theta1 <- matrix(0, n_samples, N)
  theta2 <- matrix(0, n_samples, N)
  nu1 <- numeric(n_samples)
  nu2 <- numeric(n_samples)
  delta1 <- numeric(n_samples)
  delta2 <- numeric(n_samples)
  sigma <- numeric(n_samples)
  
  # starting values
  alpha[1,] <- runif(N, -.1, .1)
  mu[1] <- runif(1, -.1, .1)
  etaAlpha[1] <- runif(1, 0, 1)
  beta[1,] <- runif(N, -.1, .1)
  phiBeta[1] <- runif(1, -.1, .1)
  etaBeta[1] <- runif(1, 0, 1)
  gamma[1,] <- runif(N, -.1, .1)
  phiGamma[1] <- runif(1, -.1, .1)
  etaGamma[1] <- runif(1, 0, 1)
  theta1[1,] <- runif(N, .1, .5)
  theta2[1,] <- runif(N, -.5, -.1)
  rho[1,] <- rdirichlet(1, c(1, 1, 1))
  z[1,,] <- rmultinom(N, 1, rho[1,])
  p[1,,] <- rho[1,]
  nu1[1] <- runif(1, 0, 1)
  nu2[1] <- runif(1, -1, 0)
  delta1[1] <- runif(1, 0, .5)
  delta2[1] <- runif(1, 0, .5)
  sigma[1] <- runif(1, 0, 1)
  
  
  ### Sampling
  for (m in 2:n_samples) {
    
    # mu
    v <- 1 / ((Ntot / sigma[m-1]) + 1 / bm)
    c <- sum(y - alpha[m-1,] - sj*beta[m-1,] - sk*gamma[m-1,] - sjk*((theta1[m-1,]^z[(m-1),1,])*(0^z[(m-1),3,])*(theta2[m-1,]^z[(m-1),2,])), na.rm = T) / 
      sigma[m-1] +  am / bm
    mu[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # alpha
    v <- 1 / ((Ji / sigma[m-1]) + 1 / etaAlpha[m-1])
    c <- rowSums(y - mu[m] - sj*beta[m-1,] - sk*gamma[m-1,] - sjk*((theta1[m-1,]^z[(m-1),1,])*(0^z[(m-1),3,])*(theta2[m-1,]^z[(m-1),2,])), na.rm = T) / 
      sigma[m-1]
    alpha[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # beta
    v <- 1 / ((Ji / sigma[m-1]) + 1 / etaBeta[m-1])
    c <- rowSums((y - mu[m] - alpha[m,] - sk*gamma[m-1,] - sjk*((theta1[m-1,]^z[(m-1),1,])*(0^z[(m-1),3,])*(theta2[m-1,]^z[(m-1),2,]))) / sj, na.rm = T) / 
      sigma[m-1] +  phiBeta[m-1] / etaBeta[m-1]
    beta[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # gamma
    v <- 1 / ((Ji / sigma[m-1]) + 1 / etaGamma[m-1])
    c <- rowSums((y - mu[m] - alpha[m,] - sj*beta[m,] - sjk*((theta1[m-1,]^z[(m-1),1,])*(0^z[(m-1),3,])*(theta2[m-1,]^z[(m-1),2,]))) / sk, na.rm = T) / 
      sigma[m-1] +  phiGamma[m-1] / etaGamma[m-1]
    gamma[m,] <- rnorm(N, mean = v*c, sd = sqrt(v))
    
    # theta1
    v <- 1 / (z[(m-1),1,] * (Ji / sigma[m-1]) + 1 / delta1[m-1])
    c <- z[(m-1),1,] * rowSums((y - mu[m]- alpha[m,] - sj*beta[m,] - sk*gamma[m,]) / sjk, na.rm = T) / 
      (sigma[m-1]) +  nu1[m-1] / delta1[m-1]
    theta1[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), lower = 0)
    
    # theta2
    v <- 1 / (z[(m-1),2,] * (Ji / sigma[m-1]) + 1 / delta2[m-1])
    c <- z[(m-1),2,] * rowSums((y - mu[m]- alpha[m,] - sj*beta[m,] - sk*gamma[m,]) / sjk, na.rm = T) / 
      (sigma[m-1]) +  nu2[m-1] / delta2[m-1]
    theta2[m,] <- rtnorm(N, mean = v*c, sd = sqrt(v), upper = 0)
    
    # z
    for (i in 1:N) {
      p[m,,i] <- .calc_rho_thiele(y[i,], sj[i,], sk[i,], sjk[i,], mu[m], alpha[m,i], beta[m,i], gamma[m,i], theta1[m,i], theta2[m,i], Ji[i], sigma[m-1], rho[m-1,])
      z[m,,i] <- rmultinom(1, 1, p[m,,i])
    }
    
    # rho
    rho[m,] <- rdirichlet(1, q + rowSums(z[m,,]))
    
    # nu1
    v <- 1 / ((N / delta1[m-1]) + 1 / b1)
    c <- sum(theta1[m,]) / delta1[m-1] +  a1 / b1
    nu1[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # nu2
    v <- 1 / ((N / delta2[m-1]) + 1 / b2)
    c <- sum(theta2[m,]) / delta2[m-1] +  a2 / b2
    nu2[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # phiBeta
    v <- 1 / ((N / etaBeta[m-1]) + 1 / bb)
    c <- sum(beta[m,]) / etaBeta[m-1] +  ab / bb
    phiBeta[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # phiGamma
    v <- 1 / ((N / etaGamma[m-1]) + 1 / bg)
    c <- sum(gamma[m,]) / etaGamma[m-1] +  ag / bg
    phiGamma[m] <- rnorm(1, mean = v*c, sd = sqrt(v))
    
    # etaAlpha
    q_ <- N / 2 + qa
    tmp <- sum(alpha[m,]^2)
    s_ <- tmp / 2 + sa
    etaAlpha[m] <- rinvgamma(1, q_, s_)
    
    # etaBeta
    q_ <- N / 2 + qb
    tmp <- sum((beta[m,] - phiBeta[m])^2)
    s_ <- tmp / 2 + sb
    etaBeta[m] <- rinvgamma(1, q_, s_)
    
    # etaGamma
    q_ <- N / 2 + qg
    tmp <- sum((gamma[m,] - phiGamma[m])^2)
    s_ <- tmp / 2 + sg
    etaGamma[m] <- rinvgamma(1, q_, s_)
    
    # delta1
    q_ <- N / 2 + q1
    tmp <- sum((theta1[m,] - nu1[m])^2)
    s_ <- tmp / 2 + s1
    delta1[m] <- rinvgamma(1, q_, s_)
    
    # delta2
    q_ <- N / 2 + q2
    tmp <- sum((theta2[m,] - nu2[m])^2)
    s_ <- tmp / 2 + s2
    delta2[m] <- rinvgamma(1, q_, s_)
    
    # sigma2
    q_ <- length(y) / 2 + qs
    tmp <- y - mu[m] - alpha[m,] - sj*beta[m,] - sk*gamma[m,] - sjk*((theta1[m,]^z[m,1,])*(0^z[m,3,])*(theta2[m,]^z[m,2,]))
    s_ <- sum(tmp^2, na.rm = T) / 2 + ss
    sigma[m] <- rinvgamma(1, q_, s_)
    
  }
  
  ### Discard burn-in samples
  alpha <- alpha[(burnin+1):n_samples,]
  beta <- beta[(burnin+1):n_samples,]
  gamma <- gamma[(burnin+1):n_samples,]
  theta1 <- theta1[(burnin+1):n_samples,]
  theta2 <- theta2[(burnin+1):n_samples,]
  z <- z[(burnin+1):n_samples,,]
  rho <- rho[(burnin+1):n_samples,]
  p <- p[(burnin+1):n_samples,,]
  
  mu <- mu[(burnin+1):n_samples]
  phiBeta <- phiBeta[(burnin+1):n_samples]
  phiGamma <- phiGamma[(burnin+1):n_samples]
  nu1 <- nu1[(burnin+1):n_samples]
  nu2 <- nu2[(burnin+1):n_samples]
  etaAlpha <- etaAlpha[(burnin+1):n_samples]
  etaBeta <- etaBeta[(burnin+1):n_samples]
  etaGamma <- etaGamma[(burnin+1):n_samples]
  delta1 <- delta1[(burnin+1):n_samples]
  delta2 <- delta2[(burnin+1):n_samples]
  sigma <- sigma[(burnin+1):n_samples]
  
  #Samples for combined theta
  theta <- z[,1,] * theta1 + z[,2,] * theta2 + z[,3,] * 0
  
  ### Return chains
  list(
    # arrays/matrices
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    p = p,
    rho = rho,
    theta = theta,
    theta1 = theta1,
    theta2 = theta2,
    z = z,
    # vectors
    delta1 = delta1,
    delta2 = delta2,
    etaAlpha = etaAlpha,
    etaBeta = etaBeta,
    etaGamma = etaGamma,
    mu = mu,
    nu1 = nu1,
    nu2 = nu2,
    phiBeta = phiBeta,
    phiGamma = phiGamma,
    sigma = sigma
  )
  
}

# log of multivariate beta function
.lmbeta <- function (alpha) {
  
  sum(lgamma(alpha)) - lgamma(sum(alpha))
  
}

# Naive Estimation and classification
naid_thiele <- function (dat, ci = .95) {
  
  df <- dat |> 
    filter(bothDiff == 1) |> 
    mutate(
      sub = dense_rank(sub),
      sj = 2 * sizeChange - 5,
      sk = 2 * angleChange - 5
    )
  
  I <- max(df$sub)	
  obs <- se <- dof <- numeric(I)
  for (i in 1:I) {
    tmp <- df[df$sub == i, ]
    w <- summary(lm(rt ~ sj*sk, data = tmp))
    dof[i] <- w$df[2]
    obs[i] <- w$coefficients[4,1]
    se[i] <- w$coefficients[4,2]
  }
  
  ll <- obs - qt(1 - (1 - ci)/2, dof) * se
  ul <- obs + qt(1 - (1 - ci)/2, dof) * se
  
  list(
    theta = obs, 
    ll = ll, 
    ul = ul
  )
  
}


# General function
baid_thiele <- function (data, priors, n_samples, burnin, chains,
                         check_convergence = FALSE, bayes_factors = FALSE,
                         h1 = c(1, 0, 0), progress = FALSE, parallel = TRUE) {
  
  if (is.data.frame(data)) {
    
    data <- .transform_data(data)
    
  }
  
  if (progress)
    message("Running MCMC for mixture_3 model with ", chains, " chains...")
  
  if (parallel) {
    
    future::plan(multisession, workers = chains)
    
    mcmc_samples_list <- furrr::future_map(
      .x = 1:chains,
      .f = \(x) {
        .mcmc_mixture_thiele(data, priors, n_samples, burnin)
      },
      .options = furrr_options(seed = TRUE),
      .progress = FALSE
    )
    
    
  } else {
    
    mcmc_samples_list <- purrr::map(
      .x = 1:chains,
      .f = \(x) {
        .mcmc_mixture_thiele(data, priors, n_samples, burnin)
      },
      .progress = progress
    )
    
  }
  
  mcmc_samples <- do.call(dplyr::bind_rows, args = list(mcmc_samples_list))
  
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
        
        bf[[toString(i)]] <- .calc_bfs(data, priors, mcmc_samples, "mixture_3", h1 = i)
        
      }
      
      
    } else {
      
      bf <- .calc_bfs(data, priors, mcmc_samples, "mixture_3", h1)
      
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
