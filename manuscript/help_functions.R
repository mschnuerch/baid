## ========================================================================= ##
##  HELP FUNCTIONS
## ========================================================================= ##

# Restructure Data
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

# requires df with columns sub, cond, y
.transform_data <- function (df) {
  
  list(
    N = length(unique(df$sub)),
    k = tapply(df$sub, df$sub, length) |> as.vector(),
    y = .make_matrix(split(df$y, df$sub)),
    x = .make_matrix(split(df$cond, df$sub)),
    Ntot = nrow(df)
  )
  
}

# Conditional Posterior rho* (Mix with 3 Components)
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

# Convergence checks
.check_convergence <- function (samps, chains) {
  
  tmp <- as.matrix(samps)
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

# Bayes factors
.calc_bfs <- function (data, priors, mcmc_samples, method, h1) {
  
  if (method == "mixture_2") {
    
    prior_prob <- beta(data$N + priors$rho[1], priors$rho[2]) / 
      beta(priors$rho[1], priors$rho[2])
    post_prob_population <- mean(mcmc_samples$pp_rho)
    post_prob_sample <- mean(mcmc_samples$pp_p)
    
    # BF sample
    bf_samp <- post_prob_sample / prior_prob
    
    # BF population
    bf_pop <- post_prob_population / prior_prob
    
    # return
    bf <- c(
      sample = bf_samp,
      population = bf_pop
    )
    
  } else if (method == "mixture_3") {
    
    # Prior prob for constraint
    par_a <- sum(priors$q * h1)
    par_b <- priors$q[!h1]
    prior_prob <- exp(.lmbeta(c(par_a + data$N, par_b)) - .lmbeta(c(par_a, par_b)))
    
    # Posterior prob population
    post_prob_population <- exp(data$N * log(rowSums(mcmc_samples$rho[,h1 * 1:3, drop = F]))) |> mean()
    
    # Posterior prob sample
    p_sum <- apply(mcmc_samples$p, c(1, 3), \(x) sum(x[h1 * c(1, 2, 3)]))
    post_prob_sample <- apply(p_sum, 1, \(x) exp(sum(log(x)))) |> mean()
    
    # BF sample
    bf_samp <- post_prob_sample / prior_prob
    
    # BF population
    bf_pop <- post_prob_population / prior_prob
    
    # return
    bf <- c(
      sample = bf_samp,
      population = bf_pop
    )
    
  } else if (method == "normal") {
    
    # prior prob
    n_reps <- 50000
    nus <- rnorm(n_reps, priors$nu[1], sqrt(priors$nu[2]))
    deltas <- rinvgamma(n_reps, .5, .5 * priors$delta2^2)
    p_reps <- data$N * pnorm(0, nus, sqrt(deltas), lower.tail = F, log.p = T)
    prior_prob <- mean(exp(p_reps))
    
    # posterior prob
    p_post <- data$N * pnorm(0, mcmc_samples$nu, sqrt(mcmc_samples$delta), 
                            lower.tail = F, log.p = T)
    post_prob <- mean(exp(p_post))
    
    # BF
    bf <- post_prob / prior_prob
    
    
  }
  
  return(bf)
}

# Generate data
gen_data <- function(sim_settings, true_model) {
  
  N <- sim_settings$N
  J <- sim_settings$J
  tMu <- sim_settings$tMu
  tEta <- sim_settings$tEta
  tRho <- sim_settings$tRho
  tNu1 <- sim_settings$tNu1
  tDelta1 <- sim_settings$tDelta1
  tNu2 <- sim_settings$tNu2
  tDelta2 <- sim_settings$tDelta2
  tSigma <- sim_settings$tSigma
  tNu <- sim_settings$tNu
  tDelta <- sim_settings$tDelta 
  
  x <- rep(c(-.5, .5), J/2) |> matrix(nrow = N, ncol = J, byrow = T)
  
  tAlpha <- rnorm(N, tMu, sqrt(tEta))
  
  if (true_model == "mixture_2") {
    
    tZ <- rbinom(N, 1, tRho)
    tTheta <- ifelse(
      tZ == 1,
      rtnorm(sum(tZ), tNu1, sqrt(tDelta1), lower = 0),
      rtnorm(sum(1 - tZ), tNu2, sqrt(tDelta2), upper = 0)
    )
    
  } else if (true_model == "mixture_3") {
    
    tZ <- rmultinom(N, 1, tRho)
    tTheta <- numeric(N)
    for (i in 1:N) {
      
      tTheta[i] <- case_when(
        tZ[1,i] == 1 ~ rtnorm(1, tNu1, sqrt(tDelta1), lower = 0),
        tZ[2,i] == 1 ~ rtnorm(1, tNu2, sqrt(tDelta2), upper = 0),
        tZ[3,i] == 1 ~ 0
      )
      
    }
    
  } else {
    
    tTheta <- rnorm(N, tNu, sqrt(tDelta))
    tZ <- ifelse(tTheta < 0, 0, 1)
    
  }
  
  y <- tAlpha  + x * tTheta + rnorm(N * J, 0, sqrt(tSigma))
  
  list(
    N = N,
    k = rep(J, N),
    y = y,
    x = x,
    Ntot = N * J,
    tZ = tZ,
    tTheta = tTheta
  )
  
}

# Get Data for Example 2 (Nadarevic & Rinnewitz, 2011)
data_te <- function () {
  
  dat <- read_csv(url("https://osf.io/jfz25/download"),
                  show_col_types = F) %>% 
    dplyr::filter(runningtrial == "TruthJudgment") %>% 
    dplyr::select(subject, statement, status, repeated, truthrating) %>% 
    dplyr::mutate(item = ifelse(status, 0.5, -0.5),
                  sub = as.integer(as.factor(subject)),
                  cond = ifelse(repeated == "No", -0.5, 0.5),
                  Y = truthrating) %>% 
    dplyr::select(sub, item, cond, Y) %>% 
    dplyr::mutate(y = (Y - 3)/2)
  
  return (dat)
  
}

# Get Data for Example 3, In-Out-Effect (Ingendahl et al., 2023, Exp. 1)
data_ioe <- function () {
  
  url <- "https://osf.io/download/vwjax/"
  filename <- "studyraw2OS.xlsx"
  httr::GET(url, httr::write_disk(filename, overwrite = TRUE))
  
  df <-readxl::read_xlsx(filename) |>
    mutate(
      wordtype = case_when(
        origin == "samec" ~ "Same Consonants",
        origin == "newc"  ~ "New Consonants",
        TRUE              ~ "Conditioned"  
      )
    ) |>
    filter(!OSid %in% c("766109203", "415721123")) |>
    filter(current == "rating") |>
    mutate(
      direction = ifelse(direction == "in", 0.5, -0.5),
      liking_norm = ((as.numeric(response))-5)/4,
      id = dense_rank(OSid)
    ) |>
    filter(wordtype == "New Consonants") |>
    dplyr::select(c("direction", "liking_norm", "id")) |>
    arrange(id)
  
  unlink(filename)
  
  dat <- df |> 
    rename(sub = id, 
           cond = direction,
           y = liking_norm) |> 
    select(sub, cond, y)
  
  return (dat)
  
}



