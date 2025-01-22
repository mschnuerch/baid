
my_theme <- theme_classic() +
  theme(
    text = element_text(size = 6, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.ticks.length = unit(.1, "cm"),
    plot.title = element_text(size = 12, face = "bold", hjust = -.25),
    legend.key.size = unit(.25, 'cm')
  ) 




plot_naive <- function(theta, ll, ul, xlim = NULL, ylim = NULL, col_vals = c("firebrick", "white", "yellow"),
                       point_size = 1) {
  
  # Combine
  df <- data.frame(
    y = theta,
    ll = ll,
    ul = ul
  ) |> 
    mutate(col = if_else(ul < 0, "a", "b"),
           col = if_else(ll > 0, "c", col))
  
  
  p <- 
    df |> 
    arrange(y) |> 
    mutate(x = 1:n()) |> 
    ggplot(aes(x)) +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = .5) +
    geom_segment(
      aes(xend = x, y = ll, yend = ul), 
      color = "grey",
      lwd = .25,
      arrow = arrow(angle = 90, ends = "both", length = unit(.05, "inches"))) +
    geom_point(
      aes(y = y, fill = col), 
      size = point_size, 
      shape = 21) +
    labs(x = "Individuals", y = "Observed Effect") +
    scale_fill_manual(
      values = col_vals
    ) +
    scale_x_continuous(breaks = c(1, nrow(df))) +
    my_theme +
    theme(legend.position = "none") +
    coord_capped_cart(bottom = "both", left = "both")
  
  if (!is.null(xlim)) 
    p <- p + xlim(xlim)
  
  if (!is.null(ylim)) 
    p <- p + ylim(ylim)
  
  p
  
}

plot_hlmm <- function (res_naive, res_mix, ylim, point_size = 1) {
  
  data.frame(
    obs = res_naive$theta, 
    est = colMeans(res_mix$mcmc_samples$theta),
    ll = apply(res_mix$mcmc_samples$theta, 2, \(x) quantile(x, probs = .025)),
    ul = apply(res_mix$mcmc_samples$theta, 2, \(x) quantile(x, probs = .975))
  ) |> 
    arrange(obs) |> 
    mutate(index = 1:n()) |> 
    ggplot(aes(x = index)) +
    geom_hline(yintercept = 0, lwd = .5, linetype = "dashed") +
    geom_ribbon(aes(ymin = ll, ymax = ul), fill = "dodgerblue2",
                alpha = .15) +
    geom_line(aes(y = obs), lwd = .75, color = "darkgrey") +
    geom_point(aes(y = est), shape = 21, size = point_size, fill = "dodgerblue2") +
    labs(x = "Individuals", y = "Estimated Effect") +
    scale_y_continuous(
      limits = ylim
    ) +
    scale_x_continuous(breaks = c(1, length(res_naive$theta))) +
    theme_classic() +
    coord_capped_cart(bottom = "both", left = "both") +
    my_theme
  
}

plot_classification <- function (res_naive, res_mix, legend_pos = c(.85, .25)) {
  
  tmp = cbind(
    res_naive$theta,
    t(apply(res_mix$mcmc_samples$p, c(2,3), mean)),
    t(apply(res_mix$mcmc_samples$p, c(2,3), \(x) quantile(x, probs = c(.025)) |> as.vector())),
    t(apply(res_mix$mcmc_samples$p, c(2,3), \(x) quantile(x, probs = c(.975)) |> as.vector()))
  ) |> 
    data.frame()
  names(tmp) = c("y_obs", "p1", "p2", "p3", "p1_l", "p2_l", "p3_l", "p1_u", "p2_u", "p3_u")
  
  ggplot(tmp, aes(x = y_obs)) +
    geom_vline(xintercept = 0, linetype = "dashed", lwd = .5) +
    geom_hline(yintercept = .5, linetype = "dashed", lwd = .5) +
    geom_ribbon(aes(ymin = p3_l, ymax = p3_u), fill = "dodgerblue4",
                alpha = .15) +
    geom_ribbon(aes(ymin = p2_l, ymax = p2_u), fill = "firebrick",
                alpha = .15) +
    geom_ribbon(aes(ymin = p1_l, ymax = p1_u), fill = "darkgreen",
                alpha = .15) +
    geom_line(aes(y = p3, color = "c"), lwd = .75) +
    geom_line(aes(y = p2, color = "b"), lwd = .75) +
    geom_line(aes(y = p1, color = "a"), lwd = .75) +
    labs(x = "Observed Effect", y = "Class Probabilities") +
    scale_color_manual(
      element_blank(),
      values = c("darkgreen", "firebrick", "dodgerblue4"),
      labels = c(expression(italic(p)[i1]), 
                 expression(italic(p)[i2]), 
                 expression(italic(p)[i3]))
    ) +
    lims(y = c(0, 1)) +
    theme_classic() +
    coord_capped_cart(bottom = "both", left = "both") +
    my_theme +
    theme(legend.position = "inside",
          legend.position.inside = legend_pos,
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent', color = NA))
  
}

plot_bfs <- function (res_mix, h1, legend_pos, ylim = NULL) {
  
  bfs <- res_mix$bf |> unlist() |> log()
  type <- rep(c("Sample", "Population"), length(res_mix$bf))
  h1 <- rep(h1, each = 2) |> factor(levels = h1)
  
  p <- data.frame(
    bfs, type, h1
  ) |> 
    ggplot(aes(x = h1, fill = type)) +
    geom_bar(aes(y = bfs), stat = "identity", position = position_dodge(),
             width = .5, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = .5) +
    scale_fill_manual(element_blank(), values = c("dimgrey", "grey90")) +
    labs(x = element_blank(), y = "log(BF) vs General Model") +
    my_theme +
    theme(legend.position = "inside",
          legend.position.inside =  legend_pos,
          # legend.text = element_text(size = 5),
          # legend.background = element_rect(color = "black", linewidth = .1)
          )
  
  if (!is.null(ylim)) {
    p +
      coord_capped_cart(bottom = "both", ylim = ylim)
  } else {
    p +
      coord_capped_cart(bottom = "both")
  }
  
}

### Deprecated

plot_comparison <- function (samples_mix, samples_norm, tZ, type = "prob") {
  
  # general settings
  agree = "green"
  unclear = "grey"
  disagree = "red"
  alpha1 = .1
  alpha2 = .4
  
  # extract p
  p_mix = colMeans(samples_mix$mcmc_samples$p)
  p_normal = colMeans(samples_norm$mcmc_samples$p)
  
  if (type == "prob") {
    
    df = data.frame(
      tZ = tZ,
      p_mix,
      p_normal
    )
    
    ggplot(df, aes(x = p_normal, y = p_mix)) +
      geom_rect(data = data.frame(xmin = 0, xmax = .25, 
                                  ymin = 0, ymax = .25),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha2, fill = agree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .75, xmax = 1, 
                                  ymin = .75, ymax = 1),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha2, fill = agree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .25, xmax = .75, 
                                  ymin = .25, ymax = .75),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha1, fill = unclear, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .25, xmax = .75, 
                                  ymin = 0, ymax = .25),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha1, fill = disagree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .25, xmax = .75, 
                                  ymin = .75, ymax = 1),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha1, fill = disagree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = 0, xmax = .25, 
                                  ymin = .25, ymax = .75),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha1, fill = disagree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .75, xmax = 1, 
                                  ymin = .25, ymax = .75),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha1, fill = disagree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = .75, xmax = 1, 
                                  ymin = 0, ymax = .25),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha2, fill = disagree, inherit.aes = FALSE) +
      geom_rect(data = data.frame(xmin = 0, xmax = .25, 
                                  ymin = .75, ymax = 1),
                mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = alpha2, fill = disagree, inherit.aes = FALSE) +
      geom_hline(
        yintercept = .25,
        lwd = .5
      ) +
      geom_hline(
        yintercept = .75,
        lwd = .5
      ) +
      geom_vline(
        xintercept = .25,
        lwd = .5
      ) +
      geom_vline(
        xintercept = .75,
        lwd = .5
      ) +
      geom_vline(
        xintercept = .5,
        lwd = .5,
        linetype = "dashed"
      ) +
      geom_hline(
        yintercept = .5,
        lwd = .5,
        linetype = "dashed"
      ) +
      geom_point(
        aes(fill = as.factor(tZ)),
        size = 4,
        shape = 21) +
      scale_fill_manual(
        "True Class",
        values = c("firebrick", "dodgerblue2"),
        labels = c("Negative", "Positive")
      ) +
      scale_x_continuous(
        "\nNormal",
        breaks = c(0, .125, .25, .5, .75, .875, 1),
        labels = c("", "Certain", "3-to-1", "Uncertain", "3-to-1", "Certain", "")
      ) +
      scale_y_continuous(
        "Mixture",
        breaks = c(0, .125, .25, .5, .75, .875, 1),
        labels = c("", "Certain", "3-to-1", "Uncertain", "3-to-1", "Certain", "")
      ) +
      my_theme +
      theme(
        axis.ticks = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
      coord_capped_cart(
        xlim = c(0, 1),
        ylim = c(0, 1),
        bottom = "both",
        left = "both"
      ) +
      ggtitle("Classification Probability (Positive)")
    
  } else if (type == "log-odds") {
    
    # transform p to lodds
    p_mix = p_mix |> qlogis()
    p_normal = p_normal |> qlogis()
    
    df = data.frame(
      tZ = tZ,
      p_mix,
      p_normal
    )
    
    ggplot(df, aes(x = p_normal, y = p_mix)) +
      geom_hline(
        yintercept = 0,
        lwd = .5,
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = 0,
        lwd = .5,
        linetype = "dashed"
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        linewidth = .75
      ) +
      geom_point(
        aes(fill = as.factor(tZ)),
        size = 4,
        shape = 21) +
      scale_fill_manual(
        "True Class",
        values = c("firebrick", "dodgerblue2"),
        labels = c("Negative", "Positive")
      ) +
      labs(x = "Normal", y = "Mixture") +
      theme_bw() +
      theme(
        text = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle("Classification Log-Odds (Positive)")
    
  }
  
  
  
}

plot_theta <- function (theta, data, ci = .90, ylim = NULL, xlim = NULL, ...) {
  
  ### Get data
  N = data$N
  J = data$J
  y = data$y
  x = data$x
  
  # Observed 
  y_obs <- numeric(N)
  for (i in 1:N) {
    tmp = tapply(y[i,], x[i,], mean)
    y_obs[i] <- as.vector(tmp[2] - tmp[1])
  }
  
  # Combine
  df <- data.frame(
    y_obs = y_obs,
    y = colMeans(theta),
    ll = apply(theta, 2, quantile, probs = (1 - ci) / 2),
    ul = apply(theta, 2, quantile, probs = 1 - (1 - ci) / 2)
  )
  
  
  p <- df |> 
    arrange(y_obs) |> 
    mutate(x = 1:n()) |> 
    ggplot(aes(x)) +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = .75) +
    geom_segment(
      aes(xend = x, y = ll, yend = ul), 
      lwd = .75,
      color = "firebrick",
      alpha = .3) +
    geom_line(
      aes(y = y_obs),
      lwd = 1, 
      color = "darkgrey"
    ) +
    geom_point(
      aes(y = y), 
      size = 3, 
      shape = 21, 
      fill = "firebrick") +
    labs(x = "Individuals", y = expression(theta)) +
    my_theme +
    coord_capped_cart(bottom = "both", left = "both")
  
  if (!is.null(xlim)) 
    p <- p + xlim(xlim)
  
  if (!is.null(ylim)) 
    p <- p + ylim(ylim)
  
  p
  
}

plot_class_prob <- function (p, data, ci = .95, ylim = NULL, xlim = NULL, ...) {
  
  ### Get data
  N = data$N
  y = data$y
  x = data$x
  
  # Observed 
  y_obs <- numeric(N)
  for (i in 1:N) {
    tmp = tapply(y[i,], x[i,], mean)
    y_obs[i] <- as.vector(tmp[2] - tmp[1])
  }
  
  # Combine
  df <- data.frame(
    y_obs = y_obs,
    y = colMeans(p),
    ll = apply(p, 2, quantile, probs = (1 - ci) / 2),
    ul = apply(p, 2, quantile, probs = 1 - (1 - ci) / 2)
  )
  
  
  p <- df |> 
    arrange(y_obs) |> 
    ggplot(aes(x = y_obs)) +
    geom_hline(yintercept = .5, linetype = "dashed", lwd = .5) +
    geom_vline(xintercept = 0, linetype = "dashed", lwd = .5) +
    # geom_segment(
    #   aes(xend = y_obs, y = ll, yend = ul), 
    #   lwd = .75,
    #   color = "dodgerblue4",
    #   alpha = .3) +
    geom_ribbon(aes(ymin = ll, ymax = ul), fill = "darkgreen",
                alpha = .15) +
    geom_point(
      aes(y = y), 
      size = 3, 
      shape = 21, 
      fill = "darkgreen") +
    labs(x = "Observed Effect", y = expression(rho)) +
    my_theme +
    coord_capped_cart(ylim = c(0, 1), bottom = "both", left = "both")
  
  if (!is.null(xlim)) 
    p <- p + xlim(xlim)
  
  if (!is.null(ylim)) 
    p <- p + ylim(ylim)
  
  p
  
}

plot_ternary <- function (grid, bw = .02) {
  
  ggtern(data = grid, aes(x = p1, y = p2, z = p3)) +
    geom_hex_tern(aes(value = density), fun = mean, 
                  binwidth = bw) +
    scale_fill_gradient("Density", high = "firebrick",low = "ghostwhite",) +
    theme_noticks() +
    theme_light() +
    scale_T_continuous(
      expression(rho[3]), 
      breaks = seq(0, 1, by = 0.2), 
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")
    ) +
    scale_L_continuous(
      expression(rho[2]), 
      breaks = seq(0, 1, by = 0.2), 
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")
    ) +
    scale_R_continuous(
      expression(rho[1]), 
      breaks = seq(0, 1, by = 0.2), 
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")
    ) +
    theme(axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 12, color = "black", hjust = 0.5),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

plot_marginal <- function (q, ylim = c(0, 3)) {
  
  ggplot(data.frame(x = 0), aes(x)) +
    stat_function(
      fun = dbeta, 
      args = list(shape1 = q[1], shape2 = sum(q[2:3])), 
      xlim = c(.001, .999), n = 500,
      lwd = 1) +
    theme_classic() +
    scale_x_continuous(
      expression(paste(rho[1], ", ", rho[2], ", ", rho[3])),
      expand = c(0,0), 
      limits = c(0, 1.01),
      breaks = seq(0, 1, .2), labels = seq(0, 1, .2)) +
    scale_y_continuous(
      "Density",
      expand = c(0, 0)
    ) +
    coord_capped_cart(ylim = ylim, bottom = "right", left = "top", gap = 0) +
    theme(
      axis.ticks.length = unit(.25, "cm"),
      axis.text = element_text(color = "black", size = 10),
      text = element_text(color = "black", size = 10)
    )
  # ggtitle(title)
  
}

