generate_nhats_like_functional_data.1012 <- function(n_strata = 20,
                                                    beta_true = beta_true,
                                                    n_psu = 5,
                                                    n_samples = 200) {
  L = dim(beta_true)[2]
  data <- data.frame(
    strata = rep(1:n_strata, each = n_psu * (n_samples / (n_strata * n_psu))),
    PSU = rep(rep(1:n_psu, each = n_samples / (n_strata * n_psu)), n_strata),
    id = 1:n_samples
  )
  data$weight <- runif(n_samples, 0.5, 1.5) * (1 + 5 * data$strata / n_strata)

  # generate exposure X
  data$X <- 3 * data$strata / n_strata + rnorm(n_samples, mean = data$PSU, sd = 2)
  # 0.5 * data$C1 + 0.3 * data$C2 + 0.2 * data$C3 + rnorm(n_samples, 0, 0.5)

  # generate outcomes Y
  Y_obs <- matrix(NA, n_samples, L)
  Y_obs <- data$X %x% t(beta_true[2,]) + rnorm(n_samples, 0, 2)
  Y_obs <- data.frame(Y_obs); colnames(Y_obs) = sub('X', 'Y', colnames(Y_obs))
  data = cbind(data, Y_obs)

  return(data)
}


