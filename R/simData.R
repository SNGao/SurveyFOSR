simData.S4 <- function(n_strata = 20,
                       n_psu = 5,
                       n_samples = 200) {
  data <- data.frame(
    strata = rep(1:n_strata, each = n_psu * (n_samples / (n_strata * n_psu))),
    psu = rep(rep(1:n_psu, each = n_samples / (n_strata * n_psu)), n_strata),
    id = 1:n_samples
  )
  data$weight <- runif(n_samples, 0.5, 1.5) * (1 + 5 * data$strata / n_strata)
  data$X <- 2 * data$strata / n_strata + 0.8 * data$strata +
    0.4 * data$psu + rnorm(n_samples, 0, 0.5)

  # generate outcomes Y
  data$Y <- 2 * data$X + rnorm(n_samples, 0, 2)

  return(data)
}


simData.S1.functional <- function(n_strata = 20,
                                  beta_true = beta_true,
                                  n_psu = 5,
                                  n_samples = 200) {
  L = dim(beta_true)[2]
  data <- data.frame(
    strata =
      # rep(1:n_strata, each = n_psu * (n_samples / (n_strata * n_psu))),
    psu = rep(rep(1:n_psu, each = n_samples / (n_strata * n_psu)), n_strata),
    id = 1:n_samples
  )
  data$weight <- runif(n_samples, 0.5, 1.5) * (1 + 5 * data$strata / n_strata)
  data$X <- 2 * data$strata / n_strata + 0.8 * data$strata +
    0.4 * data$psu + rnorm(n_samples, 0, 1)

  # generate outcomes Y (functional)
  Y_obs <- matrix(NA, n_samples, L)
  for (l in 1:L) {
    Y_obs[,l] <- beta_true[1,l] + beta_true[2,l] * data$X + rnorm(n_samples, 0, 2)
  }
  Y_obs <- data.frame(Y_obs)
  colnames(Y_obs) <- paste0("Y", 1:L)
  data <- cbind(data, Y_obs)

  return(data)
}



simData.S2.functional <- function(n_strata = 20,
                                  beta_true = beta_true,
                                  n_psu = 5,
                                  n_samples = 200) {
  # Check the dimension
  L <- dim(beta_true)[2]
  if (nrow(beta_true) < 2 || is.null(L)) {
    stop("beta_true dimensions do not match expected size.")
  }

  # Generate data frame
  data <- data.frame(
    strata = rep(1:n_strata, length.out = n_samples),
    psu = rep(1:n_psu, length.out = n_samples),
    id = 1:n_samples
  )

  # Calculate the weights and the argument X
  data$weight <- runif(n_samples, 0.5, 1.5) * (1 + 5 * (data$strata - 1) / (n_strata - 1))
  data$X <- 2 * data$strata / n_strata + 0.8 * data$strata +
    0.4 * data$psu + rnorm(n_samples, 0, 1)

  # Generate multiple Y result variables
  Y_obs <- matrix(NA, n_samples, L)
  for (l in 1:L) {
    Y_obs[, l] <- beta_true[1, l] + beta_true[2, l] * data$X + rnorm(n_samples, 0, 2)
  }
  Y_obs <- data.frame(Y_obs)
  colnames(Y_obs) <- paste0("Y", 1:L)

  # combine results
  data <- cbind(data, Y_obs)
  return(data)
}
