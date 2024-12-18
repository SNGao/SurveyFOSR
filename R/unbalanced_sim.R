rm(list = ls())
library(refund)
library(dplyr)

# source('fui.origin/lfosr3s.R')
source('R/wFUI_1.R')
source('R/lforsim_fixed_design.R')

################################################################################
## Set parameters
################################################################################
I <- 100 ## number of subjects
L <- 10 ## dimension of the functional domain
SNR_sigma <- 1 ## signal-to-noise ratio
family <- "gaussian" ## family of longitudinal functional data

## specify true fixed and random effects functions
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 2, L) # time-dependent beta-true trajectory
scenario <- 1 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")


################################################################################
## Do simulations on a local laptop (results in the manuscript are from JHPCE)
################################################################################
nsim <- 30
sim_local <- list() ## store simulation results

for(iter in 1:nsim){

  ################################################################################
  ## Generate simulated data
  ################################################################################
  set.seed(iter)
  data <- lfosrsim_fixed_design(family, I=1000, L=L,
                         beta_true = beta_true, SNR_sigma = SNR_sigma)

  ## pre-process simulated data
  Y_mat <- data[, grepl('Y.', colnames(data))]
  dat.fit <- data.frame(Y = Y_mat,
                        data[,!grepl('Y.', colnames(data))])
  formula.FUI = as.formula(paste0('Y~', 'X'))

  ################################################################################
  ## Implement different estimation methods
  ################################################################################
  ## fit the "Rao-Wu-Yue-Beaumont" bootstrap model
  ptm <- proc.time()
  mod.RWYB <- fui.weight(formula = formula.FUI,
                         analytic = FALSE,
                         argvals = seq(from = 1, to = L, by = 1),
                         parallel = TRUE,
                         weighted = TRUE,
                         boot_type = 'Rao-Wu-Yue-Beaumont',
                         num_boots = 100,
                         nknots_min = 20,
                         num_cores = 10,
                         random = FALSE,
                         data = dat.fit) # run successfully, consider next step for coverage ratio

  # mod.original <- fui.weight(formula = formula.FUI,
  #                        analytic = FALSE,
  #                        argvals = seq(from = 1, to = L, by = 1),
  #                        parallel = TRUE,
  #                        weighted = TRUE,
  #                        boot_type = 'ref_noDesign',
  #                        num_boots = 100,
  #                        nknots_min = 20,
  #                        num_cores = 10,
  #                        random = FALSE,
  #                        data = dat.fit) # run successfully, consider next step for coverage ratio

  lfosr3stime <- (proc.time() - ptm)[3]

  ################################################################################
  ## Organize simulation results
  ################################################################################
  mod.Output = mod.RWYB
  ## lfosr3s
  MISE_lfosr3s <- rowMeans((mod.Output$betaHat - beta_true)^2)
  cover_lfosr3s <- rep(NA, nrow(beta_true))
  cover_lfosr3s_pw <- matrix(FALSE, nrow(beta_true), L)
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(mod.Output$betaHat[p,]+mod.Output$qn[p]*sqrt(diag(mod.Output$betaHat.var[,,p])) > beta_true[p,])
    cover_lower <- which(mod.Output$betaHat[p,]-mod.Output$qn[p]*sqrt(diag(mod.Output$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s[p] <- length(intersect(cover_lower, cover_upper)) == L
    cover_upper_pw <- which(mod.Output$betaHat[p,]+1.96*sqrt(diag(mod.Output$betaHat.var[,,p])) > beta_true[p,])
    cover_lower_pw <- which(mod.Output$betaHat[p,]-1.96*sqrt(diag(mod.Output$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE
  }

  ## results of a single simulation
  sim_result <- list(MISE_lfosr3s = MISE_lfosr3s,
                     cover_lfosr3s = cover_lfosr3s,
                     cover_lfosr3s_pw = cover_lfosr3s_pw,
                     time_lfosr3s = lfosr3stime,
                     I = I, L = L, SNR_sigma = SNR_sigma,
                     family = family, scenario = scenario)

  sim_local[[iter]] <- sim_result
  print(iter)
}


################################################################################
## Obtain MISE, coverage, and computing time
################################################################################
## MISE
MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
colMeans(MISE_lfosr3s)

## coverage
### joint
cover_lfosr3s <- t(lapply(sim_local, '[[', 2) %>% bind_cols())
### pointwise
cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # intercept
  cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # fixed effects
}
### print results
print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))

## time
time_lfosr3s <- lapply(sim_local, '[[', 4) %>% bind_rows()
print(apply(time_lfosr3s, 2, median))
