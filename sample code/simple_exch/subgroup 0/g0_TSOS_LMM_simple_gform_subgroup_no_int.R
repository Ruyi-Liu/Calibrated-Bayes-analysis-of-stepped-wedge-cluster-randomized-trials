library(MCMCglmm)
library(dplyr) 
library(openxlsx) 

dat <- readRDS("/Users/ruyi/Desktop/TSOS/data_analysis/TSOS_imputed_once.rds")

thin_rate <- 1
mctotal <- 30000
mcburnin <- 20000
num_new_dataset <- 100

set.seed(20260913)


## -- resample whole clusters (all periods/rows), then relabel clusters 1..I
bootstrap_clusters_idx <- function(dat, cluster_var = "cluster",
                                   period_var = "period", period_f_var = "period_f") {
  cl <- as.integer(factor(dat[[cluster_var]]))
  ucl <- sort(unique(cl))
  samp <- sample(ucl, length(ucl), replace = TRUE)
  
  idx_vec <- integer(0)
  out <- vector("list", length(samp))
  for (new_id in seq_along(samp)) {
    rows <- which(cl == samp[new_id])              # original-row indices for this sampled cluster
    d <- dat[rows, , drop = FALSE]
    d[[cluster_var]] <- new_id                     # relabel to 1..I
    out[[new_id]] <- d
    idx_vec <- c(idx_vec, rows)                    # accumulate mapping to original rows IN BOOT ORDER
  }
  res <- do.call(rbind, out)
  res[[period_f_var]] <- factor(res[[period_var]], levels = levels(dat[[period_f_var]]))
  res$cp <- interaction(res[[cluster_var]], res[[period_f_var]], drop = TRUE)
  list(dat = res, index = idx_vec)
}



## build counterfactual design matrix for a constant treatment value
build_cf_design <- function(dat, fixed_formula, treat_const,
                            treat_var = "treat", period_f_var = "period_f") {
  nd <- dat
  nd[[treat_var]] <- as.numeric(treat_const)
  P <- model.matrix(reformulate(period_f_var, intercept = FALSE), data = nd)
  Jloc <- ncol(P)
  for (m in 2:(Jloc - 1)) nd[[paste0("treat_P", m)]] <- nd[[treat_var]] * P[, m]
  model.matrix(fixed_formula, data = nd)
}



## Model-robust i-ATE & c-ATE from observed Ybar_ij and predicted E[Ybar_ij|Z=z,Xij,Nij]
## Y0, Y1: n x B matrices of individual-level predictions under z=0 and z=1
ATEs_from_predictions_gform_subgroup <- function(dat, Y0, Y1,
                                                 subgroup_var = "x1",
                                                 cluster_var  = "cluster",
                                                 period_f_var = "period_f") {
  
  
  stopifnot(nrow(Y0) == nrow(dat), nrow(Y1) == nrow(dat), ncol(Y0) == ncol(Y1))
  
  S_r <- as.integer(dat[[subgroup_var]] == 0)
  
  
  i_vec <- as.integer(dat[[cluster_var]])
  j_vec <- as.integer(dat[[period_f_var]])
  Jmax  <- max(j_vec)
  
  cp_df  <- unique(data.frame(i = i_vec, j = j_vec))
  cp_df  <- cp_df[order(cp_df$i, cp_df$j), ]
  row_key <- paste(i_vec, j_vec, sep = "_")
  cp_key  <- paste(cp_df$i, cp_df$j, sep = "_")
  cp_index <- factor(row_key, levels = cp_key)
  
  cp_i   <- cp_df$i
  cp_j   <- cp_df$j
  Nij_cp <- as.numeric(table(cp_index))
  
  
  keep_mid <- (cp_j >= 2 & cp_j <= (Jmax - 1))
  
  ## subgroup cp sizes: N^{(r)}_{ij}
  N_cp <- as.numeric(rowsum(S_r, cp_index))
  
  
  one_draw <- function(y0, y1) {
    
    cp_y0 <- as.numeric(rowsum(y0 * S_r, cp_index))
    cp_y1 <- as.numeric(rowsum(y1 * S_r, cp_index))
    
    denom_i <- sum(N_cp[keep_mid])
    i_ATE <- if (denom_i > 0)
      (sum(cp_y1[keep_mid]) - sum(cp_y0[keep_mid])) / denom_i else NA_real_
    
    ###################
    
    # subgroup c-ATE with mid-period normalization by subgroup sizes
    sum_y1_by_i <- tapply(cp_y1[keep_mid], cp_i[keep_mid], sum)
    sum_y0_by_i <- tapply(cp_y0[keep_mid], cp_i[keep_mid], sum)
    sum_Nr_by_i <- tapply(N_cp [keep_mid], cp_i[keep_mid], sum)  # subgroup counts
    
    nonzero <- sum_Nr_by_i > 0
    c_ATE <- mean(
      (sum_y1_by_i[nonzero] / sum_Nr_by_i[nonzero]) -
        (sum_y0_by_i[nonzero] / sum_Nr_by_i[nonzero])
    )
    
    c(i_ATE = i_ATE, c_ATE = c_ATE)
  }
  
  B <- ncol(Y0)
  out <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("i_ATE","c_ATE")))
  for (b in seq_len(B)) out[b, ] <- one_draw(Y0[, b], Y1[, b])
  as.data.frame(out)
}



### reuse code from simulation

## ---------- Priors (weakly-informative) ----------
# Build period dummies (P1..PJ) and explicit treat-by-period (>=2 but < J)
P <- model.matrix(~ -1 + period_f, data = dat)
Jloc <- ncol(P)
for (m in 2:(Jloc - 1)) {
  dat[[paste0("treat_P", m)]] <- dat$treat * P[, m]
}

fixed_formula <- as.formula(
  paste(
    "y ~ -1 + period_f +",
    paste(paste0("treat_P", 2:(Jloc-1)), collapse = " + "),
    "+ x1 + x2 + x3 + N_ij")
)


p <- ncol(model.matrix(fixed_formula, data = dat))



prior <- list(
  B = list(mu = rep(0, p), V = diag(1e3, p)),
  G = list(
    G1 = list(V = matrix(1,1,1), nu = 0.002)   # cluster RE (α_i)
  ),
  R = list(V = 1, nu = 2)
)

fit <- MCMCglmm(
  fixed   = fixed_formula,
  random  = ~ cluster,
  family  = "gaussian",
  data    = dat,
  prior   = prior,
  nitt    = mctotal, burnin = mcburnin, thin = thin_rate,
  verbose = FALSE
)



## posterior draws for fixed effects (to make predictions once)
beta_names <- colnames(model.matrix(fixed_formula, data = dat))
Beta_draws <- fit$Sol[, beta_names, drop = FALSE]   # draws x p

## build design matrices & predict potential outcomes
X0 <- build_cf_design(dat, fixed_formula, treat_const = 0)
X1 <- build_cf_design(dat, fixed_formula, treat_const = 1)

Y0 <- X0 %*% t(Beta_draws)   # n x B
Y1 <- X1 %*% t(Beta_draws)   # n x B

## ATE draws from *predicted* outcomes (model-agnostic)
posterior_ATEs <- ATEs_from_predictions_gform_subgroup(dat, Y0, Y1)


est_c <- mean(posterior_ATEs$c_ATE, na.rm = TRUE)
est_i <- mean(posterior_ATEs$i_ATE, na.rm = TRUE)

upper_CI_cATE <- quantile(posterior_ATEs$c_ATE, 0.975, na.rm = TRUE)
lower_CI_cATE <- quantile(posterior_ATEs$c_ATE, 0.025, na.rm = TRUE)
upper_CI_iATE <- quantile(posterior_ATEs$i_ATE, 0.975, na.rm = TRUE)
lower_CI_iATE <- quantile(posterior_ATEs$i_ATE, 0.025, na.rm = TRUE)


## -------------Convergence Check-------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(posterior_ATEs$i_ATE, type = "l",
     main = "Traceplot: i-ATE",
     xlab = "MCMC iteration", ylab = "i-ATE")

plot(posterior_ATEs$c_ATE, type = "l",
     main = "Traceplot: c-ATE",
     xlab = "MCMC iteration", ylab = "c-ATE")

par(mfrow = c(1, 1))
coda::geweke.diag(posterior_ATEs$i_ATE)
coda::geweke.diag(posterior_ATEs$c_ATE)

print(mean(posterior_ATEs$c_ATE))
print(mean(posterior_ATEs$i_ATE))
## --------------------------------------------------

## ---------- Corrected variance via cluster bootstrap (unchanged logic) ----------
B <- ncol(Y0)
Variance_matrix_delta_C <- matrix(NA_real_, nrow = num_new_dataset, ncol = B)
Variance_matrix_delta_I <- matrix(NA_real_, nrow = num_new_dataset, ncol = B)
Delta_C_list <- as.numeric(posterior_ATEs$c_ATE)
Delta_I_list <- as.numeric(posterior_ATEs$i_ATE)

for (m in 1:num_new_dataset) {
  boot <- bootstrap_clusters_idx(dat)
  idx_m <- boot$index
  Y0_boot <- Y0[idx_m, , drop = FALSE]
  Y1_boot <- Y1[idx_m, , drop = FALSE]
  deltas_m <- ATEs_from_predictions_gform_subgroup(boot$dat, Y0_boot, Y1_boot)
  Variance_matrix_delta_C[m, ] <- deltas_m$c_ATE
  Variance_matrix_delta_I[m, ] <- deltas_m$i_ATE
}




row_means_E_deltaC <- apply(Variance_matrix_delta_C, 1, mean)
row_means_E_deltaI <- apply(Variance_matrix_delta_I, 1, mean)

Var_deltaC_corrected <- var(row_means_E_deltaC) + var(Delta_C_list)
Var_deltaI_corrected <- var(row_means_E_deltaI) + var(Delta_I_list)

SE_c_corrected <- sqrt(Var_deltaC_corrected)
SE_i_corrected <- sqrt(Var_deltaI_corrected)

upper_CI_cATE_corrected <- est_c + 1.96 * SE_c_corrected
lower_CI_cATE_corrected <- est_c - 1.96 * SE_c_corrected
upper_CI_iATE_corrected <- est_i + 1.96 * SE_i_corrected
lower_CI_iATE_corrected <- est_i - 1.96 * SE_i_corrected

### output results:

results_summary <- data.frame(
  Estimand = c("c-ATE", "i-ATE"),
  
  Estimate = c(est_c, est_i),
  
  CI_lower_quant = c(lower_CI_cATE, lower_CI_iATE),
  CI_upper_quant = c(upper_CI_cATE, upper_CI_iATE),
  
  CI_lower_corrected = c(lower_CI_cATE_corrected, lower_CI_iATE_corrected),
  CI_upper_corrected = c(upper_CI_cATE_corrected, upper_CI_iATE_corrected)
)







hist(posterior_ATEs$c_ATE, breaks = 30, freq = FALSE,
     main = "Posterior c-ATE with Normal Overlay", xlab = "c-ATE")
curve(dnorm(x,
            mean = mean(posterior_ATEs$c_ATE),
            sd   = sd(posterior_ATEs$c_ATE)),
      col = "red", lwd = 2, add = TRUE)


library(openxlsx)

# target folder
out_dir <- "/Users/ruyi/Desktop/TSOS/data_analysis/simple_exch/subgroup 0/results"
out_file <- file.path(out_dir, "TSOS_LMM_simple_gform_subgroup_no_int_results.xlsx")

# write Excel
write.xlsx(
  transform(results_summary,
            Estimate = round(Estimate, 3),
            CI_lower_quant = round(CI_lower_quant, 3),
            CI_upper_quant = round(CI_upper_quant, 3),
            CI_lower_corrected = round(CI_lower_corrected, 3),
            CI_upper_corrected = round(CI_upper_corrected, 3)),
  file = out_file,
  rowNames = FALSE
)


print(
  transform(results_summary,
            Estimate = round(Estimate, 3),
            CI_lower_quant = round(CI_lower_quant, 3),
            CI_upper_quant = round(CI_upper_quant, 3),
            CI_lower_corrected = round(CI_lower_corrected, 3),
            CI_upper_corrected = round(CI_upper_corrected, 3)),
  row.names = FALSE
)
