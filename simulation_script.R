# simulation: multi-reference calibration with three domains and three-level factor outcome
#
# setup:
#   sp1 contains x1 and x2
#   sp2 contains x3 and x4
#   x1 is a three-level domain variable
#   y is a three-level categorical outcome
#   final evaluation table = (3 x1 domains) by (3 y levels) by (2 precali cases)
#
# evaluation metrics follow the original script:
#   RB(%), v(*100), MSE(*100), SER, CP

#install.packages(nonprobsampling)

library(nonprobsampling)
library(sampling)
library(survey)
library(knitr)

options(scipen = 999)
options(digits = 4)

set.seed(123456)

# -----------------------------------------------------------------------------
# settings
# -----------------------------------------------------------------------------

N <- 500000
n_sim <- 100

np1 <- 12500
np2 <- 25000
nc  <- 2500

dom_levels <- c("A", "B", "C")
y_levels   <- c("Low", "Medium", "High")
precali_cases <- c("No precalibration", "Precalibration")


# change beta_x4 from -0.10 to -0.27 for extreme population
beta_x4 <- -0.10

# -----------------------------------------------------------------------------
# finite population generation
# -----------------------------------------------------------------------------

v2 <- runif(N, min = 0, max = 2)
v3 <- rexp(N, rate = 1)
v4 <- rchisq(N, df = 4)

# x1: three-level categorical variable
x1 <- factor(
  sample(dom_levels, size = N, replace = TRUE, prob = c(1 / 2, 1 / 4, 1 / 4)),
  levels = dom_levels
)

# indicators for data generation
I1_B <- as.integer(x1 == "B")
I1_C <- as.integer(x1 == "C")

# covariates
x2 <- v2 + 0.3 * I1_B + 0.3 * I1_C
x3 <- v3 + 0.2 * (I1_B + I1_C + x2)
x4 <- v4 + 0.1 * (I1_B + I1_C + x2 + x3)

# -----------------------------------------------------------------------------
# three-level factor outcome
# -----------------------------------------------------------------------------

# eta_y is linear predictor, observed outcome y should be a three-level factor
eta_y  <- -0.5 * I1_B - I1_C - x2 + x3 + x4

# add random noise
y_star <- eta_y + rnorm(N)

# cut points: 35% low, 35% medium, 30% high
y_cut <- quantile(y_star, probs = c(0.35, 0.70), names = FALSE)

# categorical y generation
y <- cut(
  y_star,
  breaks = c(-Inf, y_cut, Inf),
  labels = y_levels,
  right = TRUE
)
y <- factor(y, levels = y_levels)

# check association between covariates and high-outcome category (results are as desired)
print(cor(I1_B, as.integer(y == "High")))
print(cor(I1_C, as.integer(y == "High")))
print(cor(x2,   as.integer(y == "High")))
print(cor(x3,   as.integer(y == "High")))
print(cor(x4,   as.integer(y == "High")))

# true domain-by-outcome prevalence
truth_mat <- prop.table(table(x1 , y), margin = 1)
print(truth_mat)

# -----------------------------------------------------------------------------
# reference-sample inclusion probabilities
# -----------------------------------------------------------------------------

# sp1: Poisson sampling
# since y is now a factor, use y_star for numerical calculation
a <- x3 + 0.03 * y_star
amin <- min(a)
amax <- max(a)
cnst_sp1 <- (amax - 20 * amin) / 19
q <- cnst_sp1 + a

pi_p1 <- np1 * q / sum(q)
di1   <- 1 / pi_p1

stopifnot(all(pi_p1 > 0 & pi_p1 < 1))

# sp2: randomized systematic PPS
const2 <- 0.05
z <- const2 + x2

pi_p2 <- np2 * z / sum(z)
di2   <- 1 / pi_p2

stopifnot(all(pi_p2 > 0 & pi_p2 < 1))

# -----------------------------------------------------------------------------
# nonprobability-sample participation probabilities
# -----------------------------------------------------------------------------

eta <- 0.18 * I1_B + 0.18 * I1_C +
  0.18 * x2 - 0.27 * x3 + beta_x4 * x4

exp_beta0 <- nc / sum(exp(eta))
beta0 <- log(exp_beta0)

pi_c <- exp(eta + beta0)

stopifnot(all(pi_c > 0 & pi_c < 1))
cat("Expected nonprobability sample size:", sum(pi_c), "\n")
cat("Range of true weights", range(1/pi_c), "\n")

# -----------------------------------------------------------------------------
# finite population data frame
# -----------------------------------------------------------------------------

fp <- data.frame(
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x4 = x4,
  y  = y,
  pi_sp1 = pi_p1,
  wt_sp1 = di1,
  pi_sp2 = pi_p2,
  wt_sp2 = di2,
  wt_sc  = 1 / pi_c
)


# -----------------------------------------------------------------------------
# helper functions
# -----------------------------------------------------------------------------

# evaluate simulation results for one domain, outcome level, and precali case
# inputs:
#   est: estimated means across simulation rounds
#   se: estimated standard errors across simulation rounds
#   truth: true finite-population value
# output:
#   a named vector with RB(%), v(*100), MSE(*100), SER, and CP
evaluation <- function(est, se, truth) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]
  se  <- se[ok]
  n_valid <- length(est)

  if (n_valid < 2 || !is.finite(truth) || truth == 0) {
    return(c(
      `RB(%)` = NA_real_,
      `v(*100)` = NA_real_,
      `MSE(*100)` = NA_real_,
      SER = NA_real_,
      CP = NA_real_
    ))
  }

  c(
    `RB(%)` = mean((est - truth) / truth) * 100,
    `v(*100)` = mean(se^2) * 100,
    `MSE(*100)` = mean((est - truth)^2) * 100,
    SER = mean(se) / sd(est),
    CP = mean(est - 1.96 * se <= truth & est + 1.96 * se >= truth)
  )
}

# build the final evaluation table across all domains, outcome levels, and precalibration cases
# inputs:
#   result_array: simulation estimates stored by iter, precali_case, domain, outcome, and stat
#   truth_mat: true domain-by-outcome prevalences
#   digits: number of decimal places for rounded metrics
# output:
#   a data frame with domain, outcome, precali_case, and evaluation metrics
make_final_table <- function(result_array, truth_mat, digits = 4) {
  grid <- expand.grid(
    domain = dimnames(result_array)$domain,
    outcome = dimnames(result_array)$outcome,
    precali_case = dimnames(result_array)$precali_case
  )

  metric_mat <- matrix(
    NA_real_,
    nrow = nrow(grid),
    ncol = 5,
    dimnames = list(NULL, c("RB(%)", "v(*100)", "MSE(*100)", "SER", "CP"))
  )

  for (r in seq_len(nrow(grid))) {
    dom <- grid$domain[r]
    out <- grid$outcome[r]
    case <- grid$precali_case[r]

    metric_mat[r, ] <- evaluation(
      est = result_array[, case, dom, out, "mean"],
      se = result_array[, case, dom, out, "se"],
      truth = truth_mat[dom, out]
    )
  }

  final <- cbind(grid, as.data.frame(metric_mat, check.names = FALSE))

  num_cols <- c("RB(%)", "v(*100)", "MSE(*100)", "SER", "CP")
  final[num_cols] <- lapply(final[num_cols], round, digits = digits)

  final
}

# -----------------------------------------------------------------------------
# storage array: store the adjusted mean and SE from each simulation round
# -----------------------------------------------------------------------------
# dimensions:
#   iter by precali_case by domain by outcome by stat
result_multi <- array(
  NA_real_,
  dim = c(n_sim, length(precali_cases), length(dom_levels), length(y_levels), 2),
  dimnames = list(
    iter = seq_len(n_sim),
    precali_case = precali_cases,
    domain = dom_levels,
    outcome = y_levels,
    stat = c("mean", "se")
  )
)

# -----------------------------------------------------------------------------
# main simulation loop
# -----------------------------------------------------------------------------

for (i in seq_len(n_sim)) {

  # ----------------------------
  # nonprobability sample sc
  # ----------------------------
  sc <- fp[rbinom(N, size = 1, prob = pi_c) == 1,
           c("x1", "x2", "x3", "x4", "y", "wt_sc")]

  sc$x1 <- factor(sc$x1, levels = dom_levels)
  sc$y  <- factor(sc$y,  levels = y_levels)

  # ----------------------------
  # probability reference sample sp1: x1, x2 only
  # ----------------------------
  sp1 <- fp[rbinom(N, size = 1, prob = pi_p1) == 1,
            c("x1", "x2", "wt_sp1")]

  sp1$x1 <- factor(sp1$x1, levels = dom_levels)
  sp1$pi_sp1 <- 1 / sp1$wt_sp1

  des_sp1 <- survey::svydesign(
    ids   = ~1,
    probs = ~pi_sp1,
    data  = sp1
  )

  # ----------------------------
  # probability reference sample sp2: x3, x4 only
  # ----------------------------
  s2 <- sampling::UPrandomsystematic(fp$pi_sp2)

  sp2 <- fp[s2 == 1,
            c("x3", "x4", "wt_sp2", "pi_sp2")]

  des_sp2 <- survey::svydesign(ids = ~1, fpc = ~pi_sp2, data = sp2, pps = "brewer")

  # ----------------------------
  # multi-reference calibration: no precalibration
  # ----------------------------
  fit_multi_no_pc <- tryCatch(
    est_pw(
      data = list(sc, des_sp1, des_sp2),
      method = "multi",
      precali = FALSE,
      p_formula = list(
        ~ x1 + x2,
        ~ x3 + x4
      )
    ),
    error = function(e) {
      message("Iteration ", i, ", no precalibration failed: ", conditionMessage(e))
      NULL
    }
  )

  # ----------------------------
  # multi-reference calibration: with precalibration
  # ----------------------------
  fit_multi_pc <- tryCatch(
    est_pw(
      data = list(sc, des_sp1, des_sp2),
      method = "multi",
      precali = TRUE,
      p_formula = list(
        ~ x1 + x2,
        ~ x3 + x4
      )
    ),
    error = function(e) {
      message("Iteration ", i, ", precalibration failed: ", conditionMessage(e))
      NULL
    }
  )

  fits <- list(
    "No precalibration" = fit_multi_no_pc,
    "Precalibration" = fit_multi_pc
  )

  for (case in names(fits)) {
    if (is.null(fits[[case]])) next

    est_df <- tryCatch(
      pwmean(fits[[case]], y = "y", zcol = "x1")$estimates,
      error = function(e) {
        message("Iteration ", i, ", pwmean ", case, " failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(est_df)) next

    for (dom in dom_levels) {
      for (lv in y_levels) {
        row_i <- est_df$domain == paste("x1 =", dom) &
          est_df$category == paste("y =", lv)

        if (sum(row_i) == 1) {
          result_multi[i, case, dom, lv, "mean"] <- est_df$adjusted_mean[row_i]
          result_multi[i, case, dom, lv, "se"] <- est_df$adjusted_se[row_i]
        }
      }
    }
  }

  if (i %% 100 == 0) cat(i, " ")
}

cat("\nSimulation finished.\n")

# -----------------------------------------------------------------------------
# final evaluation table: 3 domains * 3 outcome levels * 2 precalibration cases = 18 rows
# -----------------------------------------------------------------------------

evaluation_table <- make_final_table(result_multi, truth_mat, digits = 4)

print(knitr::kable(evaluation_table, align = "c"))

