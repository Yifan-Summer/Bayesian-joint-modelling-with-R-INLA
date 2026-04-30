library(survivalROC)
library(survival)
library(dplyr)
library(ggplot2)

#ROC
ptid_df <- dataL[ , c("cd4cd8_confirm", "ptid")]
n_total <- nrow(Jointmodel$summary.fitted.values)
nL  <- nrow(dataL)
nE1 <- nrow(dataE1)
nE2 <- nrow(dataE2)

idx_E1 <- (2 * nL + 1):(2 * nL + nE1)
idx_E2 <- (2 * nL + nE1 + 1):(2 * nL + nE1 + nE2)
mean_value1  <- Jointmodel$summary.fitted.values[idx_E1, "mean"]
mean_value11 <- Jointmodel$summary.fitted.values[idx_E2, "mean"]
mean_value2  <- Jointmodel1$summary.fitted.values[idx_E1, "mean"]
mean_value22 <- Jointmodel1$summary.fitted.values[idx_E2, "mean"]


quantiles <- quantile(dataE1$Time, probs = c( 0.1, 0.2, 0.3, 0.4, 0.5, 
                                              0.6, 0.7, 0.8, 0.9), na.rm = TRUE)
results <- data.frame(Cutoff = numeric(), FP = numeric(), TP = numeric(), AUC = numeric())
for (cutoff in quantiles) {
  ROC <- survivalROC(
    Stime = dataE1$Time,
    status = dataE1$event,
    marker = -mean_value1,
    predict.time = cutoff,
    method = "KM"
  )
  results <- rbind(results, data.frame(Cutoff = cutoff, FP = ROC$FP, TP = ROC$TP, AUC = ROC$AUC))
}


quantiles <- quantile(dataE1$Time, probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                             0.6, 0.7, 0.8, 0.9), na.rm = TRUE)
results <- data.frame(Cutoff = numeric(), FP = numeric(), TP = numeric(), AUC = numeric())
for (cutoff in quantiles) {
  ROC <- survivalROC(
    Stime = dataE1$Time,
    status = dataE1$event,
    marker = -mean_value2,
    predict.time = cutoff,
    method = "KM"
  )
  results <- rbind(results, data.frame(Cutoff = cutoff, FP = ROC$FP, TP = ROC$TP, AUC = ROC$AUC))
}



quantiles <- quantile(dataE2$Time, probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                             0.6, 0.7, 0.8, 0.9), na.rm = TRUE)
results <- data.frame(Cutoff = numeric(), FP = numeric(), TP = numeric(), AUC = numeric())
for (cutoff in quantiles) {
  ROC <- survivalROC(
    Stime = dataE2$Time,
    status = dataE2$event,
    marker = -mean_value11,
    predict.time = cutoff,
    method = "KM"
  )
  results <- rbind(results, data.frame(Cutoff = cutoff, FP = ROC$FP, TP = ROC$TP, AUC = ROC$AUC))
}



quantiles <- quantile(dataE2$Time, probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                             0.6, 0.7, 0.8, 0.9), na.rm = TRUE)
results <- data.frame(Cutoff = numeric(), FP = numeric(), TP = numeric(), AUC = numeric())
for (cutoff in quantiles) {
  ROC <- survivalROC(
    Stime = dataE2$Time,
    status = dataE2$event,
    marker = -mean_value22,
    predict.time = cutoff,
    method = "KM"
  )
  results <- rbind(results, data.frame(Cutoff = cutoff, FP = ROC$FP, TP = ROC$TP, AUC = ROC$AUC))
}


#iBS

event_time <- dataS$Time
event_type <- dataS$C

times <- as.numeric(
  quantile(event_time, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE)
)


S_fun <- function(u, lambda1, lambda2, alpha1, alpha2) {
  u_pos <- pmax(u, 1e-10)
  exp(-(lambda1 * u_pos^alpha1 + lambda2 * u_pos^alpha2))
}

h1_fun <- function(u, lambda1, alpha1) {
  u_pos <- pmax(u, 1e-10)
  alpha1 * u_pos^(alpha1 - 1) * lambda1
}

h2_fun <- function(u, lambda2, alpha2) {
  u_pos <- pmax(u, 1e-10)
  alpha2 * u_pos^(alpha2 - 1) * lambda2
}

cif1_fun <- function(t, lambda1, lambda2, alpha1, alpha2) {
  if (t <= 0) return(0)
  integrate(
    function(u) {
      S_fun(u, lambda1, lambda2, alpha1, alpha2) *
        h1_fun(u, lambda1, alpha1)
    },
    lower = 0,
    upper = t,
    subdivisions = 1000L,
    rel.tol = 1e-6
  )$value
}

cif2_fun <- function(t, lambda1, lambda2, alpha1, alpha2) {
  if (t <= 0) return(0)
  integrate(
    function(u) {
      S_fun(u, lambda1, lambda2, alpha1, alpha2) *
        h2_fun(u, lambda2, alpha2)
    },
    lower = 0,
    upper = t,
    subdivisions = 1000L,
    rel.tol = 1e-6
  )$value
}


censor_ind <- as.numeric(event_type == 3)
km_cens <- survfit(Surv(time = event_time, event = censor_ind) ~ 1)

G_eval <- function(t, kmfit) {
  sf <- summary(kmfit, times = t, extend = TRUE)$surv
  if (length(sf) == 0 || is.na(sf)) return(1)
  sf
}

brier_ipcw <- function(time_point, event_time, event_type, pred_prob, cause, kmfit) {
  n <- length(event_time)
  val <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    Ti <- event_time[i]
    Di <- event_type[i]
    
    y_it <- as.numeric(Ti <= time_point && Di == cause)
    
    if (Ti <= time_point && Di != 3) {
      Gt <- G_eval(Ti, kmfit)
      w <- ifelse(Gt > 0, 1 / Gt, NA_real_)
    } else if (Ti > time_point) {
      Gt <- G_eval(time_point, kmfit)
      w <- ifelse(Gt > 0, 1 / Gt, NA_real_)
    } else {
      w <- NA_real_
    }
    
    val[i] <- w * (y_it - pred_prob[i])^2
  }
  
  mean(val, na.rm = TRUE)
}

ibs_trapz <- function(time_vec, bs_vec) {
  o <- order(time_vec)
  t <- time_vec[o]
  b <- bs_vec[o]
  area <- sum(diff(t) * (head(b, -1) + tail(b, -1)) / 2)
  area / (max(t) - min(t))
}


calc_brier_for_model <- function(model_obj, model_name, event_time, event_type, times, km_cens) {
  # 40779:44956 -> cause 1
  # 44957:49134 -> cause 2
  eta1 <- model_obj$summary.linear.predictor[40779:44956, "mean"]
  eta2 <- model_obj$summary.linear.predictor[44957:49134, "mean"]
  
  alpha1 <- model_obj$summary.hyperpar["alpha parameter for weibullsurv[3]", "mean"]
  alpha2 <- model_obj$summary.hyperpar["alpha parameter for weibullsurv[4]", "mean"]
  
  lambda1 <- exp(eta1)
  lambda2 <- exp(eta2)
  
  surv_dat <- data.frame(
    id = seq_along(event_time),
    time = event_time,
    event = event_type,
    eta1 = eta1,
    eta2 = eta2,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
  
  pred_cif1 <- sapply(times, function(tt) {
    sapply(seq_len(nrow(surv_dat)), function(i) {
      cif1_fun(
        t = tt,
        lambda1 = surv_dat$lambda1[i],
        lambda2 = surv_dat$lambda2[i],
        alpha1 = alpha1,
        alpha2 = alpha2
      )
    })
  })
  
  pred_cif2 <- sapply(times, function(tt) {
    sapply(seq_len(nrow(surv_dat)), function(i) {
      cif2_fun(
        t = tt,
        lambda1 = surv_dat$lambda1[i],
        lambda2 = surv_dat$lambda2[i],
        alpha1 = alpha1,
        alpha2 = alpha2
      )
    })
  })
  
  pred_cif1 <- as.matrix(pred_cif1)
  pred_cif2 <- as.matrix(pred_cif2)
  
  bs1 <- sapply(seq_along(times), function(j) {
    brier_ipcw(
      time_point = times[j],
      event_time = surv_dat$time,
      event_type = surv_dat$event,
      pred_prob = pred_cif1[, j],
      cause = 1,
      kmfit = km_cens
    )
  })
  
  bs2 <- sapply(seq_along(times), function(j) {
    brier_ipcw(
      time_point = times[j],
      event_time = surv_dat$time,
      event_type = surv_dat$event,
      pred_prob = pred_cif2[, j],
      cause = 2,
      kmfit = km_cens
    )
  })
  
  ibs1 <- ibs_trapz(times, bs1)
  ibs2 <- ibs_trapz(times, bs2)
  
  list(
    model_name = model_name,
    alpha1 = alpha1,
    alpha2 = alpha2,
    pred_cif1 = pred_cif1,
    pred_cif2 = pred_cif2,
    brier = data.frame(
      Time = times,
      Brier_AIDS_onset = bs1,
      Brier_preAIDS_death = bs2,
      Model = model_name
    ),
    ibs = data.frame(
      Model = model_name,
      Outcome = c("AIDS onset", "pre-AIDS death"),
      IBS = c(ibs1, ibs2)
    )
  )
}

res_gp <- calc_brier_for_model(
  model_obj = Jointmodel,
  model_name = "GP",
  event_time = event_time,
  event_type = event_type,
  times = times,
  km_cens = km_cens
)

res_gau <- calc_brier_for_model(
  model_obj = Jointmodel1,
  model_name = "Gaussian",
  event_time = event_time,
  event_type = event_type,
  times = times,
  km_cens = km_cens
)


cat("Range of GP CIF1:", range(res_gp$pred_cif1, na.rm = TRUE), "\n")
cat("Range of GP CIF2:", range(res_gp$pred_cif2, na.rm = TRUE), "\n")
cat("Range of Gaussian CIF1:", range(res_gau$pred_cif1, na.rm = TRUE), "\n")
cat("Range of Gaussian CIF2:", range(res_gau$pred_cif2, na.rm = TRUE), "\n")


brier_compare <- bind_rows(res_gp$brier, res_gau$brier)
ibs_compare <- bind_rows(res_gp$ibs, res_gau$ibs)

print(brier_compare)
print(ibs_compare)

