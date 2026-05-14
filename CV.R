library(readxl)
library(joineR)
library(INLA)
library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(POT)
library(writexl)
library(readr)
library(parallel)
library(sf)
library(survival)
library(timeROC)

dataL <- read_xlsx("final.xlsx")
dataS <-  read_xlsx("final.xlsx")
map <- read_sf("map.shp")
nb <- poly2nb(map)
head(nb)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

dataL <- dataL %>%
  group_by(ptid) %>% 
  mutate(Time = ifelse(C == 3, max(cd4cd8_confirm, na.rm = TRUE), Time)) %>%  
  ungroup()

matched_dataL <- dataL %>%
  inner_join(dataS %>% select(ptid, C), by = c("ptid", "C")) %>%  
  select(ptid, C, Time) %>%  
  distinct(ptid, .keep_all = TRUE)  

dataS <- dataS %>%
  left_join(matched_dataL, by = c("ptid", "C"), suffix = c(".x", ".y")) %>%  
  mutate(Time = ifelse(C == 3, Time.y, Time.x)) %>%  
  select(-Time.y)  

dataL$logcd4 <- log(dataL$cd4)
dataL$CD4min <- max(dataL$logcd4) - dataL$logcd4
CD4min_90th_percentile <- quantile(dataL$CD4min, 0.90)
dataL <- dataL[dataL$CD4min > CD4min_90th_percentile, ]
ptid_counts <- table(dataL$ptid)
ptids_to_keep <- names(ptid_counts[ptid_counts >= 1])
dataL <- dataL[dataL$ptid %in% ptids_to_keep, ]
existing_ptids <- dataL$ptid
dataS <- dataS[dataS$ptid %in% existing_ptids, ]
dataL$CD4_9 <- dataL$CD4min - CD4min_90th_percentile


make_ptid_folds <- function(dataS, K = 5, seed = 2026) {
  set.seed(seed)
  id_df <- dataS %>%
    distinct(ptid, C)
  id_df$fold <- NA_integer_
  for (cc in sort(unique(id_df$C))) {
    idx <- which(id_df$C == cc)
    id_df$fold[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
  }
  folds <- split(id_df$ptid, id_df$fold)
  return(folds)
}

K_cv <- 5

folds <- make_ptid_folds(
  dataS = dataS,
  K = K_cv,
  seed = 2026
)

fold_check <- data.frame(
  fold = seq_along(folds),
  n = sapply(folds, length),
  AIDS_onset = sapply(folds, function(ids) sum(dataS$ptid %in% ids & dataS$C == 1)),
  Pre_AIDS_death = sapply(folds, function(ids) sum(dataS$ptid %in% ids & dataS$C == 2)),
  censored_or_other = sapply(folds, function(ids) sum(dataS$ptid %in% ids & !(dataS$C %in% c(1, 2))))
)

print(fold_check)


fit_gp_one_fold <- function(k, dataL, dataS, folds, g) {
  cat("\n==============================\n")
  cat("Running GP fold", k, "of", length(folds), "\n")
  cat("==============================\n")
  valid_ids <- folds[[k]]
  valid_pos <- which(dataS$ptid %in% valid_ids)
  train_pos <- setdiff(seq_len(nrow(dataS)), valid_pos)
  
  dataE1 <- dataS
  dataE1$event <- ifelse(dataE1$C == 1, 1, 0)
  
  dataE2 <- dataS
  dataE2$event <- ifelse(dataE2$C == 2, 1, 0)
  
  event1_full <- dataE1$event
  event2_full <- dataE2$event
  
  TimeS_model <- as.numeric(dataS$Time)
  
  event1_mask <- event1_full
  event2_mask <- event2_full
  time1_mask <- TimeS_model
  time2_mask <- TimeS_model
  
  event1_mask[valid_pos] <- NA
  event2_mask[valid_pos] <- NA
  time1_mask[valid_pos] <- NA
  time2_mask[valid_pos] <- NA
  
  C1r1_id <- dataS$ptid
  C2r1_id <- dataS$ptid
  
  C1r1_id[valid_pos] <- NA
  C2r1_id[valid_pos] <- NA
  
  nL <- nrow(dataL)
  nS <- nrow(dataS)
  
  fixed.eff <- data.frame(
    mu = as.factor(c(rep(1, nL), rep(1, nL), rep(2, nS), rep(3, nS))),
    age_confirmL = c(dataL$age_confirm, dataL$age_confirm, rep(0, nS), rep(0, nS)),
    age_confirm1 = c(rep(0, nL), rep(0, nL), dataS$age_confirm, rep(0, nS)),
    age_confirm2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$age_confirm),
    sexL = c(dataL$sex, dataL$sex, rep(0, nS), rep(0, nS)),
    sex1 = c(rep(0, nL), rep(0, nL), dataS$sex, rep(0, nS)),
    sex2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$sex),
    marriageL = c(dataL$marriage, dataL$marriage, rep(0, nS), rep(0, nS)),
    marriage1 = c(rep(0, nL), rep(0, nL), dataS$marriage, rep(0, nS)),
    marriage2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$marriage),
    routeL = c(dataL$route, dataL$route, rep(0, nS), rep(0, nS)),
    route1 = c(rep(0, nL), rep(0, nL), dataS$route, rep(0, nS)),
    route2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$route),
    confirm_artL = c(dataL$confirm_art, dataL$confirm_art, rep(0, nS), rep(0, nS)),
    confirm_art1 = c(rep(0, nL), rep(0, nL), dataS$confirm_art, rep(0, nS)),
    confirm_art2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$confirm_art),
    tL = c(dataL$t, dataL$t, rep(0, nS), rep(0, nS)),
    jie1 = c(rep(0, nL), rep(0, nL), dataS$jie, rep(0, nS)),
    jie2 = c(rep(0, nL), rep(0, nL), rep(0, nS), dataS$jie)
  )
  
  random.eff <- list(
    timeL = c(dataL$Time, dataL$Time, rep(NA, nS), rep(NA, nS)),
    idareaS1 = c(rep(NA, nL), rep(NA, nL), dataS$idarea, rep(NA, nS)),
    idareaS2 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$idarea),
    linpredL = c(rep(NA, nL), dataL$ptid, rep(NA, nS), rep(NA, nS)),
    linpredL2 = c(rep(NA, nL), rep(-1, nL), rep(NA, nS), rep(NA, nS)),
    beta1 = c(rep(NA, nL), rep(NA, nL), dataS$ptid, rep(NA, nS)),
    beta2 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$ptid),
    Lr1 = c(dataL$ptid, dataL$ptid, rep(NA, nS), rep(NA, nS)),
    C1r1 = c(rep(NA, nL), rep(NA, nL), C1r1_id, rep(NA, nS)),
    C2r1 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), C2r1_id)
  )
  
  jointdata <- c(fixed.eff, random.eff)
  y.long <- c(dataL$CD4_9, rep(NA, nL), rep(NA, nS), rep(NA, nS))
  y.eta <- c(rep(NA, nL), rep(0, nL), rep(NA, nS), rep(NA, nS))
  
  y.survC1 <- inla.surv(
    time = c(rep(NA, nL), rep(NA, nL), time1_mask, rep(NA, nS)),
    event = c(rep(NA, nL), rep(NA, nL), event1_mask, rep(NA, nS)))
  
  y.survC2 <- inla.surv(
    time = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), time2_mask),
    event = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), event2_mask))
  
  jointdata$Y <- list(y.long, y.eta, y.survC1, y.survC2)
  
  formula.model <- Y ~
    sex1 + sex2 + sexL +
    marriage1 + marriage2 + marriageL +
    route1 + route2 + routeL +
    age_confirm1 + age_confirm2 + age_confirmL +
    confirm_art1 + confirm_art2 + confirm_artL +
    tL + jie1 + jie2 +
    f(inla.group(timeL, n = 50),
      model = "rw2",
      scale.model = TRUE,
      hyper = list(
        prec = list(prior = "pc.prec", param = c(1, 0.01)))
    ) +
    f(idareaS1,
      model = "bym2",
      graph = g,
      hyper = list(
        prec = list(prior = "pc.prec", param = c(0.5 / 0.31, 0.01)),
        phi = list(prior = "pc", param = c(0.1, 4 / 5)))
    ) +
    f(idareaS2,
      copy = "idareaS1",
      hyper = list(beta = list(fixed = FALSE))
    ) +
    f(linpredL,
      linpredL2,
      model = "iid",
      hyper = list(prec = list(initial = -6, fixed = TRUE))
    ) +
    f(Lr1, model = "iid") +
    f(C1r1, model = "iid") +
    f(C2r1, model = "iid") +
    f(beta1,
      copy = "linpredL",
      hyper = list(beta = list(fixed = FALSE))
    ) +
    f(beta2,
      copy = "linpredL",
      hyper = list(beta = list(fixed = FALSE))
    )
  

  fit <- tryCatch(
    inla(
      formula.model,
      family = c("gp", "gaussian", "weibullsurv", "weibullsurv"),
      control.family = list(
        list(control.link = list(model = "quantile", quantile = 0.6)),
        list(),
        list(),
        list()
      ),
      control.predictor = list(compute = TRUE),
      data = jointdata,
      verbose = TRUE
    ),
    error = function(e) {
      message("Fold ", k, " failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  

  idx_surv1 <- (2 * nL + 1):(2 * nL + nS)
  idx_surv2 <- (2 * nL + nS + 1):(2 * nL + 2 * nS)
  

  marker1_all <- fit$summary.fitted.values[idx_surv1, "mean"]
  marker2_all <- fit$summary.fitted.values[idx_surv2, "mean"]
  

  eta1_all <- fit$summary.linear.predictor[idx_surv1, "mean"]
  eta2_all <- fit$summary.linear.predictor[idx_surv2, "mean"]
  

  alpha1 <- fit$summary.hyperpar["alpha parameter for weibullsurv[3]", "mean"]
  alpha2 <- fit$summary.hyperpar["alpha parameter for weibullsurv[4]", "mean"]
  
  Time_eval <- as.numeric(dataS$Time)
  
  train_dat <- data.frame(
    row_id = train_pos,
    ptid = dataS$ptid[train_pos],
    Time = Time_eval[train_pos],
    C = dataS$C[train_pos],
    event1 = event1_full[train_pos],
    event2 = event2_full[train_pos],
    marker1 = marker1_all[train_pos],
    marker2 = marker2_all[train_pos],
    eta1 = eta1_all[train_pos],
    eta2 = eta2_all[train_pos],
    alpha1 = alpha1,
    alpha2 = alpha2,
    fold = k
  ) %>%
    filter(is.finite(Time), Time > 0)
  
  valid_dat <- data.frame(
    row_id = valid_pos,
    ptid = dataS$ptid[valid_pos],
    Time = Time_eval[valid_pos],
    C = dataS$C[valid_pos],
    event1 = event1_full[valid_pos],
    event2 = event2_full[valid_pos],
    marker1 = marker1_all[valid_pos],
    marker2 = marker2_all[valid_pos],
    eta1 = eta1_all[valid_pos],
    eta2 = eta2_all[valid_pos],
    alpha1 = alpha1,
    alpha2 = alpha2,
    fold = k
  ) %>%
    filter(is.finite(Time), Time > 0)
  
  fit <- NULL
  gc()
  
  return(list(
    fold = k,
    train_dat = train_dat,
    valid_dat = valid_dat
  ))
}


trapz_mean <- function(x, y) {
  ok <- !(is.na(x) | is.na(y))
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
}

make_delta <- function(C, event_code) {
  if (event_code == 1) {
    delta <- ifelse(C == 1, 1, ifelse(C == 2, 2, 0))
  } else if (event_code == 2) {
    delta <- ifelse(C == 2, 1, ifelse(C == 1, 2, 0))
  } else {
    stop("event_code must be 1 or 2")
  }
  return(delta)
}

predict_surv_cox <- function(fit, newdata, times) {
  lp <- predict(fit, newdata = newdata, type = "lp")
  bh <- basehaz(fit, centered = FALSE)
  H0 <- approx(
    x = c(0, bh$time),
    y = c(0, bh$hazard),
    xout = times,
    method = "constant",
    f = 0,
    rule = 2
  )$y
  
  S <- outer(
    exp(lp),
    H0,
    FUN = function(risk, h0) exp(-h0 * risk)
  )
  return(S)
}

calc_auc_from_risk <- function(valid_dat, risk_mat, event_code, outcome_name,
                               times, min_events = 20, min_controls = 20) {
  delta <- make_delta(valid_dat$C, event_code)
  out <- data.frame(
    time = times,
    AUC = NA_real_,
    n_event_by_t = NA_integer_,
    n_control_by_t = NA_integer_,
    outcome = outcome_name
  )
  
  for (j in seq_along(times)) {
    tt <- times[j]
    risk_j <- risk_mat[, j]
    ok <- is.finite(valid_dat$Time) &
      is.finite(delta) &
      is.finite(risk_j)
    T0 <- valid_dat$Time[ok]
    delta0 <- delta[ok]
    risk0 <- risk_j[ok]
    n_event <- sum(delta0 == 1 & T0 <= tt, na.rm = TRUE)
    n_control <- sum(T0 > tt, na.rm = TRUE)
    out$n_event_by_t[j] <- n_event
    out$n_control_by_t[j] <- n_control
    
    if (n_event < min_events || n_control < min_controls) next
    roc_obj <- tryCatch(
      timeROC(
        T = T0,
        delta = delta0,
        marker = risk0,
        cause = 1,
        weighting = "marginal",
        times = tt,
        iid = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(roc_obj)) {
      if (!is.null(roc_obj$AUC_2)) {
        auc_j <- as.numeric(roc_obj$AUC_2[1])
      } else {
        auc_j <- as.numeric(roc_obj$AUC[1])
      }
      
      if (is.finite(auc_j)) {
        out$AUC[j] <- min(max(auc_j, 0), 1)
      }
    }
  }
  
  return(out)
}


res_fold1 <- fit_gp_one_fold(
  k = 1,
  dataL = dataL,
  dataS = dataS,
  folds = folds,
  g = g
)

saveRDS(res_fold1, "cv_gp_fold_1.rds")


res_fold2 <- fit_gp_one_fold(
  k = 2,
  dataL = dataL,
  dataS = dataS,
  folds = folds,
  g = g
)

saveRDS(res_fold2, "cv_gp_fold_2.rds")


res_fold3 <- fit_gp_one_fold(
  k = 3,
  dataL = dataL,
  dataS = dataS,
  folds = folds,
  g = g
)

saveRDS(res_fold3, "cv_gp_fold_3.rds")

res_fold4 <- fit_gp_one_fold(
  k = 4,
  dataL = dataL,
  dataS = dataS,
  folds = folds,
  g = g
)

saveRDS(res_fold4, "cv_gp_fold_4.rds")

res_fold5 <- fit_gp_one_fold(
  k = 5,
  dataL = dataL,
  dataS = dataS,
  folds = folds,
  g = g
)

saveRDS(res_fold5, "cv_gp_fold_5.rds")


cv_results <- list(
  readRDS("cv_gp_fold_1.rds"),
  readRDS("cv_gp_fold_2.rds"),
  readRDS("cv_gp_fold_3.rds"),
  readRDS("cv_gp_fold_4.rds"),
  readRDS("cv_gp_fold_5.rds")
)

cv_results_ok <- cv_results[!sapply(cv_results, is.null)]

cv_valid_all <- bind_rows(
  lapply(cv_results_ok, function(x) x$valid_dat)
)

cv_train_all <- bind_rows(
  lapply(cv_results_ok, function(x) x$train_dat)
)

cv_times <- sort(unique(as.numeric(
  quantile(cv_valid_all$Time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
)))

print(cv_times)



risk_valid_event1_all <- matrix(
  NA_real_,
  nrow = nrow(cv_valid_all),
  ncol = length(cv_times)
)

risk_valid_event2_all <- matrix(
  NA_real_,
  nrow = nrow(cv_valid_all),
  ncol = length(cv_times)
)

S_valid_event1_all <- matrix(
  NA_real_,
  nrow = nrow(cv_valid_all),
  ncol = length(cv_times)
)

S_valid_event2_all <- matrix(
  NA_real_,
  nrow = nrow(cv_valid_all),
  ncol = length(cv_times)
)

for (k in sort(unique(cv_valid_all$fold))) {
  
  idx_valid <- which(cv_valid_all$fold == k)
  
  train_k <- cv_train_all %>%
    filter(fold == k)
  
  valid_k <- cv_valid_all %>%
    filter(fold == k)
  
  if (sum(train_k$event1 == 1, na.rm = TRUE) >= 10 &&
      length(unique(train_k$marker1)) > 5) {
    
    cox1 <- coxph(
      Surv(Time, event1) ~ marker1,
      data = train_k,
      x = TRUE,
      y = TRUE
    )
    
    S1 <- predict_surv_cox(
      fit = cox1,
      newdata = valid_k,
      times = cv_times
    )
    
    S_valid_event1_all[idx_valid, ] <- S1
    risk_valid_event1_all[idx_valid, ] <- 1 - S1
  }

  if (sum(train_k$event2 == 1, na.rm = TRUE) >= 10 &&
      length(unique(train_k$marker2)) > 5) {
    
    cox2 <- coxph(
      Surv(Time, event2) ~ marker2,
      data = train_k,
      x = TRUE,
      y = TRUE
    )
    
    S2 <- predict_surv_cox(
      fit = cox2,
      newdata = valid_k,
      times = cv_times
    )
    
    S_valid_event2_all[idx_valid, ] <- S2
    risk_valid_event2_all[idx_valid, ] <- 1 - S2
  }
}


## Cross-validated AUC and iAUC
auc_event1 <- calc_auc_from_risk(
  valid_dat = cv_valid_all,
  risk_mat = risk_valid_event1_all,
  event_code = 1,
  outcome_name = "AIDS onset",
  times = cv_times,
  min_events = 20,
  min_controls = 20
)

auc_event2 <- calc_auc_from_risk(
  valid_dat = cv_valid_all,
  risk_mat = risk_valid_event2_all,
  event_code = 2,
  outcome_name = "Pre-AIDS death",
  times = cv_times,
  min_events = 20,
  min_controls = 20
)

cv_gp_auc_curve <- bind_rows(auc_event1, auc_event2)

cv_gp_iauc <- cv_gp_auc_curve %>%
  group_by(outcome) %>%
  summarise(
    GP_iAUC = trapz_mean(time, AUC),
    GP_iAUC_percent = GP_iAUC * 100,
    n_valid_timepoints = sum(!is.na(AUC)),
    .groups = "drop"
  )

print(cv_gp_auc_curve)
print(cv_gp_iauc)


## Cross-validated IBS 

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

censor_ind <- as.numeric(cv_valid_all$C == 3)

km_cens_cv <- survfit(
  Surv(time = cv_valid_all$Time, event = censor_ind) ~ 1
)

G_eval <- function(t, kmfit) {
  sf <- summary(kmfit, times = t, extend = TRUE)$surv
  if (length(sf) == 0 || is.na(sf)) return(1)
  sf
}

brier_ipcw_cr <- function(time_point, event_time, event_type, pred_prob, cause, kmfit) {
  
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

cv_valid_all$lambda1 <- exp(cv_valid_all$eta1)
cv_valid_all$lambda2 <- exp(cv_valid_all$eta2)

pred_cif1_cv <- sapply(cv_times, function(tt) {
  sapply(seq_len(nrow(cv_valid_all)), function(i) {
    cif1_fun(
      t = tt,
      lambda1 = cv_valid_all$lambda1[i],
      lambda2 = cv_valid_all$lambda2[i],
      alpha1 = cv_valid_all$alpha1[i],
      alpha2 = cv_valid_all$alpha2[i]
    )
  })
})

pred_cif2_cv <- sapply(cv_times, function(tt) {
  sapply(seq_len(nrow(cv_valid_all)), function(i) {
    cif2_fun(
      t = tt,
      lambda1 = cv_valid_all$lambda1[i],
      lambda2 = cv_valid_all$lambda2[i],
      alpha1 = cv_valid_all$alpha1[i],
      alpha2 = cv_valid_all$alpha2[i]
    )
  })
})

pred_cif1_cv <- as.matrix(pred_cif1_cv)
pred_cif2_cv <- as.matrix(pred_cif2_cv)

cat("Range of CV CIF1:", range(pred_cif1_cv, na.rm = TRUE), "\n")
cat("Range of CV CIF2:", range(pred_cif2_cv, na.rm = TRUE), "\n")

brier_event1 <- sapply(seq_along(cv_times), function(j) {
  brier_ipcw_cr(
    time_point = cv_times[j],
    event_time = cv_valid_all$Time,
    event_type = cv_valid_all$C,
    pred_prob = pred_cif1_cv[, j],
    cause = 1,
    kmfit = km_cens_cv
  )
})

brier_event2 <- sapply(seq_along(cv_times), function(j) {
  brier_ipcw_cr(
    time_point = cv_times[j],
    event_time = cv_valid_all$Time,
    event_type = cv_valid_all$C,
    pred_prob = pred_cif2_cv[, j],
    cause = 2,
    kmfit = km_cens_cv
  )
})

cv_gp_ibs_curve <- data.frame(
  time = cv_times,
  Brier_AIDS_onset = brier_event1,
  Brier_Pre_AIDS_death = brier_event2
)

cv_gp_ibs <- data.frame(
  outcome = c("AIDS onset", "Pre-AIDS death"),
  GP_IBS = c(
    trapz_mean(cv_times, brier_event1),
    trapz_mean(cv_times, brier_event2)
  )
)

print(cv_gp_ibs_curve)
print(cv_gp_ibs)