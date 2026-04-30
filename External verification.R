library(joineR)
library(INLA)
library(dplyr)
library(survivalROC)

#External verification 1

data(aids, package = "joineR")

dat <- aids

dat <- dat %>%
  mutate(
    id        = as.numeric(id),
    death     = as.numeric(death),
    logCD4    = log(CD4),
    drug_bin  = ifelse(drug   == "ddI", 1, 0),
    female    = ifelse(gender == "female", 1, 0),
    prevOI_bin= ifelse(prevOI == "AIDS", 1, 0),
    AZT_bin   = ifelse(AZT    == "failure", 1, 0)
  )


dataS <- dat %>%
  arrange(id, obstime) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(id, time, death, drug_bin, female, prevOI_bin, AZT_bin)

dataL_gau <- dat %>%
  select(id, obstime, logCD4, drug_bin, female, prevOI_bin, AZT_bin)


max_logCD4 <- max(dat$logCD4, na.rm = TRUE)

dat <- dat %>%
  mutate(
    CD4min = max_logCD4 - logCD4
  )

u <- as.numeric(quantile(dat$CD4min, 0.90, na.rm = TRUE))

dataL_gp <- dat %>%
  mutate(excess = CD4min - u) %>%
  filter(excess > 0)

keep_id <- names(which(table(dataL_gp$id) >= 2))
dataL_gp <- dataL_gp %>% filter(id %in% keep_id)
dataS_gp <- dataS %>% filter(id %in% keep_id)

dataL_gau_cmp <- dataL_gau %>% filter(id %in% keep_id)
dataS_gau_cmp <- dataS %>% filter(id %in% keep_id)

nL_gp <- nrow(dataL_gp)
nS_gp <- nrow(dataS_gp)

fixed_gp <- data.frame(
  mu = factor(c(rep(1, nL_gp), rep(1, nL_gp), rep(2, nS_gp))),
  
  obstimeL = c(dataL_gp$obstime, dataL_gp$obstime, rep(0, nS_gp)),
  
  drugL    = c(dataL_gp$drug_bin,   dataL_gp$drug_bin,   rep(0, nS_gp)),
  femaleL  = c(dataL_gp$female,     dataL_gp$female,     rep(0, nS_gp)),
  prevOIL  = c(dataL_gp$prevOI_bin, dataL_gp$prevOI_bin, rep(0, nS_gp)),
  AZTL     = c(dataL_gp$AZT_bin,    dataL_gp$AZT_bin,    rep(0, nS_gp)),
  
  drugS    = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp$drug_bin),
  femaleS  = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp$female),
  prevOIS  = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp$prevOI_bin),
  AZTS     = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp$AZT_bin)
)

random_gp <- list(
  linpredL  = c(rep(NA, nL_gp), dataL_gp$id, rep(NA, nS_gp)),
  linpredL2 = c(rep(NA, nL_gp), rep(-1, nL_gp), rep(NA, nS_gp)),
  beta      = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$id)
)

joint_gp <- c(fixed_gp, random_gp)

y.long.gp <- c(dataL_gp$excess, rep(NA, nL_gp), rep(NA, nS_gp))
y.eta.gp  <- c(rep(NA, nL_gp), rep(0, nL_gp), rep(NA, nS_gp))
y.surv.gp <- inla.surv(
  time  = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$time),
  event = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$death)
)

joint_gp$Y <- list(y.long.gp, y.eta.gp, y.surv.gp)

formula_gp <- Y ~
  obstimeL +
  drugL + femaleL + prevOIL + AZTL +
  drugS + femaleS + prevOIS + AZTS +
  f(linpredL, linpredL2, model = "iid",
    hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(beta, copy = "linpredL",
    hyper = list(beta = list(fixed = FALSE)))

fit_gp <- inla(
  formula_gp,
  family = c("gp", "gaussian", "weibullsurv"),
  control.family = list(
    list(control.link = list(model = "quantile", quantile = 0.6)),
    list(),
    list()
  ),
  control.compute = list(dic = TRUE, waic = TRUE),
  data = joint_gp,
  verbose = FALSE
)

summary(fit_gp)


nL_gau <- nrow(dataL_gau_cmp)
nS_gau <- nrow(dataS_gau_cmp)

fixed_gau <- data.frame(
  mu = factor(c(rep(1, nL_gau), rep(1, nL_gau), rep(2, nS_gau))),
  
  obstimeL = c(dataL_gau_cmp$obstime, dataL_gau_cmp$obstime, rep(0, nS_gau)),
  
  drugL    = c(dataL_gau_cmp$drug_bin,   dataL_gau_cmp$drug_bin,   rep(0, nS_gau)),
  femaleL  = c(dataL_gau_cmp$female,     dataL_gau_cmp$female,     rep(0, nS_gau)),
  prevOIL  = c(dataL_gau_cmp$prevOI_bin, dataL_gau_cmp$prevOI_bin, rep(0, nS_gau)),
  AZTL     = c(dataL_gau_cmp$AZT_bin,    dataL_gau_cmp$AZT_bin,    rep(0, nS_gau)),
  
  drugS    = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_cmp$drug_bin),
  femaleS  = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_cmp$female),
  prevOIS  = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_cmp$prevOI_bin),
  AZTS     = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_cmp$AZT_bin)
)

random_gau <- list(
  linpredL  = c(rep(NA, nL_gau), dataL_gau_cmp$id, rep(NA, nS_gau)),
  linpredL2 = c(rep(NA, nL_gau), rep(-1, nL_gau), rep(NA, nS_gau)),
  beta      = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau_cmp$id)
)

joint_gau <- c(fixed_gau, random_gau)

y.long.gau <- c(dataL_gau_cmp$logCD4, rep(NA, nL_gau), rep(NA, nS_gau))
y.eta.gau  <- c(rep(NA, nL_gau), rep(0, nL_gau), rep(NA, nS_gau))
y.surv.gau <- inla.surv(
  time  = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau_cmp$time),
  event = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau_cmp$death)
)

joint_gau$Y <- list(y.long.gau, y.eta.gau, y.surv.gau)

formula_gau <- Y ~
  obstimeL +
  drugL + femaleL + prevOIL + AZTL +
  drugS + femaleS + prevOIS + AZTS +
  f(linpredL, linpredL2, model = "iid",
    hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(beta, copy = "linpredL",
    hyper = list(beta = list(fixed = FALSE)))

fit_gau <- inla(
  formula_gau,
  family = c("gaussian", "gaussian", "weibullsurv"),
  control.family = list(
    list(),
    list(),
    list()
  ),
  control.compute = list(dic = TRUE, waic = TRUE),
  data = joint_gau,
  verbose = FALSE
)

summary(fit_gau)


idx_gp_surv  <- (2 * nL_gp  + 1):(2 * nL_gp  + nS_gp)
idx_gau_surv <- (2 * nL_gau + 1):(2 * nL_gau + nS_gau)

marker_gp  <- -fit_gp$summary.fitted.values[idx_gp_surv,  "mean"]
marker_gau <- -fit_gau$summary.fitted.values[idx_gau_surv, "mean"]


roc_time <- as.numeric(quantile(dataS_gp$time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE))

roc_res <- data.frame(
  time = roc_time,
  AUC_GP = NA_real_,
  AUC_Gaussian = NA_real_
)

for (j in seq_along(roc_time)) {
  tt <- roc_time[j]
  
  roc1 <- survivalROC(
    Stime = dataS_gp$time,
    status = dataS_gp$death,
    marker = marker_gp,
    predict.time = tt,
    method = "KM"
  )
  
  roc2 <- survivalROC(
    Stime = dataS_gp$time,
    status = dataS_gp$death,
    marker = marker_gau,
    predict.time = tt,
    method = "KM"
  )
  
  roc_res$AUC_GP[j]       <- roc1$AUC
  roc_res$AUC_Gaussian[j] <- roc2$AUC
}

print(roc_res)


trapz_mean <- function(x, y) {
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
}

iAUC_GP       <- 1-trapz_mean(roc_res$time, roc_res$AUC_GP)
iAUC_Gaussian <- 1-trapz_mean(roc_res$time, roc_res$AUC_Gaussian)

cat("iAUC_GP =", iAUC_GP, "\n")
cat("iAUC_Gaussian =", iAUC_Gaussian, "\n")


cat("\n=== Model comparison ===\n")
cat("GP   DIC :", fit_gp$dic$dic, "\n")
cat("GP   WAIC:", fit_gp$waic$waic, "\n")
cat("GAU  DIC :", fit_gau$dic$dic, "\n")
cat("GAU  WAIC:", fit_gau$waic$waic, "\n")
cat("GP   iAUC:", iAUC_GP, "\n")
cat("GAU iAUC:", iAUC_Gaussian, "\n")




#External verification 2

data <- read_xlsx("Data.xlsx")

names(data) <- make.names(names(data), unique = TRUE)

data$ptid <- seq_len(nrow(data))


num_vars <- c(
  "AGE_DDVIH", "DUREE_VIH", "duration_ART",
  "Time.to.LTFU", "Time.to.death", "Time.to.Event",
  "fvcd4", "lvcd4", "CD42018"
)

for (v in intersect(num_vars, names(data))) {
  data[[v]] <- suppressWarnings(as.numeric(data[[v]]))
}

date_vars <- c(
  "dat_fcd4", "dat_lcd4", "date_vis_2012",
  "datedc", "datearv", "ddvih", "datremp", "datdv"
)

for (v in intersect(date_vars, names(data))) {
  x <- data[[v]]
  d1 <- suppressWarnings(as.Date(x, origin = "1899-12-30"))
  d2 <- suppressWarnings(as.Date(x))
  data[[v]] <- ifelse(is.na(d1), d2, d1)
  data[[v]] <- as.Date(data[[v]], origin = "1970-01-01")
}


data$sex <- case_when(
  data$SEXE == 2 ~ 1,
  data$SEXE == 1 ~ 0,
  TRUE ~ 0
)


data$age <- data$AGE_DDVIH
data$age[is.na(data$age)] <- 0


data$artdur <- data$duration_ART
data$artdur[is.na(data$artdur)] <- 0


dataS <- data

dataS$C <- 0

death_first <- !is.na(dataS$Time.to.death) &
  (is.na(dataS$Time.to.LTFU) | dataS$Time.to.death <= dataS$Time.to.LTFU)

ltfu_first <- !is.na(dataS$Time.to.LTFU) &
  (is.na(dataS$Time.to.death) | dataS$Time.to.LTFU < dataS$Time.to.death)

dataS$C[death_first] <- 1
dataS$C[ltfu_first]  <- 2

if ("Patient.outcomes" %in% names(dataS)) {
  out1 <- tolower(trimws(as.character(dataS$Patient.outcomes)))
  dataS$C[dataS$C == 0 & grepl("dead|death|deces|died|4", out1)] <- 1
  dataS$C[dataS$C == 0 & grepl("ltfu|loss|perdu|3|2", out1)] <- 2
}

if ("New.Outcome" %in% names(dataS)) {
  out2 <- suppressWarnings(as.numeric(dataS$New.Outcome))
  dataS$C[dataS$C == 0 & out2 == 2 & !is.na(dataS$Time.to.LTFU)] <- 2
  dataS$C[dataS$C == 0 & out2 == 2 & !is.na(dataS$Time.to.death)] <- 1
}


dataS$Time <- NA_real_
dataS$Time[dataS$C == 1] <- dataS$Time.to.death[dataS$C == 1]
dataS$Time[dataS$C == 2] <- dataS$Time.to.LTFU[dataS$C == 2]
dataS$Time[dataS$C == 0] <- dataS$Time.to.Event[dataS$C == 0]

fallback_time <- pmin(dataS$Time.to.Event, dataS$Time.to.death, dataS$Time.to.LTFU, na.rm = TRUE)
fallback_time[is.infinite(fallback_time)] <- NA_real_
dataS$Time[is.na(dataS$Time)] <- fallback_time[is.na(dataS$Time)]

dataS <- dataS %>%
  transmute(
    ptid,
    centre,
    sex,
    age,
    artdur,
    C,
    Time,
    Time.to.Event,
    Time.to.death,
    Time.to.LTFU
  ) %>%
  filter(!is.na(Time), Time > 0)


long_list <- list(
  data.frame(
    ptid = data$ptid,
    visit = "first",
    date = data$dat_fcd4,
    cd4 = data$fvcd4
  ),
  data.frame(
    ptid = data$ptid,
    visit = "last",
    date = data$dat_lcd4,
    cd4 = data$lvcd4
  ),
  data.frame(
    ptid = data$ptid,
    visit = "y2018",
    date = as.Date("2018-12-31"),
    cd4 = data$CD42018
  )
)

dataL <- bind_rows(long_list) %>%
  filter(!is.na(cd4), cd4 > 0)


dataL <- dataL %>%
  left_join(
    dataS %>% select(ptid, sex, age, artdur, C, Time_surv = Time),
    by = "ptid"
  )


first_date_df <- dataL %>%
  group_by(ptid) %>%
  summarise(first_date = min(date, na.rm = TRUE), .groups = "drop")

dataL <- dataL %>%
  left_join(first_date_df, by = "ptid")


dataL$Time <- ifelse(
  !is.na(dataL$date) & !is.na(dataL$first_date),
  as.numeric(difftime(dataL$date, dataL$first_date, units = "days")) / 30.4375,
  NA_real_
)


dataL$Time[is.na(dataL$Time) & dataL$visit == "first"] <- 0
dataL$Time[is.na(dataL$Time) & dataL$visit == "last"]  <- dataL$Time_surv[is.na(dataL$Time) & dataL$visit == "last"] * 0.7
dataL$Time[is.na(dataL$Time) & dataL$visit == "y2018"] <- dataL$Time_surv[is.na(dataL$Time) & dataL$visit == "y2018"] * 0.9

dataL <- dataL %>%
  filter(!is.na(Time), Time >= 0) %>%
  arrange(ptid, Time)


dataL$logcd4 <- log(dataL$cd4)


max_logcd4 <- max(dataL$logcd4, na.rm = TRUE)
dataL$CD4min <- max_logcd4 - dataL$logcd4

u <- as.numeric(quantile(dataL$CD4min, 0.60, na.rm = TRUE))

dataL_gp <- dataL %>%
  mutate(CD4_excess = CD4min - u)


keep_gp <- unique(dataL_gp$ptid[dataL_gp$CD4_excess > 0])

dataL_gp_use <- dataL_gp %>%
  filter(ptid %in% keep_gp, CD4_excess > 0)

dataS_gp_use <- dataS %>%
  filter(ptid %in% keep_gp)

dataL_gau_use <- dataL
dataS_gau_use <- dataS


dataE1_gp <- dataS_gp_use
dataE1_gp$event <- ifelse(dataE1_gp$C == 1, 1, 0)

dataE2_gp <- dataS_gp_use
dataE2_gp$event <- ifelse(dataE2_gp$C == 2, 1, 0)

dataE1_gau <- dataS_gau_use
dataE1_gau$event <- ifelse(dataE1_gau$C == 1, 1, 0)

dataE2_gau <- dataS_gau_use
dataE2_gau$event <- ifelse(dataE2_gau$C == 2, 1, 0)


gp_id_map <- data.frame(
  ptid_old = sort(unique(dataS_gp_use$ptid)),
  ptid_new = seq_along(sort(unique(dataS_gp_use$ptid)))
)

dataS_gp_use <- dataS_gp_use %>%
  left_join(gp_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataL_gp_use <- dataL_gp_use %>%
  left_join(gp_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataE1_gp <- dataE1_gp %>%
  left_join(gp_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataE2_gp <- dataE2_gp %>%
  left_join(gp_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)


gau_id_map <- data.frame(
  ptid_old = sort(unique(dataS_gau_use$ptid)),
  ptid_new = seq_along(sort(unique(dataS_gau_use$ptid)))
)

dataS_gau_use <- dataS_gau_use %>%
  left_join(gau_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataL_gau_use <- dataL_gau_use %>%
  left_join(gau_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataE1_gau <- dataE1_gau %>%
  left_join(gau_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)

dataE2_gau <- dataE2_gau %>%
  left_join(gau_id_map, by = c("ptid" = "ptid_old")) %>%
  mutate(ptid = ptid_new) %>%
  select(-ptid_new)


nL_gp <- nrow(dataL_gp_use)
nS_gp <- nrow(dataS_gp_use)
gp_values <- sort(unique(dataL_gp_use$ptid))

fixed_gp <- data.frame(
  mu = as.factor(c(rep(1, nL_gp), rep(1, nL_gp), rep(2, nS_gp), rep(3, nS_gp))),
  
  ageL = c(dataL_gp_use$age, dataL_gp_use$age, rep(0, nS_gp), rep(0, nS_gp)),
  age1 = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp_use$age, rep(0, nS_gp)),
  age2 = c(rep(0, nL_gp), rep(0, nL_gp), rep(0, nS_gp), dataS_gp_use$age),
  
  sexL = c(dataL_gp_use$sex, dataL_gp_use$sex, rep(0, nS_gp), rep(0, nS_gp)),
  sex1 = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp_use$sex, rep(0, nS_gp)),
  sex2 = c(rep(0, nL_gp), rep(0, nL_gp), rep(0, nS_gp), dataS_gp_use$sex),
  
  artdurL = c(dataL_gp_use$artdur, dataL_gp_use$artdur, rep(0, nS_gp), rep(0, nS_gp)),
  artdur1 = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp_use$artdur, rep(0, nS_gp)),
  artdur2 = c(rep(0, nL_gp), rep(0, nL_gp), rep(0, nS_gp), dataS_gp_use$artdur),
  
  tL = c(dataL_gp_use$Time, dataL_gp_use$Time, rep(0, nS_gp), rep(0, nS_gp))
)

random_gp <- list(
  linpredL  = c(rep(NA, nL_gp), dataL_gp_use$ptid, rep(NA, nS_gp), rep(NA, nS_gp)),
  linpredL2 = c(rep(NA, nL_gp), rep(-1, nL_gp), rep(NA, nS_gp), rep(NA, nS_gp)),
  beta1     = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp_use$ptid, rep(NA, nS_gp)),
  beta2     = c(rep(NA, nL_gp), rep(NA, nL_gp), rep(NA, nS_gp), dataS_gp_use$ptid)
)

joint_gp <- c(fixed_gp, random_gp)

y.long.gp <- c(dataL_gp_use$CD4_excess, rep(NA, nL_gp), rep(NA, nS_gp), rep(NA, nS_gp))
y.eta.gp  <- c(rep(NA, nL_gp), rep(0, nL_gp), rep(NA, nS_gp), rep(NA, nS_gp))

y.surv1.gp <- inla.surv(
  time  = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp_use$Time, rep(NA, nS_gp)),
  event = c(rep(NA, nL_gp), rep(NA, nL_gp), dataE1_gp$event, rep(NA, nS_gp))
)

y.surv2.gp <- inla.surv(
  time  = c(rep(NA, nL_gp), rep(NA, nL_gp), rep(NA, nS_gp), dataS_gp_use$Time),
  event = c(rep(NA, nL_gp), rep(NA, nL_gp), rep(NA, nS_gp), dataE2_gp$event)
)

joint_gp$Y <- list(y.long.gp, y.eta.gp, y.surv1.gp, y.surv2.gp)

formula_gp <- Y ~
  ageL + age1 + age2 +
  sexL + sex1 + sex2 +
  artdurL + artdur1 + artdur2 +
  tL +
  f(linpredL, linpredL2, model = "iid",
    values = gp_values,
    hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(beta1, copy = "linpredL",
    values = gp_values,
    hyper = list(beta = list(fixed = FALSE))) +
  f(beta2, copy = "linpredL",
    values = gp_values,
    hyper = list(beta = list(fixed = FALSE)))

fit_gp <- inla(
  formula_gp,
  family = c("gp", "gaussian", "weibullsurv", "weibullsurv"),
  control.family = list(
    list(control.link = list(model = "quantile", quantile = 0.6)),
    list(),
    list(),
    list()
  ),
  control.compute = list(dic = TRUE, waic = TRUE),
  data = joint_gp,
  verbose = FALSE
)


nL_gau <- nrow(dataL_gau_use)
nS_gau <- nrow(dataS_gau_use)
gau_values <- sort(unique(dataL_gau_use$ptid))

fixed_gau <- data.frame(
  mu = as.factor(c(rep(1, nL_gau), rep(1, nL_gau), rep(2, nS_gau), rep(3, nS_gau))),
  
  ageL = c(dataL_gau_use$age, dataL_gau_use$age, rep(0, nS_gau), rep(0, nS_gau)),
  age1 = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_use$age, rep(0, nS_gau)),
  age2 = c(rep(0, nL_gau), rep(0, nL_gau), rep(0, nS_gau), dataS_gau_use$age),
  
  sexL = c(dataL_gau_use$sex, dataL_gau_use$sex, rep(0, nS_gau), rep(0, nS_gau)),
  sex1 = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_use$sex, rep(0, nS_gau)),
  sex2 = c(rep(0, nL_gau), rep(0, nL_gau), rep(0, nS_gau), dataS_gau_use$sex),
  
  artdurL = c(dataL_gau_use$artdur, dataL_gau_use$artdur, rep(0, nS_gau), rep(0, nS_gau)),
  artdur1 = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau_use$artdur, rep(0, nS_gau)),
  artdur2 = c(rep(0, nL_gau), rep(0, nL_gau), rep(0, nS_gau), dataS_gau_use$artdur),
  
  tL = c(dataL_gau_use$Time, dataL_gau_use$Time, rep(0, nS_gau), rep(0, nS_gau))
)

random_gau <- list(
  linpredL  = c(rep(NA, nL_gau), dataL_gau_use$ptid, rep(NA, nS_gau), rep(NA, nS_gau)),
  linpredL2 = c(rep(NA, nL_gau), rep(-1, nL_gau), rep(NA, nS_gau), rep(NA, nS_gau)),
  beta1     = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau_use$ptid, rep(NA, nS_gau)),
  beta2     = c(rep(NA, nL_gau), rep(NA, nL_gau), rep(NA, nS_gau), dataS_gau_use$ptid)
)

joint_gau <- c(fixed_gau, random_gau)

y.long.gau <- c(dataL_gau_use$logcd4, rep(NA, nL_gau), rep(NA, nS_gau), rep(NA, nS_gau))
y.eta.gau  <- c(rep(NA, nL_gau), rep(0, nL_gau), rep(NA, nS_gau), rep(NA, nS_gau))

y.surv1.gau <- inla.surv(
  time  = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau_use$Time, rep(NA, nS_gau)),
  event = c(rep(NA, nL_gau), rep(NA, nL_gau), dataE1_gau$event, rep(NA, nS_gau))
)

y.surv2.gau <- inla.surv(
  time  = c(rep(NA, nL_gau), rep(NA, nL_gau), rep(NA, nS_gau), dataS_gau_use$Time),
  event = c(rep(NA, nL_gau), rep(NA, nL_gau), rep(NA, nS_gau), dataE2_gau$event)
)

joint_gau$Y <- list(y.long.gau, y.eta.gau, y.surv1.gau, y.surv2.gau)

formula_gau <- Y ~
  ageL + age1 + age2 +
  sexL + sex1 + sex2 +
  artdurL + artdur1 + artdur2 +
  tL +
  f(linpredL, linpredL2, model = "iid",
    values = gau_values,
    hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(beta1, copy = "linpredL",
    values = gau_values,
    hyper = list(beta = list(fixed = FALSE))) +
  f(beta2, copy = "linpredL",
    values = gau_values,
    hyper = list(beta = list(fixed = FALSE)))

fit_gau <- inla(
  formula_gau,
  family = c("gaussian", "gaussian", "weibullsurv", "weibullsurv"),
  control.family = list(
    list(),
    list(),
    list(),
    list()
  ),
  control.compute = list(dic = TRUE, waic = TRUE),
  control.inla = list(strategy = "simplified.laplace"),
  data = joint_gau,
  verbose = FALSE
)


idx_gp_surv1  <- (2 * nL_gp + 1):(2 * nL_gp + nS_gp)
idx_gp_surv2  <- (2 * nL_gp + nS_gp + 1):(2 * nL_gp + 2 * nS_gp)

idx_gau_surv1 <- (2 * nL_gau + 1):(2 * nL_gau + nS_gau)
idx_gau_surv2 <- (2 * nL_gau + nS_gau + 1):(2 * nL_gau + 2 * nS_gau)

marker_gp_1  <- fit_gp$summary.fitted.values[idx_gp_surv1, "mean"]
marker_gp_2  <- fit_gp$summary.fitted.values[idx_gp_surv2, "mean"]

marker_gau_1 <- fit_gau$summary.fitted.values[idx_gau_surv1, "mean"]
marker_gau_2 <- fit_gau$summary.fitted.values[idx_gau_surv2, "mean"]


roc_time_gp  <- as.numeric(quantile(dataS_gp_use$Time,  probs = seq(0.1, 0.9, 0.1), na.rm = TRUE))
roc_time_gau <- as.numeric(quantile(dataS_gau_use$Time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE))

calc_auc_series <- function(Stime, status, marker, times, min_events = 3) {
  out <- data.frame(time = times, AUC = NA_real_, n_event_by_t = NA_integer_)
  
  for (j in seq_along(times)) {
    tt <- times[j]
    
    n_evt <- sum(status == 1 & Stime <= tt, na.rm = TRUE)
    out$n_event_by_t[j] <- n_evt
    
    if (n_evt < min_events) {
      out$AUC[j] <- NA_real_
      next
    }
    
    roc_obj <- try(
      survivalROC(
        Stime = Stime,
        status = status,
        marker = marker,
        predict.time = tt,
        method = "KM"
      ),
      silent = TRUE
    )
    
    if (!inherits(roc_obj, "try-error")) {
      auc_j <- roc_obj$AUC
      
      if (is.finite(auc_j)) {
        out$AUC[j] <- min(max(auc_j, 0), 1)
      } else {
        out$AUC[j] <- NA_real_
      }
    }
  }
  
  out
}


roc_gp_death <- calc_auc_series(
  Stime  = dataS_gp_use$Time,
  status = dataE1_gp$event,
  marker = marker_gp_1,
  times  = roc_time_gp,
  min_events = 3
)
roc_gp_death$model <- "GP"
roc_gp_death$outcome <- "death"


roc_gau_death <- calc_auc_series(
  Stime  = dataS_gau_use$Time,
  status = dataE1_gau$event,
  marker = marker_gau_1,
  times  = roc_time_gau,
  min_events = 3
)
roc_gau_death$model <- "Gaussian"
roc_gau_death$outcome <- "death"


roc_gp_ltfu <- calc_auc_series(
  Stime  = dataS_gp_use$Time,
  status = dataE2_gp$event,
  marker = marker_gp_2,
  times  = roc_time_gp,
  min_events = 3
)
roc_gp_ltfu$model <- "GP"
roc_gp_ltfu$outcome <- "LTFU"


roc_gau_ltfu <- calc_auc_series(
  Stime  = dataS_gau_use$Time,
  status = dataE2_gau$event,
  marker = marker_gau_2,
  times  = roc_time_gau,
  min_events = 3
)
roc_gau_ltfu$model <- "Gaussian"
roc_gau_ltfu$outcome <- "LTFU"

roc_all <- bind_rows(
  roc_gp_death,
  roc_gau_death,
  roc_gp_ltfu,
  roc_gau_ltfu
)

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

iauc_res <- roc_all %>%
  group_by(model, outcome) %>%
  summarise(
    iAUC = trapz_mean(time, AUC),
    n_valid_timepoints = sum(!is.na(AUC)),
    .groups = "drop"
  ) %>%
  mutate(
    iAUC = sprintf("%.4f", iAUC)
  )


print(roc_all)
print(iauc_res)