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
library(survivalROC)
library(survival)

#Sensitivity analysis results of alternative thresholds

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
CD4min_95th_percentile <- quantile(dataL$CD4min, 0.95)
dataL <- dataL[dataL$CD4min > CD4min_95th_percentile, ]
dataL <- dataL[dataL$ptid %in% ptids_to_keep, ]
existing_ptids <- dataL$ptid
dataS <- dataS[dataS$ptid %in% existing_ptids, ]
dataL$CD4_9 <- dataL$CD4min - CD4min_95th_percentile

data1<-dataS
dataE1<-data1
dataE1$event<-dataE1$C
dataE1$event[dataE1$event!=1]<-0
dataE2<-data1
dataE2$event<-dataE2$C
dataE2$event[dataE2$event!=2]<-0
dataE2$event<-dataE2$event/2

nL<-nrow(dataL)
nS<-nrow(dataS)

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      age_confirmL=c(dataL$age_confirm,dataL$age_confirm, rep(0,nS),rep(0,nS)),
                      age_confirm1=c(rep(0,nL),rep(0,nL),dataS$age_confirm,rep(0,nS)),
                      age_confirm2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$age_confirm),
                      sexL=c(dataL$sex,dataL$sex,rep(0,nS),rep(0,nS)),
                      sex1=c(rep(0,nL),rep(0,nL),dataS$sex,rep(0,nS)),
                      sex2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$sex),
                      marriageL=c(dataL$marriage,dataL$marriage,rep(0,nS),rep(0,nS)),
                      marriage1=c(rep(0,nL),rep(0,nL),dataS$marriage,rep(0,nS)),
                      marriage2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$marriage),
                      routeL=c(dataL$route,dataL$route,rep(0,nS),rep(0,nS)),
                      route1=c(rep(0,nL),rep(0,nL),dataS$route,rep(0,nS)),
                      route2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$route),
                      confirm_artL=c(dataL$confirm_art,dataL$confirm_art,rep(0,nS),rep(0,nS)),
                      confirm_art1=c(rep(0,nL),rep(0,nL),dataS$confirm_art,rep(0,nS)),
                      confirm_art2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$confirm_art),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie),
                      ARTdyL=c(dataL$ARTdy,dataL$ARTdy,rep(0,nS),rep(0,nS)),
                      ARTdy1=c(rep(0,nL),rep(0,nL),dataS$ARTdy,rep(0,nS)),
                      ARTdy2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$ARTdy),
                      tL=c(dataL$t,dataL$t,rep(0,nS),rep(0,nS)),
                      t1=c(rep(0,nL),rep(0,nL),dataS$t,rep(0,nS)),
                      t2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$t),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie))

random.eff<-list(timeL =c(dataL$Time,dataL$Time,rep(NA,nS),rep(NA,nS)),
                 idareaS1=c(rep(NA,nL),rep(NA,nL),dataS$idarea,rep(NA,nS)),
                 idareaS2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$idarea),
                 linpredL=c(rep(NA,nL),dataL$ptid,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 beta1=c(rep(NA,nL),rep(NA,nL),dataS$ptid,rep(NA,nS)),
                 beta2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ptid),
                 Lr1 = c(dataL$ptid, dataL$ptid, rep(NA, nS), rep(NA, nS)),
                 C1r1 = c(rep(NA, nL), rep(NA, nL), dataS$ptid, rep(NA, nS)),
                 C2r1 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$ptid))

jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$CD4_9,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
dataS$Time <- as.numeric(as.factor(dataS$Time))
y.survC1 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nL),dataE1$event,rep(NA,nS)))
y.survC2 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataE2$event))
y.joint<-list(y.long,y.eta,y.survC1,y.survC2)
jointdata$Y=y.joint

formula.model2=Y~sex1+sex2+sexL+marriage1+marriage2+marriageL+
  route1+route2+routeL+age_confirm1+age_confirm2+age_confirmL+confirm_art1+confirm_art2+confirm_artL+tL+jie1+jie2+
  f(inla.group(timeL,n=50),model="rw2", scale.model = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(idareaS1, model = "bym2", graph = g, hyper = list
    (prec = list(prior = "pc.prec",param = c(0.5 / 0.31, 0.01)),
      phi = list(prior = "pc",param = c(0.1, 4 / 5))))+
  f(idareaS2, copy="idareaS1", hyper = list(beta = list(fixed = FALSE)))+
  f(linpredL, linpredL2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(Lr1, model="iid")+
  f(C1r1, model="iid")+
  f(C2r1, model="iid")+
  f(beta1, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(beta2, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))

Jointmodel2 <- inla(formula.model2, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                   control.family = list(
                     list(control.link = list(model = "quantile", quantile = 0.6)),
                     list(),
                     list(),
                     list()
                   ),
                   control.compute = list(dic=TRUE, waic=TRUE),
                   data = jointdata, verbose = TRUE
)

summary(Jointmodel2)


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
CD4min_975th_percentile <- quantile(dataL$CD4min, 0.975)
dataL <- dataL[dataL$CD4min > CD4min_975th_percentile, ]
dataL <- dataL[dataL$ptid %in% ptids_to_keep, ]
existing_ptids <- dataL$ptid
dataS <- dataS[dataS$ptid %in% existing_ptids, ]
dataL$CD4_9 <- dataL$CD4min - CD4min_975th_percentile

data1<-dataS
dataE1<-data1
dataE1$event<-dataE1$C
dataE1$event[dataE1$event!=1]<-0
dataE2<-data1
dataE2$event<-dataE2$C
dataE2$event[dataE2$event!=2]<-0
dataE2$event<-dataE2$event/2

nL<-nrow(dataL)
nS<-nrow(dataS)

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      age_confirmL=c(dataL$age_confirm,dataL$age_confirm, rep(0,nS),rep(0,nS)),
                      age_confirm1=c(rep(0,nL),rep(0,nL),dataS$age_confirm,rep(0,nS)),
                      age_confirm2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$age_confirm),
                      sexL=c(dataL$sex,dataL$sex,rep(0,nS),rep(0,nS)),
                      sex1=c(rep(0,nL),rep(0,nL),dataS$sex,rep(0,nS)),
                      sex2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$sex),
                      marriageL=c(dataL$marriage,dataL$marriage,rep(0,nS),rep(0,nS)),
                      marriage1=c(rep(0,nL),rep(0,nL),dataS$marriage,rep(0,nS)),
                      marriage2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$marriage),
                      routeL=c(dataL$route,dataL$route,rep(0,nS),rep(0,nS)),
                      route1=c(rep(0,nL),rep(0,nL),dataS$route,rep(0,nS)),
                      route2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$route),
                      confirm_artL=c(dataL$confirm_art,dataL$confirm_art,rep(0,nS),rep(0,nS)),
                      confirm_art1=c(rep(0,nL),rep(0,nL),dataS$confirm_art,rep(0,nS)),
                      confirm_art2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$confirm_art),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie),
                      ARTdyL=c(dataL$ARTdy,dataL$ARTdy,rep(0,nS),rep(0,nS)),
                      ARTdy1=c(rep(0,nL),rep(0,nL),dataS$ARTdy,rep(0,nS)),
                      ARTdy2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$ARTdy),
                      tL=c(dataL$t,dataL$t,rep(0,nS),rep(0,nS)),
                      t1=c(rep(0,nL),rep(0,nL),dataS$t,rep(0,nS)),
                      t2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$t),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie))

random.eff<-list(timeL =c(dataL$Time,dataL$Time,rep(NA,nS),rep(NA,nS)),
                 idareaS1=c(rep(NA,nL),rep(NA,nL),dataS$idarea,rep(NA,nS)),
                 idareaS2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$idarea),
                 linpredL=c(rep(NA,nL),dataL$ptid,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 beta1=c(rep(NA,nL),rep(NA,nL),dataS$ptid,rep(NA,nS)),
                 beta2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ptid),
                 Lr1 = c(dataL$ptid, dataL$ptid, rep(NA, nS), rep(NA, nS)),
                 C1r1 = c(rep(NA, nL), rep(NA, nL), dataS$ptid, rep(NA, nS)),
                 C2r1 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$ptid))

jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$CD4_9,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
dataS$Time <- as.numeric(as.factor(dataS$Time))
y.survC1 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nL),dataE1$event,rep(NA,nS)))
y.survC2 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataE2$event))
y.joint<-list(y.long,y.eta,y.survC1,y.survC2)
jointdata$Y=y.joint

formula.model3=Y~sex1+sex2+sexL+marriage1+marriage2+marriageL+
  route1+route2+routeL+age_confirm1+age_confirm2+age_confirmL+confirm_art1+confirm_art2+confirm_artL+tL+jie1+jie2+
  f(inla.group(timeL,n=50),model="rw2", scale.model = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(idareaS1, model = "bym2", graph = g, hyper = list
    (prec = list(prior = "pc.prec",param = c(0.5 / 0.31, 0.01)),
      phi = list(prior = "pc",param = c(0.1, 4 / 5))))+
  f(idareaS2, copy="idareaS1", hyper = list(beta = list(fixed = FALSE)))+
  f(linpredL, linpredL2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(Lr1, model="iid")+
  f(C1r1, model="iid")+
  f(C2r1, model="iid")+
  f(beta1, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(beta2, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))

Jointmodel3 <- inla(formula.model3, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                   control.family = list(
                     list(control.link = list(model = "quantile", quantile = 0.6)),
                     list(),
                     list(),
                     list()
                   ),
                   control.compute = list(dic=TRUE, waic=TRUE),
                   data = jointdata, verbose = TRUE
)

summary(Jointmodel3)


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
CD4min_99th_percentile <- quantile(dataL$CD4min, 0.99)
dataL <- dataL[dataL$CD4min > CD4min_99th_percentile, ]
dataL <- dataL[dataL$ptid %in% ptids_to_keep, ]
existing_ptids <- dataL$ptid
dataS <- dataS[dataS$ptid %in% existing_ptids, ]
dataL$CD4_9 <- dataL$CD4min - CD4min_99th_percentile

data1<-dataS
dataE1<-data1
dataE1$event<-dataE1$C
dataE1$event[dataE1$event!=1]<-0
dataE2<-data1
dataE2$event<-dataE2$C
dataE2$event[dataE2$event!=2]<-0
dataE2$event<-dataE2$event/2

nL<-nrow(dataL)
nS<-nrow(dataS)

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      age_confirmL=c(dataL$age_confirm,dataL$age_confirm, rep(0,nS),rep(0,nS)),
                      age_confirm1=c(rep(0,nL),rep(0,nL),dataS$age_confirm,rep(0,nS)),
                      age_confirm2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$age_confirm),
                      sexL=c(dataL$sex,dataL$sex,rep(0,nS),rep(0,nS)),
                      sex1=c(rep(0,nL),rep(0,nL),dataS$sex,rep(0,nS)),
                      sex2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$sex),
                      marriageL=c(dataL$marriage,dataL$marriage,rep(0,nS),rep(0,nS)),
                      marriage1=c(rep(0,nL),rep(0,nL),dataS$marriage,rep(0,nS)),
                      marriage2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$marriage),
                      routeL=c(dataL$route,dataL$route,rep(0,nS),rep(0,nS)),
                      route1=c(rep(0,nL),rep(0,nL),dataS$route,rep(0,nS)),
                      route2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$route),
                      confirm_artL=c(dataL$confirm_art,dataL$confirm_art,rep(0,nS),rep(0,nS)),
                      confirm_art1=c(rep(0,nL),rep(0,nL),dataS$confirm_art,rep(0,nS)),
                      confirm_art2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$confirm_art),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie),
                      ARTdyL=c(dataL$ARTdy,dataL$ARTdy,rep(0,nS),rep(0,nS)),
                      ARTdy1=c(rep(0,nL),rep(0,nL),dataS$ARTdy,rep(0,nS)),
                      ARTdy2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$ARTdy),
                      tL=c(dataL$t,dataL$t,rep(0,nS),rep(0,nS)),
                      t1=c(rep(0,nL),rep(0,nL),dataS$t,rep(0,nS)),
                      t2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$t),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie))

random.eff<-list(timeL =c(dataL$Time,dataL$Time,rep(NA,nS),rep(NA,nS)),
                 idareaS1=c(rep(NA,nL),rep(NA,nL),dataS$idarea,rep(NA,nS)),
                 idareaS2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$idarea),
                 linpredL=c(rep(NA,nL),dataL$ptid,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 beta1=c(rep(NA,nL),rep(NA,nL),dataS$ptid,rep(NA,nS)),
                 beta2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ptid),
                 Lr1 = c(dataL$ptid, dataL$ptid, rep(NA, nS), rep(NA, nS)),
                 C1r1 = c(rep(NA, nL), rep(NA, nL), dataS$ptid, rep(NA, nS)),
                 C2r1 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$ptid))

jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$CD4_9,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
dataS$Time <- as.numeric(as.factor(dataS$Time))
y.survC1 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nL),dataE1$event,rep(NA,nS)))
y.survC2 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataE2$event))
y.joint<-list(y.long,y.eta,y.survC1,y.survC2)
jointdata$Y=y.joint

formula.model4=Y~sex1+sex2+sexL+marriage1+marriage2+marriageL+
  route1+route2+routeL+age_confirm1+age_confirm2+age_confirmL+confirm_art1+confirm_art2+confirm_artL+tL+jie1+jie2+
  f(inla.group(timeL,n=50),model="rw2", scale.model = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(idareaS1, model = "bym2", graph = g, hyper = list
    (prec = list(prior = "pc.prec",param = c(0.5 / 0.31, 0.01)),
      phi = list(prior = "pc",param = c(0.1, 4 / 5))))+
  f(idareaS2, copy="idareaS1", hyper = list(beta = list(fixed = FALSE)))+
  f(linpredL, linpredL2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(Lr1, model="iid")+
  f(C1r1, model="iid")+
  f(C2r1, model="iid")+
  f(beta1, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(beta2, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))

Jointmodel4 <- inla(formula.model4, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                    control.family = list(
                      list(control.link = list(model = "quantile", quantile = 0.6)),
                      list(),
                      list(),
                      list()
                    ),
                    control.compute = list(dic=TRUE, waic=TRUE),
                    data = jointdata, verbose = TRUE
)

summary(Jointmodel4)


#Sensitivity analysis results of alternative PC priors

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
dataL <- dataL[dataL$ptid %in% ptids_to_keep, ]
existing_ptids <- dataL$ptid
dataS <- dataS[dataS$ptid %in% existing_ptids, ]
dataL$CD4_9 <- dataL$CD4min - CD4min_90th_percentile

data1<-dataS
dataE1<-data1
dataE1$event<-dataE1$C
dataE1$event[dataE1$event!=1]<-0
dataE2<-data1
dataE2$event<-dataE2$C
dataE2$event[dataE2$event!=2]<-0
dataE2$event<-dataE2$event/2

nL<-nrow(dataL)
nS<-nrow(dataS)

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      age_confirmL=c(dataL$age_confirm,dataL$age_confirm, rep(0,nS),rep(0,nS)),
                      age_confirm1=c(rep(0,nL),rep(0,nL),dataS$age_confirm,rep(0,nS)),
                      age_confirm2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$age_confirm),
                      sexL=c(dataL$sex,dataL$sex,rep(0,nS),rep(0,nS)),
                      sex1=c(rep(0,nL),rep(0,nL),dataS$sex,rep(0,nS)),
                      sex2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$sex),
                      marriageL=c(dataL$marriage,dataL$marriage,rep(0,nS),rep(0,nS)),
                      marriage1=c(rep(0,nL),rep(0,nL),dataS$marriage,rep(0,nS)),
                      marriage2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$marriage),
                      routeL=c(dataL$route,dataL$route,rep(0,nS),rep(0,nS)),
                      route1=c(rep(0,nL),rep(0,nL),dataS$route,rep(0,nS)),
                      route2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$route),
                      confirm_artL=c(dataL$confirm_art,dataL$confirm_art,rep(0,nS),rep(0,nS)),
                      confirm_art1=c(rep(0,nL),rep(0,nL),dataS$confirm_art,rep(0,nS)),
                      confirm_art2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$confirm_art),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie),
                      ARTdyL=c(dataL$ARTdy,dataL$ARTdy,rep(0,nS),rep(0,nS)),
                      ARTdy1=c(rep(0,nL),rep(0,nL),dataS$ARTdy,rep(0,nS)),
                      ARTdy2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$ARTdy),
                      tL=c(dataL$t,dataL$t,rep(0,nS),rep(0,nS)),
                      t1=c(rep(0,nL),rep(0,nL),dataS$t,rep(0,nS)),
                      t2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$t),
                      jieL=c(dataL$jie,dataL$jie,rep(0,nS),rep(0,nS)),
                      jie1=c(rep(0,nL),rep(0,nL),dataS$jie,rep(0,nS)),
                      jie2=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$jie))

random.eff<-list(timeL =c(dataL$Time,dataL$Time,rep(NA,nS),rep(NA,nS)),
                 idareaS1=c(rep(NA,nL),rep(NA,nL),dataS$idarea,rep(NA,nS)),
                 idareaS2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$idarea),
                 linpredL=c(rep(NA,nL),dataL$ptid,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 beta1=c(rep(NA,nL),rep(NA,nL),dataS$ptid,rep(NA,nS)),
                 beta2=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ptid),
                 Lr1 = c(dataL$ptid, dataL$ptid, rep(NA, nS), rep(NA, nS)),
                 C1r1 = c(rep(NA, nL), rep(NA, nL), dataS$ptid, rep(NA, nS)),
                 C2r1 = c(rep(NA, nL), rep(NA, nL), rep(NA, nS), dataS$ptid))

jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$CD4_9,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
dataS$Time <- as.numeric(as.factor(dataS$Time))
y.survC1 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nL),dataE1$event,rep(NA,nS)))
y.survC2 <- inla.surv(time = c(rep(NA, nL),rep(NA, nL),rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataE2$event))
y.joint<-list(y.long,y.eta,y.survC1,y.survC2)
jointdata$Y=y.joint

formula.model5 = Y ~ sex1 + sex2 + sexL + marriage1 + marriage2 + marriageL +
  route1 + route2 + routeL + age_confirm1 + age_confirm2 + age_confirmL +
  confirm_art1 + confirm_art2 + confirm_artL + tL + jie1 + jie2 +
  f(inla.group(timeL, n = 50), model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(idareaS1, model = "bym2", graph = g,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01)),   
      phi = list(prior = "pc", param = c(0.5, 0.5))       
    )) +
  f(idareaS2, copy = "idareaS1",
    hyper = list(beta = list(fixed = FALSE))) +
  f(linpredL, linpredL2, model = "iid",
    hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(Lr1, model = "iid") +
  f(C1r1, model = "iid") +
  f(C2r1, model = "iid") +
  f(beta1, copy = "linpredL", hyper = list(beta = list(fixed = FALSE))) +
  f(beta2, copy = "linpredL", hyper = list(beta = list(fixed = FALSE)))


Jointmodel5 <- inla(formula.model5, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                    control.family = list(
                      list(control.link = list(model = "quantile", quantile = 0.6)),
                      list(),
                      list(),
                      list()
                    ),
                    control.compute = list(dic = TRUE, waic = TRUE),
                    data = jointdata, verbose = TRUE
)

summary(Jointmodel5)


formula.model6 = Y ~ sex1 + sex2 + sexL + marriage1 + marriage2 + marriageL +
  route1 + route2 + routeL + age_confirm1 + age_confirm2 + age_confirmL + confirm_art1 + confirm_art2 + confirm_artL + 
  tL + jie1 + jie2 +
  f(inla.group(timeL, n = 50), model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.9, 0.01)))) +  # RW2 prior changed here
  f(idareaS1, model = "bym2", graph = g,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(0.5 / 0.31, 0.01)),
      phi = list(prior = "pc", param = c(0.1, 4 / 5))
    )) +
  f(idareaS2, copy = "idareaS1", hyper = list(beta = list(fixed = FALSE))) +
  f(linpredL, linpredL2, model = "iid", hyper = list(prec = list(initial = -6, fixed = TRUE))) +
  f(Lr1, model = "iid") +
  f(C1r1, model = "iid") +
  f(C2r1, model = "iid") +
  f(beta1, copy = "linpredL", hyper = list(beta = list(fixed = FALSE))) +
  f(beta2, copy = "linpredL", hyper = list(beta = list(fixed = FALSE)))

Jointmodel6 <- inla(formula.model6, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                    control.family = list(
                      list(control.link = list(model = "quantile", quantile = 0.6)),
                      list(), list(), list()),
                    control.compute = list(dic = TRUE, waic = TRUE),
                    data = jointdata, verbose = TRUE)

summary(Jointmodel6)


#Sensitivity analysis results of alternative outcome definition same to the Joint model

#Sensitivity analysis results of independent dataset
data(aids, package = "joineR")

dat <- aids
dat %>%
  count(id, name = "n") %>%
  filter(n %in% c(1, 2, 3, 4,5,6)) %>%
  count(n, name = "n_patients")

dat %>%
  filter(!is.na(CD4)) %>%
  count(id, name = "n_cd4") %>%
  summarise(
    min  = min(n_cd4),
    max  = max(n_cd4),
    mean = mean(n_cd4)
  )

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

keep_id <- names(which(table(dataL_gau$id) >= 2))
dataL_gp <- dataL_gau %>% filter(id %in% keep_id)
dataL_gau <- dataL_gp
dataS_gp <- dataS %>% filter(id %in% keep_id)
dataS_gau <- dataS_gp

max_logCD4 <- max(dataL_gp$logCD4, na.rm = TRUE)

dataL_gp <- dataL_gp %>%
  mutate(
    CD4min = max_logCD4 - logCD4
  )

u <- as.numeric(quantile(dataL_gp$CD4min, 0.90, na.rm = TRUE))

dataL_gp <- dataL_gp %>%
  mutate(excess = CD4min - u) %>%
  filter(excess > 0)
dataS_gp <- dataS_gp %>%
  filter(id %in% dataL_gp$id)

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
  AZTS     = c(rep(0, nL_gp), rep(0, nL_gp), dataS_gp$AZT_bin))

random_gp <- list(
  linpredL  = c(rep(NA, nL_gp), dataL_gp$id, rep(NA, nS_gp)),
  linpredL2 = c(rep(NA, nL_gp), rep(-1, nL_gp), rep(NA, nS_gp)),
  beta      = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$id))
joint_gp <- c(fixed_gp, random_gp)
y.long.gp <- c(dataL_gp$excess, rep(NA, nL_gp), rep(NA, nS_gp))
y.eta.gp  <- c(rep(NA, nL_gp), rep(0, nL_gp), rep(NA, nS_gp))
y.surv.gp <- inla.surv(
  time  = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$time),
  event = c(rep(NA, nL_gp), rep(NA, nL_gp), dataS_gp$death))
joint_gp$Y <- list(y.long.gp, y.eta.gp, y.surv.gp)

formula_gp <- Y ~
  drugL + femaleL + prevOIL + AZTL +
  drugS + femaleS + prevOIS + AZTS +
  f(inla.group(obstimeL, n = 50), model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
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

nL_gau <- nrow(dataL_gau)
nS_gau <- nrow(dataS_gau)

fixed_gau <- data.frame(
  mu = factor(c(rep(1, nL_gau), rep(1, nL_gau), rep(2, nS_gau))),
  obstimeL = c(dataL_gau$obstime, dataL_gau$obstime, rep(0, nS_gau)),
  drugL    = c(dataL_gau$drug_bin,   dataL_gau$drug_bin,   rep(0, nS_gau)),
  femaleL  = c(dataL_gau$female,     dataL_gau$female,     rep(0, nS_gau)),
  prevOIL  = c(dataL_gau$prevOI_bin, dataL_gau$prevOI_bin, rep(0, nS_gau)),
  AZTL     = c(dataL_gau$AZT_bin,    dataL_gau$AZT_bin,    rep(0, nS_gau)),
  drugS    = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau$drug_bin),
  femaleS  = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau$female),
  prevOIS  = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau$prevOI_bin),
  AZTS     = c(rep(0, nL_gau), rep(0, nL_gau), dataS_gau$AZT_bin))

random_gau <- list(
  linpredL  = c(rep(NA, nL_gau), dataL_gau$id, rep(NA, nS_gau)),
  linpredL2 = c(rep(NA, nL_gau), rep(-1, nL_gau), rep(NA, nS_gau)),
  beta      = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau$id))
joint_gau <- c(fixed_gau, random_gau)
y.long.gau <- c(dataL_gau$logCD4, rep(NA, nL_gau), rep(NA, nS_gau))
y.eta.gau  <- c(rep(NA, nL_gau), rep(0, nL_gau), rep(NA, nS_gau))
y.surv.gau <- inla.surv(
  time  = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau$time),
  event = c(rep(NA, nL_gau), rep(NA, nL_gau), dataS_gau$death))
joint_gau$Y <- list(y.long.gau, y.eta.gau, y.surv.gau)

formula_gau <- Y ~
  drugL + femaleL + prevOIL + AZTL +
  drugS + femaleS + prevOIS + AZTS +
  f(inla.group(obstimeL, n = 50), model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
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

marker_gp  <- fit_gp$summary.fitted.values[idx_gp_surv,  "mean"]
marker_gau <- fit_gau$summary.fitted.values[idx_gau_surv, "mean"]

stopifnot(length(marker_gp)  == nrow(dataS_gp))
stopifnot(length(marker_gau) == nrow(dataS_gau))

## GP ROC
roc_time_gp <- as.numeric(
  quantile(dataS_gp$time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
)

roc_gp <- data.frame(
  Cutoff = numeric(),
  FP = I(list()),
  TP = I(list()),
  AUC = numeric()
)

for (cutoff in roc_time_gp) {
  ROC <- survivalROC(
    Stime = dataS_gp$time,
    status = dataS_gp$death,
    marker = marker_gp,
    predict.time = cutoff,
    method = "KM"
  )
  
  roc_gp <- rbind(
    roc_gp,
    data.frame(
      Cutoff = cutoff,
      FP = I(list(ROC$FP)),
      TP = I(list(ROC$TP)),
      AUC = ROC$AUC
    )
  )
}

print(roc_gp[, c("Cutoff", "AUC")])

trapz_mean <- function(x, y) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / 
    (max(x) - min(x))
}

trapz_mean <- function(x, y) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / 
    (max(x) - min(x))
}
iAUC_gp <- trapz_mean(roc_gp$Cutoff, roc_gp$AUC)

## Gaussian ROC
roc_time_gau <- as.numeric(
  quantile(dataS_gau$time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
)

roc_gau <- data.frame(
  Cutoff = numeric(),
  FP = I(list()),
  TP = I(list()),
  AUC = numeric()
)

for (cutoff in roc_time_gau) {
  ROC <- survivalROC(
    Stime = dataS_gau$time,
    status = dataS_gau$death,
    marker = marker_gau,
    predict.time = cutoff,
    method = "KM"
  )
  
  roc_gau <- rbind(
    roc_gau,
    data.frame(
      Cutoff = cutoff,
      FP = I(list(ROC$FP)),
      TP = I(list(ROC$TP)),
      AUC = ROC$AUC
    )
  )
}

print(roc_gau[, c("Cutoff", "AUC")])
mean(roc_gau$AUC, na.rm = TRUE)
mean(roc_gp$AUC, na.rm = TRUE)


## Integrated Brier Score, IBS

get_alpha <- function(fit) {
  hp <- fit$summary.hyperpar
  idx <- grep("alpha parameter for weibullsurv", rownames(hp))
  hp[idx[1], "mean"]
}

pred_event_weibull <- function(eta, alpha, times) {
  lambda <- exp(eta)
  outer(lambda, times^alpha, function(lam, ta) {
    1 - exp(-lam * ta)
  })
}

G_eval <- function(t, kmfit) {
  sf <- summary(kmfit, times = t, extend = TRUE)$surv
  ifelse(is.na(sf) | sf <= 0, NA_real_, sf)
}

brier_ipcw_event <- function(time, status, pred_mat, times) {
  status <- as.integer(status)
  km_cens <- survfit(Surv(time, 1 - status) ~ 1)
  bs <- numeric(length(times))
  for (j in seq_along(times)) {
    tt <- times[j]
    pred <- pred_mat[, j]
    y <- as.numeric(time <= tt & status == 1)
    w <- rep(NA_real_, length(time))
    idx_event <- time <= tt & status == 1
    idx_alive <- time > tt
    G_event <- G_eval(time[idx_event], km_cens)
    G_alive <- G_eval(tt, km_cens)
    w[idx_event] <- ifelse(!is.na(G_event), 1 / G_event, NA_real_)
    w[idx_alive] <- ifelse(!is.na(G_alive), 1 / G_alive, NA_real_)
    bs[j] <- mean(w * (y - pred)^2, na.rm = TRUE)
  }
  bs
}

ibs_trapz <- function(times, bs) {
  o <- order(times)
  times <- times[o]
  bs <- bs[o]
  sum(diff(times) * (head(bs, -1) + tail(bs, -1)) / 2) /
    (max(times) - min(times))
}

eta_gp <- fit_gp$summary.linear.predictor[idx_gp_surv, "mean"]
alpha_gp <- get_alpha(fit_gp)

time_gp <- as.numeric(
  quantile(dataS_gp$time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
)

pred_gp <- pred_event_weibull(
  eta = eta_gp,
  alpha = alpha_gp,
  times = time_gp
)

bs_gp <- brier_ipcw_event(
  time = dataS_gp$time,
  status = dataS_gp$death,
  pred_mat = pred_gp,
  times = time_gp
)

IBS_gp <- ibs_trapz(time_gp, bs_gp)

bs_gp_res <- data.frame(
  Model = "GP",
  Time = time_gp,
  Brier = bs_gp
)

print(bs_gp_res)
cat("IBS_GP =", IBS_gp, "\n")

eta_gau <- fit_gau$summary.linear.predictor[idx_gau_surv, "mean"]
alpha_gau <- get_alpha(fit_gau)

time_gau <- as.numeric(
  quantile(dataS_gau$time, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
)

pred_gau <- pred_event_weibull(
  eta = eta_gau,
  alpha = alpha_gau,
  times = time_gau
)

bs_gau <- brier_ipcw_event(
  time = dataS_gau$time,
  status = dataS_gau$death,
  pred_mat = pred_gau,
  times = time_gau
)

IBS_gau <- ibs_trapz(time_gau, bs_gau)

bs_gau_res <- data.frame(
  Model = "Gaussian",
  Time = time_gau,
  Brier = bs_gau
)

print(bs_gau_res)
cat("IBS_Gaussian =", IBS_gau, "\n")

