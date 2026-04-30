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
jointdata$jie1 <- factor(jointdata$jie1)
jointdata$jie2 <- factor(jointdata$jie2)


formula.model=Y~sex1+sex2+sexL+marriage1+marriage2+marriageL+
  route1+route2+routeL+age_confirm1+age_confirm2+age_confirmL+confirm_art1+confirm_art2+confirm_artL + tL + jie1 + jie2 +
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

#GP
Jointmodel <- inla(formula.model, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                   control.family = list(
                     list(control.link = list(model = "quantile", quantile = 0.6)),
                     list(),
                     list(),
                     list()
                   ),
                   control.compute = list(dic=TRUE, waic=TRUE),
                   data = jointdata, verbose = TRUE
)

summary(Jointmodel)


#Gaussian
dataL$logcd4 <- log(dataL$cd4)
dataL$CD4_9 <- dataL$logcd4
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
jointdata$jie1 <- factor(jointdata$jie1)
jointdata$jie2 <- factor(jointdata$jie2)


formula.model=Y~sex1+sex2+sexL+marriage1+marriage2+marriageL+
  route1+route2+routeL+age_confirm1+age_confirm2+age_confirmL+confirm_art1+confirm_art2+confirm_artL + tL + jie1 + jie2 +
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


Jointmodel1 <- inla(formula.model, family = c('gaussian','gaussian', 'weibullsurv', 'weibullsurv'),
                    control.family = list(
                      list(),
                      list(),
                      list(),
                      list()
                    ),
                    control.inla = list(
                      strategy = "simplified.laplace"
                    ),
                    control.compute = list(dic=TRUE, waic=TRUE),
                    data = jointdata, verbose = TRUE
)

summary(Jointmodel1)


#Sensitivity analysis
formula.model2 = Y ~ sex1 + sex2 + sexL + marriage1 + marriage2 + marriageL +
  route1 + route2 + routeL + age_confirm1 + age_confirm2 + age_confirmL +
  confirm_art1 + confirm_art2 + confirm_artL + tL + jie1 + jie2 +
  f(inla.group(timeL, n = 50), model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(idareaS1, model = "bym2", graph = g,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01)),   # 改成更宽松的先验
      phi = list(prior = "pc", param = c(0.5, 0.5))        # 改成更均匀的先验
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

Jointmodel2 <- inla(formula.model2, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                   control.family = list(
                     list(control.link = list(model = "quantile", quantile = 0.6)),
                     list(),
                     list(),
                     list()
                   ),
                   control.compute = list(dic = TRUE, waic = TRUE),
                   data = jointdata, verbose = TRUE
)

summary(Jointmodel2)


formula.model3 = Y ~ sex1 + sex2 + sexL + marriage1 + marriage2 + marriageL +
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

Jointmodel3 <- inla(formula.model3, family = c('gp','gaussian', 'weibullsurv', 'weibullsurv'),
                   control.family = list(
                     list(control.link = list(model = "quantile", quantile = 0.6)),
                     list(), list(), list()),
                   control.compute = list(dic = TRUE, waic = TRUE),
                   data = jointdata, verbose = TRUE)

summary(Jointmodel3)
