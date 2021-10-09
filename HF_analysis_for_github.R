setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Previous_research/Harvard Forest/New")
HF<-read.csv("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Previous_research/Harvard Forest/New/HF_data.csv")
AD<-read.csv("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Previous_research/Harvard Forest/New/summary_chambers_summer_2014.csv")

library(dplyr)
library(MuMIn)
library(lme4)
library(car)
library(vegan)
library(ggplot2)
library(tidyr)
library(emmeans)
library(ggsignif)
library(piecewiseSEM)
library(DHARMa)

HF$chamber<-as.factor(HF$chamber)
HF$bag<-as.factor(HF$bag)
HF$treatment<-as.factor(HF$treatment)
HF$abundance<-rowSums(HF[,c(12:21)])
HF$richness<-apply(HF[,c(12:19)]>0,1,sum)

#total abundance and richness
x<-sum(colSums(select(HF, Acari:unknown))) #abundance 3662
y<-(colSums(select(HF, Acari:unknown))) #12 groups

#the HF dataset has two measurements for air temperature, 4 for solid temp (2 organic and 2 inorganic), will average and use the mean

AD1<- AD %>%
  mutate(avg.air.temp=rowMeans(cbind(cat1.avg, cat2.avg)), na.rm=T) %>%
  mutate(avg.soil.temp=rowMeans(cbind(csto1.avg, csto2.avg, csti1.avg, csti2.avg)), na.rm=T)

AD1<-AD1[,c("chamber","csm.avg", "avg.air.temp", "avg.soil.temp" )]
  
AD2<-merge(HF,AD1, by="chamber")

#Questions to address
## 1) Effect of invertebrate presence (treatment), air temperature, soil temperature, soil moisture, and interactions on decomposition
### Given all of the factors, need to simplify model 

#First you have to fit a linear mixed effect model with all your main factors and interactions.
## Check variance inflation (colinearity)
m1<-lmer(decomposition~csm.avg+avg.air.temp+avg.soil.temp+treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)
vif(m1) # VIF >2 for soil and air temperature, will remove soil temp 
m2<-lmer(decomposition~ csm.avg + avg.air.temp + treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)
vif(m2) # VIF <2
## now including interactions
m3<-lmer(decomposition~ csm.avg*avg.air.temp*treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)


#Then you run your dredge analysis on the model. In theory you can just write "dredge(m1)". I added R2 values to know how much variance
#the proposed model explains. Run everything (all 7 rows) to get a result
x<-dredge(m3, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  })
)

#write.csv(x, "dredge.table.csv")

#Now you get a list of models. The model that explains more variance is on the top. 
p1<-lmer(decomposition~avg.air.temp+treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
#Here we test if the model violates any assumptions of normality or homogeneity of variance
E1<-resid(p1, type="pearson")
hist(E1)
plot(resid(p1))
qqnorm(resid(p1))
simoutbin500<-simulateResiduals(fittedModel=p1, n=250)
plot(simoutbin500) # Looks good enough to continue without transformation and finish the analysis
p2<-lmer(decomposition~treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p1,p2) # avg. air temp significant
p3<-lmer(decomposition~avg.air.temp + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p1,p3) #treatment also significant
#comparisons among treatments
p1.emm.s <- emmeans(p1, "treatment")
pairs(p1.emm.s)

## 2) Effect of invertebrate presence (treatment), air temperature, soil temperature, soil moisture, and interactions on invertebrate abundance

## Check variance inflation (colinearity)
m4<-lmer(log(abundance)~ csm.avg + avg.air.temp + treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)
vif(m4) # VIF <2
## now including interactions
m5<-lmer(log(abundance)~ csm.avg*avg.air.temp*treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)

#Then you run your dredge analysis on the model. In theory you can just write "dredge(m1)". I added R2 values to know how much variance
#the proposed model explains. Run everything (all 7 rows) to get a result
y<-dredge(m5, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  })
)

#Now you get a list of models. The model that explains more variance is on the top. 
p4<-lmer(log(abundance)~avg.air.temp + treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
#Here we test if the model violates any assumptions of normality or homogeneity of variance
E2<-resid(p4, type="pearson")
hist(E2)
plot(resid(p5))
qqnorm(resid(p5))
simoutbin500<-simulateResiduals(fittedModel=p4, n=250)
plot(simoutbin500) # Looks good enough to continue without transformation and finish the analysis
p5<-lmer(log(abundance)~avg.air.temp + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p4,p5) #treatment very significant 
#comparisons among treatments
p2.emm.s <- emmeans(p4, "treatment")
pairs(p2.emm.s)
p6<-lmer(log(abundance)~ treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p4,p6) # average air temperature significant

## 3) Effect of invertebrate presence (treatment), air temperature, soil temperature, soil moisture, and interactions on invertebrate richness
#First you have to fit a linear mixed effect model with all your main factors and interactions.
## Check variance inflation (colinearity)
m6<-lmer(richness~ csm.avg + avg.air.temp + treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)
vif(m6) # VIF <2
## now including interactions
m7<-lmer(richness~ csm.avg*avg.air.temp*treatment + set + (1|chamber), data=AD2,na.action="na.fail", REML=F)

#Then you run your dredge analysis on the model. In theory you can just write "dredge(m1)". I added R2 values to know how much variance
#the proposed model explains. Run everything (all 7 rows) to get a result
z<-dredge(m7, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  })
)

#Now you get a list of models. The model that explains more variance is on the top. 
p7<-lmer(richness~avg.air.temp + csm.avg + treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
#Here we test if the model violates any assumptions of normality or homogeneity of variance
E3<-resid(p7, type="pearson")
hist(E3)
plot(resid(p7))
qqnorm(resid(p7))
simoutbin500<-simulateResiduals(fittedModel=p7, n=250)
plot(simoutbin500) # Looks good enough to continue without transformation and finish the analysis

p8<-lmer(richness~ avg.air.temp + treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p7,p8) # soil moisture marginally significant
p9<-lmer(richness~avg.air.temp + csm.avg + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p7,p9) #treatment very significant 
#comparisons among treatments
p3.emm.s <- emmeans(p7, "treatment")
pairs(p3.emm.s)
p10<-lmer(richness~ csm.avg + treatment + (1|chamber), data=AD2,na.action="na.fail", REML=F)
anova(p7,p10) # average air temperature significant

## is the community composition different by treatment and temperature as well?
## number of invertebrates counted, richness
ABND <- HF %>%
  summarise_at(vars(Acari:unknown), sum) %>%
  #rowwise() %>%
  #mutate(total=sum(c_across(Acari:unknown))) %>% #3662 invertebrates counted
  gather(taxa, count, Acari:unknown, factor_key=T)

#write.csv(ABND, "taxa.HF.csv")

#first will group by treatment and verify that communities are indeed different (treatment was effective)
AD3 <- AD2 %>%
  group_by(chamber,treatment) %>%
  summarise_at(vars(Acari:Thysanoptera), sum, na.rm=T)
AD3<-AD3[,c(-1,-2)]

AD4 <- AD2 %>%
  group_by(chamber,treatment) %>%
  summarise_at(vars(avg.air.temp, csm.avg), mean, na.rm=T)

AD4$chamber.treatment<-paste(AD4$chamber,AD4$treatment)

adonis(AD3~ treatment*avg.air.temp*csm.avg, data=AD4, perm=9999)

summary(lm(AD1$csm.avg~AD1$avg.air.temp))

#Structural equations model
mod.list=psem(lmer(log(abundance)~ avg.air.temp +  (1|chamber), data=AD2,na.action="na.fail"),
              lmer(richness~ avg.air.temp +  (1|chamber), data=AD2,na.action="na.fail"),
              lmer(decomposition~avg.air.temp + log(abundance) + richness + (1|chamber), data=AD2,na.action="na.fail"),
              richness %~~% log(abundance),
              data=AD2)

summary(mod.list)
