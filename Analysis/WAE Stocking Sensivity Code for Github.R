library(magrittr)
library(brms)
library(Matrix)
library(tidybayes)
library(gridExtra)
library(tidyverse)

############################################################
###### Send request for access to data
wae3<-read.csv("wae.final.csv",header=T)

############################################################
###### Stocked vs non-stocked model
#model
stock_mod<- brm(bf(CPUE2~Stocked*Type + (1|Waterbody2),
                   hu~Stocked*Type + (1|Waterbody2)),
                data=wae3, 
                prior = c(prior(normal(1,.5), class = "Intercept"),
                          prior(normal(-1,0.5), coef="TypeMarginal"),
                          prior(exponential(2), class="sd"),
                          prior(normal(1, .5), class = "b"),
                          prior(exponential(2), class="shape"), 
                          prior(normal(0, 2), class="b", dpar="hu")),
                family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
                #sample_prior = "only",
                file="Stocked_Type.rds",
                chains=4, iter=2000, cores=4)

#model summary and checks
print(stock_mod, digits=4)
bayes_R2(stock_mod)
pp_check(stock_mod)

#Senistivity
stock_mod.s<- brm(bf(CPUE2~Stocked*Type + (1|Waterbody2),
               hu~Stocked*Type + (1|Waterbody2)),
            data=wae3, 
            prior = c(prior(normal(1,1), class = "Intercept"),
                      prior(normal(-1,1), coef="TypeMarginal"),
                      prior(exponential(1), class="sd"),
                      prior(normal(1, 1), class = "b"),
                      prior(exponential(1), class="shape"),
                      prior(normal(0, 4), class="b", dpar="hu")),
            family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
            #sample_prior = "only",
            file="Stocked_Type_Sens.rds",
            chains=4, iter=2000, cores=4)

#model summary and checks
print(stock_mod.s, digits=4)
bayes_R2(stock_mod.s)
pp_check(stock_mod.s)

############################################################
###### Size  of stocked products * waterbody type
st <- filter(wae3, Size2!="a_no")

#model
S_T <- brm(bf(CPUE2~Size2*Type + (1|Waterbody2),
              hu~Size2*Type + (1|Waterbody2)),
           data=st, 
           prior = c(prior(normal(0,.5), class="b"),
                     prior(normal(-1,0.5), coef="TypeMarginal"),
                     prior(exponential(2), class="sd"),
                     prior(normal(1.5, 0.5), class="Intercept"),
                     prior(exponential(2), class = "shape"),
                     prior(normal(0, 2), class="b", dpar="hu")),
           family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
           #sample_prior = "only",
           file="Size_Type.rds", 
           chains=4, iter=2000, cores=4)

#model summary and check
print(S_T, digits=4)
bayes_R2(S_T)
pp_check(S_T)
plot(S_T)

#model sensitivity
S_T.s<- brm(bf(CPUE2~Size2*Type + (1|Waterbody2),
             hu~Size2*Type + (1|Waterbody2)),
          data=st, 
          prior = c(prior(normal(0,1), class="b"),
                    prior(normal(-1,1), coef="TypeMarginal"),
                    prior(exponential(1), class="sd"),
                    prior(normal(1.5, 1), class="Intercept"),
                    prior(exponential(1), class = "shape"),
                    prior(normal(0, 4), class="b", dpar="hu")),
          family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
          #sample_prior = "only",
          file="Size_Type_sens.rds", 
          chains=4, iter=2000, cores=4)

#model summary and check
print(S_T.s, digits=4)
bayes_R2(S_T.s)
pp_check(S_T.s)

############################################################
###### Management purpose
mgmt<-filter(wae3, ManagementPurpose !="None")

# model
mp<-brm(bf(CPUE2 ~ ManagementPurpose + (1|Waterbody2),
           hu ~ ManagementPurpose + (1|Waterbody2)),
        data=mgmt, 
        family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
        prior = c(prior(normal(0,.5), class="b"),
                  prior(exponential(0.1), class="sd"),
                  prior(normal(1.5, .5), class="Intercept"),
                  prior(exponential(0.5), class = "shape"), 
                  prior(normal(0,2), class="b", dpar="hu")),
        #sample_prior = "only",
        file="Management_Purpose.rds",
        chains=4, iter=2000, cores=4)

#model summary and check
print(mp, digits=4)
bayes_R2(mp)
pp_check(mp)

# model sensitivity
mp.s<-brm(bf(CPUE2 ~ ManagementPurpose + (1|Waterbody2),
           hu ~ ManagementPurpose + (1|Waterbody2)),
        data=mgmt, 
        family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
        prior = c(prior(normal(0,1), class="b"),
                  prior(exponential(0.05), class="sd"),
                  prior(normal(1.5, 1), class="Intercept"),
                  prior(exponential(0.25), class = "shape"), 
                  prior(normal(0,4), class="b", dpar="hu")),
        #sample_prior = "only",
        file="Management_Purpose_sens.rds",
        chains=4, iter=2000, cores=4)

#model summary and check
print(mp.s, digits=4)
bayes_R2(mp.s)
pp_check(mp.s)

############################################################
###### Waterbody Characteristics
#model
water_mod<-brm(bf(CPUE2 ~ Hectares_std + WaterbodyType  + (1|Waterbody2),
                  hu ~ Hectares_std + WaterbodyType   + (1|Waterbody2)),
               data=wae3, 
               family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
               prior = c(prior(normal(0,0.5), class="b"),
                         prior(exponential(0.1), class="sd"),
                         prior(normal(1, .5), class="Intercept"),
                         prior(exponential(0.5), class = "shape"), 
                         prior(normal(0,2), class="b", dpar="hu")),
               #sample_prior = "only",
               file="Waterbody_Char.rds",
               chains=4, iter=2000, cores=4)

#summary and model check
print(water_mod, digits=4)
bayes_R2(water_mod)
pp_check(water_mod)

#model sensitivity
water_mod.s<-brm(bf(CPUE2 ~ Hectares_std + WaterbodyType  + (1|Waterbody2),
                  hu ~ Hectares_std + WaterbodyType   + (1|Waterbody2)),
               data=wae3, 
               family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
               prior = c(prior(normal(0,1), class="b"),
                         prior(exponential(0.05), class="sd"),
                         prior(normal(1, 1), class="Intercept"),
                         prior(exponential(0.25), class = "shape"), 
                         prior(normal(0,4), class="b", dpar="hu")),
               #sample_prior = "only",
               file="Waterbody_Char_sens.rds",
               chains=4, iter=2000, cores=4)

#summary and model check
print(water_mod.s, digits=4)
bayes_R2(water_mod.s)
pp_check(water_mod.s)

############################################################
######Environmental variables
#model
env<-brm(bf(CPUE2~ GDD_std  + SpringPrecip_std + WSI_std + (1|Waterbody2),
            hu~ GDD_std  + SpringPrecip_std + WSI_std + (1|Waterbody2)),
         data=wae3, 
         prior = c(prior(normal(1,.5), coef="WSI_std"),
                   prior(normal(1,.5), coef="GDD_std"),
                   prior(normal(1,.5), coef="SpringPrecip_std"),
                   prior(exponential(2), class="sd"),
                   prior(normal(1, 0.5), class="Intercept"),
                   prior(exponential(2), class = "shape"), 
                   prior(normal(0,2), class="b", dpar="hu")),
         family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
         #sample_prior = "only",
         file="Env_mod.rds",
         chains=4, iter=2000, cores=4)

#model summary and checks
print(env, digits=4)
bayes_R2(env)
pp_check(env)
plot(env)
 
#model sensitivity
env.s<-brm(bf(CPUE2~ GDD_std  + SpringPrecip_std + WSI_std + (1|Waterbody2),
            hu~ GDD_std  + SpringPrecip_std + WSI_std + (1|Waterbody2)),
         data=wae3, 
         prior = c(prior(normal(1,1), coef="WSI_std"),
                   prior(normal(1,1), coef="GDD_std"),
                   prior(normal(1,1), coef="SpringPrecip_std"),
                   prior(exponential(1), class="sd"),
                   prior(normal(1, 1), class="Intercept"),
                   prior(exponential(1), class = "shape"), 
                   prior(normal(0,4), class="b", dpar="hu")),
         family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
         #sample_prior = "only",
         file="Env_mod_sens.rds",
         chains=4, iter=2000, cores=4)

#model summary and checks
print(env.s, digits=4)
bayes_R2(env.s)
pp_check(env.s)
plot(env.s)

############################################################
####### Stock-size Walleye
#filter data
adultWAE<-wae3 %>% filter(!is.na(WAE.S)) %>%   
  filter(Stocked=="N")%>% 
  mutate(WAE_adj=(WAE.S-mean(WAE.S))/sd(WAE.S)) 


#model
wae_mod<-brm(bf(CPUE2~ WAE_adj + (1|Waterbody2),
                hu~WAE_adj + (1|Waterbody2)),
             data=adultWAE, 
             prior = c(prior(normal(0,1), class="b"),
                       prior(exponential(0.5), class="sd"),
                       prior(normal(1.5, 0.5), class="Intercept"),
                       prior(exponential(0.5), class = "shape"), 
                       prior(normal(0,2), class="b", dpar="hu")),
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             file="Adult_WAE.rds",
             chains=4, iter=2000, cores=4)


#model summary and checks
print(wae_mod, digits=4)
bayes_R2(wae_mod)
pp_check(wae_mod)

#model sensitivity
wae_mod.s<-brm(bf(CPUE2~ WAE_adj + (1|Waterbody2),
                hu~WAE_adj + (1|Waterbody2)),
             data=adultWAE, 
             prior = c(prior(normal(0,2), class="b"),
                       prior(exponential(0.25), class="sd"),
                       prior(normal(1.5, 1), class="Intercept"),
                       prior(exponential(0.25), class = "shape"), 
                       prior(normal(0,4), class="b", dpar="hu")),
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             #sample_prior = "only",
             file="Adult_WAE_sens.rds",
             chains=4, iter=2000, cores=4)


#model summary and checks
print(wae_mod.s, digits=4)
bayes_R2(wae_mod.s)
pp_check(wae_mod.s)

############################################################
##### Centrarchids
# Filter data to only include waters with sampled centrarchids
cent<-wae3 %>% filter(!is.na(Centrarchids)) %>% 
  mutate(Centrarchids_std=(Centrarchids-mean(Centrarchids))/sd(Centrarchids))

#model
cent_mod<-brm(bf(CPUE2 ~ Centrarchids_std + (1|Waterbody2),
                 hu~Centrarchids_std + (1|Waterbody2)),
              data=cent, 
              prior = c(prior(normal(-1,1), class="b"),
                        prior(exponential(0.5), class="sd"),
                        prior(normal(1.5, 0.5), class="Intercept"),
                        prior(exponential(0.5), class = "shape"), 
                        prior(normal(0,2), class="b", dpar="hu")),
              family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
              #sample_prior = "only",
              file="Cent_mod.rds",
              chains=4, iter=2000, cores=4)

#model summary and checks
print(cent_mod, digits=4)
bayes_R2(cent_mod)
pp_check(cent_mod)

#model sensitivity
cent_mod.s<-brm(bf(CPUE2~ Centrarchids_std + (1|Waterbody2),
                 hu~Centrarchids_std + (1|Waterbody2)),
              data=cent, 
              prior = c(prior(normal(-1,2), class="b"),
                        prior(exponential(0.25), class="sd"),
                        prior(normal(1.5, 1), class="Intercept"),
                        prior(exponential(0.25), class = "shape"), 
                        prior(normal(0,4), class="b", dpar="hu")),
              family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
              #sample_prior = "only",
              file="Cent_mod_sens.rds",
              chains=4, iter=2000, cores=4)

#model summary and checks
print(cent_mod.s, digits=4)
bayes_R2(cent_mod.s)
pp_check(cent_mod.s)

###########################################################
##### Northern Pike
# Filter data to only include waters with sampled Northern Pike
nop<-wae3 %>% filter(!is.na(NOP.All)) %>% 
  mutate(NOP.All_std=(NOP.All-mean(NOP.All))/sd(NOP.All),
         NOP.S_std=(NOP.S-mean(NOP.S))/sd(NOP.S))

#model
nop_mod<-brm(bf(CPUE2~NOP.S_std + (1|Waterbody2),
                hu~NOP.S_std + (1|Waterbody2)),
             data=nop, 
             prior = c(prior(normal(-1,1), class="b"),
                       prior(exponential(0.5), class="sd"),
                       prior(normal(1.5, 0.5), class="Intercept"),
                       prior(exponential(0.5), class = "shape"), 
                       prior(normal(0,2), class="b", dpar="hu")),
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             #sample_prior = "only",
             file="NOP_mod.rds",
             chains=4, iter=2000, cores=4)

#model summary and checks
summary(nop_mod)
bayes_R2(nop_mod)
pp_check(nop_mod)

#model sensitivity
nop_mod.s<-brm(bf(CPUE2~NOP.S_std + (1|Waterbody2),
                hu~NOP.S_std + (1|Waterbody2)),
             data=nop, 
             prior = c(prior(normal(-1,2), class="b"),
                       prior(exponential(0.25), class="sd"),
                       prior(normal(1.5, 1), class="Intercept"),
                       prior(exponential(0.25), class = "shape"), 
                       prior(normal(0,4), class="b", dpar="hu")),
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             #sample_prior = "only",
             file="NOP_mod_sens.rds",
             chains=4, iter=2000, cores=4)

#model summary and checks
summary(nop_mod.s)
bayes_R2(nop_mod.s)
pp_check(nop_mod.s)

############################################################
##### Saugeye lakes
#filter data
sae<-read.csv("wae2.1.csv", header=T)

SAE<-sae %>% filter(Waterbody=="Campbell"|Waterbody=="Elm"|Waterbody=="Goldsmith"|
                      Waterbody=="Mina"|Waterbody=="Richmond"|Waterbody=="White") %>% 
  mutate(SppStocked = plyr::mapvalues(SpeciesStocked, from=c("None", "Walleye", "Saugeye"), 
                                      to=c("a_none", "b_Walleye", "c_Saugeye")))

#model
sae_mod<-brm(bf(CPUE2 ~ SppStocked + (1|Waterbody2),
                hu ~ SppStocked + (1|Waterbody2)),
             data=SAE, 
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             prior = c(prior(normal(0.8, 0.4), class="b"),
                       prior(exponential(0.2), class="sd"),
                       prior(normal(1, .5), class="Intercept"),
                       prior(exponential(0.5), class = "shape"), 
                       prior(normal(0,2), class="b", dpar="hu")),
             #sample_prior = "only",
             file="SAE_mod.rds",
             chains=4, iter=2000, cores=4)

#model summary and checks
print(sae_mod, digits=4)
bayes_R2(sae_mod)
pp_check(sae_mod)

#model sensitivity
sae_mod.s<-brm(bf(CPUE2 ~ SppStocked + (1|Waterbody2),
                hu ~ SppStocked + (1|Waterbody2)),
             data=SAE, 
             family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
             prior = c(prior(normal(0.8, 0.8), class="b"),
                       prior(exponential(0.1), class="sd"),
                       prior(normal(1, 1), class="Intercept"),
                       prior(exponential(0.25), class = "shape"), 
                       prior(normal(0,4), class="b", dpar="hu")),
             #sample_prior = "only",
             file="SAE_mod_sens.rds",
             chains=4, iter=2000, cores=4)

#model summary and checks
print(sae_mod.s, digits=4)
bayes_R2(sae_mod.s)
pp_check(sae_mod.s)

