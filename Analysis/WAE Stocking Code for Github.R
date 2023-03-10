library(magrittr)
library(brms)
library(Matrix)
library(tidybayes)
library(gridExtra)
library(grid)
library(ggpubr)
library(tidyverse)

############################################################
###### Send request for access to data
wae2<-read.csv("wae2.csv",header=T)
wae3<-droplevels(filter(wae2,SpeciesStocked != "Saugeye"))

############################################################
###### Stocked vs non-stocked model
#identify priors
get_prior(CPUE_adj~Stocked + (1|Waterbody2),
          data=wae3, family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"))

#model
stock_mod<- brm(bf(CPUE2~Stocked + (1|Waterbody2),
               hu~Stocked + (1|Waterbody2)),
            data=wae3, 
            prior = c(prior(normal(1,.5), class = "Intercept"),
                      prior(exponential(2), class="sd"),
                      prior(normal(1, .5), class = "b"),
                      prior(exponential(2), class="shape"), 
                      prior(normal(0, 2), class="b", dpar="hu")),
            family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
            #sample_prior = "only",
            file="Stocked_nonstocked.rds",
            chains=4, iter=2000, cores=4)

#model summary and checks
summary(stock_mod)
bayes_R2(stock_mod)
pp_check(stock_mod, ndraws = 10,type = "boxplot")+theme_classic()

#Sample posterior for magnitude of difference
posts.stock <- add_epred_draws(stock_mod, re_formula = NA,
                             newdata = stock_mod$data %>% distinct(Stocked) %>% 
                               mutate(Waterbody2 = "NSDFNSDF"),
                             allow_new_levels = TRUE, dpar=T)

#calculate means and 95% quantile-based CrI
qi<-posts.stock%>% 
  mean_qi();qi

#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

po1<-posts.stock %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate") %>% 
  as.data.frame() %>% 
  select(-Waterbody2)

po1<-bind_rows(po1, wae3)
po1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy <- tibble(Stocked=rep(c("Y","N"),2),
                estimate=c("prob","prob", "WAE","WAE"),
                value=c(0,0.4,0,40))

po1 %>% 
  ggplot(aes(x=factor(Stocked, level=c("Y", "N")), y=value))+
  geom_point(data=subset(po1, estimate=="WAE"), aes(y=CPUE2), alpha=0.3, position = position_jitter(width=0.1))+
  geom_violin(fill='gray50',trim=F, alpha=0.9)+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_discrete(breaks=c("Y","N"),labels=c("Stocked", "Not Stocked"), limits=c("Y", "N"))+
  theme_classic()+
  geom_blank(data=dummy)+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
###### Size  of stocked products * waterbody type
#get prior 
get_prior(bf(CPUE2~Size2*Type + (1|Waterbody2),
             hu~Size2*Type + (1|Waterbody2)),
          data=wae3, family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"))

#model
S_T<- brm(bf(CPUE2~Size2*Type + (1|Waterbody2),
            hu~Size2*Type + (1|Waterbody2)),
         data=wae3, 
         prior = c(prior(normal(0,.5), class="b"),
                   prior(normal(-1,0.5), coef="TypeMarginal"),
                   prior(exponential(2), class="sd"),
                   prior(normal(1, 0.5), class="Intercept"),
                   prior(exponential(2), class = "shape"),
                   prior(normal(0, 2), class="b", dpar="hu")),
         family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
         #sample_prior = "only",
         file="Size_Type.rds", 
         chains=4, iter=2000, cores=4)

#model summary and check
summary(S_T)
bayes_R2(S_T)
pp_check(S_T)
plot(S_T)

#Sample from posterior
posts_S_T <- add_epred_draws(S_T, re_formula = NA,
                          newdata = S_T$data %>% distinct(Size2, Type) %>% 
                            mutate(Waterbody2 = "NSDFNSDF"),
                          allow_new_levels = TRUE, dpar=T)

#calculate means and 95% quantile-based CrI
qi_S_T<-posts_S_T %>% 
  mean_qi();qi_S_T

#Estimate probabilities of differences among groups
posts_S_T_wide<-posts_S_T %>%  
  ungroup() %>% 
  select(.draw, Size2, Type, .epred) %>% 
  pivot_wider(names_from = c(Size2, Type), values_from=.epred) %>% 
  mutate(FrySF_C=d_fr_Consistent-c_sf_Consistent,
         FryLF_C=d_fr_Consistent-b_lf_Consistent,
         FryNo_C=d_fr_Consistent-a_no_Consistent,
         FrySF_M=d_fr_Marginal-c_sf_Marginal,
         FryLF_M=d_fr_Marginal-b_lf_Marginal,
         FryNo_M=d_fr_Marginal-a_no_Marginal,
         Fry_CM = d_fr_Consistent-d_fr_Marginal,
         SF_CM = c_sf_Consistent - c_sf_Marginal, 
         LF_CM = b_lf_Consistent - b_lf_Marginal,
         No_CM = a_no_Consistent - a_no_Marginal)
sum(posts_S_T_wide$FrySF_C>0)/4000
sum(posts_S_T_wide$FryLF_C>0)/4000
sum(posts_S_T_wide$FryNo_C>0)/4000
sum(posts_S_T_wide$FrySF_M>0)/4000
sum(posts_S_T_wide$FryLF_M>0)/4000
sum(posts_S_T_wide$FryNo_M>0)/4000
sum(posts_S_T_wide$Fry_CM>0)/4000
sum(posts_S_T_wide$SF_CM>0)/4000
sum(posts_S_T_wide$LF_CM>0)/4000
sum(posts_S_T_wide$No_CM>0)/4000


#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

type1<-posts_S_T %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate") %>% 
  as.data.frame() %>% 
  select(-Waterbody2)

sizetype<-bind_rows(type1, wae3)
sizetype %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_ST <- tibble(Size2=rep(c("d_fr", "c_sf", "b_lf", "a_no"),2),
                   Type=c(rep("Marginal",4), rep("Consistent",4)), 
                  estimate=rep(c("prob","prob", "WAE","WAE"),2),
                  value=rep(c(0,0.6,0,40),2))
sizetype %>% 
  ggplot(aes(x=factor(Size2, level=c("d_fr", "c_sf", "b_lf", "a_no")), y=value, fill=Type))+
  geom_point(data=subset(sizetype, estimate=="WAE"), aes(y=CPUE2,group=Type), alpha=0.3, position = position_jitterdodge(jitter.width = .05))+
  geom_violin(trim=F, alpha=0.9)+
  xlab("Size at stocking")+
  ylab(NULL)+
  scale_fill_grey(start=0, end=.65, name = NULL)+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_discrete(breaks=c("a_no","b_lf","c_sf", "d_fr"),
                   limits=c( "d_fr","c_sf","b_lf","a_no"),
                   labels=c("Not Stocked", "Large Fingerling", "Small Fingerling", "Fry"))+
  geom_blank(data=dummy_ST)+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
###### Management purpose
# model
mp<-brm(bf(CPUE2 ~ ManagementPurpose + (1|Waterbody2),
            hu ~ ManagementPurpose + (1|Waterbody2)),
         data=wae3, 
         family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
         prior = c(prior(normal(0,.5), class="b"),
                   prior(exponential(0.1), class="sd"),
                   prior(normal(1, .5), class="Intercept"),
                   prior(exponential(0.5), class = "shape"), 
                   prior(normal(0,2), class="b", dpar="hu")),
         #sample_prior = "only",
         file="Management_Purpose.rds",
         chains=4, iter=2000, cores=4)

#model summary and check
summary(mp)
bayes_R2(mp)
pp_check(mp)
plot(mp)

#Sample from posterior
posts_mp <- add_epred_draws(mp, re_formula = NA,
                             newdata = mp$data %>% distinct(ManagementPurpose) %>% 
                               mutate(Waterbody2 = "NSDFNSDF"),
                             allow_new_levels = TRUE, dpar=T)

#calculate means and 95% quantile-based CrI
qi_mp<-posts_mp %>% 
  mean_qi();qi_mp

#Estimate probabilities of differences among groups
posts_mp_wide<-posts_mp %>%  
  ungroup() %>% 
  select(.draw, ManagementPurpose, .epred) %>% 
  pivot_wider(names_from = ManagementPurpose, values_from=.epred) %>% 
  mutate(SI = Supplemental - Introductory,
         SM = Supplemental - Maintenance, 
         SN = Supplemental - None, 
         IM = Introductory - Maintenance, 
         IN = Introductory - None, 
         MN = Maintenance - None)
sum(posts_mp_wide$SI>0)/4000
sum(posts_mp_wide$SM>0)/4000
sum(posts_mp_wide$SN>0)/4000
sum(posts_mp_wide$IM>0)/4000
sum(posts_mp_wide$IN>0)/4000
sum(posts_mp_wide$MN>0)/4000

#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

mp1<-posts_mp %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate") %>% 
  as.data.frame() %>% 
  select(-Waterbody2)

manage1<-bind_rows(mp1, wae3)
manage1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_MP <- tibble(ManagementPurpose=c("Introductory", "Maintenance", "Supplemental", "None"),
                   estimate=c("prob","prob", "WAE","WAE"),
                   value=c(0,0.8,0,40))

manage1 %>% 
  ggplot(aes(x=factor(ManagementPurpose,level=c("Introductory", "Maintenance", "Supplemental", "None")), y=value))+
  geom_point(data=subset(manage1, estimate=="WAE"), aes(y=CPUE2), alpha=0.3, position = position_jitter(width=0.1))+
  geom_violin(fill='gray50',trim=F, alpha=0.9)+
  geom_blank(data=dummy_MP)+
  ylab(NULL)+
  xlab("Mangagement purpose")+
  scale_x_discrete(limits=c("Introductory", "Maintenance", "Supplemental", "None"))+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
###### Waterbody Characteristics
#model
water_mod<-brm(bf(CPUE2 ~ Hectares_std + WaterbodyType  + (1|Waterbody2),
           hu ~ Hectares_std + WaterbodyType   + (1|Waterbody2)),
        data=wae3, 
        family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
        prior = c(prior(normal(0,.5), class="b"),
                  prior(exponential(0.1), class="sd"),
                  prior(normal(1, .5), class="Intercept"),
                  prior(exponential(0.5), class = "shape"), 
                  prior(normal(0,2), class="b", dpar="hu")),
        #sample_prior = "only",
        file_refit = "on_change",
        file="Waterbody_Char.rds",
        chains=4, iter=2000, cores=4)

#summary and model check
summary(water_mod)
bayes_R2(water_mod)
pp_check(water_mod)
plot(water_mod)

#sample from posterior
posts_water <- add_epred_draws(water_mod, re_formula = NA,
                                 newdata = tibble(expand.grid(Hectares_std=seq(min(wae3$Hectares_std), max(wae3$Hectares_std),length.out=100),
                                                  WaterbodyType = unique(wae3$WaterbodyType))) %>% 
                                   mutate(Waterbody2 = "NSDFNSDF"),
                                 allow_new_levels = TRUE, dpar=T)

#probability 
posts_water_wide<-posts_water %>%  
  ungroup() %>% 
  select(.draw, WaterbodyType,Hectares_std, .epred) %>% 
  pivot_wider(names_from = WaterbodyType, values_from=.epred) %>% 
  mutate(diffN_I = `Natural Basin` - `Impoundment/Excavated`)
sum(posts_water_wide$diffN_I>0)/400000

#summarize mean and 95% CrI
(water_meanqi<-posts_water %>%
  select(.draw, Hectares_std, WaterbodyType, .epred, hu) %>% 
  group_by(Hectares_std, WaterbodyType) %>% 
  mean_qi(.epred, hu))

#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

w1<- water_meanqi %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate")

w2<- w1 %>%
  rename(prob.lower=hu.lower, prob.upper=hu.upper, WAE.upper=.epred.upper, WAE.lower=.epred.lower) %>% 
  pivot_longer(cols=c(prob.lower,  WAE.lower), names_to ="int_low", values_to = "int_lower")
w3<- w2 %>% 
  pivot_longer(cols=c(prob.upper,  WAE.upper), names_to ="int_up", values_to = "int_upper")

w3.1<-w3 %>% filter(estimate=="prob" & 
              str_detect(int_low, "^p")& 
              str_detect(int_up, "^p"))
w3.2<-w3 %>% filter(estimate=="WAE" & 
                      str_detect(int_low, "^W")& 
                      str_detect(int_up, "^W"))

w4<-bind_rows(w3.1, w3.2)

water1<-bind_rows(w4, wae3)
water1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_w <- tibble(Hectares=rep(c(-.50,.5),4),
                  WaterbodyType=c(rep("Natural Basin",4), rep("Impoundment/Excavated",4)), 
                  estimate=rep(c("prob","prob", "WAE","WAE"),2),
                  value=rep(c(0,0.4,0,40),2))

water1 %>% ggplot()+
  geom_point(data=subset(water1, estimate=="WAE"), aes(x=Hectares_std, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=w4, aes(y=value, x=Hectares_std, ymin=int_lower, ymax=int_upper,fill=WaterbodyType), alpha=0.5)+
  geom_line(data=w4,aes(x=Hectares_std, y=value, linetype=WaterbodyType), size=1)+
  geom_blank(data=dummy_w, aes(y=value))+
  xlab("Surface area (hectares)")+
  ylab(NULL)+
  scale_fill_manual(values=c("black", "gray"),name = NULL, labels=c("Natural Basin", "Impoundment"))+
  scale_linetype_manual(values=c(1,3), name=NULL,labels=c("Natural Basin", "Impoundment"))+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_continuous(labels=c("0", "2,000", "4,000", "6,000"), 
                     breaks=c(((0-mean(wae3$Hectares))/sd(wae3$Hectares)),
                              ((2000-mean(wae3$Hectares))/sd(wae3$Hectares)),
                              ((4000-mean(wae3$Hectares))/sd(wae3$Hectares)),
                              ((6000-mean(wae3$Hectares))/sd(wae3$Hectares))))+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
######Environmental variables
#get prior
get_prior(bf(CPUE2~ GDD_std  + SpringPrecip_std + WSI_std + (1|Waterbody2),
             hu~Type * GDD_std + WaterbodyType + SpringPrecip_std + WSI_std + (1|Waterbody2)),
          data=wae3, family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"))

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
summary(env)
bayes_R2(env)
pp_check(env)
plot(env)

#Sample posterior for magnitude of difference
#Then estimate mean and 95% qi across range of values for each variable

#GDD
posts_env_GDD <- add_epred_draws(env, re_formula = NA,
                             newdata = tibble(GDD_std=seq(min(wae3$GDD_std), max(wae3$GDD_std),length.out=100),
                                              SpringPrecip_std=rep(0,100),
                                              WSI_std = rep(0,100)) %>% 
                             mutate(Waterbody2 = "NSDFNSDF"),
                             allow_new_levels = TRUE, dpar=T)

posts_env_GDD %<>%
  select(.draw, GDD_std, SpringPrecip_std, WSI_std, .epred, hu) %>% 
  group_by(GDD_std) %>% 
  mean_qi(.epred, hu) %>% 
  rename(Predictor = GDD_std) %>% 
  mutate(Parameter = "GDD_std")

#WSI
posts_env_WSI <- add_epred_draws(env, re_formula = NA,
                             newdata = tibble(GDD_std=rep(0,100),
                                              SpringPrecip_std=rep(0,100),
                                              WSI_std = seq(min(wae3$WSI_std), max(wae3$WSI_std),length.out=100)) %>% 
                               mutate(Waterbody2 = "NSDFNSDF"),
                             allow_new_levels = TRUE, dpar=T)
posts_env_WSI %<>%
  select(.draw, GDD_std, SpringPrecip_std, WSI_std, .epred, hu) %>% 
  group_by(WSI_std) %>% 
  mean_qi(.epred, hu) %>% 
  rename(Predictor = WSI_std) %>% 
  mutate(Parameter = "WSI_std")

#Spring Precip
posts_env_SP <- add_epred_draws(env, re_formula = NA,
                                 newdata = tibble(GDD_std=rep(0,100),
                                                  SpringPrecip_std=seq(min(wae3$SpringPrecip_std), max(wae3$SpringPrecip_std),length.out=100),
                                                  WSI_std = rep(0,100)) %>% 
                                   mutate(Waterbody2 = "NSDFNSDF"),
                                 allow_new_levels = TRUE, dpar=T)

posts_env_SP %<>%
  select(.draw, GDD_std, SpringPrecip_std, WSI_std, .epred, hu) %>% 
  group_by(SpringPrecip_std) %>% 
  mean_qi(.epred, hu) %>% 
  rename(Predictor = SpringPrecip_std) %>% 
  mutate(Parameter = "SpringPrecip_std")

env_dat<-bind_rows(posts_env_GDD, posts_env_WSI,posts_env_SP)

#evaluate trends in hu parameter
posts_env_GDD %>% 
  ggplot(aes(x=GDD_std, y=hu))+
  geom_line()+
  geom_ribbon(aes(ymin=hu.lower, ymax=hu.upper), alpha=0.3)+
  theme_classic()

posts_env_WSI %>% 
  ggplot(aes(x=WSI_std, y=hu))+
  geom_line()+
  geom_ribbon(aes(ymin=hu.lower, ymax=hu.upper), alpha=0.3)+
  theme_classic()

posts_env_SP %>% 
  ggplot(aes(x=SpringPrecip_std, y=hu))+
  geom_line()+
  geom_ribbon(aes(ymin=hu.lower, ymax=hu.upper), alpha=0.3)+
  theme_classic()

#plot all three relationships and combine into one figure
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

e1<- env_dat %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate")

e2<- e1 %>%
  rename(prob.lower=hu.lower, prob.upper=hu.upper, WAE.upper=.epred.upper, WAE.lower=.epred.lower) %>% 
  pivot_longer(cols=c(prob.lower,  WAE.lower), names_to ="int_low", values_to = "int_lower") %>% 
  pivot_longer(cols=c(prob.upper,  WAE.upper), names_to ="int_up", values_to = "int_upper")

e2.1<-e2 %>% filter(estimate=="prob" & 
                        str_detect(int_low, "^p")& 
                        str_detect(int_up, "^p"))
e2.2<-e2 %>% filter(estimate=="WAE" & 
                        str_detect(int_low, "^W")& 
                        str_detect(int_up, "^W"))

e3<-bind_rows(e2.1, e2.2)
env_wae <- wae3 %>% pivot_longer(cols = c(SpringPrecip_std,GDD_std, WSI_std), names_to = "Parameter", values_to = "Predictor")

env1<-bind_rows(e3, env_wae)
env1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_env <- tibble(Parameter=c(rep("GDD_std",4),rep("WSI_std",4),rep("SpringPrecip_std",4)),
                    estimate=rep(c("prob","prob", "WAE","WAE"),3),
                    value=rep(c(0,0.3,0,40),3))

p <- "SpringPrecip_std"

p1 <- ggplot()+
  geom_point(data=subset(env1, estimate=="WAE"&is.na(.interval)&Parameter==p), aes(x=Predictor, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=subset(e3, Parameter==p), aes(y=value, x=Predictor, ymin=int_lower, ymax=int_upper), alpha=0.5)+
  geom_line(data=subset(e3, Parameter==p),aes(x=Predictor, y=value), size=1)+
  geom_blank(data=subset(dummy_env, Parameter==p), aes(y=value))+
  xlab("Spring precipitation (mm)")+
  ylab(NULL)+
  scale_x_continuous(labels=c("100", "200", "300", "400", "500"), 
                     breaks=c(((100-mean(wae3$SpringPrecip))/sd(wae3$SpringPrecip)),
                              ((200-mean(wae3$SpringPrecip))/sd(wae3$SpringPrecip)),
                              ((300-mean(wae3$SpringPrecip))/sd(wae3$SpringPrecip)),
                              ((400-mean(wae3$SpringPrecip))/sd(wae3$SpringPrecip)),
                              ((500-mean(wae3$SpringPrecip))/sd(wae3$SpringPrecip))),
                     limits=c(-1.5,5))+
  scale_fill_manual(values=c("black", "gray"),name = NULL)+
  facet_wrap(~estimate, scales="free",strip.position = "left", 
             labeller= labeller(estimate=labs), ncol=2)+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside");p1

#
dummy_env2 <- tibble(Parameter=c(rep("GDD_std",4),rep("WSI_std",4),rep("SpringPrecip_std",4)),
                    estimate=rep(c("prob","prob", "WAE","WAE"),3),
                    value=rep(c(0,0.5,0,40),3))
q <- "GDD_std"

p2 <- ggplot()+
  geom_point(data=subset(env1, estimate=="WAE"&is.na(.interval)&Parameter==q), aes(x=Predictor, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=subset(e3, Parameter==q), aes(y=value, x=Predictor, ymin=int_lower, ymax=int_upper), alpha=0.5)+
  geom_line(data=subset(e3, Parameter==q),aes(x=Predictor, y=value), size=1)+
  geom_blank(data=subset(dummy_env2, Parameter==q), aes(y=value))+
  xlab("Growing degree days (GDD)")+
  ylab(NULL)+
  scale_x_continuous(labels=c("1900", "2200","2500", "2800"), 
                     breaks=c(((1900-mean(wae3$GDD))/sd(wae3$GDD)),
                              ((2200-mean(wae3$GDD))/sd(wae3$GDD)),
                              ((2500-mean(wae3$GDD))/sd(wae3$GDD)),
                              ((2800-mean(wae3$GDD))/sd(wae3$GDD))))+
  scale_fill_manual(values=c("black", "gray"),name = NULL)+
  facet_wrap(~estimate, scales="free",strip.position = "left", 
             labeller= labeller(estimate=labs), ncol=2)+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside");p2
#
dummy_env3 <- tibble(Parameter=c(rep("GDD_std",4),rep("WSI_std",4),rep("SpringPrecip_std",4)),
                     estimate=rep(c("prob","prob", "WAE","WAE"),3),
                     value=rep(c(0,0.3,0,40),3))
r <- "WSI_std"

p3 <- ggplot()+
  geom_point(data=subset(env1, estimate=="WAE"&is.na(.interval)&Parameter==r), aes(x=Predictor, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=subset(e3, Parameter==r), aes(y=value, x=Predictor, ymin=int_lower, ymax=int_upper), alpha=0.5)+
  geom_line(data=subset(e3, Parameter==r),aes(x=Predictor, y=value), size=1)+
  geom_blank(data=subset(dummy_env3, Parameter==r), aes(y=value))+
  xlab("Winter severity index (WSI)")+
  ylab(NULL)+
  scale_x_continuous(labels=c("-2,000", "-1,500", "-1,000", "-500"), 
                     breaks=c(((-2000-mean(wae3$WSI))/sd(wae3$WSI)),
                              ((-1500-mean(wae3$WSI))/sd(wae3$WSI)),
                              ((-1000-mean(wae3$WSI))/sd(wae3$WSI)),
                              ((-500-mean(wae3$WSI))/sd(wae3$WSI))))+
  scale_fill_manual(values=c("black", "gray"),name = NULL)+
  facet_wrap(~estimate, scales="free",strip.position = "left", 
             labeller= labeller(estimate=labs), ncol=2)+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside");p3

grid.arrange(p2, p3, p1, ncol=1)
############################################################
####### Stock-size Walleye
#filter data
adultWAE<-wae3 %>% filter(!is.na(WAE.S)) %>% 
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
        #sample_prior = "only",
        file="Adult_WAE.rds",
        chains=4, iter=2000, cores=4)


#model summary and checks
summary(wae_mod)
bayes_R2(wae_mod)
pp_check(wae_mod)
plot(wae_mod)

#sample posterior
posts_wae <- add_epred_draws(wae_mod, re_formula = NA,
                              newdata = tibble(WAE_adj=seq(min(adultWAE$WAE_adj), 
                                                         max(adultWAE$WAE_adj),
                                                         length.out=100)) %>% 
                                mutate(Waterbody2 = "NSDFNSDF"),
                              allow_new_levels = TRUE, dpar=T)

#summarize means and 95% CrI
wae_meanqi <- posts_wae %>%
  select(.draw, WAE_adj, .epred, hu) %>% 
  group_by(WAE_adj) %>% 
  mean_qi(.epred, hu)


#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

st1<- wae_meanqi %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate")

st2<- st1 %>%
  rename(prob.lower=hu.lower, prob.upper=hu.upper, WAE.upper=.epred.upper, WAE.lower=.epred.lower) %>% 
  pivot_longer(cols=c(prob.lower,  WAE.lower), names_to ="int_low", values_to = "int_lower")
st3<- st2 %>% 
  pivot_longer(cols=c(prob.upper,  WAE.upper), names_to ="int_up", values_to = "int_upper")

st3.1<-st3 %>% filter(estimate=="prob" & 
                      str_detect(int_low, "^p")& 
                      str_detect(int_up, "^p"))
st3.2<-st3 %>% filter(estimate=="WAE" & 
                      str_detect(int_low, "^W")& 
                      str_detect(int_up, "^W"))

st4<-bind_rows(st3.1, st3.2)

st_wae1<-bind_rows(st4, adultWAE)
st_wae1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_wae <- tibble(WAE_adj=rep(c(-.50,.5),2),
                  estimate=c("prob","prob", "WAE","WAE"),
                  value=c(0,1,0,40))
ggplot()+
  geom_point(data=subset(st_wae1, estimate=="WAE"), aes(x=WAE_adj, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=st4, aes(y=value, x=WAE_adj, ymin=int_lower, ymax=int_upper), alpha=0.5)+
  geom_line(data=st4,aes(x=WAE_adj, y=value), size=1)+
  geom_blank(data=dummy_wae, aes(y=value))+
  xlab("Stock-length Walleye CPGN")+
  ylab(NULL)+
  scale_fill_manual(values=c("black", "gray"),name = NULL)+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_continuous(labels=c(0, 25, 50, 75, 100), 
                     breaks=c(((0-mean(adultWAE$WAE.S))/sd(adultWAE$WAE.S)),
                              ((25-mean(adultWAE$WAE.S))/sd(adultWAE$WAE.S)),
                              ((50-mean(adultWAE$WAE.S))/sd(adultWAE$WAE.S)),
                              ((75-mean(adultWAE$WAE.S))/sd(adultWAE$WAE.S)),
                              ((100-mean(adultWAE$WAE.S))/sd(adultWAE$WAE.S))))+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
##### Centrarchids
# Filter data to only include waters with sampled centrarchids
cent<-wae3 %>% filter(!is.na(Centrarchids)) %>% 
  mutate(Centrarchids_std=(Centrarchids-mean(Centrarchids))/sd(Centrarchids))

#model
cent_mod<-brm(bf(CPUE2~ Centrarchids_std + (1|Waterbody2),
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
summary(cent_mod)
bayes_R2(cent_mod)
pp_check(cent_mod)
plot(cent_mod)

#Sample posterior for magnitude of difference
#GDD
posts_cent <- add_epred_draws(cent_mod, re_formula = NA,
                                 newdata = tibble(Centrarchids_std=seq(min(cent$Centrarchids_std), 
                                                                       max(cent$Centrarchids_std),
                                                                       length.out=100)) %>% 
                                   mutate(Waterbody2 = "NSDFNSDF"),
                                 allow_new_levels = TRUE, dpar=T)

#mean and 95% quantile interval
cent_meanqi <- posts_cent %>%
  select(.draw, Centrarchids_std, .epred, hu) %>% 
  group_by(Centrarchids_std) %>% 
  mean_qi(.epred, hu)

#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

c1<- cent_meanqi %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate")

c2<- c1 %>%
  rename(prob.lower=hu.lower, prob.upper=hu.upper, WAE.upper=.epred.upper, WAE.lower=.epred.lower) %>% 
  pivot_longer(cols=c(prob.lower,  WAE.lower), names_to ="int_low", values_to = "int_lower")
c3<- c2 %>% 
  pivot_longer(cols=c(prob.upper,  WAE.upper), names_to ="int_up", values_to = "int_upper")

c3.1<-c3 %>% filter(estimate=="prob" & 
                        str_detect(int_low, "^p")& 
                        str_detect(int_up, "^p"))
c3.2<-c3 %>% filter(estimate=="WAE" & 
                        str_detect(int_low, "^W")& 
                        str_detect(int_up, "^W"))

c4<-bind_rows(c3.1, c3.2)

cent1<-bind_rows(c4, cent)
cent1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_cent <- tibble(Centrarchids_std=rep(c(-.50,.5),2),
                    estimate=c("prob","prob", "WAE","WAE"),
                    value=c(0,1,0,40))
ggplot()+
  geom_point(data=subset(cent1, estimate=="WAE"), aes(x=Centrarchids_std, y=CPUE2), alpha=0.3)+
  geom_ribbon(data=c4, aes(y=value, x=Centrarchids_std, ymin=int_lower, ymax=int_upper), alpha=0.5)+
  geom_line(data=c4,aes(x=Centrarchids_std, y=value), size=1)+
  geom_blank(data=dummy_cent,aes(y=value))+
  xlab("Stock-length Centrarchid CPTN")+
  ylab(NULL)+
  scale_fill_manual(values=c("black", "gray"),name = NULL)+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_continuous(labels=c(0,100,200,300, 400, 500), 
                     breaks=c(((0-mean(cent$Centrarchids))/sd(cent$Centrarchids)),
                              ((100-mean(cent$Centrarchids))/sd(cent$Centrarchids)),
                              ((200-mean(cent$Centrarchids))/sd(cent$Centrarchids)),
                              ((300-mean(cent$Centrarchids))/sd(cent$Centrarchids)),
                              ((400-mean(cent$Centrarchids))/sd(cent$Centrarchids)), 
                              ((500-mean(cent$Centrarchids))/sd(cent$Centrarchids))))+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

############################################################
##### Saugeye lakes
#filter data
SAE<-wae2 %>% filter(Waterbody=="Campbell"|Waterbody=="Elm"|Waterbody=="Goldsmith"|
                       Waterbody=="Mina"|Waterbody=="Richmond"|Waterbody=="White") %>% 
  mutate(SppStocked = plyr::mapvalues(SpeciesStocked, from=c("None", "Walleye", "Saugeye"), 
                                      to=c("a_none", "b_Walleye", "c_Saugeye")))

#model
sae_mod<-brm(bf(CPUE2 ~ SppStocked + (1|Waterbody2),
            hu ~ SppStocked + (1|Waterbody2)),
         data=SAE, 
         family=hurdle_gamma(link="log", link_hu = "logit",link_shape = "log"),
         prior = c(prior(normal(0.5, 0.25), class="b"),
                   prior(exponential(0.1), class="sd"),
                   prior(normal(1, .5), class="Intercept"),
                   prior(exponential(0.5), class = "shape"), 
                   prior(normal(0,2), class="b", dpar="hu")),
         #sample_prior = "only",
         file="SAE_mod.rds",
         chains=4, iter=2000, cores=4)

#model summary and checks
summary(sae_mod)
bayes_R2(sae_mod)
pp_check(sae_mod)
plot(sae_mod)

#Sample posterior for magnitude of difference
posts_sae <- add_epred_draws(sae_mod, re_formula = NA,
                              newdata = sae_mod$data %>% distinct(SppStocked) %>% 
                                mutate(Waterbody2 = "NSDFNSDF"),
                              allow_new_levels = TRUE, dpar=T)

#mean and 95% quantile interval
(sae_meanqi <- posts_sae %>%
  select(.draw, SppStocked, .epred, hu) %>% 
  group_by(SppStocked) %>% 
  mean_qi(.epred, hu))

#Probabilities
posts_sae2<-posts_sae %>%  
  ungroup() %>% 
  select(.draw, SppStocked, .epred) %>% 
  pivot_wider(names_from = SppStocked, values_from=.epred) %>% 
  mutate(S2W=c_Saugeye-b_Walleye, 
         S2N=c_Saugeye-a_none)
sum(posts_sae2$S2W>0)/4000
sum(posts_sae2$S2N>0)/4000

#double plot
labs <- c("Probability age-2 Walleye CPGN=0", "Age-2 Walleye CPGN")
names(labs) <- c("prob", "WAE")

sae1<-posts_sae %>%
  rename(prob=hu, WAE=.epred) %>% 
  pivot_longer(cols=c(prob, WAE), names_to ="estimate") %>% 
  as.data.frame() %>% 
  select(-Waterbody2)

saugeye1<-bind_rows(sae1, SAE)
saugeye1 %<>% mutate(estimate=replace_na(estimate, "WAE"))

dummy_sae <- tibble(SppStocked=c(rep("c_Saugeye",4), rep("b_Walleye",4),rep("a_none",4)),
                     estimate=rep(c("prob","prob", "WAE","WAE"),3),
                     value=rep(c(0,1,0,30),3))
saugeye1 %>% 
  ggplot(aes(x=factor(SppStocked, level=c("c_Saugeye","b_Walleye","a_none")), y=value))+
  geom_violin(fill='gray50',trim=F)+
  geom_point(data=subset(saugeye1, estimate=="WAE"), aes(y=CPUE2), alpha=0.3, position = position_jitter(width=0.1))+
  geom_blank(data=dummy_sae, aes(y=value))+
  ylab(NULL)+
  xlab("Species Stocked")+
  facet_wrap(~estimate, scales="free", strip.position = "left", 
             labeller= labeller(estimate=labs))+
  scale_x_discrete(breaks=c("c_Saugeye","b_Walleye","a_none"),
                   labels=c("Saugeye", "Walleye", "None"))+
  theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        strip.placement = "outside")

###############
# plot all hu plots together
###############
final<-grid.arrange(hu_plot_stock, hu_plot_ST, hu_plot_mp, hu_plot_water, 
             hu_plot_env, hu_plot_wae, hu_plot_cent, hu_plot_sae, 
             ncol=2, left = textGrob("Probability of age-2 Walleye CPGN=0", gp=gpar(fontsize=14), rot=90));final

ggsave(plot = final, "Probability that age-2 Walleye CPGN equals zero_letters.jpeg", width=10, height=7, dpi=600, units="in")
 