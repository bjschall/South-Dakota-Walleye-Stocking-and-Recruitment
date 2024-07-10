library(readxl)
library(magrittr)
library(dplyr)
library(Matrix)
library(tidyverse)
wae<-read_excel("WAE Stocking Analysis Final 2-2-21.xlsx", sheet="CleanData")
wae%<>%mutate(CPUE2=YC2.conv2GFP,County=as.factor(County),
              Waterbody = as.factor(Waterbody), 
              Size=as.factor(Size), 
              ManagementPurpose=as.factor(ManagementPurpose),
              WaterbodyType=as.factor(WaterbodyType),
              SpeciesStocked=as.factor(SpeciesStocked),
              Stocked=as.factor(Stocked),
              Hectares=Acres*0.404686) %<>%
  mutate(ManagementPurpose=plyr::mapvalues(ManagementPurpose,from="Winterkill", to="Introductory"))
wae<-droplevels(filter(wae, !is.na(YC2.conv2GFP)&Size!="Multiple"&Size!="Adult"&Size!="Juvenile"&Size!="Fingerling"|is.na(Size)))
wae2<-droplevels(filter(wae, !is.na(YC2.conv2GFP)))
wae2 %<>% mutate(CPUE_adj=(YC2.conv2GFP+0.1),
                 Size2=as.character(as.numeric(Size))) %>% 
  mutate(Size2=plyr::mapvalues(Size2, from=c("1","2", "3", "4"),
                         to=c("d_fr","b_lf","a_no","c_sf"))) %>% 
  mutate(Size2=as.factor(Size2),
         Waterbody2=group_indices(wae2, .dots=c("Waterbody", "County")),
         Centrarchids=BLC.S+BLG.S,
         Type=ifelse(Waterbody=="Antelope"|Waterbody=="Bitter"| Waterbody=="Blue Dog"|Waterbody=="Brant"|Waterbody=="Brush"|Waterbody=="Buffalo South"|
                       Waterbody=="Clear"&County=="Hamlin"|Waterbody=="Clear"&County=="Marshall"|Waterbody=="Dry 2"|Waterbody=="Enemy Swim"|
                       Waterbody=="Goose"|Waterbody=="Hazeldon"|Waterbody=="Horseshoe"|Waterbody=="Hwy 81 West"|Waterbody=="Kampeska"| Waterbody=="Long"|
                       Waterbody=="Lynn"|Waterbody=="Madison"|Waterbody=="Middle Lynn"|Waterbody=="Opitz"|Waterbody=="Pickerel"|Waterbody=="Piyas"|Waterbody=="Poinsett"|
                       Waterbody=="Reetz"|Waterbody=="Reid"|Waterbody=="Sinai"|Waterbody=="Thompson"|Waterbody=="Twin"&County=="Minnehaha"|Waterbody=="Waubay", "Consistent","Marginal")) %>% 
  select(-WAE.All,-WHB.All,-BLB.All, -COC.All,-WHS.All, -YEP.All,-YEP.psd,-YEP.psdp,
         -Region,-YC3.conv2GFP,-GNType,-YC3.GNType, -YC3, -COC.psd,-COC.psdp,-WHS.psd,-WHS.psdp, -CPUE3,
         -BLB.psd,-BLB.psdp,-WHB.psd,-WHB.psdp,L1,L3,L4,L5,L6,L7,L8,L9,L10, N1,N3,N4,N5,N6,N7,N8,N9,`N10+`)
#add NOP data
nop<-read.csv("Nop.csv")

nop2<-nop %>% 
  select(Year.cohort, Waterbody, County, Method, total_nop, total_nop_stock,Station) %>% 
  group_by(Year.cohort, Waterbody, County) %>% 
  summarize(total_nop=sum(total_nop), 
            total_nop_stock=sum(total_nop_stock),
            Station=mean(Station),
            Method=unique(Method)) 

waenop<-nop2 %>% select(-Station) %>% 
  right_join(wae2, by=c("County", "Waterbody", "Year.cohort"))

wae2.1<- waenop %>% mutate(NOP.All.raw=total_nop/NbofGN,
                           NOP.S.raw = total_nop_stock/NbofGN) %>% 
  mutate(NOP.All = ifelse(Method=="std exp gill net", NOP.All.raw, (10^(1.187*(log10(NOP.All.raw+1))+0.135))-1),
         NOP.S = ifelse(Method=="std exp gill net", NOP.S.raw, (10^(1.187*(log10(NOP.S.raw+1))+0.135))-1)) %>% 
  select(-NOP.All.raw, -NOP.S.raw, -Method); wae2.1
  
  
write.csv(wae2.1, "wae2.1.csv")



##
wae2.2<-wae2.1 %>% 
  filter(SpeciesStocked != "Saugeye") 
wae.final<- wae2.2 %>% ungroup() %>% 
  mutate(GDD_std=(as.numeric(GDD)-mean(as.numeric(GDD)))/sd(as.numeric(GDD)),
         SpringPrecip_std=(SpringPrecip-mean(SpringPrecip))/sd(SpringPrecip),
         Density_std=(NumStockedperAcre-mean(NumStockedperAcre))/sd(NumStockedperAcre),
         WSI_std=(WSI-mean(WSI))/sd(WSI),
         Hectares_std=(Hectares-mean(Hectares))/sd(Hectares))

write.csv(wae.final, "wae.final.csv")

