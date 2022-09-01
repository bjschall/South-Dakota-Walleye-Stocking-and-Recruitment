library(readxl)
library(FSA)
library(magrittr)
library(plyr)
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
  mutate(ManagementPurpose=mapvalues(ManagementPurpose,from="Winterkill", to="Introductory"))
wae<-droplevels(filter(wae, !is.na(YC2.conv2GFP)&Size!="Multiple"&Size!="Adult"&Size!="Juvenile"&Size!="Fingerling"|is.na(Size)))
wae2<-droplevels(filter(wae, !is.na(YC2.conv2GFP)))
wae2 %<>% mutate(CPUE_adj=(YC2.conv2GFP+0.1),
                 Size2=as.character(as.numeric(Size))) %>% 
  mutate(Size2=mapvalues(Size2, from=c("1","2", "3", "4"),
                         to=c("d_fr","b_lf","a_no","c_sf"))) %>% 
  mutate(Size2=as.factor(Size2),
         Waterbody2=group_indices(wae2, .dots=c("Waterbody", "County")),
         GDD_std=(GDD-mean(GDD))/sd(GDD),
         SpringPrecip_std=(SpringPrecip-mean(SpringPrecip))/sd(SpringPrecip),
         Density_std=(NumStockedperAcre-mean(NumStockedperAcre))/sd(NumStockedperAcre),
         WSI_std=(WSI-mean(WSI))/sd(WSI),
         Hectares_std=(Hectares-mean(Hectares))/sd(Hectares),
         Centrarchids=BLC.All+BLG.All,
         Type=ifelse(Waterbody=="Antelope"|Waterbody=="Bitter"| Waterbody=="Blue Dog"|Waterbody=="Brant"|Waterbody=="Brush"|Waterbody=="Buffalo South"|
                       Waterbody=="Clear"&County=="Hamlin"|Waterbody=="Clear"&County=="Marshall"|Waterbody=="Dry 2"|Waterbody=="Enemy Swim"|
                       Waterbody=="Goose"|Waterbody=="Hazeldon"|Waterbody=="Horseshoe"|Waterbody=="Hwy 81 West"|Waterbody=="Kampeska"| Waterbody=="Long"|
                       Waterbody=="Lynn"|Waterbody=="Madison"|Waterbody=="Middle Lynn"|Waterbody=="Opitz"|Waterbody=="Pickerel"|Waterbody=="Piyas"|Waterbody=="Poinsett"|
                       Waterbody=="Reetz"|Waterbody=="Reid"|Waterbody=="Sinai"|Waterbody=="Thompson"|Waterbody=="Twin"&County=="Minnehaha"|Waterbody=="Waubay", "Consistent","Marginal")) %>% 
  select(-WAE.All,-WHB.All,-BLB.All, -COC.All,-WHS.All, -YEP.All,-YEP.psd,-YEP.psdp,
         -Region,-YC3.conv2GFP,-YC3.GNType, -YC3, -COC.psd,-COC.psdp,-WHS.psd,-WHS.psdp, -CPUE3,
         -BLB.psd,-BLB.psdp,-WHB.psd,-WHB.psdp,L1,L3,L4,L5,L6,L7,L8,L9,L10, N1,N3,N4,N5,N6,N7,N8,N9,`N10+`)
wae3<-droplevels(filter(wae2,SpeciesStocked != "Saugeye"))
