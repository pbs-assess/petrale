# Code for getting MLE values from the Petrale Sole SS stock assessment 2024
# note that this was a two sex model and that year ranges are fixed below and specific to this assessment

library(r4ss)
library(here)
library(tidyverse)
library(DescTools)
library(ggthemes)

species <- "Petrale Sole"
stock <- "Coastwide"
scenario<-300

wd<-getwd()

out<-SS_output(dir=paste(wd,"Scenarios", scenario, sep="/"))

year<-out$timeseries$Yr
biomass<-out$timeseries$SpawnBio
catch<-out$timeseries$`sel(B):_1`#+out$timeseries$`sel(B):_2`
harvest_rate <-out$timeseries$`F:_1`#+out$timeseries$`F:_2`
recruits<-out$timeseries$Recruit_0


df <- data.frame(cbind(year, catch, biomass, harvest_rate, recruits))

VulnBFem<-rowSums(out$batage[out$batage$`Beg/Mid`=='B' & out$batage$Sex==1 & out$batage$Yr>1937,13:35]*out$ageselex[out$ageselex$Sex==1 & out$ageselex$Factor=='Asel2' & out$ageselex$Fleet==1 & out$ageselex$Yr>1937,8:30])

VulnBMal<-rowSums(out$batage[out$batage$`Beg/Mid`=='B' & out$batage$Sex==2 & out$batage$Yr>1937,13:35]*out$ageselex[out$ageselex$Sex==2 & out$ageselex$Factor=='Asel2' & out$ageselex$Fleet==1 & out$ageselex$Yr>1937,8:30])

vbiomass <-VulnBFem+VulnBMal
year <-1938:2028

dfv <- data.frame(cbind(year, vbiomass))

year<-out$recruit$Yr
rdev<-out$recruit$dev

df_rdev <- data.frame(cbind(year, rdev))

df <- left_join(df, dfv)|> left_join(df_rdev)
df$species <- species
df$stock <- stock
df$scenario <- scenario


df <- df %>% mutate(
  vbiomass_lead1 = lead(vbiomass),
  production = (vbiomass_lead1 + catch - vbiomass),
  p_by_biomass = production / vbiomass
)


saveRDS(df, paste0("data-generated/timeseries-", scenario,".rds"))

