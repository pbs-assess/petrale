# Code for getting MCMC samples from the Petrale Sole SS stock assessment 2024
# 100 random samples from the posteriors for scenario 300
# note that this was a two sex model and that year ranges are fixed below and specific to this assessment

scenario <- 300
species <- "Petrale Sole"
stock <- "Coastwide"

out<-readRDS("/mcmcresults.RData")
set.seed(34)
for(s in 1:100){

  year<-out[[1]]$timeseries$Yr
  recruits_mcmc<-biomass_mcmc<-catch_mcmc<-harvest_rate_mcmc<-matrix(NA, nrow=14000,ncol=length(year))
  for (i in 1:14000){
    biomass_mcmc[i,]<-out[[i]]$timeseries$SpawnBio
    catch_mcmc[i,]<-out[[i]]$timeseries$`sel(B):_1`#+out$timeseries$`sel(B):_2`
    harvest_rate_mcmc[i,]<-out[[i]]$timeseries$`F:_1`#+out$timeseries$`F:_2`
    recruits_mcmc[i,]<-out[[i]]$timeseries$Recruit_0
  }

  sample_num<-sample(1:14000,size=1)

  biomass<-biomass_mcmc[sample_num,]
  catch<-catch_mcmc[sample_num,]
  harvest_rate<-harvest_rate_mcmc[sample_num,]
  recruits<-recruits_mcmc[sample_num,]

  df <- data.frame(cbind(year, catch, biomass, harvest_rate, recruits))

  year <-1938:2028
  VulnBFem_mcmc<-VulnBMal_mcmc<-vbiomass_mcmc<-matrix(NA, nrow=14000, ncol=length(year))
  for(i in 1:14000){
    VulnBFem_mcmc[i,]<-rowSums(out[[i]]$batage[out[[i]]$batage$`Beg/Mid`=='B' & out[[i]]$batage$Sex==1 & out[[i]]$batage$Yr>1937,13:35]*out[[i]]$ageselex[out[[i]]$ageselex$Sex==1 & out[[i]]$ageselex$Factor=='Asel2' & out[[i]]$ageselex$Fleet==1 & out[[i]]$ageselex$Yr>1937,8:30])
    VulnBMal_mcmc[i,]<-rowSums(out[[i]]$batage[out[[i]]$batage$`Beg/Mid`=='B' & out[[i]]$batage$Sex==2 & out[[i]]$batage$Yr>1937,13:35]*out[[i]]$ageselex[out[[i]]$ageselex$Sex==2 & out[[i]]$ageselex$Factor=='Asel2' & out[[i]]$ageselex$Fleet==1 & out[[i]]$ageselex$Yr>1937,8:30])
    vbiomass_mcmc[i,]<-VulnBFem_mcmc[i,]+VulnBMal_mcmc[i,]
  }

  vbiomass<-vbiomass_mcmc[sample_num,]

  dfv <- data.frame(cbind(year, vbiomass))

  year<-out[[1]]$recruit$Yr
  rdev_mcmc<-matrix(NA, nrow=14000, ncol=length(year))
  for(i in 1:14000){rdev_mcmc[i,]<-out[[i]]$recruit$dev}

  rdev<-rdev_mcmc[sample_num,]

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

  saveRDS(df,paste0("C:/Users/fischn/Documents/300/dfs/df",s,".RData"))
}
