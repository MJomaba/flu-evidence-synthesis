library(plyr)
library(MASS)
library(reshape2)
library(triangle)
library(gplots)

#HOME="C:/Documents and Settings/Marc.Baguelin/My Documents/Dropbox/taf/CEA_paper_flu_kid"
#HOME="F:/CEA_paper_flu_kid"
HOME="/Volumes/Transcend/CEA_paper_flu_kid"
#HOME="/Users/marc/Dropbox/taf/CEA_paper_flu_kid"

USE_CLUSTER=F
SIMU<-"CEA_paper_final_simus"

n_simu=1000
horizon_eval=10
n_sample_inference=1000
n_scenario=52
discount=3.5

lnmu  =function(mean, sd) log(mean)-0.5*log(1+(sd/mean)^2)
lnsig =function(mean, sd) sqrt(log(1+(sd/mean)^2))

source(paste(HOME,"R/HypergeoConf.R",sep="/"))
source(paste(HOME,"R/Utils.R",sep="/"))

setwd(HOME)

dir_data<-paste(HOME,"data",sep="/")
dir_Rsave<-paste(HOME,"R/Rsave",sep="/")
dir_pdf<-paste(HOME,"pdf",sep="/")

dir_res<-paste(HOME,SIMU,sep="/")
dir_res_pdf<-paste(dir_res,"pdf",sep="/")
dir_res_Rsave<-paste(dir_res,"Rsave",sep="/")

if(!file.exists(dir_res_pdf)){        
  dir.create(dir_res_pdf,recursive=T)
}
if(!file.exists(dir_res_Rsave)){        
  dir.create(dir_res_Rsave,recursive=T)
}	

z_name<-"z_hyper"

name_scenario<-c(as.character(0:34),"no_vaccination","vaccination")
id_scenario_ref<-c("no_vaccination",33,"vaccination",7,34)
id_scenario_test<-(0:34)[-(c(7,33,34)+1)]
id_scenario_lic<-35:64 #scenario for Low Income Country (LIC)
name_ref<-c("no vaccination","65+ at 29.3% (as pre-2000)","historical","65+ at 70% (as post-2000)","all HR at 75%")

target<-c("S1: 0.5-4 y (LR)","S3: 50-64 y (LR)","S2: 5-16 y (LR)","S5: S1+S3","S4: S1+S2","S6: S1+S2+S3","S7: Universal")
target_lic<-c("All risk groups","5-14 y","All risk groups and 65+ y")

coverage<-c("30%","50%","15%","70%")
n_dose_max_lic<-6.5e6
n_dose_lic<-seq(0.1,1,0.1)*n_dose_max_lic

strains<-c("H1N1","H3N2","B")
years<-1995:2008
years_HARDELID<-1999:2008
years_CROMER<-2000:2007
seasons<-year_to_season(years)
seasons_wo200203<-year_to_season(years[years!=2002])
seasons_glm<-year_to_season(years_CROMER[years_CROMER!=2002])

col_scen=c("#33232322","#aa232f22","#aaaa2322","#33ee2322","#33d3ee22","#a323d322","#332aa322","#ee2a5322","#332d2322","#3323aa22","#33a32322","#aa552e22")
col_scen_rein=c("#332323ee","#aa232fee","#aaaa23ee","#33ee23ee","#33d3eeee","#a323d3ee","#332aa3ee","#ee2a53ee","#332d23ee","#3323aaee","#33a323ee","#aa552eee")

age_cut_MODEL<-c(0,1,5,15,25,45,65)
age_cut_RCGP<-c(0,5,15,45,65)
age_cut_YO<-c(0,15)#YO means Young Old, this is for computation of Rchild and Radult

age_group_name<-c("0 - 4 y","5 - 14 y","15 - 44 y","45 - 64 y","65+ y","Total")
age_group_MODEL<-c("0 - 1 y","1 - 4 y","5 - 14 y","15 - 24 y","25 - 44 y","45 - 64 y","65+ y")
age_group_YO<-c("0 - 15 y","15+ y")
age_group_RCGP<-c("0 - 4 y","5 - 14 y","15 - 44 y","45 - 64 y","65+ y")
age_group_POST_DIST<-c("0 - 14 y","15 - 64 y","65+ y")
age_group_CONTACT_MATRIX<-c("0 - 1","1 - 4","5 - 14","15 - 24","25 - 44","45 - 64","65+")
age_group_CONTACT_MATRIX<-c("0-1","1-4","5-14","15-24","25-44","45-64","65+")
age_group_CONTACT_MATRIX<-factor(age_group_CONTACT_MATRIX,levels=age_group_CONTACT_MATRIX)
risk_group<-c("low risk","high risk")

age_group_CROMER<-c("0 - 6 m","6 m - 4 y","5 - 14 y","15 - 44 y","45 - 64 y","65+ y")
age_group_deaths_CROMER<-c("0 - 14 y","15 - 64 y","65+ y")
age_group_HARDELID<-c("0 - 14 y","15 - 44 y","45 - 64 y","65+ y")
age_group_HARDELID_OBS<-c("0 - 14 y","15 - 44 y","45 - 74 y","75+ y")#obsolete

#death_QALY_loss_3_5<-c(24.07,23.90,23.17,20,14.65,7.27,22.61,23.15,22.32,19,12.75,6.52)
if(discount==3.5) death_QALY_loss<-c(24.07,24.05,23.18,20.33,14.62,7.28,22.61,23.00,22.32,18.84,12.71,6.52)
if(discount==1.5) death_QALY_loss<-c(40.77,40.62,37.84,30.40,19.10,8.40,37.06,37.61,35.19,27.16,16.09,7.41)
if(discount==0) death_QALY_loss<-c(69.96,69.45,61.64,44.47,24.14,9.46,60.93,61.62,54.96,38.23,19.75,8.24)

#death_QALY_loss_3_5_sd<-c(0.6162,0.5642,0.3863,0.16,0.1528,0.1397,0.5878,0.5543,0.3761,0.14,0.1206,0.1172)

prop_risk_MODEL<-c(2.1,5.5,9.8,8.7,9.2,18.3,45.0)

names_age_risk_group<-paste(age_group_MODEL,rep(c("LR","HR"),each=length(age_group_MODEL)))

N_AGE_GROUP<-length(age_group_MODEL)
N_AGE_GROUP_YO<-length(age_group_YO)

SCALAR_NAME<-c("Transmissibility","Basic reproduction number","Effective reproduction number")
#SCALAR_NAME<-c("Transmissibility (q)","R[0]","paste(R[e],group(\"(\",t==0,\")\"))")
SUSCEPTIBILITY_NAME<-"Susceptibility"
CONSULTATION_NAME<-"Ascertainment probability"
POSTERIOR_NAMES<-c(paste("epsilon",1:3,sep="_"),"psi","q",paste("sigma",1:3,sep="_"),"l")

SAMPLE_SIZE<-1000
PRIOR_SAMPLE_SIZE<-10000
BINOM_CONF_DATA=F
HYPERGEO_CONF_DATA=!BINOM_CONF_DATA
CI_SAMPLE_SIZE=10
CI_QUANTILE=c(0.025,0.25,0.5,0.75,0.975)
CI_NAMES<-c("low_95","low_50","median","up_50","up_95")

R_extension<-3 #expand Re plot posterior by 3%

consult_prior_param<-data.frame(AG=age_group_POST_DIST,mean=c(-4.493789,-4.117028,-2.977965),sd=c(0.2860455,0.4751615,1.331832))
suscept_prior_param<-data.frame(AG=age_group_POST_DIST,mean=c(0.688,0.529,0.523),sd=c(0.083,0.122,0.175))

#FOR_MARC: mettre ici le nouveau prior pour la transmissibilit??
transmi_prior<-data.frame(mean=c(0.156889),sd=c(0.02480252))

objo5000=function(x,y,lim=5000) {(lim-y)/as.integer(as.Date("8/10/2014","%d/%m/%Y")-as.Date(x,"%d/%m/%Y"))}

CROMER_BURDEN_DATA<-function(){
  cat("start of CROMER_BURDEN_DATA\n")
  year_names<-year_to_season(2000:2007)
  
  #influenza GP consultation for influenza A in Cromer et al. 2013 for 8 years (2000/01 to 2007/08)
  
  #read the GP fluA consultation mean per year
  temp<-read.csv(paste(dir_data,"Cromer_GP_fluA.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the GP fluA standard deviation per year
  temp<-read.csv(paste(dir_data,"Cromer_GP_fluA_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_GP_fluA<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_GP_fluA,file=paste(dir_Rsave,"tab_GP_fluA.R",sep="/"))
  cat("Data GP fluA from Cromer loaded.\n")
  
  #read the GP fluB consultation mean per year
  temp<-read.csv(paste(dir_data,"Cromer_GP_fluB.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the GP fluB standard deviation per year
  temp<-read.csv(paste(dir_data,"Cromer_GP_fluB_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_GP_fluB<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_GP_fluB,file=paste(dir_Rsave,"tab_GP_fluB.R",sep="/"))
  cat("Data GP fluB from Cromer loaded.\n")
  
  #influenza hospitalisation for influenza A in Cromer et al. 2013 for 8 years (2000/01 to 2007/08)
  
  #read the hospitalisation mean per year (LR fluA)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluA_LR.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the hospitalisation standard deviation per year (LR fluA)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluA_LR_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_hosp_fluA_LR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_hosp_fluA_LR,file=paste(dir_Rsave,"tab_hosp_fluA_LR.R",sep="/"))
  cat("Data hosp fluA low risk from Cromer loaded.\n")
  
  #read the hospitalisation mean per year (HR fluA)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluA_HR.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the hospitalisation standard deviation per year (HR fluA)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluA_HR_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_hosp_fluA_HR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_hosp_fluA_HR,file=paste(dir_Rsave,"tab_hosp_fluA_HR.R",sep="/"))
  cat("Data hosp fluA high risk from Cromer loaded.\n")
  
  #influenza hospitalisation for influenza B in Cromer et al. 2013 for 8 years (2000/01 to 2007/08)
  
  #read the hospitalisation mean per year (LR fluB)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluB_LR.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the hospitalisation standard deviation per year (LR fluB)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluB_LR_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_hosp_fluB_LR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_hosp_fluB_LR,file=paste(dir_Rsave,"tab_hosp_fluB_LR.R",sep="/"))
  cat("Data hosp fluB low risk from Cromer loaded.\n")
  
  #read the hospitalisation mean per year (HR fluB)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluB_HR.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the hospitalisation standard deviation per year (HR fluB)
  temp<-read.csv(paste(dir_data,"Cromer_hosp_fluB_HR_sd.csv",sep="/"))
  temp[,1]<-age_group_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_hosp_fluB_HR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_hosp_fluB_HR,file=paste(dir_Rsave,"tab_hosp_fluB_HR.R",sep="/"))
  cat("Data hosp fluB high risk from Cromer loaded.\n")
  
  #read the death mean per year (LR)
  temp<-read.csv(paste(dir_data,"Cromer_death_LR.csv",sep="/"))
  temp[,1]<-age_group_deaths_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the death standard deviation per year (LR)
  temp<-read.csv(paste(dir_data,"Cromer_death_LR_sd.csv",sep="/"))
  temp[,1]<-age_group_deaths_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_death_LR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_death_LR,file=paste(dir_Rsave,"tab_death_LR.R",sep="/"))
  cat("Data death low risk from Cromer loaded.\n")
  
  #read the death mean per year (HR)
  temp<-read.csv(paste(dir_data,"Cromer_death_HR.csv",sep="/"))
  temp[,1]<-age_group_deaths_CROMER
  names(temp)<-c("age_group",year_names)
  temp1<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="mean")
  
  #read the death standard deviation per year (HR)
  temp<-read.csv(paste(dir_data,"Cromer_death_HR_sd.csv",sep="/"))
  temp[,1]<-age_group_deaths_CROMER
  names(temp)<-c("age_group",year_names)
  temp2<-melt(temp,id.vars=c("age_group"),variable.name="year",value.name="sd")
  
  #merge the two tables
  tab_death_HR<-merge(temp1,temp2,by=c("age_group","year"))
  
  save(tab_death_HR,file=paste(dir_Rsave,"tab_death_HR.R",sep="/"))
  cat("Data death low risk from Cromer loaded.\n")
  
  cat("end of CROMER_BURDEN()\n\n")
  
}

MODEL_2_CROMER<-function(tab_M){
  
  #template:
  #7 age_groups (0-1,1-4,5-14,15-24,25-44,45-64,65+)
  #3 risk groups (low,high,preg)
  #HARDELID: 0-14,15-44,45-64,65+
  tmp<-tab_M
  #low risk +pregnant
  tmp[,1:7]<-tmp[,1:7]+tmp[,15:21]
  tmp<-tmp[,-c(15:21)]
  
  #group 15-24 and 25-44y
  tmp[,c(4,11)]<-tmp[,c(4,11)]+tmp[,c(5,12)]
  #split the 0-1 between 0-6m
  tmp[,c(1,8)]<-tmp[,c(1,8)]/2
  #add the 6m-1y to 1-4
  tmp[,c(2,9)]<-tmp[,c(2,9)]+tmp[,c(1,8)]
  
  tmp<-tmp[,-c(5,12)]
   
  colnames(tmp)<-paste(rep(age_group_CROMER,2),rep(risk_group,each=length(age_group_CROMER)))
  return(tmp)
  
}

MODEL_2_DEATH_CROMER<-function(tab_M){
  
  #template:
  #7 age_groups (0-1,1-4,5-14,15-24,25-44,45-64,65+)
  #3 risk groups (low,high,preg)
  #HARDELID: 0-14,15-44,45-64,65+
  tmp<-tab_M
  #low risk +pregnant
  tmp[,1:7]<-tmp[,1:7]+tmp[,15:21]
  tmp<-tmp[,-c(15:21)]
  
  #group 0-1, 1-4, and 5-14
  tmp[,c(1,8)]<-tmp[,c(1,8)]+tmp[,c(2,9)]+tmp[,c(3,10)]
  #group 15-44 and  45-64
  tmp[,c(4,11)]<-tmp[,c(4,11)]+tmp[,c(5,12)]+tmp[,c(6,13)]
  
  tmp<-tmp[,-c(2,9,3,10,5,12,6,13)]
  
  colnames(tmp)<-paste(rep(age_group_deaths_CROMER,2),rep(risk_group,each=length(age_group_deaths_CROMER)))
  return(tmp)
  
}

COMPUTE_EXPLANATORY_VARIABLES_GLM_SAMPLE<-function(){
  cat("start of COMPUTE_EXPLANATORY_VARIABLES_GLM_SAMPLE_GP_fluA()\n")
  #compute the mean incidence predicted by the model, by age (CROMER age group), risk and year (CROMER paper)
  
  tab_ply<-expand.grid(strain=strains,year=seasons_glm)
  
  cat("generate incidence samples\n")
  df_AR<-ddply(tab_ply,c("year","strain"),function(df){
    
    year<-df$year
    strain<-df$strain
    
    file_name<-"Scenario_vaccination_final_size.txt"
    tab=read.table(paste(dir_res,year,strain,"scenarii",file_name,sep="/"),h=F)
    tmp<-MODEL_2_CROMER(tab)
    
    #incidence in the low risk
    tmp2<-tmp[,1:6]
    names(tmp2)<-age_group_CROMER
    tmp2$sample<-1:nrow(tmp2)
    tmp2$risk_group<-rep("LR",nrow(tmp2))
    
    #incidence in the high risk
    tmp3<-tmp[,7:12]
    names(tmp3)<-age_group_CROMER
    tmp3$sample<-1:nrow(tmp3)
    tmp3$risk_group<-rep("HR",nrow(tmp3))
    
    return(rbind(tmp2,tmp3))
  },.progress="text")
  
  gdf_AR<-melt(df_AR,id.vars=c("year","strain","sample","risk_group"),variable.name="age_group",value.name="incidence")
  
  tab_incidence_sample<-dcast(gdf_AR,formula=sample+age_group+year~strain+risk_group,value.var="incidence")
  
  cat("generate GP samples for fluA\n")
  load(paste(dir_Rsave,"tab_GP_fluA.R",sep="/"))
  tab_GP_fluA_sample<-ddply(tab_GP_fluA,c("age_group","year"),function(df){
    tmp<-data.frame(GP_consultations_fluA=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  tab_glm_sample<-merge(tab_GP_fluA_sample,tab_incidence_sample,by=c("sample","age_group","year"))
  
  cat("generate GP samples for fluB\n")
  load(paste(dir_Rsave,"tab_GP_fluB.R",sep="/"))
  tab_GP_fluB_sample<-ddply(tab_GP_fluB,c("age_group","year"),function(df){
    tmp<-data.frame(GP_consultations_fluB=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_GP_fluB_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  cat("generate hospitalisation samples for fluA (low risk)\n")
  load(paste(dir_Rsave,"tab_hosp_fluA_LR.R",sep="/"))
  tab_hosp_fluA_LR_sample<-ddply(tab_hosp_fluA_LR,c("age_group","year"),function(df){
    tmp<-data.frame(hospitalisations_fluA_LR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_hosp_fluA_LR_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  cat("generate hospitalisation samples for fluA (high risk)\n")
  load(paste(dir_Rsave,"tab_hosp_fluA_HR.R",sep="/"))
  tab_hosp_fluA_HR_sample<-ddply(tab_hosp_fluA_HR,c("age_group","year"),function(df){
    tmp<-data.frame(hospitalisations_fluA_HR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_hosp_fluA_HR_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  cat("generate hospitalisation samples for fluB (low risk)\n")
  load(paste(dir_Rsave,"tab_hosp_fluB_LR.R",sep="/"))
  tab_hosp_fluB_LR_sample<-ddply(tab_hosp_fluB_LR,c("age_group","year"),function(df){
    tmp<-data.frame(hospitalisations_fluB_LR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_hosp_fluB_LR_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  cat("generate hospitalisation samples for fluB (high risk)\n")
  load(paste(dir_Rsave,"tab_hosp_fluB_HR.R",sep="/"))
  tab_hosp_fluB_HR_sample<-ddply(tab_hosp_fluB_HR,c("age_group","year"),function(df){
    tmp<-data.frame(hospitalisations_fluB_HR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_hosp_fluB_HR_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  cat("end of COMPUTE_EXPLANATORY_VARIABLES_GLM_SAMPLE()\n\n")
}

COMPUTE_EXPLANATORY_DEATHS_GLM_SAMPLE<-function(){
  cat("start of COMPUTE_EXPLANATORY_DEATHS_GLM_SAMPLE\n")
  #compute the mean incidence of deaths predicted by the model, by age (CROMER death age group), risk and year (CROMER paper)

  tab_ply<-expand.grid(strain=strains,year=seasons_glm)
  
  cat("generate incidence samples\n")
  df_AR<-ddply(tab_ply,c("year","strain"),function(df){
    
    year<-df$year
    strain<-df$strain
    
    file_name<-"Scenario_vaccination_final_size.txt"
    tab=read.table(paste(dir_res,year,strain,"scenarii",file_name,sep="/"),h=F)
    tmp<-MODEL_2_DEATH_CROMER(tab)
    
    #incidence in the low risk
    tmp2<-tmp[,1:3]
    names(tmp2)<-age_group_deaths_CROMER
    tmp2$sample<-1:nrow(tmp2)
    tmp2$risk_group<-rep("LR",nrow(tmp2))
    
    #incidence in the high risk
    tmp3<-tmp[,4:6]
    names(tmp3)<-age_group_deaths_CROMER
    tmp3$sample<-1:nrow(tmp3)
    tmp3$risk_group<-rep("HR",nrow(tmp3))
    
    return(rbind(tmp2,tmp3))
  },.progress="text")
  
  gdf_AR<-melt(df_AR,id.vars=c("year","strain","sample","risk_group"),variable.name="age_group",value.name="incidence")
  
  tab_incidence_sample<-dcast(gdf_AR,formula=sample+age_group+year~strain+risk_group,value.var="incidence")
  
  cat("generate  samples for death (low risk)\n")
  load(paste(dir_Rsave,"tab_death_LR.R",sep="/"))
  tab_death_LR_sample<-ddply(tab_death_LR,c("age_group","year"),function(df){
    tmp<-data.frame(death_LR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_death_LR_sample,tab_incidence_sample,by=c("sample","age_group","year"))
  
  cat("generate  samples for death (high risk)\n")
  load(paste(dir_Rsave,"tab_death_HR.R",sep="/"))
  tab_death_HR_sample<-ddply(tab_death_HR,c("age_group","year"),function(df){
    tmp<-data.frame(death_HR=rnorm_pos1(SAMPLE_SIZE,df$mean,df$sd))
    tmp$sample<-1:SAMPLE_SIZE
    return(tmp)
  },.progress="text")
  
  tab_glm_sample<-merge(tab_death_HR_sample,tab_glm_sample,by=c("sample","age_group","year"))
  save(tab_glm_sample,file=paste(dir_res_Rsave,"tab_glm_death_sample.R",sep="/"))
  
  cat("end of COMPUTE_EXPLANATORY_VARIABLES_GLM_SAMPLE()\n\n")
  
}

COMPUTE_CFR_SAMPLES<-function(outcome="hospitalisations", strain="fluA", RG="LR", mf=F,mfac){

  cat("start of COMPUTE_CFR_SAMPLES()\n") 
  
  if(outcome=="GP_consultations")
    outcome_field=paste(outcome,strain,sep="_") else
      if(outcome=="hospitalisations")
        outcome_field=paste(outcome,strain,RG,sep="_")
          if(outcome=="death")
            outcome_field=paste(outcome,RG,sep="_")
  
  if(outcome=="death")
    load(paste(dir_res_Rsave,"tab_glm_death_sample.R",sep="/")) else
  load(paste(dir_res_Rsave,"tab_glm_sample.R",sep="/"))
  
  tab_glm<-tab_glm_sample
  tab_glm$season<-season_to_season_slash(tab_glm$year)
  tab_glm<-tab_glm[,-which(names(tab_glm)=="year")]
  tab_glm$age_group<-factor(tab_glm$age_group,levels=unique(tab_glm$age_group))
  
  #round the response variable to have an integer
  tab_glm[,outcome_field]<-round(tab_glm[,outcome_field])
  
  
  #compute formula and explanatory variables
  if(mf){
    if(strain=="fluA")
    {
      tab_glm$H1N1agg=tab_glm$H1N1_LR+mfac*tab_glm$H1N1_HR
      tab_glm$H3N2agg=tab_glm$H3N2_LR+mfac*tab_glm$H3N2_HR
      expformula="H1N1agg + H3N2agg"
      names.coeff<-c("H1N1agg","H3N2agg")     
    }
    else
    {
      mfac=rep(1.5,SAMPLE_SIZE)
      tab_glm$Bagg=tab_glm$B_LR+mfac*tab_glm$B_HR
      expformula="Bagg"
      names.coeff<-"Bagg"
      #expformula=paste("B_",RG,"-1",sep="")
    }
    
  }
  if(outcome=="hospitalisations"){
    if(strain=="fluA")
    {
      expformula=paste("H1N1_",RG,"+H3N2_",RG,sep="")
      names.coeff<-c(paste("H1N1_",RG,sep=""),paste("H3N2_",RG,sep=""))     
    }
      else
      {
      expformula=paste("B_",RG,sep="")
      names.coeff<-paste("B_",RG,sep="")
      #expformula=paste("B_",RG,"-1",sep="")
      }
    }
  if(outcome=="death"){
      expformula=paste("H1N1_",RG,"+H3N2_",RG,"+B_",RG,sep="")
      names.coeff<-c(paste("H1N1_",RG,sep=""),paste("H3N2_",RG,sep=""),paste("B_",RG,sep=""))     
  }
  fmla=as.formula(paste(outcome_field,expformula,sep="~"))
  
  #for each sample and age group perform a glm.nb
  cat(paste("compute glm :",outcome_field,"~",expformula,"\n"))
  glm_sample<-dlply(tab_glm,c("age_group","sample"),function(df){
    #df=tab_glm[(tab_glm$age_group=="0 - 6 m")&(tab_glm$sample==1),]
    #reg<-try(glm.nb(data=df,hospitalisations_fluA_LR~H1N1_LR+H3N2_LR-1,link="identity"))
    
    if(sum(df[,outcome_field])==0) return(0)
    
    reg<-try(glm.nb(data=df,fmla,link="identity"))
    #reg<-try(glm(data=df,fmla,family = poisson(link = "identity")))
    
    if(length(reg)>1){
      
      cfr<-coef(reg)[-1]
      interc<-coef(reg)[1]
      ind<-which(cfr<0)
      
      while(length(ind)>0){
        
        if(length(ind)==length(cfr)){
          cat("all coef negatives!!\n")
          #print(df)
          return(NA)
        }
        
        keep_strain<-names(cfr)[-ind]
        
        #if(interc>0)
        #{
        #  fmla<-as.formula(paste(outcome_field," ~",paste(keep_strain,collapse="+")))
        #}
        #else
        #{
        #  fmla<-as.formula(paste(outcome_field," ~",paste(keep_strain,collapse="+"),"-1")) 
        #}
        
        fmla<-as.formula(paste(outcome_field," ~",paste(keep_strain,collapse="+")))
        
        reg<-try(glm.nb(data=df,fmla,link="identity"))
        #reg<-try(glm(data=df,fmla,family = poisson(link = "identity")))
        
        if(length(reg)>1){
          cfr<-coef(reg)[-1]
          interc<-coef(reg)[1]
          ind<-which(cfr<0)
        }else{
          return(NA)
        }
      }
      
      #print(reg)
      return(reg)
      
    }else{
      return(NA)
    }
  },.progress="text")
  
  save(glm_sample,file=paste(dir_res_Rsave,"glm_sample.R",sep="/"))
  load(paste(dir_res_Rsave,"glm_sample.R",sep="/"))
  
  #extract coef (cfr distribution by age and strain)
  ind<-which(!is.na(glm_sample))
  #names.coeff<-names(coef(glm_sample[[min(ind)]]))[-1]
  tmp<-ldply(glm_sample[ind],function(reg){
    cfr<-rep(0,length(names.coeff))
    names(cfr)<-names.coeff
    if(is.list(reg))
    {
      co<-coef(reg)
      cfr[names.coeff]<-co[names.coeff]
      cfr[is.na(cfr)]=0      
    }
    return(cfr)
  })
  coef.tab<-which((names(tmp)!="(Intercept)")&(names(tmp)!=".id"))
  if(length(coef.tab)>1)  
    ind<-which((tmp$"(Intercept)"<0)|apply(tmp[,coef.tab]>1,1,prod)) else
    ind<-which((tmp$"(Intercept)"<0)|(tmp[,coef.tab]>1))
  if(length(ind))
    tmp<-tmp[-ind,]
  tmp$age_group<-extract_string(tmp$.id,".",1)
  tmp$sample<-extract_string(tmp$.id,".",2)
  df_cfr<-tmp[,c("age_group","sample",names.coeff)]
  save(df_cfr,file=paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  
  #
  gdf_cfr<-melt(df_cfr,id.vars=c("age_group","sample"),variable.name="strain","value.name"="cfr")
  
  tmp<-table(gdf_cfr[,c("age_group","strain")])
  print(tmp)
  #           strain
  # age_group   H1N1 H3N2   B
  #  0 - 14 y   921  921 921
  #  15 - 44 y  951  951 951
  #  45 - 64 y  990  990 990
  #  65+ y      966  966 966
  
  cat("end of COMPUTE_CFR_SAMPLE()\n")
  
}

epidate_day<-function(year,week)
{
  return(as.Date(paste("01/01/",year,sep=""),"%d/%m/%Y")+7*(week-1)+3)
}

COMPUTE_MF_GP_ANTON<-function()
{
  cat("COMPUTE_MF_GP_ANTON()\n")
  
  tab_mf_GP_sample=rnorm(n=SAMPLE_SIZE,mean=1.51,sd=0.18)
  
  save(tab_mf_GP_sample,file=paste(dir_Rsave,"tab_mf_GP_sample.R",sep="/"))
  cat("COMPUTE_MF_GP_ANTON()\n")
}

COMPARE_LABBASE_RCGP<-function()
{
  cat("COMPARE_LABBASE_RCGP\n")
  tab<-read.csv(paste(dir_data,"LabBaseRespiratory_Oct2011.csv",sep="/"))
 
  tab$midWeek<-epidate_day(year=tab$year,week=tab$week)
  
  set2plot=(tab$organism=="INFLUENZA B")&(tab$ageGp==1)
  plot(tab$midWeek[set2plot],tab$count[set2plot],type="h")
  
}

COMPUTE_HOSP<-function()
{
  cat("start of COMPUTE_HOSP()\n")
  
  #load mf sample
  load(paste(dir_Rsave,"tab_mf_GP_sample.R",sep="/"))
  
  tab_ply<-expand.grid(strain=strains,scenario=name_scenario)
  
  cat("generate cumulated incidence samples\n")
  df_AR<-ddply(tab_ply,c("scenario","strain"),function(df){
    
    strain<-df$strain
    scenario<-df$scenario
    file_name<-paste("Scenario",scenario,"final_size.txt",sep="_")
    
    all_years<-ldply(seasons,function(season){
      
      tab=read.table(paste(dir_res,season,strain,"scenarii",file_name,sep="/"),h=F)
      tmp<-MODEL_2_CROMER(tab)
      #sum low and high risk with multplicative factor
      tmp2<-tmp[,1:6]+tab_mf_GP_sample*tmp[,7:12]
      names(tmp2)<-age_group_CROMER
      tmp2$sample<-1:nrow(tmp2)
      tmp2$season<-rep(season,nrow(tmp2))
      return(tmp2)
      
    })
    
    #cumulated<-ddply(all_years,'sample',function(df) colSums(df[,age_group_CROMER]))
    
    return(all_years)
  },.progress="text")
  
  save(df_AR,file=paste(dir_res_Rsave,"df_AR.R",sep="/"))
  
  cat("end of COMPUTE_HOSP()\n")
}

GENERATE_PARAMETERS_SCENARIOS<-function()
{
  cat("start of GENERATE_PARAMETERS_SCENARIOS()\n")
  
  #temp<-matrix(sample(years,horizon_eval*n_simu,replace=T),ncol=10)
  parameters_MC<-data.frame(set=1:n_simu)
  
  parameters_MC<-ddply(parameters_MC,'set',function(df)
  {
    res<-paste(sample(years,1),sample.int(n_sample_inference,1),sep="-")
    for(i in 2:horizon_eval)
      res<-paste(res,paste(sample(years,1),sample.int(n_sample_inference,1),sep="-"),sep=",")
    return(res)
  },.progress="text")
  names(parameters_MC)<-c("set","epi_sample")
  
  parameters_MC<-data.frame(parameters_MC,cost.TIV=rtriangle(n_simu,a=12,b=20,c=15.55))   #mean price =15.85
  parameters_MC<-data.frame(parameters_MC,cost.LAIV=rtriangle(n_simu,a=12+6,b=20+6,c=15.55+6))   #mean price =15.85
  parameters_MC<-data.frame(parameters_MC,cost.adjuvanted=rtriangle(n_simu,a=12+6,b=20+6,c=15.55+6))   #mean price =15.85
  #parameters_MC<-data.frame(parameters_MC,febrile.case=rtriangle(n_simu,a=0.25,b=0.51,c=0.37))
  parameters_MC<-data.frame(parameters_MC,febrile.case=rtriangle(n_simu,a=0.309,b=0.513,c=0.396))
  parameters_MC<-data.frame(parameters_MC,ARI.case=rtriangle(n_simu,a=0.546,b=0.733,c=0.656)-parameters_MC$febrile.case) #TO CHANGE
  
  QALY.loss.AJ<-read.csv(paste(dir_data,"AJ_QALY_list.csv",sep="/"))
  parameters_MC<-data.frame(parameters_MC,QALY.loss.case=sample(QALY.loss.AJ$v,n_simu,replace=TRUE))
  
  parameters_MC<-data.frame(parameters_MC,QALY.loss.hospitalisation=rnorm(n_simu,mean=0.018,sd=0.0018))
  parameters_MC<-data.frame(parameters_MC,QALY.loss.ARI=rnorm(n_simu,mean= 0.37,sd=0.03061224)/365)
  parameters_MC<-data.frame(parameters_MC,hosp.cost=sapply(n_simu, rlnorm, meanlog=lnmu(839,192.1),sdlog=lnsig(839,192.1)))
  parameters_MC<-data.frame(parameters_MC,GP.cost=sapply(n_simu, rlnorm, meanlog=lnmu(37,8.4),sdlog=lnsig(37,8.4)))

  save(parameters_MC,file=paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  cat("end of GENERATE_PARAMETERS_SCENARIOS()\n")  
}

COMPUTE_DOSES_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_DOSES_SCENARIOS_TABLE()\n\n")
  
  #load coverage figures for each years and scenarios
  tab_coverage<-read.csv(paste(dir_res,"vaccine_coverage_scenarios.csv",sep="/"),header=T)
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  #load the population sizes
  AGP<-c(0.021,0.055,0.098,0.087,0.092,0.183,0.45) #proportion in high risk groups
  pop_by_seasons<-matrix(unlist(lapply(seasons,function(s){
    tab=read.table(paste(dir_res,s,"H3N2","age_groups_model.txt",sep="/"),h=F)})),ncol=14)
  
  pop_by_seasons<-t(rbind(pop_by_seasons*(1-AGP),pop_by_seasons*AGP)) #split by low and high risk 
  
  rownames(pop_by_seasons)=seasons
  colnames(pop_by_seasons)=names_age_risk_group
  
  save(pop_by_seasons,file=paste(dir_res_Rsave,"pop_by_seasons.R",sep="/"))  
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  cov_sst_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-df$scenario
    
    epi_tab<-matrix(as.integer(unlist(strsplit(strsplit(parameters_MC$epi_sample[set],",")[[1]],"-"))),nrow=2)
    
    #    tmp<-aaply(epi_tab,2,function(a){return(colSums(inc_yss_array[year=year_to_season(a[1]),,scenario,a[2],]))})
    tmp<-aaply(epi_tab,2,function(a){
      coverage=tab_coverage[(tab_coverage$year==year_to_season(a[1]))&(tab_coverage$scenario==scenario),3:16]
      doses=coverage*pop_by_seasons[year_to_season(a[1]),]
      doses.LAIV=sum(doses*c(0,1,1,rep(0,11)))
      doses.TIV=sum(doses*c(1,0,0,rep(1,11)))
      tmp2=c(doses.TIV,doses.LAIV)
      names(tmp2)=c("TIV","LAIV")
      return(tmp2)
    })
    return(apply(tmp,2,mean))
    
  },.progress="text")
  
  save(cov_sst_simus,file=paste(dir_res_Rsave,"cov_sst_simus.R",sep="/"))
  
  cat("end of COMPUTE_DOSES_SCENARIOS_TABLE()\n\n")
}

RETRIEVE_INCIDENCE_SAMPLE<-function(){
  cat("start of RETRIEVE_INCIDENCE_SAMPLE\n")
  #compute the incidence of cases predicted by the model, by age (CROMER death age group), risk and year (CROMER paper)
  #and scenarios, then save it in an array for further calculations
  
  tab_ply<-expand.grid(strain=strains,year=seasons,scenario=1:(n_scenario+2))
  #tab_ply<-expand.grid(strain=strains,year=seasons,scenario=1:5)
  
  name_scenario_file=paste("Scenario_",c("vaccination","no_vaccination",as.character(0:(n_scenario-1))),"_final_size.txt",sep="")
  
  cat("generate incidence samples\n")
  inc_yss_array<-daply(tab_ply,c("year","strain","scenario"),function(df){
    
    year<-df$year
    strain<-df$strain
    scenario<-df$scenario
    
    file_name<-name_scenario_file[scenario]
    tab=read.table(paste(dir_res,year,strain,"scenarii",file_name,sep="/"),h=F)
    tmp<-MODEL_2_CROMER(tab)
    
    return(as.matrix(tmp))
  },.progress="text")
    
  save(inc_yss_array,file=paste(dir_res_Rsave,"inc_yss_array.R",sep="/"))
  
  cat("end of RETRIEVE_INCIDENCE_SAMPLE()\n\n")
  
}

COMPUTE_INCIDENCE_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_INCIDENCE_SCENARIOS\n")
  
  #load the results of the inference and simulations follodwing methhod escribed in PLoS Med paper 
  load(paste(dir_res_Rsave,"inc_yss_array.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  inc_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-df$scenario
    
    epi_tab<-matrix(as.integer(unlist(strsplit(strsplit(parameters_MC$epi_sample[set],",")[[1]],"-"))),nrow=2)
    
    tmp<-aaply(epi_tab,2,function(a){return(colSums(inc_yss_array[year=year_to_season(a[1]),,scenario,a[2],]))})    
    
    return(colMeans(tmp))
  },.progress="text")
  
  save(inc_ss_simus,file=paste(dir_res_Rsave,"inc_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_INCIDENCE_SCENARIOS()\n\n")
}

COMPUTE_CASES_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_CASES_SCENARIOS_TABLE\n")

  #load parameters of the sets
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  #load the different incidences by set, scenarios and RAG
  load(paste(dir_res_Rsave,"inc_ss_simus.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  cases_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp=inc_ss_simus[set,scenario=scenario,]*parameters_MC$febrile.case[set]
    return(tmp)
  },.progress="text")
  
  save(cases_ss_simus,file=paste(dir_res_Rsave,"cases_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_CASES_SCENARIOS_TABLE()\n\n")
}

COMPUTE_QALY_LOSS_CASES_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_QALY_LOSS_CASES_SCENARIOS_TABLE\n")
  
  #load parameters of the sets
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  #load the different incidences by set, scenarios and RAG
  load(paste(dir_res_Rsave,"cases_ss_simus.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  QALY_cases_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp=cases_ss_simus[set,scenario=scenario,]*parameters_MC$QALY.loss.case[set]
    return(tmp)
  },.progress="text")
  
  save(QALY_cases_ss_simus,file=paste(dir_res_Rsave,"QALY_cases_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_QALY_LOSS_CASES_SCENARIOS_TABLE()\n\n")
}

COMPUTE_HOSP_RISKS<-function()
{
  cat("start of COMPUTE_HOSP_RISKS\n")
  
  COMPUTE_CFR_SAMPLES(outcome="hospitalisations",strain="fluA",RG="LR",mf=F)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_A_LR=df_cfr
  
  COMPUTE_CFR_SAMPLES(outcome="hospitalisations",strain="fluB",RG="LR",mf=F)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_B_LR=df_cfr
  
  COMPUTE_CFR_SAMPLES(outcome="hospitalisations",strain="fluA",RG="HR",mf=F)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_A_HR=df_cfr
  
  COMPUTE_CFR_SAMPLES(outcome="hospitalisations",strain="fluB",RG="HR",mf=F)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_B_HR=df_cfr
  
  #risk for B
  tmp<-merge(df_cfr_B_HR,df_cfr_B_LR)
  tmp<-melt(tmp, id=c("sample","age_group"))
  tmp$variable<-sapply(tmp$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp<-dcast(tmp,formula=sample~age_group+variable,value.var="value")
  tmp<-tmp[rowSums(is.na(tmp))==0,] #remove the lines with NA in them
  tmp<-tmp[sample(length(tmp[,1]),n_simu,replace=TRUE),]
  tab_risk_hosp_B<-as.matrix(tmp[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_hosp_B,file=paste(dir_res_Rsave,"tab_risk_hosp_B.R",sep="/"))
  
  #risk for H1N1
  tmp<-merge(df_cfr_A_HR,df_cfr_A_LR)
  tmp<-melt(tmp, id=c("sample","age_group"))
  tmp2=tmp[(tmp$variable=="H1N1_HR")|(tmp$variable=="H1N1_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tab_risk_hosp_H1N1<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_hosp_H1N1,file=paste(dir_res_Rsave,"tab_risk_hosp_H1N1.R",sep="/"))
  
  #risk for H3N2
  tmp2=tmp[(tmp$variable=="H3N2_HR")|(tmp$variable=="H3N2_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tab_risk_hosp_H3N2<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_hosp_H3N2,file=paste(dir_res_Rsave,"tab_risk_hosp_H3N2.R",sep="/"))
  
  cat("end of COMPUTE_HOSP_RISKS()\n\n")
}

COMPUTE_GP_RISKS<-function()
{
  cat("start of COMPUTE_GP_RISKS\n")
  
  load(paste(dir_Rsave,"tab_mf_GP_sample.R",sep="/"))
  mfac=tab_mf_GP_sample
  
  COMPUTE_CFR_SAMPLES(outcome="GP_consultations",strain="fluA",mf=T,mfac=mfac)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_A=df_cfr
  
  COMPUTE_CFR_SAMPLES(outcome="GP_consultations",strain="fluB",mf=T,,mfac=mfac)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_B=df_cfr
  
  #risk for B
  df_cfr_B$B_HR=mfac[as.integer(df_cfr_B$sample)]*df_cfr_B$Bagg
  df_cfr_B$B_LR=df_cfr_B$Bagg
  df_cfr_B=df_cfr_B[,-which(names(df_cfr_B)=="Bagg")]
  tmp<-melt(df_cfr_B, id=c("sample","age_group"))
  tmp$variable<-sapply(tmp$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp<-dcast(tmp,formula=sample~age_group+variable,value.var="value")
  tmp<-tmp[rowSums(is.na(tmp))==0,] #remove the lines with NA in them
  tmp<-tmp[sample(length(tmp[,1]),n_simu,replace=TRUE),]
  tab_risk_GP_B<-as.matrix(tmp[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_GP_B,file=paste(dir_res_Rsave,"tab_risk_GP_B.R",sep="/"))
  
  #divide in LR and HR for the A strains
  df_cfr_A$H1N1_HR=mfac[as.integer(df_cfr_A$sample)]*df_cfr_A$H1N1agg
  df_cfr_A$H1N1_LR=df_cfr_A$H1N1agg
  df_cfr_A$H3N2_HR=mfac[as.integer(df_cfr_A$sample)]*df_cfr_A$H3N2agg
  df_cfr_A$H3N2_LR=df_cfr_A$H3N2agg
  df_cfr_A=df_cfr_A[,-which((names(df_cfr_A)=="H1N1agg")|(names(df_cfr_A)=="H3N2agg"))]
  tmp<-melt(df_cfr_A, id=c("sample","age_group"))
  
  #risk for H1N1  
  tmp2=tmp[(tmp$variable=="H1N1_HR")|(tmp$variable=="H1N1_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tab_risk_GP_H1N1<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_GP_H1N1,file=paste(dir_res_Rsave,"tab_risk_GP_H1N1.R",sep="/"))
  
  #risk for H3N2
  tmp2=tmp[(tmp$variable=="H3N2_HR")|(tmp$variable=="H3N2_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tab_risk_GP_H3N2<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_GP_H3N2,file=paste(dir_res_Rsave,"tab_risk_GP_H3N2.R",sep="/"))
  
  cat("end of COMPUTE_GP_RISKS()\n\n")
}

COMPUTE_DEATH_RISKS<-function()
{
  cat("start of COMPUTE_DEATH_RISKS\n")
    
  COMPUTE_CFR_SAMPLES(outcome="death", RG="LR", mf=F,mfac)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_LR=df_cfr
  
  COMPUTE_CFR_SAMPLES(outcome="death", RG="HR", mf=F,mfac)
  load(paste(dir_res_Rsave,"df_cfr.R",sep="/"))
  df_cfr_HR=df_cfr
  
  tmp<-merge(df_cfr_LR,df_cfr_HR)
  tmp<-melt(tmp, id=c("sample","age_group"))
  
  #risk for H1N1
  tmp2=tmp[(tmp$variable=="H1N1_HR")|(tmp$variable=="H1N1_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tmp2$"0 - 6 m_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"6 m - 4 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"5 - 14 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"15 - 44 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"45 - 64 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"0 - 6 m_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"6 m - 4 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"5 - 14 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"15 - 44 y_HR"<-tmp2$"15 - 64 y_HR"
  tmp2$"45 - 64 y_HR"<-tmp2$"15 - 64 y_HR"
  tab_risk_death_H1N1<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_death_H1N1,file=paste(dir_res_Rsave,"tab_risk_death_H1N1.R",sep="/"))
  
  #risk for H3N2
  tmp2=tmp[(tmp$variable=="H3N2_HR")|(tmp$variable=="H3N2_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tmp2$"0 - 6 m_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"6 m - 4 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"5 - 14 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"15 - 44 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"45 - 64 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"0 - 6 m_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"6 m - 4 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"5 - 14 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"15 - 44 y_HR"<-tmp2$"15 - 64 y_HR"
  tmp2$"45 - 64 y_HR"<-tmp2$"15 - 64 y_HR"
  tab_risk_death_H3N2<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_death_H3N2,file=paste(dir_res_Rsave,"tab_risk_death_H3N2.R",sep="/"))
  
  #risk for B
  tmp2=tmp[(tmp$variable=="B_HR")|(tmp$variable=="B_LR"),]
  tmp2$variable<-sapply(tmp2$variable,function(str){strsplit(as.character(str),"_")[[1]][2]})
  tmp2<-dcast(tmp2,formula=sample~age_group+variable,value.var="value")
  tmp2<-tmp2[rowSums(is.na(tmp2))==0,] #remove the lines with NA in them
  tmp2<-tmp2[sample(length(tmp2[,1]),n_simu,replace=TRUE),]
  tmp2$"0 - 6 m_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"6 m - 4 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"5 - 14 y_LR"<-tmp2$"0 - 14 y_LR"
  tmp2$"15 - 44 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"45 - 64 y_LR"<-tmp2$"15 - 64 y_LR"
  tmp2$"0 - 6 m_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"6 m - 4 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"5 - 14 y_HR"<-tmp2$"0 - 14 y_HR"
  tmp2$"15 - 44 y_HR"<-tmp2$"15 - 64 y_HR"
  tmp2$"45 - 64 y_HR"<-tmp2$"15 - 64 y_HR"
  tab_risk_death_B<-as.matrix(tmp2[,c("0 - 6 m_LR","6 m - 4 y_LR","5 - 14 y_LR","15 - 44 y_LR","45 - 64 y_LR","65+ y_LR","0 - 6 m_HR","6 m - 4 y_HR","5 - 14 y_HR","15 - 44 y_HR","45 - 64 y_HR","65+ y_HR")])
  save(tab_risk_death_B,file=paste(dir_res_Rsave,"tab_risk_death_B.R",sep="/"))
  
  cat("end of COMPUTE_DEATH_RISKS()\n\n")
} 

COMPUTE_HOSPITALISATION_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_HOSPITALISATION_SCENARIOS\n")
  
  #load the results of the inference and simulations follodwing methhod escribed in PLoS Med paper 
  load(paste(dir_res_Rsave,"inc_yss_array.R",sep="/"))
  
  #load risks
  load(paste(dir_res_Rsave,"tab_risk_hosp_B.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_hosp_H1N1.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_hosp_H3N2.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  hosp_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    epi_tab<-matrix(as.integer(unlist(strsplit(strsplit(parameters_MC$epi_sample[set],",")[[1]],"-"))),nrow=2)
    
    tmp<-aaply(epi_tab,2,function(a){
      hosp_H3N2=inc_yss_array[year=year_to_season(a[1]),strain="H3N2",scenario,a[2],]*tab_risk_hosp_H3N2[set,]
      hosp_H1N1=inc_yss_array[year=year_to_season(a[1]),strain="H1N1",scenario,a[2],]*tab_risk_hosp_H1N1[set,]
      hosp_B=inc_yss_array[year=year_to_season(a[1]),strain="B",scenario,a[2],]*tab_risk_hosp_B[set,]
      #return(colSums(inc_yss_array[year=year_to_season(a[1]),,scenario,a[2],]))
      return(hosp_H3N2+hosp_H1N1+hosp_B)
      })    
    
    return(colMeans(tmp))
  },.progress="text")
  
  save(hosp_ss_simus,file=paste(dir_res_Rsave,"hosp_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_HOSPITALISATION_SCENARIOS()\n\n")
}

COMPUTE_GP_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_GP_SCENARIOS\n")
  
  #load the results of the inference and simulations follodwing methhod escribed in PLoS Med paper 
  load(paste(dir_res_Rsave,"inc_yss_array.R",sep="/"))
  
  #load risks
  load(paste(dir_res_Rsave,"tab_risk_GP_B.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_GP_H1N1.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_GP_H3N2.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  GP_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    epi_tab<-matrix(as.integer(unlist(strsplit(strsplit(parameters_MC$epi_sample[set],",")[[1]],"-"))),nrow=2)
    
    tmp<-aaply(epi_tab,2,function(a){
      GP_H3N2=inc_yss_array[year=year_to_season(a[1]),strain="H3N2",scenario,a[2],]*tab_risk_GP_H3N2[set,]
      GP_H1N1=inc_yss_array[year=year_to_season(a[1]),strain="H1N1",scenario,a[2],]*tab_risk_GP_H1N1[set,]
      GP_B=inc_yss_array[year=year_to_season(a[1]),strain="B",scenario,a[2],]*tab_risk_GP_B[set,]
      #return(colSums(inc_yss_array[year=year_to_season(a[1]),,scenario,a[2],]))
      return(GP_H3N2+GP_H1N1+GP_B)
    })    
    
    return(colMeans(tmp))
  },.progress="text")
  
  save(GP_ss_simus,file=paste(dir_res_Rsave,"GP_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_GP_SCENARIOS()\n\n")
}

COMPUTE_DEATH_SCENARIOS_TABLE<-function(scenarios){
  cat("start of COMPUTE_DEATH_SCENARIOS\n")
  
  #load the results of the inference and simulations follodwing methhod escribed in PLoS Med paper 
  load(paste(dir_res_Rsave,"inc_yss_array.R",sep="/"))
  
  #load risks
  load(paste(dir_res_Rsave,"tab_risk_death_B.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_death_H1N1.R",sep="/"))
  load(paste(dir_res_Rsave,"tab_risk_death_H3N2.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  death_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    epi_tab<-matrix(as.integer(unlist(strsplit(strsplit(parameters_MC$epi_sample[set],",")[[1]],"-"))),nrow=2)
    
    tmp<-aaply(epi_tab,2,function(a){
      death_H3N2=inc_yss_array[year=year_to_season(a[1]),strain="H3N2",scenario,a[2],]*tab_risk_death_H3N2[set,]
      death_H1N1=inc_yss_array[year=year_to_season(a[1]),strain="H1N1",scenario,a[2],]*tab_risk_death_H1N1[set,]
      death_B=inc_yss_array[year=year_to_season(a[1]),strain="B",scenario,a[2],]*tab_risk_death_B[set,]
      #return(colSums(inc_yss_array[year=year_to_season(a[1]),,scenario,a[2],]))
      return(death_H3N2+death_H1N1+death_B)
    })    
    
    return(colMeans(tmp))
  },.progress="text")
  
  save(death_ss_simus,file=paste(dir_res_Rsave,"death_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_DEATH_SCENARIOS()\n\n")
}

COMPUTE_QALY_DEATH_LOSS_TABLE<-function(scenarios){
  cat("start of COMPUTE_QALY_DEATH_LOSS_TABLE\n")
  
  #load the results of the inference and simulations follodwing method escribed in PLoS Med paper 
  load(paste(dir_res_Rsave,"death_ss_simus.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  death_QALY_loss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp<-death_ss_simus[set,scenario,]*death_QALY_loss     
    
    return(tmp)
  },.progress="text")
  
  save(death_QALY_loss_simus,file=paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  
  cat("end of COMPUTE_QALY_DEATH_LOSS_TABLE()\n\n")
}

COMPUTE_QALY_NON_DEATH_LOSS_TABLE<-function(scenarios){
  cat("start of COMPUTE_QALY_NON_DEATH_LOSS_TABLE\n")
  
  #load the hospitalisations and sympt cases
  load(paste(dir_res_Rsave,"hosp_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"cases_ss_simus.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  non_death_QALY_loss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp<-hosp_ss_simus[set,scenario,]*parameters_MC$QALY.loss.hospitalisation[set]+cases_ss_simus[set,scenario,]*parameters_MC$QALY.loss.case[set]     
    
    return(tmp)
  },.progress="text")
  
  save(non_death_QALY_loss_simus,file=paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  
  cat("end of COMPUTE_QALY_NON_DEATH_LOSS_TABLE()\n\n")
}

COMPUTE_QALY_ARI_LOSS_TABLE<-function(scenarios){
  cat("start of COMPUTE_QALY_ARI_LOSS_TABLE\n")
  
  load(paste(dir_res_Rsave,"cases_ss_simus.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  ARI_QALY_loss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp<-cases_ss_simus[set,scenario,]/parameters_MC$febrile.case[set]*parameters_MC$ARI.case[set]*parameters_MC$QALY.loss.ARI[set]     
    
    return(tmp)
  },.progress="text")
  
  save(ARI_QALY_loss_simus,file=paste(dir_res_Rsave,"ARI_QALY_loss_simus.R",sep="/"))
  
  cat("end of COMPUTE_QALY_ARI_LOSS_TABLE()\n\n")
}

COMPUTE_PROGRAMME_COST<-function(scenarios){
  cat("start of COMPUTE_PROGRAMME_COST\n")
  
  #load parameters of the sets
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  #load the number of doses
  load(paste(dir_res_Rsave,"cov_sst_simus.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  programme_cost_ss_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    #tmp=cov_sst_simus[set,scenario,1]*parameters_MC$cost.TIV[set]+cov_sst_simus[set,scenario,2]*parameters_MC$cost.LAIV[set]
    tmp=cov_sst_simus[set,scenario,1]*parameters_MC$cost.TIV[set]+cov_sst_simus[set,scenario,2]*parameters_MC$cost.TIV[set]
    return(tmp)
  },.progress="text")
  
  save(programme_cost_ss_simus,file=paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  
  cat("end of COMPUTE_PROGRAMME_COST()\n\n")
}

COMPUTE_HEALTH_CARE_COST<-function(scenarios){
  cat("start of COMPUTE_HEALTH_CARE_COST\n")
  
  #load the hospitalisations and sympt cases
  load(paste(dir_res_Rsave,"hosp_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"GP_ss_simus.R",sep="/"))
  
  #load parameters
  load(paste(dir_res_Rsave,"parameters_MC.R",sep="/"))
  
  tab_ply<-expand.grid(set=1:n_simu,scenario=scenarios)
  
  healthcare_cost_simus<-daply(tab_ply,c("set","scenario"),function(df){
    
    set<-df$set
    scenario<-as.character(df$scenario)
    
    tmp<-hosp_ss_simus[set,scenario,]*parameters_MC$hosp.cost[set]+GP_ss_simus[set,scenario,]*parameters_MC$GP.cost[set]     
    
    return(tmp)
  },.progress="text")
  
  save(healthcare_cost_simus,file=paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  
  cat("end of COMPUTE_HEALTH_CARE_ss_COST()\n\n")
}

SCATTER_PLOT<-function(add,incremental_on="10",scn_ext=11:17){
  
  load(paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  
  total_cost_ss<-programme_cost_ss_simus+apply(healthcare_cost_simus,c(1,2),sum)
  #diff.cov=(programme_cost_ss_simus+hc_ss-programme_cost_ss_simus[,10]-hc_ss[,10])/1e6
  diff.cov=(programme_cost_ss_simus-programme_cost_ss_simus[,"10"])/1e6
  diff.cov=(total_cost_ss-total_cost_ss[,incremental_on])/1e6
  min.cov.ext=min(diff.cov[,as.character(scn_ext)])
  max.cov.ext=max(diff.cov[,as.character(scn_ext)])
  
  flat.inc=apply(death_QALY_loss_simus+non_death_QALY_loss_simus,c(1,2),sum)
  diff.inc=flat.inc[,incremental_on]-flat.inc
  min.inc.ext=min(diff.inc[,as.character(scn_ext)])
  max.inc.ext=max(diff.inc[,as.character(scn_ext)])
  
  
  #plot(NULL,ylim=c(min.cov.ext,max.cov.ext),xlim=c(min.inc.ext,max.inc.ext))
  tiff(paste(dir_res_pdf,"CE_plan.tiff",sep="/"),height=16,width=16,units="cm",res=600,compression="lzw")
  #png(paste(dir_res_pdf,"CE_plan.png",sep="/"),height=16,width=16,units="cm",res=200)
  #cairo_pdf(paste(dir_res_pdf,"CE_plan.pdf",sep="/"),height=7,width=7)
  plot(NULL,ylim=c(-10,max.cov.ext),xlim=c(min.inc.ext,80000),ylab="Incremental cost in ??millions/year",xlab="QALY gained")
  
  abline(h=0,col=grey(.7))
  abline(v=0,col=grey(.7))
  abline(a=0,b=0.020000,col="#eeaa44ee")
  abline(a=0,b=0.030000,col="#eeaa44ee",lty=2)
  
  i=3
  for(scenario in as.character(scn_ext))
  {
    dominated=scn_ext[c(2,4)]
    if(scenario %in% dominated)
      pch.t=1 else pch.t=19
    points(rep(mean(diff.inc[,scenario]),100),rep(mean(diff.cov[,scenario]),100),col=col_scen_rein[i],pch=pch.t,cex=1.5,lwd=2)
    #points(diff.inc[,scenario],diff.cov[,scenario],col=col_scen[scenario],pch=19)
    #points(rep(mean(diff.inc[,scenario]),100),rep(mean(diff.cov[,scenario]),100),col=col_scen_rein[i],pch=19)
    #lines(c(0,mean(diff.inc[,scenario])),c(0,mean(diff.cov[,scenario])),col="grey")
    if(add=="contour")
    {
      d2d<-kde2d(diff.inc[,scenario],diff.cov[,scenario],n=100)
      level.cuts<-c(0.9)
      #contour(d2d$x,d2d$y,1-d2d$z/max(d2d$z),levels=level.cuts,col=col_scen_rein[scenario],add=TRUE)
      contour(d2d$x,d2d$y,1-d2d$z/max(d2d$z),levels=level.cuts,col=col_scen_rein[i],add=TRUE)
    }
    i=i+1
  }
  
  path.scenarios=as.character(scn_ext[-c(2,4)])
  arrows(0,0,mean(diff.inc[,path.scenarios[1]]),mean(diff.cov[,path.scenarios[1]]),col="#55555555",lwd=2,length=0.1,angle=20)
  for(i in 2:5)
    arrows(mean(diff.inc[,path.scenarios[i-1]]),mean(diff.cov[,path.scenarios[i-1]]),mean(diff.inc[,path.scenarios[i]]),mean(diff.cov[,path.scenarios[i]]),col="#55555555",lwd=2,length=0.1,angle=20)  
  #arrows(0,0,mean(diff.inc[,3]),mean(diff.cov[,3]))
  
  #arrows(mean(diff.inc[,2]),mean(diff.cov[,2]),mean(diff.inc[,3]),mean(diff.cov[,3]),col=grey(.5),lwd=2,lty=2)
  #arrows(mean(diff.inc[,path.scenarios[i-1]]),mean(diff.cov[,path.scenarios[i-1]]),mean(diff.inc[,path.scenarios[i]]),mean(diff.cov[,path.scenarios[i]]),col=grey(.5),lwd=2)
  
  points(0,0,col="#eeaa99ee",pch=15,cex=1.5)
  
  legend("bottomright",legend=c("2-4 years","50-64 years","5-16 years","2-4 & 50-64 years", "2-16 years", "2-16 & 50-64", "2-64 years"),col=col_scen_rein[3:9],pch=19, cex=.8)

  dev.off()  
}

QALY_BD_PLOT<-function()
{
  load(paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"ARI_QALY_loss_simus.R",sep="/"))
  
  risk_pop_cromer<-c(305340,2691769,5961073,19843631,10423303,4632786,6550,145442,647655,1974225,2334718,3790462)
  
  diff_ND_QALY<-non_death_QALY_loss_simus[,"10",]-non_death_QALY_loss_simus[,"15",]
  diff_ARI_QALY<-ARI_QALY_loss_simus[,"10",]-ARI_QALY_loss_simus[,"15",]
  diff_D_QALY<-death_QALY_loss_simus[,"10",]-death_QALY_loss_simus[,"15",]
  
  tiff(paste(dir_res_pdf,"QALD_distribution.tiff",sep="/"),height=22,width=17.9,units="cm",res=600,compression="lzw")
  #png(paste(dir_res_pdf,"QALD_distribution.png",sep="/"),height=27,width=22,units="cm",res=200)
  
  par(mfrow=c(3,1),cex=1.5,mar=c(2,5,2,2))
  mp<-barplot2(colMeans(diff_D_QALY)/risk_pop_cromer*365,ci.l=apply(diff_D_QALY,2,quantile,probs=.05)/risk_pop_cromer*365,ci.u=apply(diff_D_QALY,2,quantile,probs=.95)/risk_pop_cromer*365,ylim=c(0,.7),ylab="Death-associated\n QALD gain/year/person",col=c(rep("hotpink",6),rep("steelblue",6)), axisnames = FALSE,plot.grid = T,plot.ci=T,cex.names=5)
  text(1.2,0.64,labels="Low risk",cex=1.3)
  text(8.7,0.64,labels="High risk",cex=1.3)
  #axis(2,cex=.7)
  axis(1, tick=F, at = mp, labels = rep(age_group_CROMER,2), cex.axis = .45)

  #abline(v=7.3,col=grey(.7))
  mp<-barplot2(colMeans(diff_ND_QALY)/risk_pop_cromer*365,ci.l=apply(diff_ND_QALY,2,quantile,probs=.05)/risk_pop_cromer*365,ci.u=apply(diff_ND_QALY,2,quantile,probs=.95)/risk_pop_cromer*365,ylim=c(0,.7),ylab="Non fatal ILI-associated\n QALD gain/year/person",col=c(rep("hotpink",6),rep("steelblue",6)), axisnames = FALSE,plot.grid = T,plot.ci=T)
  text(1.2,0.64,labels="Low risk",cex=1.3)
  text(8.7,0.64,labels="High risk",cex=1.3)
  #axis(2,cex=.7)
  axis(1, tick=F, at = mp, labels = rep(age_group_CROMER,2), cex.axis = .45)
  mp<-barplot2(colMeans(diff_ARI_QALY)/risk_pop_cromer*365,ci.l=apply(diff_ARI_QALY,2,quantile,probs=.05)/risk_pop_cromer*365,ci.u=apply(diff_ARI_QALY,2,quantile,probs=.95)/risk_pop_cromer*365,ylim=c(0,.07),ylab="Non fatal ARI-associated\n QALD gain/year/person",col=c(rep("hotpink",6),rep("steelblue",6)), axisnames = FALSE,plot.grid = T,plot.ci=T)
  text(1.2,0.064,labels="Low risk",cex=1.3)
  text(8.7,0.064,labels="High risk",cex=1.3)
  #axis(2,cex=.7)
  axis(1, tick=F, at = mp, labels = rep(age_group_CROMER,2), cex.axis = .45)
  dev.off()
  
  load(paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  
  total_cost_ss<-programme_cost_ss_simus+apply(healthcare_cost_simus,c(1,2),sum)
  #t<-diff_ND_QALY[,3]
  t<-diff_ND_QALY[,3]+diff_ND_QALY[,2]*3/4.5+diff_ND_QALY[,4]*2/30
  c<-total_cost_ss[,15]-total_cost_ss[,10]
  d<-t*25000-c
  
  print(quantile(d,c(.5,0.05,0.95)))
  
  h<-d/25000/(cov_sst_simus[,15,2]+cov_sst_simus[,15,1]-cov_sst_simus[,10,2]-cov_sst_simus[,10,1])
  print(365*quantile(h,c(.5,0.05,0.95)))
}  

OUTPUT_TABLE<-function()
{
  output_scenario=c("10","2","11","12","13","14","15","16","17")
  
  final_table=NULL
  
  load(paste(dir_res_Rsave,"cases_ss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(cases_ss_simus,output_scenario)/1000)
  
  load(paste(dir_res_Rsave,"GP_ss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(GP_ss_simus,output_scenario))
  
  load(paste(dir_res_Rsave,"hosp_ss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(hosp_ss_simus,output_scenario))
  
  load(paste(dir_res_Rsave,"death_ss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(death_ss_simus,output_scenario))
  
  load(paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(non_death_QALY_loss_simus,output_scenario))
  
  load(paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(death_QALY_loss_simus,output_scenario))
  
  final_table=rbind(final_table,COL_LHM_SCENARIOS_TABLE(non_death_QALY_loss_simus+death_QALY_loss_simus,output_scenario))

  load(paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  
  pc_tab=COL_LHM_SCENARIOS_TABLE(programme_cost_ss_simus,output_scenario)/1000
  hc_tab=COL_LHM_SCENARIOS_TABLE(healthcare_cost_simus,output_scenario)/1000
  final_table=rbind(final_table,pc_tab,hc_tab,pc_tab+hc_tab)
   
  write.csv(final_table,file=paste(dir_res_Rsave,"final_outcome.csv",sep="/"))
}

COL_LHM_SCENARIOS_TABLE<-function(distribution_simus,output_scenario)
{
  tmp<-apply(distribution_simus,c(1,2),sum)
  tmp_mean<-apply(tmp,2,mean)[output_scenario]
  tmp_low<-apply(tmp,2,function(x){quantile(x,0.025)})[output_scenario]
  tmp_high<-apply(tmp,2,function(x){quantile(x,0.975)})[output_scenario]
  return(rbind(tmp_mean,tmp_low,tmp_high))
}

SENSITIVITY_COVERAGE<-function()
{
  dominated=c(2,4)
  s.names=c("2-4 years","50-64 years","5-16 years","2-4 & 50-64 years", "2-16 years", "2-16 & 50-64", "2-64 years")
  s.names.short=c("2-4","50-64","5-16","2-4 & 50-64", "2-16", "2-16 & 50-64", "2-64")
  load(paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  total_cost=apply(healthcare_cost_simus,c(1,2),sum)+programme_cost_ss_simus
  
  load(paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  total_QALY=apply(death_QALY_loss_simus+non_death_QALY_loss_simus,c(1,2),sum)
  
  scenario.inc.coverage=matrix(rep("0",28),ncol=4)
  scenario.inc.coverage[1,]=c("18","3","11","25")
  for(i in 1:6)
    scenario.inc.coverage[1+i,]=as.character(as.integer(scenario.inc.coverage[1,])+i)  
  
  tiff(paste(dir_res_pdf,"Sensitivity_to_coverage.tiff",sep="/"),height=10.16,width=19.05,units="cm",res=600,compression="lzw")
  #png(paste(dir_res_pdf,"Sensitivity_to_coverage.png",sep="/"),height=16,width=30,units="cm",res=200)
  par(mfrow=c(1,2))
  plot(NULL,xlim=c(10,100),ylim=c(0,1000),axes=F,ylab="Incremental net benefit over the current programme in ?M",xlab="Coverage of extension",cex.lab=0.7)
  i=1
  coverage=c(15,30,50,70)
  NB=rep(0,4)
  for(i in 1:7)
  {
    j=1
    for(scenario in scenario.inc.coverage[i,])
    {
      net.benefice<-25000*(total_QALY[,"10"]-total_QALY[,scenario])+(total_cost[,"10"]-total_cost[,scenario])    
      NB[j]=mean(net.benefice)/1e6
      j=j+1
    }
    if(i %in% dominated)
    {
      lt=2
      lc=1
    } else
    {
      lt=1
      lc=19
    }
    lines(coverage,NB,col=col_scen_rein[i+2],lwd=2,lty=lt)
    points(coverage,NB,lwd=2,col=col_scen_rein[i+2],pch=lc)
  }
  axis(1, tick=T, at = coverage, labels = c("15%","30%","50%","70%"), cex.axis = 0.7)
  axis(2)
  legend("topright",legend=s.names,col=col_scen_rein[3:9],pch=c(19,1,19,1,19,19,19), cex=.5, bty="n")
  
  NB=matrix(rep(0,3*7),ncol=3)
  j=1
  for(scenario in scenario.inc.coverage[,3])
  {
    net.benefice<-25000*(total_QALY[,"10"]-total_QALY[,scenario])+(total_cost[,"10"]-total_cost[,scenario])    
    NB[j,1]=mean(net.benefice)/1e6
    NB[j,c(2,3)]=quantile(net.benefice,c(0.05,0.95))/1e6
    j=j+1
  }
  
  mp<-barplot2(NB[,1],ci.l=NB[,2],ci.u=NB[,3],col=c(col_scen_rein[3],'white',col_scen_rein[5],'white',col_scen_rein[7:9]),border=col_scen_rein[3:9],axisnames = FALSE,plot.grid = T,plot.ci=T)
  axis(1, tick=F, at=mp, labels = s.names.short, cex.axis = .4)
  
  
  dev.off()
}

EVAL_TRANSITION<-function(scenario1="10",scenario2="15")
{
  load(paste(dir_res_Rsave,"programme_cost_ss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"healthcare_cost_simus.R",sep="/"))
  total_cost=apply(healthcare_cost_simus,c(1,2),sum)+programme_cost_ss_simus
  
  load(paste(dir_res_Rsave,"non_death_QALY_loss_simus.R",sep="/"))
  load(paste(dir_res_Rsave,"death_QALY_loss_simus.R",sep="/"))
  total_QALY=apply(death_QALY_loss_simus+non_death_QALY_loss_simus,c(1,2),sum)
  
  ICER<-mean(total_cost[,scenario2]-total_cost[,scenario1])/mean(total_QALY[,scenario1]-total_QALY[,scenario2])
  ICER.dis<-(total_cost[,scenario2]-total_cost[,scenario1])/(total_QALY[,scenario1]-total_QALY[,scenario2])
  
  prop.20k<-sum(ICER.dis<20000)/length(ICER.dis)
  prop.30k<-sum(ICER.dis<30000)/length(ICER.dis)
  
  net.benefice<-25000*(total_QALY[,scenario1]-total_QALY[,scenario2])+(total_cost[,scenario1]-total_cost[,scenario2])
  
  res<-data.frame(transition=paste(scenario1,"->",scenario2),ICER=ICER,p20k=prop.20k,p30k=prop.30k,NB=t(quantile(net.benefice/1e6,c(0.5,0.025,0.975))))
  
  return(res)
}

TABLE_TRANSITION<-function()
{
  path.scenarios=c("10","11","13","15","16","17")
  
  tab_res=NULL
  
  for(i in 2:length(path.scenarios))
    tab_res<-rbind(tab_res,EVAL_TRANSITION(path.scenarios[i-1],path.scenarios[i]))
  
  write.csv(tab_res,file=paste(dir_res_Rsave,"table_transition.csv",sep="/"))
  
  return(tab_res)  
}

TABLE_DISCOUNT<-function()
{
  path.scenarios=c("11","12","13","14","15","16","17")
  
  tab_res=NULL
  
  for(i in 1:length(path.scenarios))
    tab_res<-rbind(tab_res,EVAL_TRANSITION("10",path.scenarios[i]))
  
  write.csv(tab_res,file=paste(dir_res_Rsave,"/table_discount_",discount,".csv",sep=""))
  
  return(tab_res)
  
}

TABLE_RISK<-function()
{
  path.scenarios=c("41","42","43","44","45","46","47")
  
  tab_res=NULL
  
  for(i in 1:length(path.scenarios))
    tab_res<-rbind(tab_res,EVAL_TRANSITION("40",path.scenarios[i]))
  
  write.csv(tab_res,file=paste(dir_res_Rsave,"table_risk_75.csv",sep="/"))
  
  return(tab_res)
  
}

main<-function()
{
  #COMPARE_LABBASE_RCGP()
  #CROMER_BURDEN_DATA()
  #COMPUTE_EXPLANATORY_VARIABLES_GLM_SAMPLE()
  #COMPUTE_EXPLANATORY_DEATHS_GLM_SAMPLE()
  #COMPUTE_CFR_SAMPLES(outcome="hospitalisations",strain="fluA",RG="LR",mf=F)
  
  #COMPUTE_MF_GP_ANTON()
  
  #GENERATE_PARAMETERS_SCENARIOS()
  #RETRIEVE_INCIDENCE_SAMPLE()
  
  working.scenarios=1:54
  
  #COMPUTE_INCIDENCE_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_DOSES_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_CASES_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_QALY_LOSS_CASES_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_QALY_ARI_LOSS_TABLE(working.scenarios)
  #COMPUTE_HOSP_RISKS()
  #COMPUTE_GP_RISKS()
  #COMPUTE_DEATH_RISKS()
  #COMPUTE_HOSPITALISATION_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_GP_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_DEATH_SCENARIOS_TABLE(working.scenarios)
  #COMPUTE_QALY_DEATH_LOSS_TABLE(working.scenarios)
  #COMPUTE_QALY_NON_DEATH_LOSS_TABLE(working.scenarios)
  #COMPUTE_PROGRAMME_COST(working.scenarios)
  #COMPUTE_HEALTH_CARE_COST(working.scenarios)
  
  #SCATTER_PLOT(add="contour",incremental_on="10",scn_ext=11:17)
  #QALY_BD_PLOT()
  #SENSITIVITY_COVERAGE()
  #OUTPUT_TABLE()
  
  #TABLE_TRANSITION()
  #TABLE_DISCOUNT()
  #TABLE_RISK()

}

main()