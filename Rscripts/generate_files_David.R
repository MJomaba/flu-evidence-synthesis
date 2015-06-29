HOME="/Volumes/Transcend/CEA_paper_flu_kid/R"

age_group_MODEL<-c("0 - 1 y","1 - 4 y","5 - 14 y","15 - 24 y","25 - 44 y","45 - 64 y","65+ y")
names_age_risk_group<-paste(age_group_MODEL,rep(c("LR","HR"),each=length(age_group_MODEL)))
season.fnames<-c("199596","199697","199798","199899","199900","200001","200102","200203","200304","200405","200506","200607","200708","200809")

load.data<-function()
{
  cov.figures<-read.csv(paste(input.files.path,"coverage_figures.csv",sep=""))
  uptake.figures<-read.csv(paste(input.files.path,"uptake_by_month.csv",sep=""))
  
  year.match<-matrix(rep(TRUE,14*3),ncol=3)
  year.match<-read.csv(paste(input.files.path,"match_of_vaccine.csv",sep=""))
  colnames(year.match)<-c("Season","H3N2","H1N1","B")
  
  save(cov.figures,uptake.figures,year.match,file=paste(HOME,"Rsave/vaccine_env.R",sep="/"))
}
 
generate.vaccine.calendar.file<-function(coverage.scenarios, efficacy,name.file="vaccine_calendar.txt")
{
    load(file=paste(HOME,"Rsave/vaccine_env.R",sep="/"))
    
    n_scenarios=length(coverage.scenarios[,1])
    
    i<-which(uptake.figures$Year==year)
    
    cov.array.year.scen <- matrix(rep(0,16*((n_scenarios)*14)),ncol=16)
    colnames(cov.array.year.scen)<-c("year","scenario",names_age_risk_group)
    
    vacc.cal<-matrix(rep(0,123*21),ncol=21)
    new.vacc.cal<-matrix(rep(0,123*21),ncol=21)
    
    #load the coverage figures
    #cov.array.year.scen[(i-1)*(n_scenario+2)+n_s,]<-c(season.fnames[i],n_s,LR.cov/100,HR.cov/100)
    
    for(s in 1:length(coverage.scenarios[,1]))
    {
    LR.cov<-as.numeric(coverage.scenarios[s,2:8])
    HR.cov<-as.numeric(coverage.scenarios[s,9:15])
    
    #Octobre
    for(j in 1:31)
      new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake.figures[i,2],6),uptake.figures[i,3],rep(uptake.figures[i,2],6),uptake.figures[i,3],rep(uptake.figures[i,2],6),uptake.figures[i,3])/(31*100)  
    #Novembre
    for(j in 32:61)
      new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake.figures[i,4]-uptake.figures[i,2],6),uptake.figures[i,5]-uptake.figures[i,3],rep(uptake.figures[i,4]-uptake.figures[i,2],6),uptake.figures[i,5]-uptake.figures[i,3],rep(uptake.figures[i,4]-uptake.figures[i,2],6),uptake.figures[i,5]-uptake.figures[i,3])/(30*100)
    #Decembre
    for(j in 62:92)
      new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake.figures[i,6]-uptake.figures[i,4],6),uptake.figures[i,7]-uptake.figures[i,5],rep(uptake.figures[i,6]-uptake.figures[i,4],6),uptake.figures[i,7]-uptake.figures[i,5],rep(uptake.figures[i,6]-uptake.figures[i,4],6),uptake.figures[i,7]-uptake.figures[i,5])/(31*100)
    #Janvier
    for(j in 93:123)
      new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(100-uptake.figures[i,6],6),100-uptake.figures[i,7],rep(100-uptake.figures[i,6],6),100-uptake.figures[i,7],rep(100-uptake.figures[i,6],6),100-uptake.figures[i,7])/(31*100)
    
    if(s==1)
      app.flag=FALSE else
        app.flag=TRUE
    
      write(paste("#",s," ",coverage.scenarios$name.scenario[s],sep=""),file=name.file,append=app.flag)
      write(paste("#",s," Vaccine efficacy",sep=""),file=name.file,append=TRUE)
      
      write(paste(as.character(efficacy[s,]),collapse=' '),file=name.file,append=TRUE)
      write(paste("#",s," Vaccine coverage",sep=""),file=name.file,append=TRUE)
      write.table(signif(new.vacc.cal,4),file=name.file,row.names = FALSE,col.names = FALSE,append=TRUE)
    }
    
    #write the table containing the final coverage for each scenario, year and RAG
    #write.csv(cov.array.year.scen,file=paste(path.files,"vaccine_coverage_scenarios.csv",sep=""),row.names = FALSE)

    
}

name.scenario<-c("scenario1","scenario2")
cov<-t(matrix(c(rep(.3,6),.7,rep(.3,6),.7,rep(.5,6),.7,rep(.5,6),.7),nrow=14))

coverage.scenarios<-data.frame(cbind(name.scenario,cov))
eff.ex1=c(rep(.7,6),.5)
eff.ex2=c(rep(.5,6),.5)
eff.ex=rbind(eff.ex1,eff.ex2)

generate.vaccine.calendar.file(coverage.scenarios,efficacy=eff.ex,name.file="vaccine_calendar_test.txt")

