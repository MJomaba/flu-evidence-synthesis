vaccine.calendar <- function( coverage, efficacy, uptake )
{
    new.vacc.cal<-matrix(rep(0,123*21),ncol=21)

    LR.cov<-as.numeric(coverage[1:7])
    HR.cov<-as.numeric(coverage[8:14])

    #October
    for(j in 1:31)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[2],6),uptake[3],rep(uptake[2],6),uptake[3],rep(uptake[2],6),uptake[3])/(31*100)  
    #November
    for(j in 32:61)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3],rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3],rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3])/(30*100)
    #December
    for(j in 62:92)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5],rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5],rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5])/(31*100)
    #Januari
    for(j in 93:123)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(100-uptake[6],6),100-uptake[7],rep(100-uptake[6],6),100-uptake[7],rep(100-uptake[6],6),100-uptake[7])/(31*100)

    cal <- list( "efficacy"=efficacy, "calendar"=new.vacc.cal )
    cal
}
