###
### Calculate vaccine impact under high pop-at-risk scenario
###
library(maptools)
library(plyr)
setwd("..") ##Assumes code being run from "code" directory

##UN Pop data
un.pop=read.csv(file="data/un_pop_estimates.csv",header=T,as.is=T)
##Estimated pop at risk
c.df=read.csv("data/JEV_2020_pop_at_risk.csv",header=T,as.is=T)

###Use rice/wetland habitat 
c.df$pop.lo.pct=c.df$Pct20_rw  

##Country
args = commandArgs(trailingOnly=TRUE)
cname=args[1]
nruns=1000

##Scenario data
vacc.df=read.csv(file="data/coverage_campaign.csv",header=T,as.is=T)
vacc.df=vacc.df[vacc.df$country_code==cname,]


##Pop data
if(cname=="TWN"){
  tdat=read.csv(file="data/vacc_coverage_foi_by_year_age.csv")
  tdat=tdat[which(tdat$country_code=="TWN"&tdat$year<2020),]
  tdat$pop=tdat$pop*1000
  yrs=unique(tdat$year)

  tdat=tdat[,c(1:5,7:8)]
  names(tdat)[grep("year",names(tdat))]="yr"
  tdat$pop.pct=NA
  for(ti in yrs){
    tdat$pop.pct[tdat$yr==ti]=tdat$pop[tdat$yr==ti]/sum(tdat$pop[tdat$yr==ti])
  }
  tdat$pop.risk=tdat$pop*c.df$pop.lo.pct[c.df$ISO==cname]
  pop.df=tdat  
}else{
  cpop=un.pop[which(un.pop$country_code==cname),]
  yr.cols=grep("X",names(cpop))
  yrs=as.numeric(substring(names(cpop)[yr.cols],first=2))
  yr.cols=yr.cols[which(yrs>1949&yrs<2020)]
  yrs=yrs[which(yrs>1949&yrs<2020)] 
  
  pop.df1=cpop[,1:5]
  for(i in 1:length(yrs)){
    pop.df1$yr=yrs[i]
    pop.df1$pop=cpop[,yr.cols[i]]
    pop.df1$pop.pct=pop.df1$pop/sum(pop.df1$pop)
    ##Calculate at-risk portion
    pop.df1$pop.risk = pop.df1$pop*c.df$pop.lo.pct[c.df$ISO==cname]
    #pop.df$pop.hi.pct=
    if(i==1){
      pop.df=pop.df1
    }else{
      pop.df=rbind(pop.df,pop.df1)
    }
  }
}

##
## Parameter sets
##
param.df=read.csv(file="results/stochastic_template_params.csv",header=T,as.is=T)
symp.runs=param.df$symptomatic_probability
death.runs=param.df$mortality_probability
foi.runs=param.df[,which(names(param.df)==paste0("foi.",cname))]


##Start calculations
##
rout.l=list()
vacc.l=list()
novacc.l=list()
for(iter in 1:(nruns)){
    ###Calc infections - No vacc
    vimc.foi=foi.runs[iter]
    for(j in 1:length(yrs)){
      j.df=pop.df[pop.df$yr==yrs[j],]
      j.df$pr.inf=(1-exp(-vimc.foi))*exp(-vimc.foi*(j.df$age_from))
      j.df$infs=j.df$pr.inf*j.df$pop.risk
      if(j==1){
        novacc.df=j.df
      }else{
        novacc.df=rbind(novacc.df,j.df)
      }
    }
    novacc.df$run_id=iter
    
    ###Calc infections - Vaccine

    #First year of routine vaccination
    vmin.r=min(vacc.df$year[which(vacc.df$coverage>0&vacc.df$activity_type=="routine")],2020)
    
    ##Routine scenarios
    rout.df=pop.df
    rout.df$run_id=iter
    rout.df$susc.pct=NA
    rout.df$susc.pct[rout.df$yr<vmin.r]=exp(-vimc.foi*rout.df$age_from[rout.df$yr<vmin.r])
    if(vmin.r<2020){
      for(iyr in vmin.r:2019){
        yr.cover=vacc.df$coverage[vacc.df$year==iyr&vacc.df$activity_type=="routine"]
        yr.cover=ifelse(length(yr.cover)==0,0,yr.cover)
        ##Assume routine vaccination occurs at end of age 0
        rout.df$susc.pct[rout.df$yr==iyr&rout.df$age_from==0]=(1-yr.cover)
        rout.df$susc.pct[rout.df$yr==iyr&rout.df$age_from>0]=rout.df$susc.pct[rout.df$yr==(iyr-1)][-101]*exp(-vimc.foi)
      }      
    }

    ##Calculate infections
    rout.df$infs=rout.df$pop.risk*rout.df$susc.pct*(1-exp(-vimc.foi))

    ####
    ####Campaign
    ####
    vmin.c=min(vacc.df$year[which(vacc.df$coverage>0)],2020)
    
    ##Calc pct susceptible using s(a) from FOI
    camp.df=pop.df
    camp.df$run_id=iter
    camp.df$susc.pct=NA
    camp.df$susc.pct[camp.df$yr<vmin.c]=exp(-vimc.foi*camp.df$age_from[camp.df$yr<vmin.c])

    ####
    ##Add vaccination
    ####
    camp.df$susc.pop=camp.df$pop.risk*camp.df$susc.pct
    ##
    ## Campaigns
    ##
    campaign.yrs=NA
    if(any(grepl("campaign",vacc.df$activity_type))){
      ##Which years have campaigns
      cids=which(vacc.df$target>0&vacc.df$coverage>0)
      if(length(cids)>0){
        campaign.yrs=vacc.df$year[cids] 
        age.lo=vacc.df$age_first[cids]
        age.hi=vacc.df$age_last[cids]
        target.sizes=rep(NA,length(campaign.yrs))
        doses=target.sizes
        cover=rep(NA,length(doses))
        for(xx in 1:length(cids)){
          target.sizes[xx]=sum(camp.df$pop[camp.df$age_from %in% age.lo[xx]:age.hi[xx]&camp.df$yr==campaign.yrs[xx]])    
          doses[xx]=as.numeric(vacc.df$target[cids[xx]])*as.numeric(vacc.df$coverage[cids[xx]])
          cover[xx]=min(1,doses[xx]/target.sizes[xx])
        }
      }
    }
    if(vmin.c<2020){
      for(iyr in vmin.c:2019){
        yr.cover=vacc.df$coverage[vacc.df$year==iyr&vacc.df$activity_type=="routine"]
        yr.cover=ifelse(length(yr.cover)==0,0,yr.cover)
        camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from==0]=(1-yr.cover)
        ##Baseline susc from previous yr routine
        camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from>0]=camp.df$susc.pct[camp.df$yr==(iyr-1)][-101]*exp(-vimc.foi)
        if(iyr %in% campaign.yrs){
          c.ind=which(campaign.yrs==iyr)
          for(cctr in c.ind){
            ##Apply doses indiscriminantly to appropriate age class
            ##Assume that there is no targeting by previous vaccine status or whether actually live in area at risk
            camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from %in% age.lo[cctr]:age.hi[cctr]]=camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from %in% age.lo[cctr]:age.hi[cctr]]*(1-cover[cctr])
            ##Make sure doesn't go below 0        
            #camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from %in% age.lo[c.ind]:age.hi[c.ind]]=ifelse(camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from %in% age.lo[c.ind]:age.hi[c.ind]]<0,0,camp.df$susc.pct[camp.df$yr==iyr&camp.df$age_from %in% age.lo[c.ind]:age.hi[c.ind]])
            
          }
        }
      }      
    }
    ##Calculate infections
    camp.df$infs=camp.df$pop.risk*camp.df$susc.pct*(1-exp(-vimc.foi))
    
 
    ###CAses and Deaths
    rout.df$cases=rbinom(length(rout.df$infs),round(rout.df$infs,0),symp.runs[iter])
    camp.df$cases=rbinom(length(camp.df$infs),round(camp.df$infs,0),symp.runs[iter])
    novacc.df$cases=rbinom(length(novacc.df$infs),round(novacc.df$infs,0),symp.runs[iter])
    
    rout.df$deaths=rbinom(length(rout.df$infs),rout.df$cases,prob=death.runs[iter])
    camp.df$deaths=rbinom(length(camp.df$infs),camp.df$cases,prob=death.runs[iter])
    novacc.df$deaths=rbinom(length(novacc.df$infs),novacc.df$cases,prob=death.runs[iter])

    novacc.l[[iter]]=novacc.df
    rout.l[[iter]]=rout.df
    vacc.l[[iter]]=camp.df

    if(iter%%100==0) print(iter)
}
rout.dat=rbind.fill(rout.l)
vacc.dat=rbind.fill(vacc.l)
novacc.dat=rbind.fill(novacc.l)

###
### Calcuate DALYs
###
seyll=read.csv(file="data/SEYLL_GlobalDALYmethods_2000_2015.csv",header=T,as.is=T)
##merge seyll with datasets
rout.dat=merge(rout.dat,seyll,by.x="age_from",by.y="Age",all.x=T)
vacc.dat=merge(vacc.dat,seyll,by.x="age_from",by.y="Age",all.x=T)
novacc.dat=merge(novacc.dat,seyll,by.x="age_from",by.y="Age",all.x=T)

rout.dat$yll=rout.dat$deaths*rout.dat$SEYLL
rout.dat$yld=rout.dat$cases*0.133*(21/365)+rout.dat$cases*0.3*0.542+rout.dat$cases*0.2*0.203
vacc.dat$yll=vacc.dat$deaths*vacc.dat$SEYLL
vacc.dat$yld=vacc.dat$cases*0.133*(21/365)+vacc.dat$cases*0.3*0.542+vacc.dat$cases*0.2*0.203
novacc.dat$yll=novacc.dat$deaths*novacc.dat$SEYLL
novacc.dat$yld=novacc.dat$cases*0.133*(21/365)+novacc.dat$cases*0.3*0.542+novacc.dat$cases*0.2*0.203

rout.dat$dalys=rout.dat$yll+rout.dat$yld
vacc.dat$dalys=vacc.dat$yll+vacc.dat$yld
novacc.dat$dalys=novacc.dat$yll+novacc.dat$yld

##retrieve only necessary columns
vimc.f=function(dat){
  ###Remove pre-2010 - Can 
  dat=dat[dat$yr>2009,]
  #dat$disease="JE"
  dat=dat[,c("run_id","yr","age_from","country_code","country","pop","infs","deaths","cases","dalys","yll","yld")]
  names(dat)=c("run_id","year","age","country","country_name","cohort_size","infections","deaths","cases","dalys","yll","yld")
  dat=dat[order(dat$run_id,dat$year,dat$age),]
  return(dat)
}

rout.dat=vimc.f(rout.dat)
novacc.dat=vimc.f(novacc.dat)
vacc.dat=vimc.f(vacc.dat)


write.table(novacc.dat,paste("results/country_runs/je_no_vaccination_",cname,"_hi.csv",sep=""),sep=",",row.names = F,col.names = T)
write.table(rout.dat,paste("results/country_runs/je_routine-only_",cname,"_hi.csv",sep=""),sep=",",row.names = F,col.names = T)
write.table(vacc.dat,paste("results/country_runs/je_vaccination_",cname,"_hi.csv",sep=""),sep=",",row.names = F,col.names = T)
