##
## Script reads in age-specific incidence data for each country and estimates FOI
##
library(rstan)
library(rjson)
setwd("results")
ccdf=read.csv(file="../data/country_codes.csv")
ccdf=ccdf[order(ccdf$NAME),]

inf.l=fromJSON(file="../data/incidence_age_data_2020.json")
##Remove any NULL objects
inf.l=inf.l[!unlist(lapply(inf.l,is.null))]

#pct pop in each age class
pop.dat=read.csv(file="../data/un_pop_estimates.csv",header=T,as.is=T)

##Vaccination coverage dataset
vacc.dat=read.csv(file="../data/vacc_coverage_foi_by_year_age.csv")
### Loop through and get FOI estimate for each dataset
### - 
##Country
foi.dist=list()
vacc.est=list()
exp.cases=list()
warn=list()

for(j in 1:length(inf.l)){
  jdat=inf.l[[j]]
  if(is.null(jdat)){
    next()
  }
  cname=ifelse(jdat$ctry=="CAM","KHM",jdat$ctry)
  ##Use pop distrubtion from final year of study
  if(cname!="TWN"){
    jpop=pop.dat[pop.dat$country_code==cname,]   
    study.yr=jdat$year.e
    study.cn=paste0('X',study.yr)
    pop10=jpop[1:nrow(jpop),which(names(jpop)==study.cn)]
    pop10.pct=pop10/sum(pop10)  
  }else{
    jpop=vacc.dat[vacc.dat$country_code==cname,]
    study.yr=jdat$year.e
    pop10=jpop$pop[jpop$year==study.yr]*1000
    pop10.pct=pop10/sum(pop10)      
  }

  pop.pct.agec=rep(NA,length(jdat$age.c))
  for(i in 1:length(jdat$age.c)){
    sind=jdat$age.c[i]+1 #Account for indexing
    eind=jdat$upper.c[i]+1
    pop.pct.agec[i]=sum(pop10.pct[sind:eind])
  }
  pop.pct.agec=pop.pct.agec/sum(pop.pct.agec) #Pct of study pop (not total pop) in age class 
  
  #Run model with no vaccination coverage if dataset is pre-vaccine
  if(!jdat$vacc){
  #if(!cname %in% c("CHN","THA","VNM")){
      
      model.data=list(num_ages=length(jdat$age.c),
                      age_range=length(jdat$age.c[1]:jdat$upper.c[length(jdat$upper.c)]),
                      lower_bound=jdat$age.c,
                      upper_bound=jdat$upper.c,
                      cases_age=round(jdat$infs.a,0),
                      age_distr=pop.pct.agec)
    
      mfit<- stan(file="../code/constant_foi_age_groups.stan",
                  data=model.data,
                  #init=init1,
                  #pars=c("betaA"),
                  chains=4,
                  iter=10000,
                  refresh=1000,
                  control=list(adapt_delta=0.99,max_treedepth=20))
      foi.dist[[j]]=rstan::extract(mfit,permuted=T)$lambdaA
      exp.cases[[j]]=rstan::extract(mfit,permuted=T)$exp_cases_age
  }else{
    ##Need to estimate vaccination coverage in addition to FOI
    vacc.j=vacc.dat[vacc.dat$country_code==cname&vacc.dat$year==study.yr,]
    vacc.j$cover_pop_yr=vacc.j$pop*vacc.j$cover_pct_tot_yr
    vacc.pct.agec=rep(NA,length(jdat$age.c))
    ##Get % covered in each age class
    for(i in 1:length(jdat$age.c)){
      sind=jdat$age.c[i]+1 #Account for indexing
      eind=jdat$upper.c[i]+1
      vacc.pct.agec[i]=sum(vacc.j$cover_pop_yr[sind:eind])/sum(vacc.j$pop[sind:eind])
    }
    ##Correct vaccine coverage for efficacy
    ve=0.993
    phiV=5
    vacc.pct.agec=vacc.pct.agec #*ve
    ##Set priors for vaccine coverage
    vacc.pct.agec=ifelse(vacc.pct.agec==0,0.005,vacc.pct.agec)
    vacc.alphas=phiV*vacc.pct.agec
    vacc.betas=phiV*(1-vacc.pct.agec)

    model.data=list(num_ages=length(jdat$age.c),
                    age_range=length(jdat$age.c[1]:jdat$upper.c[length(jdat$upper.c)]),
                    lower_bound=jdat$age.c,
                    upper_bound=jdat$upper.c,
                    cases_age=round(jdat$infs.a,0),
                    age_distr=pop.pct.agec,
                    ve=0.993,
                    alpha_pr=vacc.alphas,
                    beta_pr=vacc.betas)
                    #vacc_cov=c(0.9,0.9,0.9,0.2,0.1,0.1)) #rev(seq(0.1,0.8,length=6)))
                    #vacc_cov=rep(0.50,length(jdat$age.c)))
    
    mfit<- stan(file="../code/constant_foi_age_groups_vaccination_est.stan",
                data=model.data,
                #init=init1,
                #pars=c("betaA"),
                chains=4,
                iter=10000,
                refresh=1000,
                control=list(adapt_delta=0.99,max_treedepth=15)) ##16, 24, and 27 had rhat>1.01 at d=.95 and tree=10
    foi.dist[[j]]=rstan::extract(mfit,permuted=T)$lambdaA
    vacc.est[[j]]=rstan::extract(mfit,permuted=T)$vacc_cov
    exp.cases[[j]]=rstan::extract(mfit,permuted=T)$exp_cases_age
  }
  saveRDS(mfit,file=paste0("stan_fit_",cname,"_",j,".rds"))
  print(j)
  print(warnings())
  warn[[j]]=warnings()
  assign("last.warning", NULL, envir = baseenv())
}
save(foi.dist,vacc.est,exp.cases,file="JEV_foi_vaccination_stan_estimates.rda")


# ###
# ### Analysis/figures of FOI results
# ###
# library(data.table)
# library(tidyverse)
# sim_ind=1:length(foi.dist[[1]])
# foi_df=rbindlist(lapply(foi.dist,function(x) (as.data.frame(t(x),col.names=sim_ind))))
# foi_df=cbind(unlist(lapply(inf.l,function(x) x$ctry)),unlist(lapply(inf.l,function(x) x$year.s)),
#              unlist(lapply(inf.l,function(x) x$year.e)),foi_df)
# names(foi_df)=c("ctry","start_yr","end_yr",paste("sim",sim_ind,sep="_"))
# foi_df$model_num=1:nrow(foi_df)
# 
# foi_df_l=pivot_longer(foi_df,cols=starts_with("sim"),names_to="sim",
#                       values_to="foi")
# ###Save FOI table for estimate of infections
# write_csv(foi_df_l,path="foi_stan_fits.csv",col_names = T)
# 
# pdf(file="figures/foi_estimates_violin.pdf",width=8,height=5.5)
# ggplot(foi_df_l,aes(x=ctry,y=foi))+geom_violin(fill="lightblue")+
#   theme(axis.text.x=element_text(angle=60,hjust=1,size=14),axis.text.y=element_text(size=14),axis.title=element_text(size=16),
#         panel.background = element_rect(fill="white",color="grey50"))+
#   #theme_bw(base_size = 16)+
#   labs(x="Country",y="Force of infection (FOI)")+
#   scale_x_discrete(breaks=ccdf$ISO3,labels=ccdf$NAME) #,limits=ccdf$ISO3[order(ccdf$NAME)])
# dev.off()
# 
# mean(foi_df_l$foi)
# quantile(foi_df_l$foi,c(0.025,0.5,0.975))
# 
# foi_sum_l=foi_df_l %>% group_by(ctry) %>%
#   summarize(med=quantile(foi,0.5),
#                        lo=quantile(foi,0.025),
#                        hi=quantile(foi,0.975))
# 
# ##Also group by study (model_num)
# foi_model_sum=foi_df_l %>% group_by(ctry,model_num,start_yr,end_yr) %>%
#   summarize(mean=mean(foi),
#             med=quantile(foi,0.5),
#             lo=quantile(foi,0.025),
#             hi=quantile(foi,0.975))
# 
# ###Look at vaccination coverage estimates
# vsets=which(unlist(lapply(vacc.est,function(x) !is.null(x))))
# 
# ##vaccine coverage age groups
# age_l=list()
# age_l[[1]]=c("<1","1-4","5-14","15-19","20-34",expression("">=35))
# age_l[[2]]=c(expression(""<=15),">15")
# age_l[[3]]=c("<1","1-14",expression("">=15))
# age_l[[4]]=c("<1","1-14",expression("">=15))
# age_l[[5]]=c("<1","1-4","5-14",expression("">=15))
# age_l[[7]]=c("0-6","7-12","13-20","21-40","41-59",expression("">=60))
# age_l[[6]]=c("0-4","5-9","10-14",expression("">=15))
# age_l[[8]]=c("0-4","5-14",expression("">=15))
# age_l[[9]]=c("<1",1:15)
# age_l[[10]]=c("<1","1-5","6-10","11-15")
# age_l[[11]]=c("<1","1-5","6-10","11-15",">15")
# age_l[[12]]=c("0-4","5-9","10-11")
# age_l[[13]]=c("<1","1-4","5-9","10-15")
# age_l[[14]]=c("0-4","5-14","15-24","25-59",expression("">=60))
# age_l[[15]]=c("0-12",expression("">=13))
# age_l[[16]]=c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59",expression("">=60))
# age_l[[17]]=c("0-14",expression("">=15))
#   age_l[[18]]=c("0-14",expression("">=15))
#   age_l[[19]]=c("0-14",expression("">=15))
#   age_l[[20]]=c("0-1","2-4","5-9","10-14","15-17",expression("">=18))
#   age_l[[21]]=c("0-5","6-10","11-17","18-25","26-40","41-50","51-59",expression("">=60))
#   age_l[[22]]=c("<1","1-5","6-10","11-15",expression("">15))
# age_l[[23]]=c("0-15","16-39",expression("">=40))
# age_l[[24]]=c("0-5","6-14","15-24",expression("">=25))
# age_l[[25]]=c("0-4","5-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79",expression("">=80))
# age_l[[26]]=c("0-9","10-19","20-29","30-39","40-49","50-59",expression("">=60))
# age_l[[27]]=c("<1","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59",expression("">=60))
# age_l[[28]]=c("0-9","10-19","20-29","30-39","40-49","50-59",expression("">=60))
# age_l[[29]]=c("<4","4-7","8-11","12-15","16-19","20-23","24-27","28-31","32-35","36-39",
#               "40-43","44-47","48-51","52-55","56-59","60-63","64-67","68-71","72-75","76-79",expression("">=80))
# 
# foi_model_sum$vacc=ifelse(foi_model_sum$model_num %in% vsets,T,F)
# write_csv(foi_model_sum,path="foi_estimates_study_stats.csv")
# 
# vacc_l=list()
# for(vi in vsets){
#   v.inf=inf.l[[vi]]
#   v.est=vacc.est[[vi]]
#   cname=ifelse(v.inf$ctry=="CAM","KHM",v.inf$ctry)
#   ##Use pop distrubtion from final year of study
#   if(cname!="TWN"){
#     vpop=pop.dat[pop.dat$country_code==cname,]   
#     study.yr=v.inf$year.e
#     study.cn=paste0('X',study.yr)
#     pop10=vpop[1:nrow(vpop),which(names(vpop)==study.cn)]
#     pop10.pct=pop10/sum(pop10)  
#   }else{
#     vpop=vacc.dat[vacc.dat$country_code==cname,]
#     study.yr=v.inf$year.e
#     pop10=vpop$pop[vpop$year==study.yr]*1000
#     pop10.pct=pop10/sum(pop10)      
#   }
#   
#   pop.pct.agec=rep(NA,length(v.inf$age.c))
#   for(i in 1:length(v.inf$age.c)){
#     sind=v.inf$age.c[i]+1 #Account for indexing
#     eind=v.inf$upper.c[i]+1
#     pop.pct.agec[i]=sum(pop10.pct[sind:eind])
#   }
#   pop.pct.agec=pop.pct.agec/sum(pop.pct.agec) #Pct of study pop (not total pop) in age class 
#   ##Need to estimate vaccination coverage in addition to FOI
#   vacc.j=vacc.dat[vacc.dat$country_code==cname&vacc.dat$year==study.yr,]
#   vacc.j$cover_pop_yr=vacc.j$pop*vacc.j$cover_pct_tot_yr
#   vacc.pct.agec=rep(NA,length(v.inf$age.c))
#   ##Get % covered in each age class
#   for(i in 1:length(v.inf$age.c)){
#     sind=v.inf$age.c[i]+1 #Account for indexing
#     eind=v.inf$upper.c[i]+1
#     vacc.pct.agec[i]=sum(vacc.j$cover_pop_yr[sind:eind])/sum(vacc.j$pop[sind:eind])
#   }
#   ##Correct vaccine coverage for efficacy
#   ve=0.993
#   phiV=5
#   vacc.pct.agec=vacc.pct.agec*ve
#   ##Set priors for vaccine coverage
#   vacc.pct.agec=ifelse(vacc.pct.agec==0,0.005,vacc.pct.agec)
#   vacc.alphas=phiV*vacc.pct.agec
#   vacc.betas=phiV*(1-vacc.pct.agec)
#   
#   model.data=list(num_ages=length(v.inf$age.c),
#                   age_range=length(v.inf$age.c[1]:v.inf$upper.c[length(v.inf$upper.c)]),
#                   lower_bound=v.inf$age.c,
#                   upper_bound=v.inf$upper.c,
#                   cases_age=round(v.inf$infs.a,0),
#                   age_distr=pop.pct.agec,
#                   ve=ve,
#                   alpha_pr=vacc.alphas,
#                   beta_pr=vacc.betas)
#   
#   colnames(v.est)=paste0("age_",1:ncol(v.est))
#   v_est_l=pivot_longer(as.data.frame(v.est),cols=starts_with("age_"),names_to="age_class",
#                        values_to="cov")
#   v_est_l$dist="Posterior"
#   va.l=list()
#   for(aind in 1:length(vacc.alphas)){
#     va.pr=rbeta(10000,vacc.alphas[aind],vacc.betas[aind]) 
#     va.l[[aind]]=data.frame(age_class=paste0("age_",aind),cov=va.pr,dist="Prior")
#   }
#   v_prior=rbindlist(va.l)
#   
#   v_dist=bind_rows(v_est_l,v_prior)
#   v_dist$dist=factor(v_dist$dist,levels=c("Prior","Posterior"))
#   v_dist$age_class=as.integer(gsub("age_","",v_dist$age_class))
# 
#   vplot=ggplot(v_dist,aes(x=as.factor(age_class),y=cov,fill=dist))+geom_boxplot()+
#     theme(axis.text.x=element_text(angle=30,hjust=1,size=14),axis.text.y=element_text(size=14),axis.title=element_text(size=16),
#           panel.background = element_rect(fill="white",color="grey50"))+
#     labs(x="Age class",y="Vaccine coverage",fill="Distribution")+
#     scale_x_discrete(labels=age_l[[vi]])
#   ggsave(vplot,filename=paste0("figures/vaccine_coverage_",v.inf$ctry,"_",vi,"_upd.pdf"),width=8,height=6.5,units="in")
#   vacc_l[[vi]]=v_dist
# }
# 
# ### Plots of expected vs. observed cases
# 
# for(cc in 1:length(inf.l)){
#   cases.c=inf.l[[cc]]$infs.a
#   c.est=exp.cases[[cc]]
#   colnames(c.est)=paste0("age_",1:ncol(c.est))
# 
#   c_est_l=pivot_longer(as.data.frame(c.est),cols=starts_with("age_"),names_to="age_class",
#                        values_to="cases")
#   c_sum=data.frame(age_class="sum",cases=rowSums(c.est))
#   c_est_l=bind_rows(c_est_l,c_sum)
#   cases.df=data.frame(cases=c(cases.c,sum(cases.c)),age=unique(c_est_l$age_class))
#   cplot=ggplot(c_est_l,aes(x=as.factor(age_class),y=cases))+geom_boxplot()+
#     theme(axis.text.x=element_text(angle=30,hjust=1,size=14),axis.text.y=element_text(size=14),axis.title=element_text(size=16),
#           panel.background = element_rect(fill="white",color="grey50"))+
#     geom_point(cases.df,mapping=aes(x=as.factor(age),y=cases),pch=18,color="red",size=3)+
#     labs(x="Age class",y="JE cases")+ scale_x_discrete(labels=c(age_l[[cc]],"Total"))
#   ggsave(cplot,filename=paste0("figures/estimated_cases/estimateed_cases_",inf.l[[cc]]$ctry,"_",cc,".pdf"),width=8,height=6.5,units="in")
#   #vacc_l[[vi]]=v_dist
#   print(cc)
#   print(paste0("Reported cases: ",sum(cases.c)))
#   print(paste0("Estimated cases: ", quantile(c_sum$cases,c(0.025,.5,.975))))
# }













# ### Save foi distributions for each country
# cfoi=unique(unlist(lapply(inf.l,function(x) x$ctry)))
# for(ic in cfoi){
#   ic.foi=which(foi_df$ctry==ic) #which(lapply(inf.l,function(x) x$ctry)==ic)
#   foi_post=foi_df[ic.foi,]
#   foi_post=unlist(foi_post[,-1])
#   save(foi_post,file=paste("~/jev/output/",ic,"_foi_mfit_rev_2019.rda",sep=""))
# }

#save(foi.dist,file=paste("~/jev/output/",cname,"_foi_mfit_upd.rda",sep=""))



