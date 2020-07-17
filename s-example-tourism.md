Data: Tourism flow
==================

Consider the tourism flow dataset:
<a href="http://www3.weforum.org/docs/WEF_TTCR19_data_for_download.xlsx" class="uri">http://www3.weforum.org/docs/WEF_TTCR19_data_for_download.xlsx</a>

    rm(list=ls())
    library(ggplot2)
    library(reshape2)
    dat=read.table("tourism-data.csv",header = TRUE,sep=";",dec=",",fileEncoding = "utf8",na.strings = "na")
    dat=dat[-(2:4),]
    var.name=dat$var.name
    names(var.name)=dat$varid
    dat=dat[,-(1:2)]
    dat=t(dat)
    colnames(dat)=names(var.name)
    dat=data.frame(dat)
    dat=na.omit(dat)
    (p=ncol(dat))

    ## [1] 65

    (n=nrow(dat))

    ## [1] 52

Common Functions
----------------

    post.baytauf=function(x){
      x=unlist(x[[1]])
      if(any(names(x)=="1")){x=x[names(x)=="1"]}else{x=0}
      return(x)
    }
    post.baytaue=function(x,v){
      x=unlist(x)
      x=x[names(x)=="1"]
      x[is.na(x)]=0
      
      return(sum(x*v))
    }
    dirprocanalysis=function(gammas) mixdir(gammas,n_latent = 10,select_latent = FALSE,alpha = 100)

    S=1000

Conventional Prior
------------------

    res.conventional <- GibbsBvs(formula = y ~ .,data = dat,n.iter = S,init.model = "Null",n.burnin =S*0.1,time.test= FALSE)

    ## In this dataset n<p and unitary Bayes factors are used for models with k>n.
    ## Info. . . .
    ## Most complex model has 65 covariates
    ## From those 1 are fixed and we should select from the remaining 64 
    ## ttindgdp, ttindempl, ttindgdpshare, ttindemplshare, constpermdays, constpermcost, startbusdays, startbuscost, corptaxrate, labortaxrate, proftaxrate, othtaxrate, terrorismincidence, homicidert, physdens, sanituse, wateruse, hospbed, hivprev, malariapc, enrol1net, enrol2gr, femlabor, netuserpct, bbsubpc, mobsubpc, mobbbsubpc, mobilenetworkcoverage, ttgovexp, compttdata, tmttdata, cntrybrandidx, visareq, openbilair, ftawto, taxaircharge, hotelcpi, ppp, fuelprice, pm25, envtreaty, waterstrs, thrspecies, forch, wasterwater, fssepi, dairseatkm, fairseatkm, deppop, airptdens...
    ## The problem has a total of 1.844674e+19 competing models
    ## Of these, 1100 are sampled with replacement
    ## Then, 1000 are kept and used to construct the summaries
    ## Working on the problem...please wait.

    renc = exp(res.conventional$modelslogBF[, ncol(res.conventional$modelslogBF)])
    gammas = res.conventional$modelslogBF[, -ncol(res.conventional$modelslogBF)]
    conv.emptau=apply(gammas,2,mean)
    conv.rentau=apply(gammas,2,function(x) sum(x*renc)/sum(renc))
    res.conv <- dirprocanalysis(gammas)
    v=res.conv$lambda
    baytauf=unlist(lapply(res.conv$category_prob,function(x) post.baytauf(x)))
    baytaue=unlist(lapply(res.conv$category_prob,function(x) post.baytaue(x,v)))
    names(baytauf)=names(baytaue)=colnames(gammas)
    ii=match(names(conv.emptau),names(baytauf))
    conv.baytauf=baytauf[ii]
    ii=match(names(conv.emptau),names(baytaue))
    conv.baytaue=baytaue[ii]

Non-local Prior
---------------

    res.nl <- modelSelection(y = dat$y,x = dat[,-1],scale=T, center=T,burnin = S*0.1,niter = S*1.1)

    ## Using default prior for continuous outcomes priorCoef=momprior(tau=0.348), priorVar=igprior(.01,.01)
    ## Greedy searching posterior mode... Done.
    ## Running Gibbs sampler.......... Done.

    renc = exp(res.nl$postProb)
    gammas = res.nl$postSample
    nl.emptau=apply(gammas,2,mean)
    nl.rentau=apply(gammas,2,function(x) sum(x*renc)/sum(renc))
    res.conv <- dirprocanalysis(gammas)
    v=res.conv$lambda
    baytauf=unlist(lapply(res.conv$category_prob,function(x) post.baytauf(x)))
    baytaue=unlist(lapply(res.conv$category_prob,function(x) post.baytaue(x,v)))
    names(baytauf)=names(baytaue)=colnames(gammas)
    ii=match(names(nl.emptau),names(baytauf))
    nl.baytauf=baytauf[ii]
    ii=match(names(nl.emptau),names(baytaue))
    nl.baytaue=baytaue[ii]

Median probability models
-------------------------

    res=data.frame(conv.emptau,conv.rentau,conv.baytaue,conv.baytauf,nl.emptau,nl.rentau,nl.baytaue,nl.baytauf)
    ii=match(rownames(res),names(var.name))
    row.names(res)=var.name[ii]
    med.mod=res>0.5
    res=res[rowSums(med.mod)>0,]
    med.mod=med.mod[rowSums(med.mod)>0,]
    oo=order(rowSums(med.mod),decreasing = TRUE)
    med.mod=med.mod[oo,]
    res=res[oo,]
    rtab=med.mod
    for(i in 1:nrow(rtab)) for(j in 1:ncol(rtab)) rtab[i,j]=paste(c("","*")[1+med.mod[i,j]],"(",c("",round(res[i,j],2))[1+med.mod[i,j]],")",sep="")
    rtab[rtab=="()"]=""
    colnames(rtab)=rep(c("empestincprob","renestincprob","bayestincprob","bayfestincprob"),2)
    rtab

    ##                                                               empestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.97)"    
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 ""           
    ## purchasing power parity                                       ""           
    ## natural tourism digital demand                                ""           
    ## openness of bilateral air service agreements                  ""           
    ## active mobile broadband internet subscriptions/100 population ""           
    ## total known species                                           ""           
    ##                                                               renestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(1)"       
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 "*(0.94)"    
    ## purchasing power parity                                       "*(0.95)"    
    ## natural tourism digital demand                                "*(0.71)"    
    ## openness of bilateral air service agreements                  "*(0.95)"    
    ## active mobile broadband internet subscriptions/100 population "*(0.82)"    
    ## total known species                                           "*(0.71)"    
    ##                                                               bayestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.95)"    
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 ""           
    ## purchasing power parity                                       ""           
    ## natural tourism digital demand                                ""           
    ## openness of bilateral air service agreements                  ""           
    ## active mobile broadband internet subscriptions/100 population ""           
    ## total known species                                           ""           
    ##                                                               bayfestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.99)"     
    ## number of operating airlines                                  "*(1)"        
    ## individuals using internet, %                                 "*(0.76)"     
    ## purchasing power parity                                       "*(1)"        
    ## natural tourism digital demand                                ""            
    ## openness of bilateral air service agreements                  "*(1)"        
    ## active mobile broadband internet subscriptions/100 population ""            
    ## total known species                                           ""            
    ##                                                               empestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.96)"    
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 "*(0.68)"    
    ## purchasing power parity                                       "*(0.61)"    
    ## natural tourism digital demand                                "*(0.61)"    
    ## openness of bilateral air service agreements                  ""           
    ## active mobile broadband internet subscriptions/100 population ""           
    ## total known species                                           ""           
    ##                                                               renestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(1)"       
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 ""           
    ## purchasing power parity                                       ""           
    ## natural tourism digital demand                                ""           
    ## openness of bilateral air service agreements                  ""           
    ## active mobile broadband internet subscriptions/100 population ""           
    ## total known species                                           ""           
    ##                                                               bayestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.87)"    
    ## number of operating airlines                                  "*(1)"       
    ## individuals using internet, %                                 "*(0.61)"    
    ## purchasing power parity                                       "*(0.55)"    
    ## natural tourism digital demand                                "*(0.55)"    
    ## openness of bilateral air service agreements                  ""           
    ## active mobile broadband internet subscriptions/100 population ""           
    ## total known species                                           ""           
    ##                                                               bayfestincprob
    ## timeliness of providing monthly/quarterly t&t data            "*(0.95)"     
    ## number of operating airlines                                  "*(1)"        
    ## individuals using internet, %                                 "*(0.79)"     
    ## purchasing power parity                                       "*(1)"        
    ## natural tourism digital demand                                "*(0.75)"     
    ## openness of bilateral air service agreements                  "*(1)"        
    ## active mobile broadband internet subscriptions/100 population "*(0.56)"     
    ## total known species                                           ""

    # xtable::xtable(rtab)
