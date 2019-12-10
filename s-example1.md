Set up
======

    rm(list=ls()) 
    set.seed(11) 
    library(BayesVarSel) 

    ## Loading required package: MASS

    ## Loading required package: mvtnorm

    ## Loading required package: parallel

    library(mixdir) 
    p=200 
    n=50 
    nsim.gibbs=10000 
    m1=3 
    signaltonoise=3 

Suppose we have 200 covariates and 50 iid observations from a regression
model where intercept and 197 coefficients are zero while the first 3
covariates are 3 and errors with variance 1.

    x=cbind(1,array(rnorm(n*p,0,1),dim=c(n,p))) 
    betas=c(0,c(1,1,1)*signaltonoise,rep(0,p-3)) 
    colnames(x)=c("I",paste("x",1:p,sep="")) 
    true.mod=(betas>0)*1 
    names(true.mod)=colnames(x) 
    y=x%*%betas+rnorm(n,0,1) 
    simdat=data.frame(y,x) 

Usual analysis with Gibbs sampling on model space.
==================================================

Suppose to explore the model space with 10^{4} steps using Gibbs
Sampling, starting from the full and assuming a constant prior on model
space (only to perform Gibb Sampling). BF factors are calculated using
the conventional prior.

    res.GibbsBvs<- GibbsBvs(formula= y ~ ., data=simdat, prior.betas="gZellner", 
                            prior.models="Constant", n.iter=nsim.gibbs, init.model="Full", n.burnin=100, 
                            time.test = FALSE) 

    ## Info. . . .
    ## Most complex model has 202 covariates
    ## From those 1 is fixed and we should select from the remaining 201 
    ## I, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52, x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x92, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192, x193, x194, x195, x196, x197, x198, x199, x200
    ## The problem has a total of 3.213876e+60 competing models
    ## Of these, 10100 are sampled with replacement
    ## Then, 10000 are kept and used to construct the summaries
    ## Working on the problem...please wait.

    summary(res.GibbsBvs) 

    ## 
    ## Call:
    ## GibbsBvs(formula = y ~ ., data = simdat, prior.betas = "gZellner", 
    ##     prior.models = "Constant", n.iter = nsim.gibbs, init.model = "Full", 
    ##     n.burnin = 100, time.test = FALSE)
    ## 
    ## Inclusion Probabilities:
    ##      Incl.prob. HPM MPM
    ## I        0.1246        
    ## x1       0.3587   *    
    ## x2        0.375   *    
    ## x3       0.3702   *    
    ## x4       0.2439        
    ## x5       0.2412        
    ## x6       0.2269        
    ## x7        0.233        
    ## x8        0.247   *    
    ## x9       0.2371        
    ## x10      0.2268        
    ## x11      0.2359        
    ## x12      0.2393        
    ## x13      0.2356        
    ## x14       0.245        
    ## x15      0.2485        
    ## x16      0.2486        
    ## x17       0.228        
    ## x18      0.2386        
    ## x19      0.2331        
    ## x20      0.2454        
    ## x21      0.2307   *    
    ## x22      0.2467        
    ## x23      0.2245        
    ## x24      0.2344        
    ## x25      0.2423        
    ## x26      0.2348        
    ## x27      0.2484   *    
    ## x28      0.2347        
    ## x29      0.2345        
    ## x30      0.2307        
    ## x31      0.2427        
    ## x32      0.2474        
    ## x33      0.2296        
    ## x34      0.2448        
    ## x35      0.2262        
    ## x36      0.2384   *    
    ## x37      0.2228        
    ## x38      0.2379        
    ## x39      0.2407        
    ## x40      0.2385        
    ## x41      0.2423   *    
    ## x42      0.2457        
    ## x43      0.2421        
    ## x44      0.2387        
    ## x45      0.2416        
    ## x46      0.2482        
    ## x47      0.2456        
    ## x48      0.2357        
    ## x49      0.2384        
    ## x50      0.2398        
    ## x51      0.2449        
    ## x52      0.2285        
    ## x53       0.243        
    ## x54       0.246        
    ## x55      0.2347        
    ## x56      0.2266   *    
    ## x57      0.2458        
    ## x58      0.2413        
    ## x59      0.2475   *    
    ## x60      0.2418        
    ## x61      0.2307        
    ## x62      0.2348        
    ## x63        0.24        
    ## x64       0.234        
    ## x65      0.2465        
    ## x66      0.2477        
    ## x67      0.2253        
    ## x68      0.2329        
    ## x69      0.2353        
    ## x70      0.2273        
    ## x71      0.2317        
    ## x72      0.2428        
    ## x73      0.2361        
    ## x74      0.2256        
    ## x75      0.2359        
    ## x76       0.237        
    ## x77      0.2433        
    ## x78      0.2431        
    ## x79      0.2382        
    ## x80      0.2404        
    ## x81       0.229        
    ## x82      0.2292        
    ## x83      0.2413        
    ## x84      0.2312        
    ## x85      0.2351        
    ## x86      0.2424        
    ## x87      0.2233        
    ## x88      0.2338        
    ## x89      0.2319        
    ## x90      0.2456        
    ## x91      0.2385   *    
    ## x92      0.2368        
    ## x93      0.2174        
    ## x94      0.2294        
    ## x95      0.2257        
    ## x96      0.2373        
    ## x97      0.2338        
    ## x98      0.2365        
    ## x99        0.25        
    ## x100     0.2398        
    ## x101     0.2449        
    ## x102     0.2429   *    
    ## x103     0.2473        
    ## x104     0.2394        
    ## x105     0.2259        
    ## x106     0.2407        
    ## x107      0.241        
    ## x108     0.2383        
    ## x109     0.2409   *    
    ## x110     0.2553        
    ## x111     0.2334        
    ## x112     0.2447        
    ## x113     0.2367        
    ## x114     0.2422        
    ## x115      0.235        
    ## x116     0.2493        
    ## x117     0.2332        
    ## x118     0.2439        
    ## x119     0.2405        
    ## x120     0.2335        
    ## x121     0.2454   *    
    ## x122     0.2283        
    ## x123     0.2301   *    
    ## x124     0.2322        
    ## x125     0.2381        
    ## x126     0.2396        
    ## x127     0.2414        
    ## x128     0.2385        
    ## x129     0.2505        
    ## x130     0.2325        
    ## x131     0.2413        
    ## x132     0.2362        
    ## x133     0.2391        
    ## x134     0.2287        
    ## x135     0.2273        
    ## x136     0.2272        
    ## x137     0.2334        
    ## x138     0.2379        
    ## x139     0.2381   *    
    ## x140     0.2288        
    ## x141     0.2252        
    ## x142     0.2375        
    ## x143     0.2467        
    ## x144     0.2366   *    
    ## x145     0.2368        
    ## x146     0.2491        
    ## x147     0.2308        
    ## x148     0.2395        
    ## x149     0.2424   *    
    ## x150     0.2432        
    ## x151     0.2381        
    ## x152     0.2382        
    ## x153      0.252        
    ## x154     0.2405        
    ## x155     0.2422        
    ## x156     0.2415        
    ## x157     0.2498        
    ## x158     0.2394   *    
    ## x159     0.2413        
    ## x160     0.2433   *    
    ## x161     0.2578        
    ## x162     0.2382        
    ## x163     0.2557        
    ## x164     0.2523   *    
    ## x165     0.2543        
    ## x166     0.2345        
    ## x167     0.2398        
    ## x168     0.2386        
    ## x169     0.2396        
    ## x170     0.2406        
    ## x171     0.2419        
    ## x172     0.2357        
    ## x173     0.2485        
    ## x174     0.2399        
    ## x175     0.2287        
    ## x176     0.2475        
    ## x177     0.2448        
    ## x178     0.2344        
    ## x179     0.2451   *    
    ## x180     0.2466        
    ## x181     0.2304        
    ## x182     0.2429        
    ## x183     0.2275        
    ## x184     0.2214        
    ## x185     0.2355        
    ## x186     0.2378        
    ## x187     0.2292        
    ## x188     0.2233        
    ## x189     0.2374        
    ## x190     0.2268        
    ## x191     0.2405        
    ## x192     0.2483   *    
    ## x193     0.2414        
    ## x194     0.2394        
    ## x195     0.2412        
    ## x196     0.2401        
    ## x197     0.2272        
    ## x198     0.2263        
    ## x199      0.229   *    
    ## x200     0.2362   *    
    ## ---
    ## Code: HPM stands for Highest posterior Probability Model and
    ##  MPM for Median Probability Model.
    ##  Results are estimates based on the visited models.

From exploration steps we need at least the set of visited models in all
performed Gibbs steps (BFs are ignored). These are our observations.

    xx=res.GibbsBvs$modelslogBF 
    dim(xx) 

    ## [1] 10000   202

    ww=xx[,ncol(xx)] 
    xx=xx[,-ncol(xx)] 
    image(xx[,1:10]) 

![](s-example1_files/figure-markdown_strict/unnamed-chunk-4-1.svg)

This matrix is the of 0/1s that provides evidence for the true model
into a part of the model space.

A refinement using DPMP
=======================

Not all models into the model space can be visited, only those with
*p* − *n* − 1 covariates (149 in this case). The full model space is
given by the hyper contingency table made by crossing *p* 0/1
categorical variables with 2 levels and thus by 2<sup>*p*</sup> cells
(1.60693810^{60} in this case). On the probability distribution of such
contingency table we assume a Dirichlet process prior (which is the
prior on model space) and the data are the samples from Gibbs step. This
model is proposed here Dunson and Xing (2009).

    mdat=data.frame(xx) 
    for(i in 1:ncol(mdat)) mdat[,i]=factor(mdat[,i]) 

Posterior of cells probabilities (and thus models) is obtained using the
variational algorithm detailed in Ahlmann-Eltze and Yau (2018)
(otherwise Gibbs sampling would be used as originally proposed in Dunson
and Xing (2009)). We assume that the latent space has dimension 2: the
space with covariates with cells with high probabilities (in which it is
supposed to lie the true model) and the set of cells with low
probabilities (in which there is not the true model)

    res <- mixdir(mdat, n_latent=2,max_iter = 1000) 
    cat(res$converged,"\n") 

    ## TRUE

Given one of the latent spaces, each covariate has a probability to be 1
or 0 (i.e. included or not included in the model). Let’s analyse each
latent space separately:

The posterior probability of each covariate for the first latent space
is (first 3 covariates belong to this latent space with high
probability):

    pls1=unlist(lapply(res$category_prob,function(x) x[[1]][2])) 
    barplot(pls1[-1],col=c(rep(2,m1),rep(1,p-m1)),las=2,cex.axis=0.5,ylim=c(0,1)) 

![](s-example1_files/figure-markdown_strict/unnamed-chunk-7-1.svg)

The posterior probability of each covariate for the second latent space
(which represents the null model):

    pls2=unlist(lapply(res$category_prob,function(x) x[[2]][2])) 
    barplot(pls2[-1],col=c(rep(2,m1),rep(1,p-m1)),las=2,cex.axis=0.5,ylim=c(0,1)) 

![](s-example1_files/figure-markdown_strict/unnamed-chunk-8-1.svg)

Comparison among the two analysis
=================================

Let’s compare the posterior inclusion probability with just Gibbs
sampling with the posterior probability using the DPMP:

    exdetail=paste("nGiibsBVS=",nsim.gibbs,
                   " signaltonoise=",signaltonoise,sep="") 
    par(mfrow=c(2,1)) 
    barplot(pls1,col=true.mod+1, 
            main=exdetail, 
            ylab="Post. Prob. DMMP",ylim=c(0,1)) 
    barplot(res.GibbsBvs$inclprob,col=true.mod+1,
            ylab="Post. Inclusion Prob.",ylim=c(0,1)) 

![](s-example1_files/figure-markdown_strict/unnamed-chunk-9-1.svg)

    par(mfrow=c(1,1)) 
    plot(res.GibbsBvs$inclprob,pls1,col=true.mod+1,cex=1,pch=19, 
         xlab="Post. Inclusion Prob.",ylab="Post. Prob. DMMP",main=exdetail,ylim=c(0,1),xlim=c(0,1)) 
    points(res.GibbsBvs$inclprob[true.mod==1],pls1[true.mod==1],col=2,cex=1,pch=23) 
    abline(0,1,lty=2) 

![](s-example1_files/figure-markdown_strict/unnamed-chunk-10-1.svg)

References
==========

Ahlmann-Eltze, Constantin, and Christopher Yau. 2018. “MixDir: Scalable
Bayesian Clustering for High-Dimensional Categorical Data.” In *2018
Ieee 5th International Conference on Data Science and Advanced Analytics
(Dsaa)*, 526–39. IEEE.

Dunson, David B, and Chuanhua Xing. 2009. “Nonparametric Bayes Modeling
of Multivariate Categorical Data.” *Journal of the American Statistical
Association* 104 (487). Taylor & Francis: 1042–51.
