library(MendelianRandomization,lib.loc = .libPaths()[2])
library(cause,lib.loc = .libPaths()[2])
library(foreach,lib.loc = .libPaths()[2])
library(doParallel,lib.loc = .libPaths()[2])
numCores <- detectCores()
registerDoParallel(numCores)

# load functions

source("/gpfs/home/yz17m/statistical genetics/simulation-functions.R")

# simulate 10000 independent SNPs
M = 40000

# define the sample size for GWAS
N = list(50000, 80000)

# define number of replicates
N.rep = 400

# mixtures
mix.lst = list(
          list(p1=0.014,p2=0.006),
          list(p1=0.01,p2=0.01),
          list(p1=0.006,p2=0.014)
)

# theta
theta.lst = as.list(c(0.2,0,-0.2))

# variances
var.lst = list(
   list(tilde.sigma2.x=5e-5,sigma2.u=0,mu.y=0,tilde.sigma2.y=5e-5),
   list(tilde.sigma2.x=4.1e-5,sigma2.u=1e-4,mu.y=0,tilde.sigma2.y=4.1e-5)
)

# define the loop
n.N = length(N)
n.mix =length(mix.lst)
n.var = length(var.lst)
n.theta =length(theta.lst)

conmix.array.est = array(NA,dim=c(n.N,n.theta,n.mix,n.var,N.rep))
conmix.array.pval = array(NA,dim=c(n.N,n.theta,n.mix,n.var,N.rep))

cause.array.est = array(NA,dim=c(n.N,n.theta,n.mix,n.var,N.rep))
cause.array.zero = array(NA,dim=c(n.N,n.theta,n.mix,n.var,N.rep))

options(warn=-1)

for (h in 1:n.N){
  for (i in 1:n.theta){
   for (j in 1:n.mix){
    for (k in 1:n.var){
      
      seed = i*10000+j*1000+k*100
      p1 = mix.lst[[j]]$p1
      p2 = mix.lst[[j]]$p2
      
      tilde.sigma2.x = var.lst[[k]]$tilde.sigma2.x
      sigma2.u = var.lst[[k]]$sigma2.u
      mu.y = var.lst[[k]]$mu.y
      tilde.sigma2.y = var.lst[[k]]$tilde.sigma2.y
      
      theta = theta.lst[[i]]
      
      Nsize = N[[h]]
      
      results = foreach(l=1:N.rep,.combine = rbind, 
                        .packages=c('MendelianRandomization','cause'))%dopar%{
        seed = seed+l*10
        sim.mixG = sim.MR.mixG(M=M,p1=p1,p2=p2,seed = seed)
        
        sim.params = sim.MR.params(sim.mixG=sim.mixG,
                                   tilde.sigma2.x=tilde.sigma2.x,
                                   sigma2.u=sigma2.u,
                                   mu.y=mu.y,tilde.sigma2.y=tilde.sigma2.y,seed=seed+1)
        
        sim.beta = sim.MR.largeN(N=Nsize,sim.params=sim.params,theta=theta,seed=seed+2)
        
        conmix.rslt = mr_conmix(mr_input(bx = sim.beta$beta.jx, bxse = sim.beta$se.jx,
                                         by = sim.beta$beta.jy, byse = sim.beta$se.jy), psi = 3, CIMin = -1, CIMax = 5, CIStep = 0.01)
        
        conmix.est= conmix.rslt$Estimate
        conmix.pval = conmix.rslt$Pvalue
        # conmix.valid = length(conmix.rslt$Valid)
        
        cause.X = data.frame(snp=sim.beta$G.id, beta_hat_1=sim.beta$beta.jx, seb1=sim.beta$se.jx, beta_hat_2=sim.beta$beta.jy, seb2=sim.beta$se.jy)
        set.seed(seed+3)
        if (M <= 1e5){ind.replace=TRUE} else {ind.replace=FALSE}
        cause.varlist <- with(cause.X, sample(snp, size=100000, replace=ind.replace))
        cause.params <- est_cause_params(cause.X, cause.varlist)
        cause.rslt <- cause(X=cause.X, param_ests = cause.params,pval_thresh = 0.05)
        cause.rslt.summ <- summary(cause.rslt)
        
        cause.est = cause.rslt.summ$quants[[2]][1,1]
        cause.ci = c(cause.rslt.summ$quants[[2]][2,1], cause.rslt.summ$quants[[2]][3,1])
        cause.zero = (cause.ci[1]<=0 & cause.ci[2]>=0)
     
        
        return(c(conmix.est,conmix.pval,cause.est,cause.zero))
        
      }
      
      conmix.array.est[h,i,j,k,]= results[,1]
      conmix.array.pval[h,i,j,k,]= results[,2]
      cause.array.est[h,i,j,k,]= results[,3]
      cause.array.zero[h,i,j,k,]= results[,4]
    }
   }
  }
}


save(list = ls(all=TRUE), file = "StatGen-Sim-study-complete-part1.RData")
