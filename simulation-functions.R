library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores)

# ---------- functions -----------#

# simulate the mixtures of IVs
sim.MR.mixG = function(M,p1,p2,p3=0.01,p4=0.97,seed=1126){
  set.seed(seed)
  
  my_prob = c(p1,p2,p3,p4)
  Z_label <- rmultinom(n=1,size=M,prob=my_prob)
  
  G = data.frame(G.id=1:M,G.label=sample(unlist(lapply(1:4,function(x)rep(x,Z_label[x,1])))))
  return(G)
}


# simulate the parameter alpha, gamma, phi based on the mixture

sim.MR.params = function(sim.mixG,sigma2.x=5e-5,tilde.sigma2.x,sigma2.u,mu.y,sigma2.y=5e-5,tilde.sigma2.y,seed=1127){
  set.seed(seed)
  M=nrow(sim.mixG)
  Z_label = table(sim.mixG$G.label)
  
  gamma = rep(0,M);phi = rep(0,M); alpha=rep(0,M)
  gamma[sim.mixG$G.label==1] = rnorm(n=Z_label[1],mean=0,sd=sqrt(sigma2.x))
  gamma[sim.mixG$G.label==2] = rnorm(n=Z_label[2], mean=0,sd=sqrt(tilde.sigma2.x))
  if (sigma2.u>0){
    phi[sim.mixG$G.label==2] = rnorm(n=Z_label[2], mean=0,sd=sqrt(sigma2.u)) 
  }
  alpha[sim.mixG$G.label==2] = rnorm(n=Z_label[2],mean=mu.y,sd=sqrt(tilde.sigma2.y))
  alpha[sim.mixG$G.label==3] = rnorm(n=Z_label[3],mean=mu.y,sd=sqrt(sigma2.y))
  
  return(list(sim.mixG=sim.mixG,gamma=gamma,phi=phi,alpha=alpha))
}

# simulate the summary statistics when N large
sim.MR.largeN = function(N,sim.params,theta.ux=0.3, theta.uy=0.3,theta,seed=1128){
  set.seed(seed)
  n.x = N
  n.y = N/2
  sim.mixG  = sim.params$sim.mixG
  M = nrow(sim.mixG)
  
  gamma=sim.params$gamma
  phi=sim.params$phi
  alpha = sim.params$alpha
  
  beta.jx = gamma + theta.ux*phi + rnorm(n=M,mean=0,sd=sqrt(1/n.x))
  beta.jy = alpha + theta.uy*phi + theta*beta.jx + rnorm(n=M,mean=0,sd=sqrt(1/n.y))
  
  se.jx = rep(sqrt(1/n.x),M) 
  se.jy = rep(sqrt(1/n.y),M) 
  
  return(data.frame(G.id=sim.mixG$G.id, G.label = sim.mixG$G.label, beta.jx=beta.jx,beta.jy=beta.jy,
                    se.jx = se.jx, se.jy = se.jy))
}




