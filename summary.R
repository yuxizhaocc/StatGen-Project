# generate output

library(SixSigma)
library(tidyverse)
library(abind)

#--------- summary of 800 replicates for report ----------#

# load the output from the two jobs

load("StatGen-Sim-study-complete-part1.RData")

cause.array.est.1=cause.array.est
cause.array.zero.1=cause.array.zero
conmix.array.est.1=conmix.array.est
conmix.array.pval.1=conmix.array.pval

load("StatGen-Sim-study-complete-part2.RData")

cause.array.est.2=cause.array.est
cause.array.zero.2=cause.array.zero
conmix.array.est.2=conmix.array.est
conmix.array.pval.2=conmix.array.pval

# combine the two output
cause.array.est.combine=abind(cause.array.est.1,cause.array.est.2,along=5)
cause.array.zero.combine=abind(cause.array.zero.1,cause.array.zero.2,along=5)
conmix.array.est.combine=abind(conmix.array.est.1,conmix.array.est.2,along=5)
conmix.array.pval.combine=abind(conmix.array.pval.1,conmix.array.pval.2,along=5)

# generate summary statistics
summ.conmix.mse = array(NA,dim=c(n.N,n.theta,n.mix,n.var))
summ.conmix.error1 = array(NA,dim=c(n.N,n.theta,n.mix,n.var))
summ.conmix.power = array(NA,dim=c(n.N,n.theta,n.mix,n.var))

summ.cause.mse = array(NA,dim=c(n.N,n.theta,n.mix,n.var))
summ.cause.error1 = array(NA,dim=c(n.N,n.theta,n.mix,n.var))
summ.cause.power = array(NA,dim=c(n.N,n.theta,n.mix,n.var))

for (h in 1:n.N){
  for (i in 1:n.theta){
    for (j in 1:n.mix){
      for (k in 1:n.var){
        
        theta = theta.lst[[i]]
        summ.conmix.mse[h,i,j,k] = mean((conmix.array.est.combine[h,i,j,k,]-theta)^2)
        summ.cause.mse[h,i,j,k] = mean((cause.array.est.combine[h,i,j,k,]-theta)^2)
        
        if (theta!=0){
          summ.conmix.power[h,i,j,k] = mean(conmix.array.pval.combine[h,i,j,k,]<0.05)
          summ.cause.power[h,i,j,k] = mean(1-cause.array.zero.combine[h,i,j,k,])
        }
        
        if (theta==0){
          summ.conmix.error1[h,i,j,k] = 1-mean(conmix.array.pval.combine[h,i,j,k,]>=0.05)
          summ.cause.error1[h,i,j,k] = 1-mean(cause.array.zero.combine[h,i,j,k,])
        }
      }
    }
  }
}

output <- expand.grid(N=gl(n.N,1,labels=c(N)),
                      theta=gl(n.theta,1,labels=c(theta.lst)),
                      mix=gl(n.mix,1,labels=c(1:n.mix)),
                      vars=gl(n.var,1,labels=c(1:n.var))
)

output2 <- data.frame(output,conmix.mse=round(c(summ.conmix.mse),4),
                      cause.mse=round(c(summ.cause.mse),4),
                      conmix.error1=c(summ.conmix.error1),
                      cause.error1=c(summ.cause.error1),
                      conmix.power=c(summ.conmix.power),
                      cause.power=c(summ.cause.power)
)
output2%>%arrange(N,theta,mix,vars)

#        N theta mix vars conmix.mse cause.mse conmix.error1 cause.error1 conmix.power cause.power
# 1  50000   0.2   1    1     0.0017    0.0054            NA           NA            1     0.54875
# 2  50000   0.2   1    2     0.0018    0.0056            NA           NA            1     0.74750
# 3  50000   0.2   2    1     0.0017    0.0069            NA           NA            1     0.48375
# 4  50000   0.2   2    2     0.0019    0.0093            NA           NA            1     0.81375
# 5  50000   0.2   3    1     0.0017    0.0073            NA           NA            1     0.42625
# 6  50000   0.2   3    2     0.0021    0.0153            NA           NA            1     0.85750
# 7  50000     0   1    1     0.0001    0.0044       0.06625      0.02125           NA          NA
# 8  50000     0   1    2     0.0001    0.0069       0.07500      0.03750           NA          NA
# 9  50000     0   2    1     0.0001    0.0049       0.07375      0.01625           NA          NA
# 10 50000     0   2    2     0.0002    0.0116       0.08875      0.06625           NA          NA
# 11 50000     0   3    1     0.0001    0.0060       0.08000      0.01500           NA          NA
# 12 50000     0   3    2     0.0002    0.0187       0.11500      0.17125           NA          NA
# 13 50000  -0.2   1    1     0.0016    0.0055            NA           NA            1     0.53375
# 14 50000  -0.2   1    2     0.0015    0.0101            NA           NA            1     0.30375
# 15 50000  -0.2   2    1     0.0017    0.0065            NA           NA            1     0.45125
# 16 50000  -0.2   2    2     0.0014    0.0173            NA           NA            1     0.10625
# 17 50000  -0.2   3    1     0.0017    0.0068            NA           NA            1     0.43000
# 18 50000  -0.2   3    2     0.0014    0.0253            NA           NA            1     0.05250
# 19 80000   0.2   1    1     0.0015    0.0021            NA           NA            1     0.93500
# 20 80000   0.2   1    2     0.0017    0.0030            NA           NA            1     0.98875
# 21 80000   0.2   2    1     0.0016    0.0028            NA           NA            1     0.85375
# 22 80000   0.2   2    2     0.0019    0.0068            NA           NA            1     0.99000
# 23 80000   0.2   3    1     0.0016    0.0032            NA           NA            1     0.80375
# 24 80000   0.2   3    2     0.0021    0.0139            NA           NA            1     0.99250
# 25 80000     0   1    1     0.0001    0.0019       0.06500      0.01625           NA          NA
# 26 80000     0   1    2     0.0001    0.0038       0.08000      0.05500           NA          NA
# 27 80000     0   2    1     0.0001    0.0023       0.07500      0.01500           NA          NA
# 28 80000     0   2    2     0.0002    0.0079       0.10250      0.13375           NA          NA
# 29 80000     0   3    1     0.0001    0.0030       0.08875      0.01250           NA          NA
# 30 80000     0   3    2     0.0002    0.0159       0.13000      0.32875           NA          NA
# 31 80000  -0.2   1    1     0.0015    0.0022            NA           NA            1     0.93125
# 32 80000  -0.2   1    2     0.0013    0.0056            NA           NA            1     0.70625
# 33 80000  -0.2   2    1     0.0016    0.0029            NA           NA            1     0.85500
# 34 80000  -0.2   2    2     0.0013    0.0113            NA           NA            1     0.32125
# 35 80000  -0.2   3    1     0.0016    0.0033            NA           NA            1     0.79625
# 36 80000  -0.2   3    2     0.0012    0.0196            NA           NA            1     0.12250