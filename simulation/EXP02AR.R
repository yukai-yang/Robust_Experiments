# bivariate AR model experiment
library(parallel)
source("utilfuncs.R")

# number of repetitions
iM = 2000
# settings for parallel computation
ncores = 5

# bivariate
n = 2

# VAR model
ar = matrix( c(0.7, 0.3, 0, 0.7), 2, 2)
ma = matrix( 0, 2, 2)
model = list(ar=ar,ma=ma)
mv = chol(matrix(c(1, .2, .2, 1), 2, 2))

# alpha
iTT = c(500, 1000)
alpha = c(0.01, 0.05, 0.1)
outliers = c(0, 5, 10, 50, 100)
kappa = c(0, 1, 2, 4)

# simulation
set.seed(1)
seeds = rnorm(iM)*100000000

rep_func <- function(im){
  set.seed(seeds[im])
  
  results = array(0, dim=c(3, length(iTT), length(alpha), length(outliers), length(kappa), n*n))
  for(tter in seq_along(iTT)){
    iT = iTT[tter]
    
    for(ater in seq_along(alpha)){
      aa = alpha[ater]
      
      for(oter in seq_along(outliers)){
        outl = outliers[oter]
        
        # zeta, outliers
        zeta = matrix(0, iT, n)
        tmp = aa * iT
        tmpp = seq(1,iT,ceiling(iT/tmp))
        zeta[tmpp,] = outl
        subsamp = rep(1, iT)
        subsamp[tmpp] = 0
        
        innov = matrix(rnorm(iT * n), iT, n) %*% mv
        dat = sim_varma(model, iT, n, innov, zeta)
        
        for(kter in seq_along(kappa)){
          est0 = est_var(dat$mx0, model=model, subsamp, kappa=kappa[kter])
          est1 = est_var(dat$mx1, model=model, subsamp, kappa=kappa[kter])
          est2 = est_var(dat$mx2, model=model, subsamp, kappa=kappa[kter])
          
          results[1,tter,ater,oter,kter,] = c(est0$ar)
          results[2,tter,ater,oter,kter,] = c(est1$ar)
          results[3,tter,ater,oter,kter,] = c(est2$ar)
        }
      }
    }
  }
  
  print(im)
  return(results)
}

# iM
# contamination (0, 1, 2)
# iT, alpha, outlier, kappa
# results = array(0, dim=c(3, length(iTT), length(alpha), length(outliers), length(kappa), n*n))

sink("output/log_VAR.txt")

cat("Experiment starts!\n")

ptm <- proc.time()

lres = mclapply(1:iM, rep_func, mc.cores = ncores)

ret = array(0, dim=c(iM, 3, length(iTT), length(alpha), length(outliers), length(kappa), n*n))
for(mter in 1:iM){
  ret[mter,,,,,,] = lres[[mter]]
}

proc.time() - ptm

cat("Done!\n")

save(ret, file="output/EXP02AR.Rdata")

sink()

# load("output/EXP02AR.Rdata")
# ret[,1,1,1,1,1,1]
