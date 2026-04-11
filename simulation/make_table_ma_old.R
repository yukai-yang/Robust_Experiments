library(kableExtra)
# VAR case
load("output/EXP02MAold.Rdata")

# true parameters
ma = matrix( c(0.7, 0.3, 0, 0.7), 2, 2)
ma = c(ma)

atables = array(0, dim=c(3, 2, 3, 5, 4))
for(cter in 1:3){
  for(tter in 1:2){
    for(ater in 1:3){
      for(oter in 1:5){
        for(kter in 1:4){
          tmp = NULL;
          for(est in 1:4){
            tmp = c(tmp,ret[,cter,tter,ater,oter,kter,est] - ma[est])
          }
          atables[cter, tter, ater, oter, kter] = sqrt(mean(tmp**2))
        }
      }
    }
  }
}


atables[1, 1, 1:3, 5, 1:4]

# contamination, no, AO, IO
cont = 1
# 0.01, 0.05, 0.1
alpha = 3
# 0, 1, 2, 4
kappa = 1
# 500, 1000
iT = 1
# 0, 5, 10, 50, 100
outlier = 5;

tables = NULL; iT = 1
for(cter in 1:2){
  subtab = NULL
  for(oter in 2:5){
    subtab = rbind(subtab,atables[cter,iT,1:3, oter,1:4])
  }
  tables = cbind(tables, subtab)
}
tab1 = tables

tables = NULL; iT = 2
for(cter in 1:2){
  subtab = NULL
  for(oter in 2:5){
    subtab = rbind(subtab,atables[cter,iT,1:3, oter,1:4])
  }
  tables = cbind(tables, subtab)
}
tab2 = tables
tables = data.frame(tab1, tab2)
tables = cbind(b=c(5,"","",10,"","",50,"","",100,"",""), a=c(1,5,10), tables)

kbl(tables, digits=4, booktabs=T, linesep="") %>%
  add_header_above(c(" ", "%", rep(c("k=0", "1", "2","4", "0", "1","2","4"),2))) %>%
  add_header_above(c(" "=2, "no"=4, "AO"=4, "no"=4, "AO"=4)) %>%
  add_header_above(c("z", "a", "T=500"=8, "T=1000"=8)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  kable_styling(font_size = 12)
