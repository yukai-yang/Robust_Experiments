library(kableExtra)
source("utilfuncs.R")

# VAR case
load("output/EXP02MA.Rdata")

# true parameters
ma = matrix( c(0.7, 0.3, 0, 0.7), 2, 2)
ma = c(ma)

atables = array(0, dim=c(2, 2, 3, 5, 4))
btables = array(0, dim=c(2, 2, 3, 5, 4))
for(cter in 1:2) for(tter in 1:2) for(ater in 1:3) for(oter in 1:5) for(kter in 1:4){
  tmp = NULL; for(est in 1:4){
    tmp = cbind(tmp,ret[,cter,tter,ater,oter,kter,est] - ma[est])
  }
  # total bias
  atables[cter, tter, ater, oter, kter] = sqrt(sum(colMeans(tmp)**2))
  # total RMSE
  btables[cter, tter, ater, oter, kter] = sqrt(sum(colMeans(tmp**2)))
}


atables[1, 1, 1:3, 5, 1:4]
btables[1, 1, 1:3, 5, 1:4]

# total bias
#tmptab = atables
# RMSE
tmptab = btables

# total bias, change atables
tab1 = make_table(tmptab, iT=1, cont=1:2, outlier=2:5, alpha=1:3,kappa=1:4)
tab2 = make_table(tmptab, iT=2, cont=1:2, outlier=2:5, alpha=1:3,kappa=1:4)
cbind(tab1,tab2)

tables = data.frame(tab1, tab2)
tables = cbind(b=c(5,"","",10,"","",50,"","",100,"",""), a=c(1,5,10), tables)

kbl(tables, digits=4, booktabs=T, linesep="") %>%
  add_header_above(c(" ","%", rep(c("k=0","2","5","9"),4))) %>%
  add_header_above(c(" "=2, "no"=4, "AO"=4, "no"=4, "AO"=4)) %>%
  add_header_above(c("z", "a", "T=500"=8, "T=1000"=8)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  kable_styling(font_size = 10)

kbl(tables, format="latex", digits=4, linesep="") %>%
  add_header_above(c(" ","%", rep(c("k=0","2","5","9"),4))) %>%
  add_header_above(c(" "=2, "no"=4, "AO"=4, "no"=4, "AO"=4)) %>%
  add_header_above(c("z", "a", "T=500"=8, "T=1000"=8)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  kable_styling(font_size = 10)
