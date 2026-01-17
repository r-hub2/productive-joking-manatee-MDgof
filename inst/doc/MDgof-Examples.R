## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(error=TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(MDgof)
B=200 
st=function() knitr::knit_exit()
examples=MDgof::examples.mdgof.vignette

## -----------------------------------------------------------------------------
ReRunExamples=FALSE

## -----------------------------------------------------------------------------
set.seed(123)

## ----ex1----------------------------------------------------------------------
# cumulative distribution function under the null hypothesis
pnull=function(x) { 
  f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), x)
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
# function that generates data under the null hypothesis
rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
x=rnull()

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex1a"]]=gof_test(x, pnull, rnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex1a"]]

## -----------------------------------------------------------------------------
dnull=function(x) {
  f=function(x) mvtnorm::dmvnorm(x)
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex1b"]]=gof_test(x, pnull, rnull, dnull=dnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex1b"]]

## ----ex2----------------------------------------------------------------------
pnull=function(x, p) {
  f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), 
      x, mean=p[1:2],  
      sigma=matrix(c(p[3], p[5], p[5], p[4]), 2, 2))
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
rnull=function(p) mvtnorm::rmvnorm(200, 
      mean=p[1:2],
      sigma=matrix(c(p[3], p[5], p[5], p[4]), 2, 2))
dnull=function(x, p) {
  f=function(x) mvtnorm::dmvnorm(x, 
      mean=p[1:2],
      sigma=matrix(c(p[3], p[5], p[5], p[4]), 2, 2))
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
phat=function(x) 
   c(apply(x, 2, mean), c(apply(x, 2, var), cov(x[,1],x[,2])))
x=rnull(c(1, 2, 4, 9,0.6))
phat(x)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex2a"]]=gof_test(x, pnull, rnull, phat, dnull=dnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex2a"]]

## -----------------------------------------------------------------------------
y=mvtnorm::rmvt(200, sigma=diag(2), df=5)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex2b"]]=gof_test(y, pnull, rnull,
#     phat, dnull=dnull,  B=B)

## -----------------------------------------------------------------------------
examples[["ex2b"]]

## ----ex3----------------------------------------------------------------------
TSextra=list(Continuous=TRUE, 
             WithEstimation=TRUE, 
             Withpvalue=TRUE)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex3a"]]=gof_test(x, pnull, rnull, phat,
#          TS=newTS, TSextra=TSextra,
#          B=B)

## -----------------------------------------------------------------------------
examples[["ex3a"]]

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex3b"]]=gof_test(y, pnull, rnull, phat,
#          TS=newTS, TSextra=TSextra,
#          B=B)

## -----------------------------------------------------------------------------
examples[["ex3b"]]

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex3c"]]=gof_test_adjusted_pvalue(x, pnull,
#                 rnull, dnull=dnull, phat=phat,
#                 B=c(B,B))

## -----------------------------------------------------------------------------
a=examples[["ex3c"]]
for(i in seq_along(a)) 
    print(paste(names(a)[i],": ", round(a[i],4)))


## ----ex4----------------------------------------------------------------------
pnull=function(x) {
  f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), x)
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
dnull=function(x) {
  f=function(x) mvtnorm::dmvnorm(x)
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
ralt=function(s) mvtnorm::rmvnorm(100, c(0, 0),
                sigma=matrix(c(1, s, s, 1), 2, 2))

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex4a"]]=gof_power(pnull, rnull, ralt, c(0, 0.8),
#           dnull=dnull, B=B, maxProcessor=1)

## -----------------------------------------------------------------------------
examples[["ex4a"]]

## -----------------------------------------------------------------------------
TSextra=list(Continuous=TRUE, 
             WithEstimation=FALSE, 
             Withpvalue=TRUE)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex4b"]]=gof_power(pnull, rnull, ralt, c(0, 0.8), TS=newTS,
#           TSextra=TSextra, B=500,
#           With.p.value=TRUE)

## -----------------------------------------------------------------------------
examples[["ex4b"]]

## ----ex5----------------------------------------------------------------------
pnull=function(x) {
  f=function(x) 
    sum(dbinom(0:x[1], 5, 0.5)*pbinom(x[2], 3, 0.5+0:x[1]/20))
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
rnull=function() {
   x=rbinom(1000, 5, 0.5)
   y=rbinom(1000, 3, 0.5+x/20)
   MDgof::sq2rec(table(x, y))
}
x=rnull()
x

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex5"]]=gof_test(x, pnull, rnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex5"]]

## ----ex6----------------------------------------------------------------------
rnull_cont=function(p) {
  x=rbeta(250, p, p)
  y=rbeta(250, x, x)
  cbind(x=x, y=y)
}  
dta_cont=rnull_cont(2)
ggplot2::ggplot(data=as.data.frame(dta_cont), ggplot2::aes(x,y)) + 
  ggplot2::geom_point()

## -----------------------------------------------------------------------------
phat_cont=function(x) {
  n=nrow(x)
  sx=sum( log(x[,1]*(1-x[,1])) )/2/n
  num=function(a) digamma(2*a)-digamma(a)+sx
  denom=function(a) 2*trigamma(2*a)-trigamma(a)
  anew=2
  repeat {
    aold=anew
    anew=aold-num(aold)/denom(aold)
    if(abs(aold-anew)<0.001) break
  }
  anew
}
phat_cont(dta_cont)

## -----------------------------------------------------------------------------
pnull=function(x, p) {
   f=function(x) {
     g=function(u) dbeta(u, p, p)*pbeta(x[2], u, u)  
     integrate(g, 0, x[1])$value
   }
   if(!is.matrix(x)) return(f(x))
   apply(x, 1, f)  
}

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex6a"]]=gof_test(dta_cont, pnull, rnull_cont,
#          phat=phat_cont,
#          Ranges=matrix(c(0, 1,0, 1), 2, 2),
#          maxProcessor = 1, B=B)

## -----------------------------------------------------------------------------
examples[["ex6a"]]

## -----------------------------------------------------------------------------
rnull_disc=function(p) {
  pnull=function(x, p) {
     f=function(x) {
       g=function(u) dbeta(u, p, p)*pbeta(x[2], u, u)  
       integrate(g, 0, x[1])$value
     }
     if(!is.matrix(x)) return(f(x))
     apply(x, 1, f)  
  }
  nb=c(50, 30)
  xvals=1:nb[1]/nb[1]
  yvals=1:nb[2]/nb[2]
  z=cbind(rep(xvals, nb[2]), rep(yvals, each=nb[1]), 0)
  prob=c(t(MDgof::p2dC(z, pnull, p)))
  z[, 3]=rbinom(prod(nb), 1e6, prob)
  z
}
dta_disc=rnull_disc(2)
dta_disc[1:10, ]

## -----------------------------------------------------------------------------
phat_disc=function(x) {
  nb=c(50, 30)
  n=sum(x[,3]) #sample size
  valsx=unique(x[,1])-1/nb[1]/2#x values, center of bins
  x=tapply(x[,3], x[,1], sum) # counts for x alone
  sx=sum(x*log(valsx*(1-valsx)))/2/n
  num=function(a) digamma(2*a)-digamma(a)+sx
  denom=function(a) 2*trigamma(2*a)-trigamma(a)
  anew=2
  repeat {
    aold=anew
    anew=aold-num(aold)/denom(aold)
    if(abs(aold-anew)<0.00001) break
  }
  anew
}
p=phat_disc(dta_disc)
p

## ----disc---------------------------------------------------------------------
set.seed(111)
x=rbeta(1e6, 2, 2)
y=rbeta(1e6, x, x)
dta_cont=cbind(x=x, y=y)
dta_disc=MDgof::discretize(dta_cont, 
                  Range=matrix(c(0,1,0,1),2,2),
                  nbins=c(50, 30))
head(dta_disc)

## ----A, eval=ReRunExamples----------------------------------------------------
# examples[["ex6b"]]=gof_test(dta_disc, pnull, rnull_disc, phat=phat_disc, B=B)

## -----------------------------------------------------------------------------
examples[["ex6b"]]

## -----------------------------------------------------------------------------
TSextra=list(Continuous=FALSE, 
             WithEstimation=TRUE, 
             Withpvalue=FALSE)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex6c"]]=gof_test(dta_disc, pnull, rnull_disc, phat_disc, TS=newTS, TSextra=TSextra,  B=B)

## -----------------------------------------------------------------------------
examples[["ex6c"]]

## ----ex7----------------------------------------------------------------------
pnull=function(x) {
  f=function(x) 
    sum(dbinom(0:x[1], 5, 0.5)*pbinom(x[2], 3, 0.5+0:x[1]/20))
  if(!is.matrix(x)) return(f(x))
  apply(x, 1, f)
}
rnull=function() {
   x=rbinom(1000, 5, 0.5)
   y=rbinom(1000, 3, 0.5+x/20)
   MDgof::sq2rec(table(x, y))
}
ralt=function(p) {
   x=rbinom(1000, 5, p)
   y=rbinom(1000, 3, 0.5+x/20)
   MDgof::sq2rec(table(x, y))
}
x=ralt(0.475)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex7a"]]=gof_test(x, pnull, rnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex7a"]]

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex7b"]]= gof_power(pnull, rnull, ralt, c(0.5, 0.475),B=200)

## -----------------------------------------------------------------------------
examples[["ex7b"]]

## -----------------------------------------------------------------------------
TSextra=list(Continuous=FALSE, 
             WithEstimation=FALSE, 
             Withpvalue=TRUE)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex7c"]]=gof_test(x, pnull, rnull,
#          TS=newTS, TSextra=TSextra,
#          B=B)

## -----------------------------------------------------------------------------
examples[["ex7c"]]

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex7d"]]=gof_power(pnull, rnull, ralt, c(0.5, 0.475),
#          TS=newTS, TSextra=TSextra,
#          With.p.value = TRUE, B=B)

## -----------------------------------------------------------------------------
examples[["ex7d"]]

## ----ex8----------------------------------------------------------------------
TSextra=list(Continuous=TRUE, 
             WithEstimation=TRUE, 
             Withpvalue=TRUE)

## ----eval=ReRunExamples-------------------------------------------------------
# a=run.studies(study=1:2, # do first two
#             Continuous=TRUE,
#             WithEstimation = TRUE,
#             TS=newTS,
#             TSextra=TSextra,
#             With.p.value=TRUE,
#             B=B)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex8"]]= run.studies(1:2,
#             Continuous=TRUE, WithEstimation = FALSE,
#             param_alt=cbind(c(0, 0.4), c(0, 0.4)),
#             alpha=0.1, nsample=500, B=B)

## -----------------------------------------------------------------------------
examples[["ex8"]][["Null"]][1:2, ]
examples[["ex8"]][["Pow"]][1:2, ]

