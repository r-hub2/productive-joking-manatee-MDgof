## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(error=TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(MDgof)
B=200 
st=function() knitr::knit_exit()
examples=MDgof::hybrid.mdgof.vignette

## -----------------------------------------------------------------------------
ReRunExamples=FALSE

## -----------------------------------------------------------------------------
set.seed(123)

## ----ex1----------------------------------------------------------------------
rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
x=rnull()

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex1a"]]=hybrid_test(x, rnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex1a"]]

## ----ex2----------------------------------------------------------------------
rnull=function(p) mvtnorm::rmvnorm(200, 
      mean=p[1:2],
      sigma=matrix(c(p[3], p[5], p[5], p[4]), 2, 2))
phat=function(x) 
   c(apply(x, 2, mean), c(apply(x, 2, var), cov(x[,1],x[,2])))
x=rnull(c(1, 2, 4, 9,0.6))
phat(x)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex2a"]]=hybrid_test(x, rnull, phat, B=B)

## -----------------------------------------------------------------------------
examples[["ex2a"]]

## -----------------------------------------------------------------------------
y=mvtnorm::rmvt(300, sigma=diag(2), df=5)

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex2b"]]=hybrid_test(y, rnull, phat, B=B)

## -----------------------------------------------------------------------------
examples[["ex2b"]]

## -----------------------------------------------------------------------------
TSextra=list(which="statistic", nbins=cbind(c(3,3), c(3,4)))

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex3a"]]=hybrid_test(x, rnull, phat,
#          TS=MD2sample::chiTS.cont, TSextra=TSextra,
#          B=B)

## -----------------------------------------------------------------------------
examples[["ex3a"]]

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex3b"]]=hybrid_test(y, rnull, phat,
#          TS=MD2sample::chiTS.cont, TSextra=TSextra,
#          B=B)

## -----------------------------------------------------------------------------
examples[["ex3b"]]

## ----ex4----------------------------------------------------------------------
rnull=function() mvtnorm::rmvnorm(100, c(0, 0))
ralt=function(s) mvtnorm::rmvnorm(100, c(0, 0),
                sigma=matrix(c(1, s, s, 1), 2, 2))

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex4a"]]=hybrid_power(rnull, ralt,
#       c(0, 0.8), B=B, maxProcessor=1)

## -----------------------------------------------------------------------------
examples[["ex4a"]]

## -----------------------------------------------------------------------------
TSextra=list(which="pvalue", nbins=c(3,4))

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex4b"]]=hybrid_power(rnull, ralt, c(0, 0.8), TS=MD2sample::chiTS.cont,
#           TSextra=TSextra, B=500,
#           With.p.value=TRUE)

## -----------------------------------------------------------------------------
examples[["ex4b"]]

## ----ex5----------------------------------------------------------------------
rnull=function() {
   x=rbinom(1000, 5, 0.5)
   y=rbinom(1000, 3, 0.5+x/20)
   MDgof::sq2rec(table(x, y))
}
x=rnull()

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex5"]]=hybrid_test(x, rnull, B=B)

## -----------------------------------------------------------------------------
examples[["ex5"]]

## -----------------------------------------------------------------------------
TSextra=list(which="statistic")

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex6c"]]=hybrid_test(x, rnull,  TS=MD2sample::chiTS.disc, TSextra=TSextra,  B=B)

## -----------------------------------------------------------------------------
examples[["ex6c"]]

## ----ex7----------------------------------------------------------------------
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
# examples[["ex7b"]]= hybrid_power(rnull, ralt, c(0.5, 0.475), B=200)

## -----------------------------------------------------------------------------
examples[["ex7b"]]

## -----------------------------------------------------------------------------
TSextra=list(which="pvalue")

## ----eval=ReRunExamples-------------------------------------------------------
# examples[["ex7d"]]=hybrid_power(rnull, ralt, c(0.5, 0.475), TS=MD2sample::chiTS.disc, TSextra=TSextra,
#          With.p.value = TRUE, B=B)

## -----------------------------------------------------------------------------
examples[["ex7d"]]

## -----------------------------------------------------------------------------
saveRDS(examples, "hybrid.mdgof.vignette.rds")

