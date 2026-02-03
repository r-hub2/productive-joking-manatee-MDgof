## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error=TRUE,
                      message = FALSE, warning = FALSE)
library(MDgof)

## ----message=TRUE-------------------------------------------------------------
MDgof::draw_case("uniform.uniform-diagonal-n", 
                  Continuous=TRUE,  
                  WithEstimation=FALSE, 
                  Dim=2)

## ----message=TRUE-------------------------------------------------------------
MDgof::draw_case("uniform.uniform-diagonal-n", 
                  Continuous=TRUE,  
                  WithEstimation=FALSE, 
                  Dim=2, AltOnly=TRUE,
                 palt=0.6)

## ----message=TRUE-------------------------------------------------------------
MDgof::draw_case(3, Dim=5, AltOnly=TRUE, Dms=c(2,4))

## -----------------------------------------------------------------------------
tmp=MDgof::case.studies(1, F, F)
x=tmp$rnull()
x

## -----------------------------------------------------------------------------
plot(x, xlim=c(0,1), ylim=c(0, 1), xlab="x", ylab="y", type="n")
for(i in 0:6) {
  segments(i/5,0,i/5,1)
  segments(0/5,i/5,1,i/5)
}
for(i in 0:4) {
  for(j in 0:4)
    text(j/5+1/10, i/5+1/10, x[5*i+j+1,3])
}

