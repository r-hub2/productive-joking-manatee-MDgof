#' Change Marginals of 2D Data
#' 
#' This function creates routines to modify the marginals
#' 
#' @param which method for modifying the data sets.
#' @param side which dimension to modify
#' @param nsample =250 sample size
#' @param null_param parameters under null hypothesis
#' @keywords internal
#' @return a list of functions
#' @export
change.marginals=function(which, side, nsample=250, null_param) {
   if(which=="normal") {
     cdf=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::pnorm(x[,1], 0, p[1])
       if(side%in%c(0,2)) x[,2]=stats::pnorm(x[,2], 0, p[2])    
       x
     }      
     dens=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::dnorm(x[,1], 0, p[1])
       if(side%in%c(0,2)) x[,2]=stats::dnorm(x[,2], 0, p[2])    
       x
     }      
     transf=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::qnorm(x[,1], 0, p[1])
       if(side%in%c(0,2)) x[,2]=stats::qnorm(x[,2], 0, p[2])    
       x
     }   
     Range=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
     if(side==1) 
       Range=matrix(c(-Inf, Inf, 0, 1), 2, 2)
     if(side==2) 
       Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
   }
   if(which=="exponential") {
     cdf=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::pexp(x[,1], p[1])
       if(side%in%c(0,2)) x[,2]=stats::pexp(x[,2], p[2])    
       x
     }      
     dens=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::dexp(x[,1], p[1])
       if(side%in%c(0,2)) x[,2]=stats::dexp(x[,2], p[2])    
       x
     }      
     transf=function(x, p) {
       if(!is.matrix(x)) x=rbind(x)
       if(side%in%c(0,1)) x[,1]=stats::qexp(x[,1], p[1])
       if(side%in%c(0,2)) x[,2]=stats::qexp(x[,2], p[2])    
       x
     }   
     Range=matrix(c(0, Inf, 0, Inf), 2, 2)
     if(side==1) 
       Range=matrix(c(0, Inf, 0, 1), 2, 2)
     if(side==2) 
       Range=matrix(c(0, 1, 0, Inf), 2, 2)
   }
   if(which=="beta22") {
    cdf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::pbeta(x[,1], p[1], p[1])
      if(side%in%c(0,2)) x[,2]=stats::pbeta(x[,2], p[2], p[2])    
      x
    }      
    dens=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::dbeta(x[,1], p[1], p[1])
      if(side%in%c(0,2)) x[,2]=stats::dbeta(x[,2], p[2], p[2])    
      x
    }      
    transf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::qbeta(x[,1], p[1], p[1])
      if(side%in%c(0,2)) x[,2]=stats::qbeta(x[,2], p[2], p[2])    
      x
    }   
    Range=matrix(c(0, 1, 0, 1), 2, 2)
  }
  if(which=="dblexp") {
    cutoff=2
    cdf=function(x, p) {
      K=1/2/(1-exp(-p*cutoff))
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) 
        x[,1]=ifelse(x[,1]<0, exp(p[1]*x[,1])-exp(-p[1]*cutoff), 
                     2-exp(-p[1]*cutoff)-exp(-p[1]*x[,1]))*K[1]
      if(side%in%c(0,2))
        x[,2]=ifelse(x[,2]<0, exp(p[2]*x[,2])-exp(-p[2]*cutoff), 
                     2-exp(-p[2]*cutoff)-exp(-p[1]*x[,2]))*K[2]
      x
    }      
    dens=function(x, p) {
      K=1/2/(1-exp(-p*cutoff))
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=p[1]*exp(-p[1]*abs(x[,1]))*K[1]
      if(side%in%c(0,2)) x[,2]=p[2]*exp(-p[2]*abs(x[,2]))*K[2] 
      x
    }      
    transf=function(x, p) {
      K=1/2/(1-exp(-p*cutoff))
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) 
        x[,1]=ifelse(x[,1]<0.5, log(x[,1]/K[1]+exp(-p[1]*cutoff)), 
                     -log(2-x[,1]/K[1]-exp(-p[1]*cutoff)))/p[1]
      if(side%in%c(0,2)) 
        x[,2]=ifelse(x[,2]<0.5, log(x[,2]/K[2]+exp(-p[2]*cutoff)), 
                     -log(2-x[,2]/K[2]-exp(-p[2]*cutoff)))/p[2]
      x
    }   
    Range=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
    if(side==1) 
      Range=matrix(c(-Inf, Inf, 0, 1), 2, 2)
    if(side==2) 
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
  }
  if(which=="linear") {
    cdf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=p[1]*x[,1]^2+(1-p[1])*x[,1]
      if(side%in%c(0,2)) x[,2]=p[2]*x[,2]^2+(1-p[2])*x[,2]
      x
    }      
    dens=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=2*p[1]*x[,1]+1-p[1]    
      if(side%in%c(0,2)) x[,2]=2*p[2]*x[,2]+1-p[2]    
      x
    }      
    transf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=(p[1]-1+sqrt((1-p[1])^2+4*p[1]*x[,1]))/2/p[1]
      if(side%in%c(0,2)) x[,2]=(p[2]-1+sqrt((1-p[2])^2+4*p[2]*x[,2]))/2/p[2]
      x
    }   
    Range=matrix(c(0, 1, 0, 1), 2, 2)
  }
  if(which=="truncexp") {
    cdf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::pexp(x[,1], p[1])/stats::pexp(1, p[1])
      if(side%in%c(0,2)) x[,2]=stats::pexp(x[,2], p[2])/stats::pexp(1, p[2])    
      x
    }      
    dens=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::dexp(x[,1], p[1])/stats::pexp(1, p[1])
      if(side%in%c(0,2)) x[,2]=stats::dexp(x[,2], p[2])/stats::pexp(1, p[2])    
      x
    }      
    transf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::qexp(x[,1]*stats::pexp(1, p[1]), p[1])
      if(side%in%c(0,2)) x[,2]=stats::qexp(x[,2]*stats::pexp(1, p[2]), p[2])    
      x
    }   
    Range=matrix(c(0, 1, 0, 1), 2, 2)
  }
  if(which=="betaa1") {
    cdf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::pbeta(x[,1], p[1], 1)
      if(side%in%c(0,2)) x[,2]=stats::pbeta(x[,2], p[2], 1)    
      x
    }      
    dens=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::dbeta(x[,1], p[1], 1)
      if(side%in%c(0,2)) x[,2]=stats::dbeta(x[,2], p[2], 1)    
      x
    }      
    transf=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      if(side%in%c(0,1)) x[,1]=stats::qbeta(x[,1], p[1], 1)
      if(side%in%c(0,2)) x[,2]=stats::qbeta(x[,2], p[2], 1)    
      x
    }
    Range=matrix(c(0, 1, 0, 1), 2, 2)
  }
  if(which%in% c("normal", "exponential", "beta22", "dblexp",
                 "linear", "truncexp", "betaa1"))
      return(list(cdf=cdf, dens=dens, transf=transf, Range=Range))
  if(which=="beta-nor-1.beta-nor-mean.marginal") {
      dens1=function(x) stats::dbeta(x, 2, 2)
      dens2=function(x, u) stats::dnorm(x, u)  
      cdf2=function(x, u) stats::pnorm(x, u) 
      rnull=function() {
        x=stats::rbeta(nsample, 2, 2)
        y=stats::rnorm(nsample, x)
        cbind(x, y)
      }
      ralt=function(p=1) {
        if(p<0) p=null_param[which]
        x=stats::rbeta(nsample, p, p)
        y=stats::rnorm(nsample, x)
        cbind(x, y)
      }
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)    
  }
  if(which=="beta-nor-1.beta-nor-stats::sd.marginal") {
      dens1=function(x) stats::dbeta(x, 2, 2)
      dens2=function(x, u) stats::dnorm(x, 0, u)  
      cdf2=function(x, u) stats::pnorm(x, 0, u) 
      rnull=function() {
        x=stats::rbeta(nsample, 2, 2)
        y=stats::rnorm(nsample, 0, x)
        cbind(x, y)
      }
      ralt=function(p=1) {
        if(p<0) p=null_param[which]
        x=stats::rbeta(nsample, p, p)
        y=stats::rnorm(nsample, 0, x)
        cbind(x, y)
      }
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
  }    
  if(which=="beta-beta-2.beta-beta-a.marginal") {
      dens1=function(x) stats::dbeta(x, 2, 2)
      dens2=function(x, u) stats::dbeta(x, u+1, u+1)  
      cdf2=function(x, u) stats::pbeta(x, u+1, u+1) 
      rnull=function() {
        x=stats::rbeta(nsample, 2, 2)
        y=stats::rbeta(nsample, x+1, x+1)
        cbind(x, y)
      }
      ralt=function(p=1) {
        if(p<0) p=null_param[which]
        x=stats::rbeta(nsample, p, p)
        y=stats::rbeta(nsample, x+1, x+1)
        cbind(x, y)
      }
      Range=matrix(c(0, 1, 0, 1), 2, 2)
  }
  if(which=="beta05.normal.marginal") {
    h=function(x) 2*x
    dens1=function(x) {
      p=1-2*stats::pbeta(0.01, 0.5, 0.5)
      ifelse(x<0.01|x>0.99, 0,
             stats::dbeta(x, 1/2, 1/2)/p)
    }
    dens2=function(x, u) stats::dnorm(x, h(u))  
    cdf2=function(x, u) stats::pnorm(x, h(u)) 
    rnull=function() {
      n=round(2*nsample)
      x=stats::rbeta(n, 1/2, 1/2)
      x=x[x>0.01&x<0.99][1:nsample]
      y=stats::rnorm(nsample, h(x))
      cbind(x, y)
    }
    ralt=function(p=1) {
      if(p<0) p=null_param[which]
      n=round(2*nsample)
      x=stats::rbeta(n, 1/2, 1/2)
      x=x[x>0.01&x<0.99][1:nsample]
      y=stats::rnorm(nsample, h(x), p)
      cbind(x, y)
    }
    Range=matrix(c(0.01, 0.99, -Inf, Inf), 2, 2)
  }
  if(which%in%c("beta-nor-1.beta-nor-mean.marginal",
                "beta-nor-1.beta-nor-stats::sd.marginal",
                "beta-beta-2.beta-beta-a.marginal",
                "beta05.normal.marginal")) 
    return(list(cdf2=cdf2, dens1=dens1, dens2=dens2,
              rnull=rnull, ralt= ralt, Range=Range))       
   
}  
