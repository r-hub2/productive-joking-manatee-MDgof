#' Create various case studies with parameter estimation
#'
#' This function creates the functions needed to run the various case studies that include parameter estimation.
#' 
#' @param which name of the case study.
#' @param nsample =250, sample size.
#' @param ReturnCaseNames =FALSE, just return names of case studies?
#' @return a list of functions
#' @export
case.studies.est=function(which, nsample=250, ReturnCaseNames = FALSE) {
  cases=c("normal-sigma.t", "Frank.Frank", "Clayton.Clayton", 
          "Gumbel.Gumbel", "Galambos.Galambos", 
          "HuslerReiss.HuslerReiss", "Joe.Joe",
          "mix.uniform.Frank", "mix.Clayton.Frank",
          "mix.Clayton.Clayton", "mix.uniform.Joe",
          "mix.uniform.Plackett", "Clayton.marginal.normal-sigma",
          "Gumbel.marginal.normal-sigma", "Joe.marginal.exponential",
          "normal-mean.t.marginal", "Clayton.marginal.normal.marginal",
          "Joe.marginal.exponential.marginal",
          "Plackett.marginal.beta22.marginal",
          "Frank.marginal.dblexp.marginal",
          "Frank.marginal.linear.marginal",
          "Joe.marginal.truncexp.marginal",
          "Clayton.marginal.betaa1.marginal",
          "beta22.normal-mean.marginal",
          "beta22.normal-stats::sd.marginal",
          "beta22.lognormal-mean.marginal",
          "beta22.lognormal-stats::sd.marginal",
          "beta22.exponential.marginal",
          "exp1.normal-stats::sd.marginal",
          "beta05.gamma.marginal"
          )
  if (ReturnCaseNames) return(cases)
  if (is.numeric(which))  which = cases[which]
  Range=matrix(c(0, 1, 0, 1), 2, 2)

# mle for mixing distributions  
  mlemix=function(x, f, g) {
    fx=f(x)
    gx=g(x)
    psi=fx-gx
    aold=0.5
    k=0
    repeat {
      k=k+1
      tmp=psi/(aold*fx+(1-aold)*gx)
      anew=aold+sum(tmp)/sum(tmp^2)
      if(anew<0) return(0.1)
      if(anew>1) return(0.9)
      if(abs(aold-anew)<1e-3) break
      aold=anew
      if(k>100) break
    }
    anew
  }
# add a diagonal stripe to data  
  add_stripe=function(n, x, a) {
    if(n==0) return(x)
    m=nrow(x)
    Range = range(x[,1])
    y1=stats::runif(n, Range[1], Range[2])
    u=(y1-min(y1))/(max(y1)-min(y1))
    y2=y1+(1-a[1]*(u-0.5)^2)*stats::runif(n, -a[2], a[2])
    y2[y2<Range[1]]=Range[1]-y2[y2<Range[1]]
    y2[y2>Range[2]]=2*Range[2]-y2[y2>Range[2]]
    if(a[3]==1) 
      y2= Range[1]+diff(Range)*(1-(y2-Range[1])/diff(Range))
    return(rbind(x[1:(m-n), ], cbind(y1, y2))) 
  }
  null_param=MDgof::power_studies_cont_est_results[[3]][,1]
# case studies with equal marginals  
  if(which=="normal-sigma.t") {
    return(list(
      pnull=function(x, s=0) {
        f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), x, sigma=matrix(c(1, s, s, 1), 2, 2))
        if(!is.matrix(x)) return(f(x))
        apply(x, 1, f)      
      },
      dnull=function(x, s=0) {
        f=function(x) mvtnorm::dmvnorm(x, sigma=matrix(c(1, s, s, 1), 2, 2))
        if(!is.matrix(x)) return(f(x))
        apply(x, 1, f)      
      },
      rnull=function(p=0) {
        round(mvtnorm::rmvnorm(nsample, sigma=matrix(c(1, p, p, 1), 2, 2)), 5)     
      }, 
      ralt=function(df) {
         if(df<0) {
           p=null_param[which]
           return(round(mvtnorm::rmvnorm(nsample, sigma=matrix(c(1, p, p, 1), 2, 2)), 5))
         }
         round(mvtnorm::rmvt(nsample, sigma=diag(2), df=df), 5)
      },
      param_alt=c(0.5, 8),
      phat=function(x) stats::cor(x)[1,2],
      Range=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
    ))
  }

# Same Copula with different parameters
  a=strsplit(which, "\\.")[[1]]
  if(length(a==2) & a[1]==a[2]) {
    family=strsplit(which, "\\.")[[1]][1]
    if(family=="Frank") p=c(2, 0.75)
    if(family=="Clayton") p=c(2.5, 0.6)
    if(family=="Gumbel") p=c(1.5, 0.5)
    if(family=="Galambos") p=c(1.5, 0.8)
    if(family=="HuslerReiss") p=c(1, 0.4)
    if(family=="Joe") p=c(2, 0.5)    
    return(list(
        pnull=function(x, p) {
          x[x==0]=1e-10
          x[x==1]=1-1e-10
          cop0=gen.cop(family, p)
          copula::pCopula(x, cop0) 
        },  
        dnull=function(x, p) {
          x[x==0]=1e-10
          x[x==1]=1-1e-10
          cop0=gen.cop(family, p)
          copula::dCopula(x, cop0)
        }, 
        rnull=function(p) {
          cop0=gen.cop(family, p)
          round(copula::rCopula(nsample, cop0), 5)
        }, 
        ralt=function(s) {
          if(s<0) {
            p=null_param[which]
            cop0=gen.cop(family, p)
            return(round(copula::rCopula(nsample, cop0), 5))   
          }
          cop0=gen.cop(family, p[1])
          x=copula::rCopula(nsample, cop0)
          round(add_stripe(round(s*nsample), x, c(1, 0.1, 0)), 5)
        },        
        phat=function(x) {
           cop0=gen.cop(family, p[1])
           copula::fitCopula(cop0, copula::pobs(x))@estimate
        },
        param_alt=p, Range=Range
      ))
  }
# Various mixtures, with mixing ratio estimated
  if(startsWith(which, "mix")) {
    family1=strsplit(which, "\\.")[[1]][2]
    family2=strsplit(which, "\\.")[[1]][3]
    if(family1=="uniform" & family2=="Frank") 
        pr=c(0.5, NA, 2, 15)
    if(family1=="Clayton" & family2=="Frank") 
        pr=c(0.5, 2.5, 2, 15)
    if(family1=="Clayton" & family2=="Clayton") 
      pr=c(0.75, 2, 5, 35)
    if(family1=="uniform" & family2=="Joe") 
      pr=c(0.5, NA, 5, 25)
    if(family1=="uniform" & family2=="Plackett") 
      pr=c(0.5, NA, 0.9, 0.1)
    
    return(list(
      pnull=function(x, p) {
        cop1=gen.cop(family1, pr[2])
        cop2=gen.cop(family2, pr[3])
        p[1]*copula::pCopula(x, cop1)+
        (1-p[1])*copula::pCopula(x, cop2)
      },  
      dnull=function(x, p) {
        cop1=gen.cop(family1, pr[2])
        cop2=gen.cop(family2, pr[3])
        p[1]*copula::dCopula(x, cop1)+
          (1-p[1])*copula::dCopula(x, cop2)
      }, 
      rnull=function(p) {
        cop1=gen.cop(family1, pr[2])
        cop2=gen.cop(family2, pr[3])
        m=round(p[1]*nsample)
        rbind(copula::rCopula(m, cop1),
              copula::rCopula(nsample-m, cop2))
      },      
      ralt=function(p) {
        if(p<0) {
          p=null_param[which]
          cop1=gen.cop(family1, pr[2])
          cop2=gen.cop(family2, pr[3])
          m=round(p[1]*nsample)
          return(rbind(copula::rCopula(m, cop1),
                copula::rCopula(nsample-m, cop2)))
        }
        cop1=gen.cop(family1, pr[2])
        cop2=gen.cop(family2, p)
        m=round(0.5*nsample)
        rbind(copula::rCopula(m, cop1),
              copula::rCopula(nsample-m, cop2))
      },        
      phat=function(x) {
        cop1=gen.cop(family1, pr[2])
        cop2=gen.cop(family2, pr[3])
        f=function(x) copula::dCopula(x, cop1)
        g=function(x) copula::dCopula(x, cop2)
        mlemix(x, f, g)
      },
      param_alt=pr[c(1,4)], Range=Range
    ))
  }
# Copula with non uniform marginals
  a=strsplit(which, "\\.")[[1]]
  if(length(a)==3 & a[2]=="marginal") {
    family=a[1]
    if(a[1]=="Clayton" & a[3]=="normal-sigma") {
       pr=c(2, 0.25)
       cdf=function(x, p) cbind(stats::pnorm(x[,1], 0, p[1]), x[,2])
       dens=function(x, p) stats::dnorm(x[,1], 0, p[1])
       transf=function(x, p) cbind(stats::qnorm(x[,1], 0, p[1]), x[,2])
       transf_alt=function(x, p) cbind(stats::qnorm(x[,1], p[1]), x[,2])             
       phat=function(x) stats::sd(x[, 1])
       Range=matrix(c(-Inf, Inf, 0, 1), 2, 2)
    }
    if(a[1]=="Gumbel" & a[3]=="normal-sigma") {
      pr=c(2, 0.25)
      cdf=function(x, p) cbind(stats::pnorm(x[,1], 0, p[1]), stats::pnorm(x[,2], 0, p[2]))
      dens=function(x, p) stats::dnorm(x[,1], 0, p[1])*stats::dnorm(x[,2], 0, p[2])
      transf=function(x, p) cbind(stats::qnorm(x[,1], 0, p[1]), stats::qnorm(x[,2], 0, p[2]))
      transf_alt=function(x, p) cbind(stats::qnorm(x[,1], p), stats::qnorm(x[,2], p))      
      phat=function(x) apply(x, 2, stats::sd)
      Range=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
    }
    if(a[1]=="Joe" & a[3]=="exponential") {
      pr=c(4, 1.35)
      cdf=function(x, p) cbind(stats::pexp(x[,1], p[1]), stats::pexp(x[,2], p[2]))
      dens=function(x, p) stats::dexp(x[,1], p[1])*stats::dexp(x[,2], p[2])
      transf=function(x, p) cbind(stats::qexp(x[,1], p[1]), stats::qexp(x[,2], p[2]))
      transf_alt=function(x, p) cbind(stats::qgamma(x[,1], p, 1), stats::qgamma(x[,2], p, 1))      
      phat=function(x) 1/apply(x, 2, mean)
      Range=matrix(c(0, Inf, 0, Inf), 2, 2)
    }
    return(list(
      pnull=function(x, p) {
        if(!is.matrix(x)) x=rbind(x)
        cop0=gen.cop(family, pr[1])
        copula::pCopula(cdf(x, p), cop0) 
      },  
      dnull=function(x, p) {
        if(!is.matrix(x)) x=rbind(x)
        cop0=gen.cop(family, pr[1])
        copula::dCopula(cdf(x, p), cop0)*dens(x, p)
      }, 
      rnull=function(p) {
        cop0=gen.cop(family, pr[1])
        x=copula::rCopula(nsample, cop0)
        x=transf(x, p)
        round(x, 5)
      }, 
      ralt=function(p) {
        if(p<0) {
          cop0=gen.cop(family, pr[1])
          x=copula::rCopula(nsample, cop0)
          x=transf(x, c(1,1))
          return(round(x, 5))
        }
        cop0=gen.cop(family, pr[1])
        x=copula::rCopula(nsample, cop0)
        x=transf_alt(x, p)
        round(x, 5)
      },        
      phat=phat,
      param_alt=pr, Range=Range
    ))
  }  
  
# case studies with unequal marginals  
  if(which=="normal-mean.t.marginal") {
    return(list(
       pnull=function(x, p=0) {
          f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), x, mean=c(p,0))
          if(!is.matrix(x)) return(f(x))
          apply(x, 1, f)      
      },
      dnull=function(x, p=0) {
          f=function(x) mvtnorm::dmvnorm(x, mean=c(p,0))
          if(!is.matrix(x)) return(f(x))
          apply(x, 1, f)      
      },
      rnull=function(p=0) round(mvtnorm::rmvnorm(nsample, c(p,0)), 5),     
      ralt=function(df) {
        if(df<0) {
           p=null_param[which]
           return(round(mvtnorm::rmvnorm(nsample, c(p, 0)), 5))     
        }
        round(mvtnorm::rmvt(nsample, sigma=diag(2), df=df), 5)
      },  
      param_alt=c(1, 8),
      phat=function(x) mean(x[,1]),
      Range=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
    ))
  }

# Copula with non uniform marginals
  a=strsplit(which, "\\.")[[1]]
  if(length(a)==4 & a[2]=="marginal"& a[4]=="marginal") {
  family=a[1]
  if(a[1]=="Clayton" & a[3]=="normal") 
     pr=c(2, 4, 1, 1, 1.2, 1.2, 0)
  if(a[1]=="Joe" & a[3]=="exponential") 
     pr=c(2, 4, 1, 1, 1.2, 1.2, 0)
  if(a[1]=="Plackett" & a[3]=="beta22") 
     pr=c(0.9, 0.5, 2, 2, 2.75, 2.75, 0)
  if(a[1]=="Frank" & a[3]=="dblexp") 
     pr=c(2, 4, 1, 1, 1.35, 1.35, 0)
  if(a[1]=="Frank" & a[3]=="linear") 
     pr=c(2, 10, 0.5, 0.5, 0.9, 0.9, 1)
  if(a[1]=="Joe" & a[3]=="truncexp") 
     pr=c(2, 5, 1, 1, 2, 2, 1)
  if(a[1]=="Clayton" & a[3]=="betaa1") 
     pr=c(2, 5, 2, 2, 2.5, 2.5, 2)  
  names(pr)=c("Cupola null", "Cupola alt",
              "Transform x null", "Transform y null",
              "Transform x alt", "Transform y alt",
              "side")
  tmp=change.marginals(a[3], pr[7], nsample, null_param)
  return(list(
    pnull=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      cop0=gen.cop(family, p)
      copula::pCopula(tmp$cdf(x, pr[3:4]), cop0) 
    },  
    dnull=function(x, p) {
      if(!is.matrix(x)) x=rbind(x)
      cop0=gen.cop(family, p)
      copula::dCopula(tmp$cdf(x, pr[3:4]), cop0)*tmp$dens(x, pr[3:4]) 
    }, 
    rnull=function(p) {
      cop0=gen.cop(family, p)
      x=copula::rCopula(nsample, cop0)
      x=tmp$transf(x, pr[3:4])
      round(x, 5)
    }, 
    ralt=function(p) {
      if(p<0) {
        p=null_param[which]
        cop0=gen.cop(family, p)
        x=copula::rCopula(nsample, cop0)
        x=tmp$transf(x, pr[3:4])
        return(round(x, 5))   
      }
      cop0=gen.cop(family, pr[1])
      x=copula::rCopula(nsample, cop0)
      x=tmp$transf(x, c(p,p))
      round(x, 5)
    },        
    phat=function(x) {
      cop0=gen.cop(family, pr[1])
      res=try(copula::fitCopula(cop0, copula::pobs(x))@estimate)
      if(inherits(res, "try-error")) res=pr[1]
      res  
    },
    param_alt=pr[c(1, 5)], Range=tmp$Range
   ))
  }
# Conditional distributions with non uniform marginals
  a=strsplit(which, "\\.")[[1]]
  if(a[[1]]%in%c("beta05", "beta22", "exp1")) {
    dens1=function(x) stats::dbeta(x, 2, 2)
    if(a[[1]]=="beta22" & a[[2]]=="normal-mean") {
      pr=c(2, 0, 3.5)
      dens2=function(x, u, p=pr[1]) stats::dnorm(x, p, u)  
      cdf2=function(x, u, p=pr[1]) stats::pnorm(x, p, u) 
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) stats::rnorm(nsample, p, x)
      phat=function(x) sum(x[,2]/x[,1]^2)/sum(1/x[,1]^2)
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
    }
    if(a[[1]]=="beta22" & a[[2]]=="normal-stats::sd") {
      pr=c(2, 1, 3)
      dens2=function(x, u, p=pr[1]) stats::dnorm(x, u, p)  
      cdf2=function(x, u, p=pr[1]) stats::pnorm(x, u, p) 
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) stats::rnorm(nsample, x, p)
      phat=function(x) sqrt(mean( (x[,1]-x[,2])^2 ))
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
    }
    if(a[[1]]=="beta22" & a[[2]]=="lognormal-mean") {
      pr=c(2, 0, 3)
      h=function(x) (x+1)/10
      dens2=function(x, u, p=pr[1]) stats::dnorm(log(x), p, h(u))/x  
      cdf2=function(x, u, p=pr[1]) stats::pnorm(log(x), p, h(u)) 
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) exp(stats::rnorm(nsample, p, h(x)))
      phat=function(x) sum(log(x[,2])/h(x[,1])^2)/sum(1/h(x[,1])^2)
      Range=matrix(c(0, 1, 0, Inf), 2, 2)
    }
    if(a[[1]]=="beta22" & a[[2]]=="lognormal-stats::sd") {
      pr=c(2, 1, 3)
      h=function(x) (x+1)/100
      dens2=function(x, u, p=pr[1]) stats::dnorm(log(x), h(u), p)/x  
      cdf2=function(x, u, p=pr[1]) stats::pnorm(log(x), h(u), p) 
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) exp(stats::rnorm(nsample, h(x), p))
      phat=function(x) sqrt(mean( (h(x[,1])-log(x[,2]))^2 ))
      Range=matrix(c(0, 1, 0, Inf), 2, 2)
    } 
    if(a[[1]]=="beta22" & a[[2]]=="exponential") {
      pr=c(2, 1, 3.5)
      h=function(x) x+1
      dens2=function(x, u, p=pr[1]) stats::dexp(x, h(u)*p)  
      cdf2=function(x, u, p=pr[1]) stats::pexp(x, h(u)*p)
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) stats::rexp(nsample, h(x)*pr[1])
      phat=function(x) 1/mean(h(x[,1])*x[,2])
      Range=matrix(c(0, 1, 0, Inf), 2, 2)
    }  
    if(a[[1]]=="exp1" & a[[2]]=="normal-stats::sd") {
      pr=c(1, 1, 3)
      h=function(x) x+1
      dens1=function(x) stats::dexp(x, pr[1])/stats::pexp(1, pr[1])
      dens2=function(x, u, p=pr[1]) stats::dnorm(x, h(u), p)  
      cdf2=function(x, u, p=pr[1]) stats::pnorm(x, h(u), p) 
      r1=function(a=pr[1]) {
         x=NULL
         repeat {
           x=c(x, stats::rexp(nsample, a))
           x=x[x<1]
           if(length(x)>=nsample) break
         }   
         x[1:nsample]
      }   
      r2=function(x, p=pr[1]) stats::rnorm(nsample, h(x), p)
      phat=function(x) sqrt(mean( ( h(x[,1])-x[,2])^2 ))
      Range=matrix(c(0, 1, -Inf, Inf), 2, 2)
    }     
    if(a[[1]]=="beta05" & a[[2]]=="gamma") {
      pr=c(1/2, 1, 0.3)
      h=function(x) x+1
      dens1=function(x) 
         ifelse(x<1e-5|x>1-1e-5, stats::dbeta(1e-5, pr[1], pr[1]), 
                                 stats::dbeta(x, pr[1], pr[1]))
      dens2=function(x, u, p=pr[1]) stats::dgamma(x, h(u), p)  
      cdf2=function(x, u, p=pr[1]) stats::pgamma(x, h(u), p) 
      r1=function(a=pr[1]) stats::rbeta(nsample, a, a)
      r2=function(x, p=pr[1]) stats::rgamma(nsample, h(x), p)
      phat=function(x) sum(h(x[,1]))/sum(x[,2])
      Range=matrix(c(0, 1, 0, Inf), 2, 2)
    }
    return(list( 
      pnull=function(x, p=pr[2]) {
        f=function(x, p) {
          g=function(u) dens1(u)*cdf2(x[2], u, p)  
          stats::integrate(g, 0, x[1])$value
        }
        if(!is.matrix(x)) x=rbind(x)
        out=rep(0, nrow(x))
        for(i in 1:nrow(x)) 
          out[i]=f(x[i, ], p)
        out
      },
      dnull=function(x, p=pr[2]) {
        if(!is.matrix(x)) x=rbind(x)
        dens1(x[,1])*dens2(x[,2], x[,1], p)
      },
      rnull=function(p=pr[2]) {
        x=r1(pr[1])
        y=r2(x, p)
        cbind(x, y)
      },
      ralt=function(l) {
        if(l<0) {
          p=null_param[which]
          x=r1(pr[1])
          y=r2(x, p)
          return(cbind(x, y))  
        }
        x=r1(l)
        y=r2(x, pr[2])
        cbind(x, y)
      },
      phat=phat,
      param_alt=pr[c(1,3)], Range=Range
    ))
  }  
  
}
