#' Create various case studies for continuous data without parameter estimation
#'
#' This function creates the functions needed to run the various case studies.
#' 
#' @param which name of the case study.
#' @param nsample =250, sample size.
#' @param ReturnCaseNames =FALSE, just return names of case studies?
#' @return a list of functions
#' @export
case.studies.cont=function(which, nsample=250, ReturnCaseNames = FALSE) {
  cases=c("uniform.uniform-diagonal-n", "uniform.uniform-diagonal-b", "normal-ind.normal-cor",
          "t-ind.t-cor", "uniform.Frank", "uniform.Clayton",
          "uniform.Gumbel", "uniform.Galambos", "uniform.HuslerReiss",
          "uniform.Joe", "mix.Clayton.Gumbel", "mix.uniform.Frank",
          "mix.Clayton.Frank", "mix.Frank.Gumbel", "mix.Frank.Joe",
          "normal.normal-shift-one.marginal", "normal.normal-shift-both.marginal",
          "normal.normal-stretch-one.marginal", "normal.normal-stretch-both.marginal",
          "normal.stretch-shift-one.marginal", "normal.stretch-shift-two.marginal",
          "uniform.rotated-uniform.marginal", "uniform.beta-one.marginal",
          "uniform.beta-two.marginal", "uni-exp-1.uni-exp-l.marginal",
          "exp-exp-1.exp-exp-l.marginal","beta-nor-1.beta-nor-mean.marginal",
          "beta-nor-1.beta-nor-sd.marginal",
          "beta-beta-2.beta-beta-a.marginal",
          "beta05.normal.marginal"
          )
  if (ReturnCaseNames) return(cases)
  if (is.numeric(which))  which = cases[which]
  Range01=matrix(c(0, 1, 0, 1), 2, 2)
  Range0Inf=matrix(c(0, Inf, 0, Inf), 2, 2)
  RangeInfInf=matrix(c(-Inf, Inf, -Inf, Inf), 2, 2)
  null_param=MDgof::power_studies_cont_results[[3]][,1]
  if(startsWith(which, "uniform")) {
     pnull=function(x) {
          if(!is.matrix(x)) x=rbind(x)
          d=ncol(x)
          y=rep(0, nrow(x))
          for(i in seq_along(y)) {
             y[i]=stats::punif(x[i,1])*stats::punif(x[i,2])
             if(d>2) 
             for(j in 3:d) y[i]=y[i]*stats::punif(x[i,j])
          }  
          y
      }
      dnull=function(x) {
          if(!is.matrix(x)) x=rbind(x)
          rep(1, nrow(x))
      }
      rnull=function() round(matrix(stats::runif(2*nsample), ncol=2), 5)     
  }
  if(startsWith(which, "normal")) {
    pnull=function(x) {
      f=function(x) mvtnorm::pmvnorm(rep(-Inf, length(x)), x)
      if(!is.matrix(x)) return(f(x))
      apply(x, 1, f)
    }
    dnull=function(x) {
      f=function(x) mvtnorm::dmvnorm(x)
      if(!is.matrix(x)) return(f(x))
      apply(x, 1, f)
    }
    rnull=function() round(mvtnorm::rmvnorm(nsample, c(0, 0)), 5)     
  } 
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
# 1:  
  if(which=="uniform.uniform-diagonal-n") {
    return(list(
      pnull=pnull, dnull=dnull, rnull=rnull,    
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        x=round(matrix(stats::runif(2*nsample), ncol=2), 5)
        add_stripe(round(alpha*nsample), x, c(1, 0.2, 0))
      },
      param_alt=c(0, 0.3),
      Range=Range01
    ))
  }  
# 2:   
  if(which=="uniform.uniform-diagonal-b") {
    return(list(
      pnull = pnull, rnull=rnull, dnull=dnull,
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        x=round(matrix(stats::runif(2*nsample), ncol=2), 5)         
        add_stripe(round(alpha*nsample), x, c(5, 0.3, 0))
      },
      param_alt=c(0, 0.3),
      Range=Range01
    ))
  }
# 3: 
  if(which=="normal-ind.normal-cor") {
    return(list(
      pnull=pnull, dnull=dnull, 
      rnull = rnull, 
      ralt=function(s) {
         if(s<0) s=null_param[which]
         round(mvtnorm::rmvnorm(nsample, sigma=matrix(c(1, s, s, 1), 2, 2)), 5)
      },
      param_alt=c(0, 0.35),
      Range=RangeInfInf
    ))
  }
# 4:  
  if(which=="t-ind.t-cor") {
    return(list(
      pnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        y=rep(0, nrow(x))
        for(i in seq_along(y))
          y[i]=mvtnorm::pmvt(rep(-Inf, ncol(x)), x[i, ],df=5)
        y
      },
      dnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        k=nrow(x)
        y=rep(0, k)
        for(i in 1:k)
          y[i]=mvtnorm::dmvt(x[i,], df=5, log=FALSE) 
        y
      },
      rnull=function() 
        round(mvtnorm::rmvt(nsample, sigma=diag(2), df=5), 5),
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvt(nsample, sigma=matrix(c(1, s, s, 1), 2, 2), 
                            df=5), 5)
      },
      param_alt=c(0, 0.35),
      Range=RangeInfInf
    ))
  }  
# Case Studies where distribution under the null hypothesis is uniform 
# 5-10: 
  if(startsWith(which, "uniform") & !(endsWith(which, "marginal"))) {
    family=strsplit(which , "\\.")[[1]][2]
    ralt=function(s) {
        if(s==0) return(rnull())
        cop=gen.cop(family, s)
        round(copula::rCopula(nsample, cop), 5)
    }
    if(family=="Frank") p=c(0, 2)
    if(family=="Joe") p=c(0, 1.55)
    if(family=="Clayton") p=c(0, 0.6)
    if(family=="Gumbel") p=c(1, 1.25)
    if(family=="Galambos") p=c(0, 0.6)
    if(family=="HuslerReiss") p=c(0, 0.95)
    return(list(
        pnull=pnull, dnull=dnull, rnull = rnull, 
        ralt=ralt, param_alt=p, Range=Range01
      ))
  }
# Case Studies of mixtures 
# 11:15:  
  if(startsWith(which, "mix")) {
    if(which=="mix.Clayton.Gumbel") { 
      cop1=gen.cop("Clayton", p=8)
      cop2=gen.cop("Gumbel", p=1.175)
      p=c(0.5, 0.25)
    }
    if(which=="mix.uniform.Frank") { 
      cop1=gen.cop("indep", p=2)
      cop2=gen.cop("Frank", p=20)
      p=c(0.5, 0.25)
    }
    if(which=="mix.Clayton.Frank") { 
      cop1=gen.cop("Clayton", p=8)
      cop2=gen.cop("Frank", p=8)
      p=c(0.5, 0.04)
    }
    if(which=="mix.Frank.Gumbel") { 
      cop1=gen.cop("Frank", p=20)
      cop2=gen.cop("Gumbel", p=2)
      p=c(0.5, 0.15)
    }
    if(which=="mix.Frank.Joe") { 
      cop1=gen.cop("Frank", p=20)
      cop2=gen.cop("Joe", p=2)
      p=c(0.5, 0.15)
    }
    return(list(
        pnull=function(x) {
          wts=c(1/2, 1/2)
          mxc=copula::mixCopula(list(cop1, cop2), w=wts)
          copula::pCopula(x, mxc)
        },
        dnull=function(x) {
          wts=c(1/2, 1/2)
          mxc=copula::mixCopula(list(cop1, cop2), w=wts)
          copula::dCopula(x, mxc)
        },       
        rnull=function() {
          wts=c(1/2, 1/2)
          mxc=copula::mixCopula(list(cop1, cop2), w=wts)
          round(copula::rCopula(nsample, mxc), 5)
        },
        ralt=function(a) {
          if(a<0) a=null_param[which]
          wts=c(a, 1-a)
          mxc=copula::mixCopula(list(cop1, cop2), w=wts)
          round(copula::rCopula(nsample, mxc), 5)
        },
        param_alt=p,
        Range=Range01
      ))
  }
# Studies with unequal marginals
# 16:  
  if(which=="normal.normal-shift-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, mean = c(0, s)), 5)
      },  
      param_alt=c(0, 0.35),
      Range=RangeInfInf
    ))
  }  
# 17:  
  if(which=="normal.normal-shift-both.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, mean = c(s, s)), 5)
      },
      param_alt=c(0, 0.25),
      Range=RangeInfInf
    ))
  }  
# 18:  
  if(which=="normal.normal-stretch-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, sigma = matrix(c(s, 0, 0, 1), 2, 2)), 5)
      },
      param_alt=c(1, 1.75),
      Range=RangeInfInf
    ))
  }  
# 19:  
  if(which=="normal.normal-stretch-both.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, sigma = matrix(c(s, 0, 0, s), 2, 2)), 5)
      },
      param_alt=c(1, 1.5),
      Range=RangeInfInf
    ))
  }  
# 20:  
  if(which=="normal.stretch-shift-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, 
                    mean=c(s, 0), sigma = matrix(c(s+1, 0, 0, 1), 2, 2)), 5)
      },
      param_alt=c(0, 0.3),
      Range=RangeInfInf
    ))
  }  
# 21:  
  if(which=="normal.stretch-shift-two.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, 
                    mean=c(s, s),sigma = matrix(c(1+s, 0, 0, 1+s), 2, 2)), 5)
      },
      param_alt=c(0, 0.15),
      Range=RangeInfInf
    ))
  } 
# 22:  
  if(which=="uniform.rotated-uniform.marginal") {
      return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        matrix(stats::runif(2 * nsample), ncol = 2) %*% 
          matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), 2, 2)
      },
      param_alt=c(0, 0.125),
      Range=Range01
    ))
  }
# 23:  
  if(which=="uniform.beta-one.marginal") {
    return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(a) {
        if(a<0) a=null_param[which]
        matrix(stats::rbeta(2 * nsample, 1, a), ncol = 2)
      },
      param_alt=c(1, 1.25),
      Range=Range01
    ))
  }
# 24:  
  if(which=="uniform.beta-two.marginal") {
    return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(a) {
        if(a<0) a=null_param[which]
        matrix(stats::rbeta(2 * nsample, a, a), ncol = 2)
      },
      param_alt=c(1, 1.35),
      Range=Range01
    ))
  }
# 25:  
  if(which=="uni-exp-1.uni-exp-l.marginal") {
    return(list(
      pnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        x[x[,2]==0, 2]=1e-6
        x[,1]+exp(-(x[,1]+1)*x[,2])/x[,2]-exp(-x[,2])/x[,2]
      },  
      dnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        (x[,1]+1)*exp(-(x[,1]+1)*x[,2])
      },  
      rnull=function() {
        x=stats::runif(nsample)
        y=stats::rexp(nsample, x+1)
        cbind(x, y)
      },
      ralt=function(l) {
        if(l<0) l=null_param[which]
        x=stats::runif(nsample)
        y=stats::rexp(nsample, l*(x+1))
        cbind(x, y)
      },
      param_alt=c(1, 1.35),
      Range=matrix(c(0, 1, 0, Inf),2,2)
    ))
  }
# 26:  
  if(which=="exp-exp-1.exp-exp-l.marginal") {
    return(list(
      pnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        x[is.infinite(x[,2]), 2]=1e6
        1-exp(-x[,1])+exp(-(x[,2]+1)*x[,1]-x[,2])/(x[,2]+1)-exp(-x[,2])/(x[,2]+1)
      },  
      dnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        (x[,1]+1)*exp(-(x[,1]+1)*x[,2]-x[,1])
      },  
      rnull=function() {
        x=stats::rexp(nsample, 1)
        y=stats::rexp(nsample, x+1)
        cbind(x, y)
      },
      ralt=function(l) {
        if(l<0) l=null_param[which]
        x=stats::rexp(nsample, 1)
        y=stats::rexp(nsample, l*(x+1))
        cbind(x, y)
      },
      param_alt=c(1, 1.35),
      Range=Range0Inf
    ))
  }
# 27-30:  
  if(which%in% c("beta-nor-1.beta-nor-mean.marginal",
                 "beta-nor-1.beta-nor-sd.marginal",
                 "beta-beta-2.beta-beta-a.marginal",
                 "beta05.normal.marginal")) {
    tmp=change.marginals(which, 0, nsample, null_param)
    if(which=="beta-nor-1.beta-nor-mean.marginal") {
      pr=c(2, 3.25)
      Ra=matrix(c(0,1,-Inf, Inf),2,2)
    }
    if(which=="beta-nor-1.beta-nor-sd.marginal") {
      pr=c(2, 3.25)
      Ra=matrix(c(0,1,-Inf, Inf),2,2)
    }
    if(which=="beta-beta-2.beta-beta-a.marginal") {
      pr=c(2, 3.25)
      Ra=Range01
    }
    if(which=="beta05.normal.marginal") {
      pr=c(1, 1.4)
      Ra=matrix(c(0,1,-Inf, Inf),2,2)
    }  
    return(list(
      pnull=function(x) {
       f=function(x) {
        g=function(u) tmp$dens1(u)*tmp$cdf2(x[2], u)  
        stats::integrate(g, tmp$Range[1,1], x[1])$value
       }
       if(!is.matrix(x)) x=rbind(x)
       out=rep(0, nrow(x))
       for(i in 1:nrow(x)) 
        out[i]=f(x[i, ])
       out
      },
      dnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        tmp$dens1(x[,1])*tmp$dens2(x[,2], x[,1])
      },
      rnull=tmp$rnull,
      ralt=tmp$ralt,
      param_alt=pr,
      Range=matrix(c(0.01, 0.99, -Inf, Inf), 2, 2)
    ))
  }
  
}
