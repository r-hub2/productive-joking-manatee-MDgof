#' Create various case studies for continuous data in 5 dimensions without parameter estimation
#'
#' This function creates the functions needed to run the various case studies.
#' 
#' @param which name of the case study.
#' @param nsample =250, sample size.
#' @param ReturnCaseNames =FALSE, just return names of case studies?
#' @return a list of functions
#' @export
case.studies.cont.D5=function(which, nsample=250, ReturnCaseNames = FALSE) {
  cases=c("uniform.diagonal.n", "uniform.diagonal.b", "normal.sigA",
          "normal.sigB", "t.sig", "uniform.Frank", "uniform.Clayton",
          "uniform.Gumbel", "uniform.Joe", "mix.uniform.Clayton",
          "mix.Clayton.Gumbel", "mix.uniform.Frank",
          "mix.Clayton.Frank", "mix.Frank.Gumbel", "mix.Frank.Joe",
          "normal.shift-one.marginal", "normal.shift-all.marginal",
          "normal.stretch-one.marginal", "normal.stretch-all.marginal",
          "normal.stretch-shift-one.marginal", "normal.stretch-shift-all.marginal",
          "uniform.rotate.marginal", "uniform.beta-one.marginal",
          "uniform.beta-two.marginal", "uniform.all-normal.marginal",
          "beta22.all-normal.marginal", "uniform.all-beta.marginal",
          "uniform.all-exp.marginal", "uniform.all-gamma.marginal",
          "uniform.different.marginal"
          )
  if (ReturnCaseNames) return(cases)
  if (is.numeric(which))  which = cases[which]
  Range=matrix(rep(c(0,1), 5), 2, 5)
  null_param=MDgof::power_studies_cont_D5_results[[3]][,1]
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
      rnull=function() round(matrix(stats::runif(5*nsample), ncol=5), 5)     
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
    rnull=function() round(mvtnorm::rmvnorm(nsample, rep(0, 5)), 5)     
  } 
  add_stripe=function(n, x, a) {
    if(n==0) return(x)
    m=nrow(x)
    Range = range(x[,1])
    y=matrix(0, n, 5)
    y[,1]=stats::runif(n, Range[1], Range[2])    
    u=(y[,1]-min(y[,1]))/(max(y[,1])-min(y[,1]))    
    for(i in 2:5) {
      y[,i] = y[,1]+(1-a[1]*(u-0.5)^2)*stats::runif(n, -a[2], a[2])
      y[y[,i] < Range[1], i]=Range[1]-y[y[,i]<Range[1], i]
      y[y[,i] > Range[2], i]=2*Range[2]-y[y[,i]>Range[2], i]
      if(a[3]==1) 
        y[,i] = Range[1]+diff(Range)*(1-(y[,i]-Range[1])/diff(Range))
    }
    return(rbind(x[1:(m-n), ], y)) 
  }
# 1:  
  if(which=="uniform.diagonal.n") {
    return(list(
      pnull=pnull, dnull=dnull, rnull=rnull,    
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        x=round(matrix(stats::runif(5*nsample), ncol=5), 5)
        add_stripe(round(alpha*nsample), x, c(1, 0.2, 0))
      },
      param_alt=c(0, 0.15),
      Range=Range
    ))
  }  
# 2:   
  if(which=="uniform.diagonal.b") {
    return(list(
      pnull = pnull, rnull=rnull, dnull=dnull,
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        x=round(matrix(stats::runif(5*nsample), ncol=5), 5)         
        add_stripe(round(alpha*nsample), x, c(5, 0.3, 0))
      },
      param_alt=c(0, 0.15),
      Range=Range
    ))
  }
# 3: 
  if(which=="normal.sigA") {
    return(list(
      pnull=pnull, dnull=dnull, 
      rnull = rnull, 
      ralt=function(s) {
         if(s<0) s=null_param[which]
         S=(1-s)*diag(5)+matrix(s, 5, 5)
         round(mvtnorm::rmvnorm(nsample, sigma=S), 5)
      },
      param_alt=c(0, 0.15),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }
# 4: 
  if(which=="normal.sigB") {
    return(list(
      pnull=pnull, dnull=dnull, 
      rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=(1-s)*diag(5)+matrix(s, 5, 5)
        S[2,3]=S[3,2]=S[2,4]=S[4,2]=S[4,3]=S[3,4]=0
        round(mvtnorm::rmvnorm(nsample, sigma=S), 5)
      },
      param_alt=c(0, 0.15),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 5:  
  if(which=="t.sig") {
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
        round(mvtnorm::rmvt(nsample, sigma=diag(5), df=5), 5),
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=(1-s)*diag(5)+matrix(s, 5, 5)
        round(mvtnorm::rmvt(nsample, sigma=S, df=5), 5)
      },
      param_alt=c(0, 0.15),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# Case Studies where distribution under the null hypothesis is uniform 
# 6-9: 
  if(startsWith(which, "uniform") & !(endsWith(which, "marginal"))) {
    family=strsplit(which , "\\.")[[1]][2]
    ralt=function(s) {
        if(s==0) return(rnull())
        cop=gen.cop(family, s, 5)
        round(copula::rCopula(nsample, cop), 5)
    }
    if(family=="Frank") p=c(0, 1)
    if(family=="Clayton") p=c(0, 0.6)
    if(family=="Gumbel") p=c(0, 1.1)
    if(family=="Joe") p=c(0, 1.15)    
    return(list(
        pnull=pnull, dnull=dnull, rnull = rnull, 
        ralt=ralt, param_alt=p, Range=Range
      ))
  }
# Case Studies of mixtures 
# 10:15:  
  if(startsWith(which, "mix")) {
    if(which=="mix.uniform.Clayton") { 
      cop1=gen.cop("indep", 5)
      cop2=gen.cop("Clayton", p=8, 5)
      p=c(0.5, 0.25)
    }    
    if(which=="mix.Clayton.Gumbel") { 
      cop1=gen.cop("Clayton", p=8, 5)
      cop2=gen.cop("Gumbel", p=1.175, 5)
      p=c(0.5, 0.25)
    }
    if(which=="mix.uniform.Frank") { 
      cop1=gen.cop("indep", 5)
      cop2=gen.cop("Frank", p=20, 5)
      p=c(0.5, 0.25)
    }
    if(which=="mix.Clayton.Frank") { 
      cop1=gen.cop("Clayton", p=8, 5)
      cop2=gen.cop("Frank", p=20, 5)
      p=c(0.5, 0.04)
    }
    if(which=="mix.Frank.Gumbel") { 
      cop1=gen.cop("Frank", p=20, 5)
      cop2=gen.cop("Gumbel", p=2, 5)
      p=c(0.5, 0.15)
    }
    if(which=="mix.Frank.Joe") { 
      cop1=gen.cop("Frank", p=20, 5)
      cop2=gen.cop("Joe", p=2, 5)
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
        Range=Range
      ))
  }
# Studies with unequal marginals
# 16:  
  if(which=="normal.shift-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, mean = c(s, rep(s, 4))), 5)
      },  
      param_alt=c(0, 0.35),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 17:  
  if(which=="normal.shift-all.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        round(mvtnorm::rmvnorm(nsample, mean = rep(s, 5)), 5)
      },
      param_alt=c(0, 0.25),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 18:  
  if(which=="normal.stretch-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=diag(5)
        S[1,1]=s
        round(mvtnorm::rmvnorm(nsample, sigma = S), 5)
      },
      param_alt=c(1, 1.75),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 19:  
  if(which=="normal.stretch-all.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=s*diag(5)
        round(mvtnorm::rmvnorm(nsample, sigma = S), 5)
      },
      param_alt=c(1, 1.5),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 20:  
  if(which=="normal.stretch-shift-one.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=diag(5)
        S[1,1]=1+s
        round(mvtnorm::rmvnorm(nsample, mean=c(s, rep(0, 4)), sigma = S), 5)
      },
      param_alt=c(0, 0.3),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  }  
# 21:  
  if(which=="normal.stretch-shift-all.marginal") {
    return(list(
      pnull=pnull, dnull=dnull, rnull = rnull, 
      ralt=function(s) {
        if(s<0) s=null_param[which]
        S=(1+s)*diag(5)
        round(mvtnorm::rmvnorm(nsample, mean=rep(s, 5),sigma = S), 5)
      },
      param_alt=c(0, 0.3),
      Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
    ))
  } 
# 22:  
  if(which=="uniform.rotate.marginal") {
      return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(alpha) {
        if(alpha<0) alpha=null_param[which]
        x=matrix(stats::runif(5 * nsample), ncol = 5) 
        A=matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), 2, 2)
        x[,1:2]=x[,1:2]%*%A
        x[,3:4]=x[,3:4]%*%A
      },
      param_alt=c(0, 0.125),
      Range=Range
    ))
  }
# 23:  
  if(which=="uniform.beta-one.marginal") {
    return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(a) {
        if(a<0) a=null_param[which]
        matrix(stats::rbeta(5 * nsample, 1, a), ncol = 5)
      },
      param_alt=c(1, 1.25),
      Range=Range
    ))
  }
# 24:  
  if(which=="uniform.beta-two.marginal") {
    return(list(
      pnull=pnull,  rnull=rnull, dnull=dnull,
      ralt = function(a) {
        if(a<0) a=null_param[which]
        matrix(stats::rbeta(5 * nsample, a, a), ncol = 5)
      },
      param_alt=c(1, 1.35),
      Range=Range
    ))
  }
# 25-30: transform marginals
# 25:
  if(which%in%c("uniform.all-normal.marginal")) {
    H=function(x) prod(x)
    rH=function() stats::runif(5*nsample)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pnorm(x),
           f2=function(x) stats::pnorm(x),
           f3=function(x) stats::pnorm(x, 0, 2),
           f4=function(x) stats::pnorm(x, 1),
           f5=function(x) stats::pnorm(x, 0, 2)
    )
    trans=list(f1=function(x) stats::qnorm(x),
             f2=function(x) stats::qnorm(x),
             f3=function(x) stats::qnorm(x, 0, 2),
             f4=function(x) stats::qnorm(x, 1),
             f5=function(x) stats::qnorm(x, 0, 2)
    )  
    p=c(1, 1.45)
    Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)    
  }
# 26:  
  if(which%in%c("beta22.all-normal.marginal")) {
    H=function(x) prod(stats::pbeta(x, 2, 2))
    rH=function() stats::rbeta(5*nsample, 2, 2)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pnorm(x),
             f2=function(x) stats::pnorm(x),
             f3=function(x) stats::pnorm(x, 0, 2),
             f4=function(x) stats::pnorm(x, 1),
             f5=function(x) stats::pnorm(x, 0, 2)
    )
    trans=list(f1=function(x) stats::qnorm(x),
               f2=function(x) stats::qnorm(x),
               f3=function(x) stats::qnorm(x, 0, 2),
               f4=function(x) stats::qnorm(x, 1),
               f5=function(x) stats::qnorm(x, 0, 2)
    )
    p=c(2, 3)
    Range=matrix(rep(c(-Inf, Inf), 5), 2, 5)
  }
# 27:  
  if(which%in%c("uniform.all-beta.marginal")) {
    H=function(x) prod(x)
    rH=function() stats::runif(5*nsample)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pbeta(x, 1, 2),
             f2=function(x) stats::pbeta(x, 2, 2),
             f3=function(x) stats::pbeta(x, 0.5, 1),
             f4=function(x) stats::pbeta(x, 3, 2),
             f5=function(x) stats::pbeta(x, 1, 0.5)
    )
    trans=list(f1=function(x) stats::qbeta(x, 1, 2),
               f2=function(x) stats::qbeta(x, 2, 2),
               f3=function(x) stats::qbeta(x, 0.5, 1),
               f4=function(x) stats::qbeta(x, 3, 2),
               f5=function(x) stats::qbeta(x, 1, 0.5)
    )
    p=c(1, 1.25)
    Range=matrix(rep(c(0, 1), 5), 2, 5)
  }
# 28:  
  if(which%in%c("uniform.all-exp.marginal")) {
    H=function(x) prod(x)
    rH=function() stats::runif(5*nsample)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pexp(x, 1),
             f2=function(x) stats::pexp(x, 2),
             f3=function(x) stats::pexp(x, 2),
             f4=function(x) stats::pexp(x, 3),
             f5=function(x) stats::pexp(x, 1)
    )
    trans=list(f1=function(x) stats::qexp(x, 1),
               f2=function(x) stats::qexp(x, 2),
               f3=function(x) stats::qexp(x, 2),
               f4=function(x) stats::qexp(x, 3),
               f5=function(x) stats::qexp(x, 1)
    )
    p=c(1, 1.45)
    Range=matrix(rep(c(0, Inf), 5), 2, 5)
  } 
# 29:  
  if(which%in%c("uniform.all-gamma.marginal")) {
    H=function(x) prod(x)
    rH=function() stats::runif(5*nsample)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pgamma(x, 2, 2),
             f2=function(x) stats::pgamma(x, 2, 3),
             f3=function(x) stats::pgamma(x, 3, 2),
             f4=function(x) stats::pgamma(x, 3, 3),
             f5=function(x) stats::pgamma(x, 3, 4)
    )
    trans=list(f1=function(x) stats::qgamma(x, 2, 2),
               f2=function(x) stats::qgamma(x, 2, 3),
               f3=function(x) stats::qgamma(x, 3, 2),
               f4=function(x) stats::qgamma(x, 3, 3),
               f5=function(x) stats::qgamma(x, 3, 4)
    )
    p=c(1, 0.7)
    Range=matrix(rep(c(0, Inf), 5), 2, 5)
  }
# 30:  
  if(which%in%c("uniform.different.marginal")) {
    H=function(x) prod(x)
    rH=function() stats::runif(5*nsample)
    rHalt=function(a) stats::rbeta(5*nsample, a, a)
    cdf=list(f1=function(x) stats::pgamma(x, 2, 2),
             f2=function(x) stats::pnorm(x),
             f3=function(x) stats::punif(x),
             f4=function(x) stats::punif(x, 1, 3),
             f5=function(x) stats::pbeta(x, 0.5, 1)
    )
    trans=list(f1=function(x) stats::qgamma(x, 2, 2),
               f2=function(x) stats::qnorm(x),
               f3=function(x) stats::qunif(x),
               f4=function(x) stats::qunif(x, 1, 3),
               f5=function(x) stats::qbeta(x, 0.5, 1)
    )
    p=c(1, 0.7)
    Range=matrix(c(0, Inf, -Inf, Inf, 0, 1, 1, 3, 0, 1), 2, 5)
  }  
  if(which%in%c("uniform.all-normal.marginal",
                "beta22.all-normal.marginal",
                "uniform.all-beta.marginal",
                "uniform.all-exp.marginal",
                "uniform.all-gamma.marginal",
                "uniform.different.marginal")) {
    return(list(
      pnull=function(x) {
        if(!is.matrix(x)) x=rbind(x)
        for(i in 1:5) x[,i]=cdf[[i]](x[,i])
        apply(x, 1, H)
      },
      rnull=function() {
        x=matrix(rH(), ncol=5)  
        for(i in 1:5) x[,i]=trans[[i]](x[,i])
        x
      },
      ralt=function(a) {
        x=matrix(rHalt(a), ncol=5)  
        for(i in 1:5) x[,i]=trans[[i]](x[,i])
        x
      },
      param_alt=p,
      Range=Range
    ))
  }
}
