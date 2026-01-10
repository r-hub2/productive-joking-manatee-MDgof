#' Chi-square test for 2D data
#' 
#' This function does the chi square goodness-of-fit test for continuous data in two dimensions.
#' 
#' @param dta  a matrix of numbers.
#' @param pnull function to calculate expected counts.
#' @param phat =function(x) -99, function to estimate parameters of pnull.
#' @param Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2), a 2x2 matrix with lower and upper bounds
#' @param nbins =c(5,5) number of bins in x and y direction
#' @param minexpcount =5 minimum counts required per bin
#' @param SuppressMessages =FALSE, should info be shown?
#' @return a matrix with statistics, p values and degree of freedoms
#' @export
chi_cont_test=function(dta, pnull, phat=function(x) -99, 
                      Ranges=matrix(c(-Inf, Inf, -Inf, Inf),2,2),
                      nbins=c(5,5), minexpcount=5, SuppressMessages =TRUE) {
    if(ncol(dta)!=2) {message("Test is for two dimensional data only!");return(NULL)}
    if(abs(phat(dta)[1]+99)<1e-5) WithEst=FALSE
    else {
       p = phat(dta)
       WithEst=TRUE
    }
    pfun=function(x1,x2,y1,y2) {
       if(WithEst)
         pnull(c(x2,y2), p)-pnull(c(x1,y2), p)-pnull(c(x2,y1), p)+pnull(c(x1,y1), p)
      else pnull(c(x2,y2))-pnull(c(x1,y2))-pnull(c(x2,y1))+pnull(c(x1,y1))
    }
    out=matrix(0, 2, 3)
    rownames(out) = c("ES", "EP")
    colnames(out) = c("Statistic", "p.value", "df")
    for(k in 1:2) {
      grd=as.list(1:2)
      for(i in 1:2) {
        if(k==1) {
          grd[[i]]=seq(min(dta[,i]), max(dta[,i]), length=nbins[i]+1)
        }  
        else  {
          grd[[i]]=quantile(dta[,i], 0:nbins[i]/nbins[i])
        }  
        grd[[i]][c(1,nbins[i]+1)]=Ranges[, i]
      } 
      A=0*dta
      for(i in 1:2) {
         A[,i]=cut(dta[,i], grd[[i]], labels=1:nbins[i],
                   include.lowest=TRUE)
      }
      O=matrix(0, nbins[1], nbins[2])
      for(i in 1:nrow(dta)) {
        O[A[i,1], A[i,2]]=O[A[i,1], A[i,2]]+1
      }
      E=O
      for(i in 1:nbins[1]) {
        for(j in 1:nbins[2]) {
            E[i,j] = pfun(grd[[1]][i], grd[[1]][i+1],grd[[2]][j], grd[[2]][j+1])
        }
      }
      E=sum(O)*E
      step=1
      repeat {
        if(min(E,na.rm=TRUE)>minexpcount) break
        lowestE = c(1, 1)
        tmp = max(E, na.rm = TRUE)
        for(i in 1:nbins[1])
         for(j in 1:nbins[2]) {
           if(!is.na(E[i,j]) && E[i,j]<tmp) {
             tmp=E[i,j]
             lowestE=c(i,j)
           }
        }
        I=lowestE[1]+c(-1,0,1)*step
        I=I[1<=I&I<=nbins[2]]
        J=lowestE[2]+c(-1,0,1)*step
        J=J[1<=J&J<=nbins[1]]
        I=I[I<=nbins[1]]
        J=J[J<=nbins[2]]
        tmp=max(E,na.rm = TRUE)
        nearestE=c(-1,-1)
        for(i in I) {
         for(j in J) {
           if(i==lowestE[1]&&j==lowestE[2]) next
           if(!is.na(E[i,j])&&E[i,j]<tmp) {
             tmp=E[i,j]
             nearestE = c(i,j)
           }
         }
        }  
        if(step>min(nbins)) {
           if(!SuppressMessages)
               message("Not enough data for chi square test! Decrease nbins.")
           return(NULL)
        }  
        if(nearestE[1]<0) {step=step+1;next}
        O[lowestE[1], lowestE[2]] = O[lowestE[1], lowestE[2]] + 
                                   O[nearestE[1], nearestE[2]]
        O[nearestE[1], nearestE[2]] = NA
        E[lowestE[1], lowestE[2]] = E[lowestE[1], lowestE[2]] + 
                                   E[nearestE[1], nearestE[2]]
        E[nearestE[1], nearestE[2]] = NA
        step=1
      }
      chi=sum(O[!is.na(O)]^2/E[!is.na(E)])-sum(O[!is.na(O)])
      df=length(O[!is.na(O)])-1-ifelse(WithEst, length(p), 0)
      if(df<1) {
        if(!SuppressMessages)
          message("degrees of freedom is not positive")
        out[k, ] = c(chi, NA, df) 
      }
      else out[k, ] = c(chi, 1-stats::pchisq(chi, df), df) 
    }
    out
}

#' Chi-square test for discrete 2D data
#' 
#' This function does the chi square goodness-of-fit test for discrete data in two dimensions.
#' 
#' @param dta  a matrix of numbers.
#' @param pnull distribution function to calculate expected counts.
#' @param dnull density to calculate expected counts.
#' @param phat =function(x) -99, function to estimate parameters of pnull.
#' @param minexpcount =5 minimum counts required per bin
#' @param SuppressMessages =TRUE, should info be shown?
#' @return a vector with statistic, p value and degree of freedom
#' @export
chi_disc_test=function(dta, pnull, dnull, phat=function(x) -99,
               minexpcount=5, SuppressMessages =FALSE) {
  if(abs(phat(dta)[1]+99)<1e-5) WithEst=FALSE
  else {
    p = phat(dta)
    WithEst=TRUE
  }
  if(missing(pnull)) {
     WithDensity=TRUE
     if(length(formals(dnull))==1) 
        h=function(x) dnull(x)
     else h=function(x) dnull(x, phat(dta))
  }
  else {
     WithDensity=FALSE
     if(length(formals(pnull))==1) 
        h=function(x) pnull(x)
     else h=function(x) pnull(x, phat(dta))
  }
  out=c(0, 0, 0)
  names(out) = c("Statistic", "p.value", "df")
  valsx=sort(unique(dta[,1]))
  valsy=sort(unique(dta[,2]), decreasing = TRUE)
  nbins=c(length(valsx), length(valsy))
  O=matrix(0, nbins[2], nbins[1])
  rownames(O)=valsy
  colnames(O)=valsx
  for(i in 1:nbins[2]) {
    for(j in 1:nbins[1]) {
      O[i, j]=dta[dta[,1]==valsx[j]&dta[,2]==valsy[i],3]
    }  
  }
  Etmp=0*O
  for(i in nbins[2]:1) {
    for(j in 1:nbins[1]) {
      tmp=c(valsx[j], valsy[i])
      Etmp[i,j] = h(tmp)
    }      
  }
  E=Etmp
  if(!WithDensity) {
    Etmp=cbind(0, Etmp)
    Etmp=rbind(Etmp, 0)
    for(i in nbins[2]:1) {
      for(j in 1:nbins[1]) {
        E[i,j] = Etmp[i,j+1]-Etmp[i, j]-Etmp[i+1, j+1]+Etmp[i+1,j]
      }      
    }
  }  
  E=sum(O)*E/sum(E)
  step=1
  repeat {
    if(min(E,na.rm=TRUE)>minexpcount) break
    lowestE = c(1, 1)
    tmp = max(E, na.rm = TRUE)
    for(i in 1:nbins[2])
      for(j in 1:nbins[1]) {
        if(!is.na(E[i,j]) && E[i,j]<tmp) {
          tmp=E[i,j]
          lowestE=c(i,j)
        }
      }
    I=lowestE[1]+c(-1,0,1)*step
    I=I[1<=I&I<=nbins[2]]
    J=lowestE[2]+c(-1,0,1)*step
    J=J[1<=J&J<=nbins[1]]
    tmp=max(E,na.rm = TRUE)
    nearestE=c(-1,-1)
    for(i in I) {
      for(j in J) {
        if(i==lowestE[1]&&j==lowestE[2]) next
        if(!is.na(E[i,j])&&E[i,j]<tmp) {
          tmp=E[i,j]
          nearestE = c(i,j)
        }
      }
    }   
    if(step>min(nbins)) {
      if(!SuppressMessages)
        message("Not enough data for chi square test! Decrease nbins.")
      return(NULL)
    } 
    if(nearestE[1]<0) {step=step+1;next}
    O[lowestE[1], lowestE[2]] = O[lowestE[1], lowestE[2]] + 
      O[nearestE[1], nearestE[2]]
    O[nearestE[1], nearestE[2]] = NA
    E[lowestE[1], lowestE[2]] = E[lowestE[1], lowestE[2]] + 
      E[nearestE[1], nearestE[2]]
    E[nearestE[1], nearestE[2]] = NA
    step=1
  }
  chi=sum(O[!is.na(O)]^2/E[!is.na(E)])-sum(O[!is.na(O)])
  df=length(O[!is.na(O)])-1-ifelse(WithEst, length(p), 0)
  if(df<1) {
    if(!SuppressMessages)
      message("degrees of freedom is not positive")
    out = c(chi, NA, df) 
  }
  else out = c(chi, 1-stats::pchisq(chi, df), df) 
  out
}


