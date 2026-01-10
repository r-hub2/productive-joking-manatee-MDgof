#' Ripley's K function test
#' 
#' this function calculates the test statistic of Ripley's K function test
#' 
#' @param x matrix with data
#' @return a number (test statistic)
#' @export
RipleyK <- function(x) {
  X=spatstat.geom::ppp(x = x[,1], y = x[,2], 
                  window = spatstat.geom::owin(c(0,1), c(0,1)))
  Kvals <- spatstat.explore::Kest(X, correction = "border")
  L_obs <- sqrt(Kvals$border / pi)
  L_theo <- Kvals$r 
  r <- Kvals$r
  valid <- (r >= 0.05) & (r <= 0.4)
  T_stat <- max(abs(L_obs[valid] - L_theo[valid]))
  return(T_stat)
}
