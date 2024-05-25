pepnorm<-function(data, method="ZpepQuad", robust = TRUE, centered = FALSE, standard = FALSE) {
  x=rownames(data)
  zW=mkZpepMx(x)
  x=gZp(method, zW)
  ynew <- gWE(ymat = data, x = x, robust, centered, standard)
}
gWE <- function(ymat, x, robust, centered, standard) {
  x <- cbind(1, x)
  iter.max <- 10
  n <- nrow(x)
  out <- apply(ymat, 2, function(y) {
    # initial fit
    fit <- lm.fit(x = x, y = y)
    if (robust) {
      for(i in 1:iter.max) {
        RSS <- sum(fit$residuals^2)
        w <- (4 + 1)/(4 + (n/RSS)*fit$residuals^2)
        fit <- lm.wfit(x = x, y = y, w = w)
      }
    }
    if (standard) {
      o <- order(fit$fitted.values)
      ro <- order(o)
      
      n.bin <- 10
      n.probe <- length(y)
      n.ppb <- floor(n.probe/n.bin)
      rem <- n.probe - n.ppb*n.bin
      bin <- factor(c(rep(1:n.bin, each = n.ppb), rep(n.bin, rem)))
      
      split.resid <- split(fit$residuals[o], bin)
      scaled.resid <- unlist(lapply(split.resid, function(x) x/sd(x)))
      fit$residuals <- scaled.resid[ro]
    }
    if (!centered) {
      fit$residuals <- fit$residuals + fit$coefficients[1]
    }
    return(fit$residuals)
  })
  return(out)
}


gZp = function(method, X)
{
  
  if(method == "ZpepQuad")
    X <- cbind(X, X^2)
  
  as.matrix(X)
}