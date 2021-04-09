#' Model prediction based on a fitted HierFabs object.
#'
#' Similar to the usual predict methods, this function returns predictions from a fitted \code{'HierFabs'} object.
#' @export
#' @param object Fitted \code{'HierFabs'} model object.
#' @param new.G Matrix of new values for genes.
#' @param new.y Vector of new values for response. 
#' @param new.E An optional new values for environment variables.
#' @param new.s An optional new values for status variables. Used in cox model prediction only. 
#' @param allpath allpath = T will output all the predictions on the solution path. allpath = FALSE will only output the optimal solution selected in the 'HierFabs' object.
#' @param \dots Not used. Other arguments to predict.
#'
#' @return A list.
#' \itemize{
#'   \item pred - The linear predict.
#'   \item mse - Mean square error for linear model.
#'   \item c.index - C index for cox model.
#' }
#' @seealso \code{\link{HierFabs}}

predict.HierFabs <- function(object, new.G, new.y, new.E, new.s, allpath=FALSE, ...) {
  if (missing(new.E) && object$ge)
    stop("Gene-Environment interactions need environment inputs.\n")
  n  = nrow(new.G)
  pG = ncol(new.G)
  
  if (allpath == FALSE) {
    id = object$opt
  } else {
    id = c(1:object$iter)
  }
  eta = matrix(rep( object$intercept[id], n ), n, length(id), byrow = TRUE)
  if (object$ge) {
    pE = ncol(new.E)
    eta = eta + new.G %*% object$theta[1:pG,id,drop = FALSE]
    eta = eta + new.E %*% object$theta[(pG+1):(pG+pE),id,drop = FALSE]
    loc = pG + pE + 1
    for (i in 1:pG) {
      eta = eta + (new.G[, i] * new.E) %*% object$theta[loc:(loc+pE-1),id,drop = FALSE]
      loc = loc + pE
    }
  } else {
    eta = eta + new.G %*% object$theta[1:pG,id,drop = FALSE]
    loc = pG + 1
    if (object$diagonal) {
      for (i in 1:pG) {
        eta = eta + (new.G[, i:pG] * new.G[,i]) %*% object$theta[loc:(loc+pG-i),id,drop = FALSE]
        loc = loc + pG - i + 1
      }
    } else {
      for (i in 1:(pG-1)) {
        eta = eta + (new.G[, (i+1):pG] * new.G[,i]) %*% object$theta[loc:(loc+pG-i-1),id,drop = FALSE]
        loc = loc + pG - i
      }
    }
  }
  
  if(object$model == 'gaussian') {
    mse = apply(eta, 2, function(x) mean((new.y - x)^2))
  }
  if(object$model == 'cox') {
    c.index = apply(eta, 2, function(etai) {
      y.order = order(new.y)
      etai = etai[y.order]
      new.s = new.s[y.order]
      n = length(new.y)
      
      sum2 = new.s %*% seq(n-1,0,-1)
      sum1 = 0
      for (i in 1:(n-1)) {
        for (k in (i+1):n) {
          sum1 = sum1 + new.s[i]*(etai[i]>etai[k])
        }
      }
      return(sum1/sum2)
    })
  }
  
  if (allpath == FALSE) eta = as.vector(eta)
  
  val = list(pred = eta)
  if(object$model == 'gaussian') val$mse = mse
  if(object$model == 'cox') val$c.index = c.index
  return(val)
}


