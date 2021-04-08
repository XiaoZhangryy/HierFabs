#' Model prediction based on a fitted HierFabs object.
#'
#' Similar to the usual predict methods, this function returns predictions from a fitted \code{'HierFabs'} object.
#' @export
#' @param object Fitted \code{'HierFabs'} model object.
#' @param new.G Matrix of new values for genes.
#' @param new.E An optional new values for environment variables.
#' @param allpath allpath = T will output all the predictions on the solution path. allpath = FALSE will only output the optimal solution selected in the 'HierFabs' object.
#' @param \dots Not used. Other arguments to predict.
#'
#' @return Predictions for each observation.
#' @seealso \code{\link{HierFabs}}

predict.HierFabs <- function(object, new.G, new.E, allpath=FALSE, ...) {
    if (missing(new.E) && object$ge)
        stop("Gene-Environment interactions need environment inputs.\n")
    n  = nrow(new.G)
    pG = ncol(new.G)

    eta = matrix(rep( object$intercept, n ), n, object$iter, byrow = TRUE)
    if (object$ge) {
        pE = ncol(new.E)
        eta = eta + new.G %*% object$theta[1:pG,]
        eta = eta + new.E %*% object$theta[(pG+1):(pG+pE),]
        loc = pG + pE + 1
        for (i in 1:pG) {
            eta = eta + (new.G[, i] * new.E) %*% object$theta[loc:(loc+pE-1),]
            loc = loc + pE
        }
    } else {
        eta = eta + new.G %*% object$theta[1:pG,]
        loc = pG + 1
        if (object$diagonal) {
            for (i in 1:pG) {
                eta = eta + (new.G[, i:pG] * new.G[,i]) %*% object$theta[loc:(loc+pG-i),,drop = FALSE]
                loc = loc + pG - i + 1
            }
        } else {
            for (i in 1:(pG-1)) {
                eta = eta + (new.G[, (i+1):pG] * new.G[,i]) %*% object$theta[loc:(loc+pG-i-1),,drop = FALSE]
                loc = loc + pG - i
            }
        }
    }
    if (allpath == FALSE) {
        eta = eta[, object$opt]
    }
    return(eta)
}
