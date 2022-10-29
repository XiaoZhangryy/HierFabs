#' A hierarchical Forward and Backward Stagewise (HierFabs) algorithm for identifying hierachical interaction of genomics data.
#'
#' @useDynLib HierFabs, .registration = TRUE
#' @export
#' @param G Gene matrix, each row is an observation vector.
#' @param y Response variable. For logistic regression model, y takes value at 1 and -1. 
#' @param E An optional environment matrix. If Z is given, the interactions between environment and gene are of interest. Otherwise, the gene-gene interactions are of interest.
#' @param weight An optional weights. Default is 1 for each observation.
#' @param model A character string representing one of the built-in models. 'gaussian' for linear model and 'cox' for cox model.
#' @param back An indicator of whether to take backward steps.
#' @param stoping An indicator of whether to stop iteration when lambda is less than lambda.min.
#' @param eps Step size. Default is 0.01.
#' @param xi A tolerate to improve program stability. Default is 10^-6.
#' @param iter Maximum number of outer-loop iterations allowed. Default is 3000.
#' @param lambda.min Smallest value for lambda. \code{NULL} by default. 
#' @param lambda.ratio Smallest value for lambda, as a fraction of lambda.max. Works when lambda.min is \code{NULL}.
#' @param hier Whether to enforce strong or weak heredity. Default is 'strong'.
#' @param max_s Limit the maximum number of variables in the model. When exceed this limit, program will early stopped.
#' @param diagonal An indicator of whether to include "pure" quadratic terms. Work when gene-gene interactions are of interest.
#' @param status A censoring indicator.
#' @param gamma A tuning parameter in EBIC.
#' @param tau parameter for quantile regression.
#' @param criteria The criteria used to select the optimal solution.
#'
#' @return A list.
#' \itemize{
#'   \item theta - The coefficients of covariates, each column is a solution.
#'   \item beta - The optimal solution.
#'   \item lambda - Lambda sequence.
#'   \item direction - Update indicator. 1 means take a forward step, -1 means take a backward step.
#'   \item iter - Number of iterations.
#'   \item EBIC - EBIC for each solution.
#'   \item loss - loss for each solution.
#'   \item df - Number of nonzero coefficients.
#'   \item opt - Position of the optimal tuning based on EBIC.
#'   \item intercept - The intercept term, which appearance is due to standardization.
#' }
#' @seealso \code{\link{predict.HierFabs}}, \code{\link{print.HierFabs}}
#'
#' @examples
#' set.seed(0)
#' n = 500
#' p = 100
#' x = matrix(rnorm(n*p),n,p)
#' eta = x[,1:4] %*% rep(1,4) + 3*x[,1]*x[,2] + 3*x[,1]*x[,4]
#' y =  eta + rnorm(n)
#' xtest = matrix(rnorm(n*p),n,p)
#' eta.test = xtest[,1:4] %*% rep(1,4) + 3*xtest[,1]*xtest[,2] + 3*xtest[,1]*xtest[,4]
#' ytest =  eta.test + rnorm(n)
#' fit.gg.strong = HierFabs(x, y)
#' y.pred.gg.s = predict(fit.gg.strong, xtest, ytest)
#' y.pred.gg.s$mse
#' print(fit.gg.strong)
#'
#' ## Weak hierarchy
#' fit.gg.weak = HierFabs(x, y, hier="weak")
#' y.pred.gg.w = predict(fit.gg.weak, xtest, ytest)
#' y.pred.gg.w$mse
#' print(fit.gg.weak)
#'
#' ## Cox model with Gene-Environment interactions
#' pz = 10
#' z = matrix(rnorm(n*pz),n,pz)
#' eta.ge = x[,1:4] %*% rep(1,4) + z[,1] + z[,2] + 3*x[,1]*z[,1] + 3*x[,2]*z[,2]
#' err = log(rweibull(n, shape = 1, scale = 1))
#' y0 = exp(-eta.ge + err)
#' cens = quantile(y0, 0.9)
#' y.ge = pmin(y0, cens)
#' status = 1 * (y0<=cens)
#' ztest = matrix(rnorm(n*pz),n,pz)
#' eta.ge.test = rowSums(xtest[,1:4]) + ztest[,1] + ztest[,2]
#' eta.ge.test = eta.ge.test + 3*xtest[,1]*ztest[,1] + 3*xtest[,2]*ztest[,2]
#' err.test = log(rweibull(n, shape = 1, scale = 1))
#' y0.test = exp(-eta.ge.test + err.test)
#' cens = quantile(y0.test, 0.9)
#' y.ge.test = pmin(y0.test, cens)
#' status.test = 1 * (y0.test<=cens)
#' fit.ge.strong = HierFabs(x, y.ge, z, model="cox", status=status)
#' y.pred.ge.s = predict(fit.ge.strong, xtest, y.ge.test, ztest, status.test)
#' y.pred.ge.s$c.index
#' print(fit.ge.strong)
#'
#' ## Weak hierarchy
#' fit.ge.weak = HierFabs(x, y.ge, z, model="cox", status=status, hier="weak")
#' y.pred.ge.w = predict(fit.ge.weak, xtest, y.ge.test, ztest, status.test)
#' y.pred.ge.w$c.index
#' print(fit.ge.weak)

HierFabs = function(G, y, E, weight = NULL, model = c("gaussian", "cox", "quantile", "logistic"), back = TRUE,
  stoping = TRUE, eps = 0.01, xi = 10^-6, iter = 3000, lambda.min = NULL, lambda.ratio = NULL, 
  hier = c("strong", "weak"), max_s = NULL, diagonal = FALSE, status = NULL, gamma = NULL, 
  tau = NUL, criteria = c("EBIC", "BIC"))
{
  n  = nrow(G)
  px = ncol(G)

  if (missing(E)) {
    q = as.integer(px*(px+1)/2 + px*diagonal)
    pz = 0
    ge = F
    E  = 1.0
  } else {
    pz = ncol(E)
    q  = px+pz+px*pz
    ge = T
    if(is.null(colnames(E))) {
      colnames(E) = paste("E", 1:pz, sep = "")
    }
  }

  if (is.null(weight))         weight = rep(1, q)
  if (is.null(status))         status = rep(1, n)
  if (is.null(gamma)) gamma = 1-log(n)/(2*log(q))
  if(is.null(colnames(G))) colnames(G) = paste("G", 1:px, sep = "")


  param     = c(n, px, q, 0, 0, pz, 0, 0)
  hier = match.arg(hier)
  criteria = match.arg(criteria)
  model     = match.arg(model)
  meany     = mean(y)
  if (model != "logistic") {
    y = y - meany
  }
  if (is.null(lambda.min)) {
    if (is.null(lambda.ratio)) lambda.ratio = {if (n > q) 1e-4 else .02}
  } else {
    param[8] = 1
    lambda.ratio = lambda.min
  }
  

  if (is.null(max_s))  {
    tmp   = min(n, q)
    max_s = as.integer(tmp/log(tmp))
    #if (tmp < 200)
    #  max_s = as.integer(tmp/log(tmp))
    #else
    #  max_s = as.integer(sqrt(tmp))
  }
  if (model == "cox") {
    y.order = order(y)
    G       = G[y.order, ]
    y       = y[y.order]
    status  = status[y.order]
    if (!missing(E)) E = E[y.order, ]
  }

  if (model != "quantile") tau = 0

  fit <- .Call("Hierarchy_Fabs",
               as.numeric(G),
               as.numeric(E),
               as.numeric(y),
               as.numeric(weight),
               as.character(model),
               as.numeric(eps),
               as.numeric(lambda.ratio),
               as.numeric(xi),
               as.integer(back),
               as.integer(stoping),
               as.integer(iter),
               as.character(hier),
               as.integer(param),
               as.integer(max_s),
               as.numeric(meany),
               as.integer(diagonal),
               as.integer(status),
               as.integer(ge),
               as.numeric(gamma),
               as.numeric(tau) )

  theta <- sparseMatrix(fit$index_i, fit$index_j, x = fit$beta, dims = c(q, fit$iter), index1 = FALSE)
  if (criteria== "BIC") {
    fit$df2 = fit$df/2 + Matrix::colSums(theta[1:(px+pz),]!=0)/2
    if (model == "cox") {
      fit$bic = 2*fit$loss+fit$df2*log(n)
    } else {
      fit$bic = 2*n*fit$loss+fit$df2*log(n)
    }
  }
  opt      = which.min(fit$bic)
  opttheta = theta[,opt,drop=FALSE]

  val = list(theta     = theta,
             beta      = opttheta,
             lambda    = fit$lambda,
             direction = fit$direction,
             iter      = fit$iter,
             EBIC       = fit$bic,
             loss      = fit$loss,
             df        = fit$df,
             opt       = opt,
             intercept = fit$intercept,
             ge        = ge,
             diagonal  = diagonal,
             G.names   = colnames(G),
             model     = model
             )
  if(ge) val$E.names   = colnames(E)
  class(val) <- "HierFabs"
  return(val)
}


