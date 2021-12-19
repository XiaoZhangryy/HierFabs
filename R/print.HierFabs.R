#' Result summary of a fitted HierFabs x.
#'
#' Similar to the usual print methods, this function summarize results from a fitted \code{'HierFabs'} object.
#' @export
#' @param x Fitted \code{'HierFabs'} model object..
#' @param digits  The number of significant digits for the coefficient estimates.
#' @param \dots Not used. Other arguments to print.
#'
#' @return No value is returned.
#' @importFrom methods as
#' @seealso \code{\link{HierFabs}}

print.HierFabs <- function(x, digits = max(getOption("digits")-2,3), ...) {
  pG = length(x$G.names)

  active = which(x$beta!=0)
  beta = sapply(x$beta, round, digits)
  if (x$ge) {
    pE = length(x$E.names)
    main.G = active[active<=pG]
    main.E = active[active<=(pG+pE)]
    inter  = setdiff(active, main.E) - (pE+pG)
    main.E = setdiff(main.E, main.G) - pG

    parient.G = (inter %/% pE) + 1
    parient.E = inter %% pE

    Total.E = union(main.E, parient.E)
    Total.G = union(main.G, parient.G)

    coe  = matrix(0, length(Total.G)+1, length(Total.E)+1)
    colnames(coe) = c("main effect", x$E.names[Total.E])
    rownames(coe) = c("main effect", x$G.names[Total.G])

    # main effect of G
    if(length(Total.G)>0) {
      id = sapply(main.G, function(x) which(Total.G==x))+1
      coe[id,1] = beta[main.G]
    }

    # main effect of E
    if(length(Total.E)>0) {
      id = sapply(main.E, function(x) which(Total.E==x))+1
      coe[1,id] = beta[main.E+pG]
    }

    # interactions
    for (i in 1:length(inter)) {
      id.E = which(Total.E == parient.E[i]) + 1
      id.G = which(Total.G == parient.G[i]) + 1
      coe[id.G, id.E] = beta[inter[i]+pE+pG]
    }
  } else {
    main.effect = active[active<=pG]
    interaction = active[active>pG] - pG
    n.inter = length(interaction)

    if(n.inter>0){
      parients = matrix(NA, 2, n.inter)
      c = 1
      for (i in 1:pG) {
        l.i = pG-i+x$diagonal
        if(interaction[c] > l.i) {
          interaction[c:n.inter] = interaction[c:n.inter] - l.i
          next
        } else {
          while (interaction[c] <= l.i) {
            if(i %in% main.effect) {
              parients[1,c] = i
              parients[2,c] = i+interaction[c]-x$diagonal
            } else {
              parients[2,c] = i
              parients[1,c] = i+interaction[c]-x$diagonal
            }
            c = c+1
            if(c > n.inter) break
          }
          interaction[c:n.inter] = interaction[c:n.inter] - l.i
          if(c > n.inter) break
        }
      }
      
      interaction = active[active>pG]
      cols = unique(parients[2,])
      coe  = matrix(0, length(main.effect), length(cols)+1)
      rownames(coe) = x$G.names[main.effect]
      colnames(coe) = c("main effect",x$G.names[cols])
      
      coe[,1] = beta[main.effect]
      loc = matrix(NA, nrow(parients), ncol(parients))
      loc[1,] = sapply(parients[1,], function(t) which(main.effect==t))
      loc[2,] = sapply(parients[2,], function(t) which(cols==t)+1)
      for (i in 1:ncol(loc)) {
        coe[loc[1,i], loc[2,i]] = beta[interaction[i]]
      }
    } else {
      coe  = matrix(0, length(main.effect), 1)
      coe[,1] = beta[main.effect]
      rownames(coe) = x$G.names[main.effect]
      colnames(coe) = c("main effect")
    }

    
  }
  # coe = as(coe, "sparseMatrix")
  # class(coe) <- "HierFabs"
  # return(coe)
  x$coef = as(coe, "sparseMatrix")
  x$coef
}
