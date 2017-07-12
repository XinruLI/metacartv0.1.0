#' A function to compute Q-between and residual heterogeneity
#'
#' @param yi the effect sizes
#' @param vi the sampling variances
#' @param mods the subgrouping moderator
#' @return Q-between and the residual heterogeneity under mixed effects model assumption
#' @keywords internal
rebetQ<- function(yi, vi, mods){
  wts = 1/vi
  wy = wts*yi
  wy2 = wts * yi^2
  Q <- tapply(wy2, mods, sum) - tapply(wy, mods, function(x) (sum(x))^2)/tapply(wts, mods, sum)
  df <- tapply(wy, mods, length)-1
  C <- tapply(wts, mods, sum) - tapply(wts, mods, function(x) sum(x^2))/ tapply(wts, mods, sum)
  tau2 <- (sum(Q) - sum(df))/sum(C)
  tau2 <- max(0, tau2)
  wstar = 1/(vi+tau2)
  wystar = wstar*yi
  wy2star = wstar*yi^2
  Qstar <- tapply(wy2star, mods, sum) - tapply(wystar, mods, function(x) (sum(x))^2)/tapply(wstar, mods, sum)
  Qstar.total <- sum(wy2star) - (sum(wystar))^2/sum(wstar)
  Qbet <- Qstar.total - sum(Qstar)
  return(c(Qbet, tau2))

}


#' A function to subgroup newdata according to the fitted tree
#'
#' @param x the fitted RE meta-tree results
#' @param newdata the new data set
#' @return a matrix consists of node lables
#' @keywords internal
REmrt.prednode <- function(x, newdata){
  tt <- terms(x$data)
  ms <- model.frame(delete.response(tt), newdata)
  oms <- model.frame(delete.response(tt), x$data)
  tree <- x[["tree"]]
  # if (any(sapply(ms, class) != sapply(oms, class)))
  #   stop("The type of the variables do not match")
  if(nrow(tree) < 2) {
    pred.node <- rep(1, nrow(ms))
  } else {
    tnode <- rep(1, nrow(ms))
    nodes <- tnode
    for (i in 1:(nrow(tree) - 1)){
      tinx <- which(tnode == tree[i+1, "pleaf"])
      tempm <- ms[tree[i+1, "mod"]]
      if(sapply(tempm, is.numeric) == TRUE) {
        tnode[tinx] <- ifelse(tempm[tinx,1] <= x[["cpt"]][[i]], 2*i, 2*i+1)
      } else {
        tnode[tinx] <- ifelse(tempm[tinx,1] %in% oms[,tree[i+1, "mod"]],
                              ifelse(tempm[tinx,1] %in% x[["cpt"]][[i]], 2*i, 2*i+1),
                              NA)
      }
      nodes <- cbind(nodes, tnode)
    }
    pred.node <- nodes
  }
  pred.node
}


#' A function to grow the tree
#'
#' @param mf the data.frame to grow the tree
#' @param newdata the new data set
#' @return a matrix consists of node lables
#' @keywords internal
REmrt.fit1<- function(mf, maxL){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(node = rep(1, nrow(mf)))
  res <- data.frame(Qb = rebetQ(y, vi, mods = nodemark)[1],
                    tau2 = rebetQ(y, vi, mods = nodemark)[2],
                    split = NA, mod = NA, pleaf = NA)
  for (i in 1:maxL){
    Dev<- -Inf
    cnode <- nodemark[ ,i]
    # Check if the parent leaves are all smaller than the minsplit
    # and only split nodes with subjects more than the minsplit
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) #[len.node >= minsplit]
    #if (length(nodes) == 0) break

    for (j in 1:length(nodes)){
      pleaf.inx <- cnode == as.numeric(nodes[j])
      for (k in 1:nmod){
        if(sapply(mods[mods.names[k]], is.factor) == TRUE) {
          chosenmod <- mods[pleaf.inx, k]
          mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
          cmod.ordinal <- mod.order[as.character(chosenmod)]
          cpoints <- sort(mod.order)
          if (length(cpoints) >= 2) {
            for (g in 1:(length(cpoints)-1)) {
              cnode.test <- cnode
              cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
              temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
              if (temp[1] > Dev) {
                Dev <- temp[1]
                msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                   split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tcpt <- names(mod.order[mod.order <= cpoints[g]])
                tnode <- cnode.test
              }
            }
          }
        } else{
          chosenmod <- mods[pleaf.inx, k]
          cpoints <- sort(unique(chosenmod))
          if (length(cpoints) >= 2) {
            for (g in 1:(length(cpoints)-1)) {
              cnode.test <- cnode
              cnode.test[pleaf.inx] <- ifelse( chosenmod <= cpoints[g], 2*i, 2*i+1)
              temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
              if (temp[1] > Dev) {
                Dev <- temp[1]
                tcpt <- cpoints[g]
                msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2], split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tnode <- cnode.test
              }
            }
          }
        }

      }
    }

    new.node <- tnode
    nodemark <- cbind(nodemark, new.node)
    res <- rbind(res, tres)
    cpt[[i]] <- tcpt
  }
  list(tree = res, node.split = nodemark, cpt = cpt, data = mf)

}
