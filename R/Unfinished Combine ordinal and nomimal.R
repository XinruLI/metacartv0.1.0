chosenmod <- mods[pleaf.inx, k]
if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
  mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
  cmod.ordinal <- mod.order[as.character(chosenmod)]
  } else {
  mod.order <- unique(chosenmod)
  cmod.ordinal <- chosenmod}
cpoints <- sort(mod.order)
  if (length(cpoints) >= 2) {
    for (g in 1:(length(cpoints)-1)) {
      cnode.test <- cnode
      cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
      if (min(table(cnode.test[pleaf.inx])) <= minbucket) {
        Dev.new <- -Inf
      } else {
        temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
        Dev.new <- temp[1]
      }
      if (Dev.new > Dev) {
        Dev <- temp[1]
        if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
          tcpt <- names(mod.order[mod.order <= cpoints[g]])
          msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
        } else {
          tcpt <- cpoints[g]
          msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
          }
        tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                           split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
        tnode <- cnode.test
      }
    }
  }
