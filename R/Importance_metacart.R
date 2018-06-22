Importance <- function(x){
  # Calculate variable importance for an REmrt tree
  # Argument:
  # x: a FEmrt object or a REmrt object
  # Return:
  # a vector of the variable importance which is computed base on
  # the contribution of the increase of between-subgroup Q-statistic
  
  if (inherits(x, "REmrt")) {
    if (nrow(x$tree) < 2) stop("no moderator was detected")
    delta.Q <- diff(x$tree$Qb)
    sort(tapply(delta.Q, x$tree$mod[-1], sum), decreasing = T)
  } else {
    if (inherits(x, "FEmrt")) {
      x$tree$variable.importance[names(x$tree$variable.importance) %in% x$moderators]
    } else {
      stop("argument should be a REmrt object")
    }
    
  }
  }

################## An Example########################
library(metacart)
fit <- REmrt(efk ~ m1 + m2 + m3 + m4 + m5, vi = vark, data = SimData)
fit2 <- FEmrt(efk ~ m1 + m2 + m3 + m4 + m5, vi = vark, data = SimData)
Importance(fit)
Importance(fit2)

