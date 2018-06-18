Importance <- function(x){
  # Calculate variable importance for an REmrt tree
  if (!inherits(REtree, "REmrt")) stop("argument should be a REmrt object")
  if (nrow(x$tree) < 2) stop("no moderator was detected")
  # compute the between-subgroups Q to which each split contributes
  delta.Q <- diff(x$tree$Qb)
  # Return a vector of the variable importance
  sort(tapply(delta.Q, x$tree$mod[-1], sum), decreasing = T)
  
}