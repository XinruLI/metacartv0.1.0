#' Prune a tree
#'
#' Prune an initial rpart tree by "c-standard-error" rule.
#' @inheritParams prune.rpart
#' @param tree A initial tree fitted by rpart, needs to an rpart object.
#' @param c A scalar to prune the  tree by selecting the tree with minum cross-validation error plus the standard error multiplied by c.
#' @param ... Additional arguments passed to prune.rpart().
#' @return The pruned tree
#' @importFrom rpart prune.rpart
#' @importFrom methods is
#' @keywords internal
treepruner <- function(tree, c, ...){
  prune <- NULL
  if (!is(tree, "rpart")) stop("The pruned tree should be an rpart object")
  else {
    if (!is.numeric(c) | length(c) != 1 | c < 0) {
      stop("The pruning parameter c should be a positive constant number")
    } else {
      tree <- tree
      c <- c
      mindex <- which.min(tree$cptable[,4])  # find the row of the minimum x-error
      cp.minse <- tree$cptable[mindex,4] + c*tree$cptable[mindex,5]  # the minimum x-error + c*SE
      cp.row <- min(which(tree$cptable[,4]<= cp.minse))  # find the smallest tree within the minimum x-error + c*SE
      cp.take <- tree$cptable[cp.row, 1]  # get the cp value for the smallest tree
      prune(tree, cp=cp.take, ...)  # prune the tree
    }
  }
}


#' Fixed effect meta-tree
#'
#' A function to fit FE meta-trees to meta-analytic data.
#' The model is assuming fixed effect within subgroups and between subgroups.
#' The tree growing process is equivalent to the approach described in Li et al.(2017) using fixed effect weights in rpart (Therneau, Atkinson & Ripley)
#' @name FEmrt
#' @aliases FEmrt
#' @inheritParams rpart
#' @param formula A formula, with a response variable (usually the effect size) and the interested moderators but no interaction terms.
#' @param vi The column name of the sampling variance in the data set.
#' @param data A data frame of a meta-analytic data set, including the effect sizes, sampling variance, and the potential moderators.
#' @param subset optional expression that selects only a subset of the rows of the data.
#' @param c A non-negative scalar.The pruning parameter to prune the initial tree by "c-standard-error" rule.
#' @param control the rpart.control object that passes to rpart
#' @param ... Additional arguments passed to rpart().
#' @return If no moderator effect is detected, the function will return a list including the following objects:
#' @return n: The total number of the studies
#' @return Q: The Q-statistics for the heterogeneity test
#' @return df: The degree of freedoms of the heterogeneity test
#' @return pval.Q: The p-value for the heterogeneity test
#' @return g: The overall effect size for all studies
#' @return se: The standard error of the overall effect size
#' @return zval: The test statistic of the overall effect sie
#' @return pval: The p-value for the test statistic of the overall effect size
#' @return ci.lb: The lower bound of the confidence interval for the overall effect size
#' @return ci.ub: The upper bound of the confidence interval for the overall effect size
#' @return call: The matched call
#' @return If  moderator effect(s) is(are) detected, the function will return a list including the following objects:
#' @return tree: An rpart object. The pruned tree that represents the moderator effects and interaction effects between them.
#' @return n: The number of the studies in each subgroup
#' @return Qb: The between-subgroups Q-statistic
#' @return df: The degree of freedoms of the between-subgroups Q test
#' @return pval.Qb: The p-value of the between-subgroups Q test
#' @return Qw: The within-subgroup Q-statistic in each subgroup
#' @return g: The subgroup overall effect sizes, based on Hedges'g
#' @return se: The standard errors of subgroup overall effect sizes
#' @return zval: The test statistics of the subgroup overall effect sizes
#' @return pval: The p-value for the test statistics of the subgroup overall effect sizes
#' @return ci.lb: The lower bounds of the confidence intervals
#' @return ci.ub: The upper bounds of the confidence intervals
#' @return call: The matched call
#' @examples data(dat.BCT2015)
#' FEtree <- FEmrt(g ~ T1 + T2+ T4 +T25, vi = vi, data = dat.BCT2015, c=0.5)
#' print(FEtree)
#' summary(FEtree)
#' plot(FEtree)
#' @references Dusseldorp, E., van Genugten, L., van Buuren, S., Verheijden, M. W., & van Empe-
#' len, P. (2014). Combinations of techniques that effectively change health behavior: Evidence from meta-cart analysis. Health Psychology, 33 (12), 1530-1540. doi:
#'      10.1037/hea0000018.
#' @references Li, X., Dusseldorp, E., & Meulman, J. J. (2016). Meta-cart: A tool to identify interactions
#' between moderators in meta-analysis. British Journal of Mathematical and Statistical Psychology. In press. doi: 10.1111/bmsp.12088.
#' @importFrom rpart rpart.control
#' @importFrom rpart rpart
#' @importFrom rpart prune.rpart
#' @importFrom stats model.response
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom rpart prune.rpart
#' @export
FEmrt <- function(formula, data, vi, subset, c = 1,
                  control = rpart.control(xval=10, minbucket=5, minsplit=10, cp=0.0001),
                  ...) {
  Call <- match.call()
  wts.metacart <- NULL
  indx <- match(c("formula", "data", "vi", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  if (indx[3] == 0L)
    stop("The sampling variances need to be specified")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(temp)
  m$wts.metacart <- c(t(1/m["(vi)"]))
  tree <- rpart(formula, weights = wts.metacart, data = m, control = control, ...)
  prunedtree <- treepruner(tree, c*sqrt(mean(m$wts.metacart)))
  tree$cptable[ ,5] <- tree$cptable[ ,5]*sqrt(mean(m$wts.metacart))
  if (length(unique(prunedtree$where)) < 2) {
    warning("no moderator effect was detected")
    y <- model.response(m)
    v <- c(t(m["(vi)"]))
    n <- length(y)
    g <- sum(y/v)/sum(1/v)
    Q <- sum((y-g)^2/v)
    df <- n - 1
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/v))
    zval <- g/se
    pval <- pnorm(abs(zval), lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se

    res <- list(n = n ,  Q = Q,
                df = df, pval.Q = pval.Q, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, call = Call, cv.res = tree$cptable)
  } else {
    y <- model.response(m)
    v <- c(t(m["(vi)"]))
    treeframe <- prunedtree$frame
    n <- treeframe[treeframe$var == "<leaf>", 2]
    Qw <- treeframe[treeframe$var == "<leaf>", 4]
    g <- treeframe[treeframe$var == "<leaf>", 5]
    Qb <- treeframe[1,4] - sum(Qw)
    df <- length(unique(prunedtree$where))-1
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(v, prunedtree$where, function(x) sqrt(1/sum(1/x)))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    mod.names <- unique(prunedtree$frame$var[prunedtree$frame$var != "<leaf>"])
    res <- list(tree =  prunedtree, n = n, moderators =  mod.names, Qb = Qb, df = df, pval.Qb = pval.Qb,
                Qw = Qw, g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
                ci.ub = ci.ub, call = Call, cv.res = tree$cptable)


  }
  class(res) <- "FEmrt"
  res
}
