#' A function to fit random effects RE meta-tree
#'
#' @name REmrt
#' @aliases REmrt
#' @param formula A formula, with a response variable (usually the effect size) and the interested moderators but no interaction terms.
#' @param vi The column name of the sampling variance in the data set.
#' @param data A data frame of a meta-analytic data set, including the effect sizes, sampling variance, and the potential moderators.
#' @param c A non-negative scalar.The pruning parameter to prune the initial tree by "c-standard-error" rule.
#' @param maxL the maximum number of splits
#' @param minsplit the minimal number of studies in a parent node to be split
#' @param delQ the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @param n.fold the number of folds to perform the cross-validation
#' @param ... Additional arguments to be passed.
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
#' @return tree: A data frame that represents the tree, with the Q-between and the tau^2 after each split.
#' @return n: The number of the studies in each subgroup
#' @return moderators: the names of identified moderators
#' @return Qb: The between-subgroups Q-statistic
#' @return tau2: The estimate for the residual heterogeneity
#' @return df: The degree of freedoms of the between-subgroups Q test
#' @return pval.Qb: The p-value of the between-subgroups Q test
#' @return g: The subgroup overall effect sizes, based on Hedges'g
#' @return se: The standard errors of subgroup overall effect sizes
#' @return zval: The test statistics of the subgroup overall effect sizes
#' @return pval: The p-value for the test statistics of the subgroup overall effect sizes
#' @return ci.lb: The lower bounds of the confidence intervals
#' @return ci.ub: The upper bounds of the confidence intervals
#' @return call: The matched call
#' @return cv.res: The cross-validation table
#' @return data: the data set subgrouped by the fitted tree
#' @importFrom stats terms model.response
#'@examples data(dat.BCT2015)
#' REtree <- REmrt(g ~ T1 + T2+ T4 +T25, vi = vi, data = dat.BCT2015, c = 0)
#' print(REtree)
#' @export

REmrt <- function(formula, data, vi, c = 1, maxL = 10L, minsplit = 2L, delQ = 0.001, n.fold = 10, ...){
  Call <- match.call()
  indx <- match(c("formula", "data", "vi"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  if (indx[3] == 0L)
    stop("The sampling variances need to be specified")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  cv.res <- REmrt.xvalid(mf, maxL = maxL, n.fold = n.fold)
  mindex <- which.min(cv.res[, 1])
  cp.minse <- cv.res[mindex,1] + c*cv.res[mindex,2]
  cp.row <- min(which(cv.res[,1]<= cp.minse))
  if (cp.row == 1) {
    warning("no moderator effect was detected")
    y <- model.response(mf)
    vi <- c(t(mf["(vi)"]))
    wts <- 1/vi
    wy <- wts*y
    wy2 <- wts * y^2
    n <- length(y)
    Q <- sum(wy2) - (sum(wy))^2/sum(wts)
    df <- nrow(mf) - 1
    C <- sum(wts) - sum(wts^2)/sum(wts)
    tau2 <- max(0, (Q-df)/C)
    vi.star <- vi + tau2
    g <- sum(y/vi.star)/sum(1/vi.star)
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/vi.star))
    zval <- g/se
    pval <- pnorm(abs(zval), lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se

    res <- list(n = n ,  Q = Q,
                df = df, pval.Q = pval.Q, tau2 = tau2, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, call = Call, data = mf, cv.res = cv.res)
    res.f <- res
  } else{
    y <- model.response(mf)
    vi <- c(t(mf["(vi)"]))
    res <- REmrt.fit1(mf, maxL = cp.row - 1, minsplit = minsplit, delQ = delQ)
    depth <- nrow(res$tree)
    tau2 <- res$tree$tau2[depth]
    vi.star <- vi + tau2
    subnodes <- c(t(res$node.split[ncol(res$node.split)]))
    wy.star <- y/vi.star
    n <- tapply(y, subnodes, length)
    g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
    df <- depth - 1
    Qb <- res$tree$Qb[depth]
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(vi.star, subnodes, function(x) sqrt(1/sum(1/x)))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    mod.names <- unique(res$tree$mod[!is.na(res$tree$mod)])
    mf$term.node <- subnodes
    res.f<- list(tree =  res$tree, n = n, moderators =  mod.names, Qb = Qb, tau2 = tau2, df = df, pval.Qb = pval.Qb,
         g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
         ci.ub = ci.ub, call = Call, cv.res = cv.res, data = mf, cpt = res$cpt)
   }
 class(res.f) <- "REmrt"
 res.f

}
