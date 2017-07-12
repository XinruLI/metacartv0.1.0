#' A simulated meta-analytic data set
#'
#' Data simuated from a true model with a three-way interaction between three
#' moderators: m1, m2 and m3. If the values of the three moderators are all "B"s
#' the ture effect size will be 0.80. Otherwise, the true effect size is 0.
#' @name  SimData
#' @docType data
#'
#' @usage data(SimData)
#'
#' @keywords simulated datasets
#' @format A data frame of 80 studies with 6 moderators
#' \itemize{
#'   \item efk: The effect size of each study
#'   \item vark: The sampling variance of each study
#'   \item m1 to m5: Five randomly generated moderators. m1 and m2 have two levels (A and B),
#'   whereas m3, m4 and m5 have three levels (A, B and C)
#'   }
"SimData"



#' A subset of meta-analytic data set in Michie et al. (2015)
#'
#' Data are heal psychological interventions that used at least one of the
#' motivation-enchancing behavior change techniques T1, T2, T4, and T25.
#'
#' @name  datphase1
#' @docType data
#'
#' @usage data(datphase1)
#'
#' @keywords data
#' @format A data frame of 106 studies with 26 moderators
#' \itemize{
#'   \item g: The effect size of each study
#'   \item vark: The sampling variance of each study
#'   \item T1 to T26: moderators. Behavior change techniques.
#'   }
"datphase1"
