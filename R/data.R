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
#' The complete data consist of 101 studies reporting 122 interventions targeted at physical activity and healthy eating.
#' In this subset of the data, the interventions that include at least one
#' of the motivation-enhancing BCTs were selected (N = 106).
#'
#' @name  dat.BCT2015
#' @docType data
#'
#' @usage data(dat.BCT2015)
#'
#' @keywords data
#' @format A data frame of 106 interventions with five motivation-enhancing BCTs
#' \itemize{
#'   \item study: The name of the intervention.
#'   \item g: The effect size of each intervention.
#'   \item vi: The sampling variance of each study.
#'   \item T1: Indicating whether the BCT1 "Provide information about
#'   behavior-health link" was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T2: Indicating whether the BCT2 "Provide information on consequences"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T3: Indicating whether the BCT3 "Provide information about other's approval"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T4: Indicating whether the BCT4 "Prompt intention formation"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T25: Indicating whether the BCT25 " Motivational interviewing"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'
#'
#'   }
"dat.BCT2015"
