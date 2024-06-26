% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmaEmpCov.R
\name{ctmaEmpCov}
\alias{ctmaEmpCov}
\title{ctmaEmpCov}
\usage{
ctmaEmpCov(
  targetVariables = NULL,
  recodeVariables = c(),
  combineVariables = c(),
  combineVariablesNames = c(),
  missingVariables = c(),
  n.latent = NULL,
  Tpoints = NULL,
  sampleSize = NULL,
  pairwiseN = NULL,
  empcov = NULL
)
}
\arguments{
\item{targetVariables}{(col-/row-) number or names of the target variables}

\item{recodeVariables}{(col-/row-) number or names of the target variables require inverse coding}

\item{combineVariables}{list of vectors, which put together the targeted variables that should be used for composite variables}

\item{combineVariablesNames}{new names for combined variables - not really important}

\item{missingVariables}{missing variables}

\item{n.latent}{number of (latent) variables - actually it is the number of all variables}

\item{Tpoints}{number of time points.}

\item{sampleSize}{sample size}

\item{pairwiseN}{matrix of same dimensions as empcov containing possible pairwiseN.}

\item{empcov}{empirical correlation matrix}
}
\value{
returns a list with two elements. The first element (results$r) contains the adapted correlation matrix, and
the second element (results$pairwiseNNew) an adapted version of a matrix of pairwise N if pariwiseN was provided for
the original correlation matrix supplied.
}
\description{
changes a full covariance matrix by selecting target variables, recoding them, combining them (compute the
mean of two or more variables), and by adding rows/columns with NA if focal variables are not available.
}
\examples{
source17 <- c()
delta_t17 <- c(12)
sampleSize17 <- 440
empcov17 <- matrix(
  c( 1.00, -0.60, -0.36,  0.20,  0.62, -0.47, -0.18,  0.20,
    -0.60,  1.00,  0.55, -0.38, -0.43,  0.52,  0.27, -0.21,
    -0.36,  0.55,  1.00, -0.47, -0.26,  0.37,  0.51, -0.28,
     0.20, -0.38, -0.47,  1.00,  0.15, -0.28, -0.35,  0.56,
     0.62, -0.43, -0.26,  0.15,  1.00, -0.63, -0.30,  0.27,
    -0.47,  0.52,  0.37, -0.28, -0.63,  1.00,  0.55, -0.37,
    -0.18,  0.27,  0.51, -0.35, -0.30,  0.55,  1.00, -0.51,
     0.20, -0.21, -0.28,  0.56,  0.27, -0.37, -0.51,  1.00),
 nrow=8, ncol=8)
moderator17 <- c(3, 2)
rownames(empcov17) <- colnames(empcov17) <-
  c("Workload_1", "Exhaustion_1", "Cynicism_1", "Values_1",
    "Workload_2", "Exhaustion_2", "Cynicism_2", "Values_2")
targetVariables17 <-
  c("Workload_1", "Exhaustion_1", "Cynicism_1",
    "Workload_2", "Exhaustion_2", "Cynicism_2")
recodeVariables17 <- c("Workload_1", "Workload_2")
combineVariables17 <- list("Workload_1", c("Exhaustion_1", "Cynicism_1"),
                           "Workload_2", c("Exhaustion_2", "Cynicism_2"))
combineVariablesNames17 <- c("Demands_1",  "Burnout_1",
                             "Demands_2",  "Burnout_2")
missingVariables17 <- c();
results17 <- ctmaEmpCov(targetVariables = targetVariables17,
                        recodeVariables = recodeVariables17,
                        combineVariables = combineVariables17,
                        combineVariablesNames = combineVariablesNames17,
                        missingVariables = missingVariables17,
                        n.latent = 2, sampleSize = sampleSize17,
                        Tpoints = 2, empcov = empcov17)
empcov17 <- results17$r

}
