% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmaOptimizeInit.R
\name{ctmaOptimizeInit}
\alias{ctmaOptimizeInit}
\title{ctmaOptimizeInit}
\usage{
ctmaOptimizeInit(
  primaryStudies = NULL,
  activeDirectory = NULL,
  problemStudy = NULL,
  reFits = NULL,
  finishsamples = NULL,
  n.latent = NULL,
  coresToUse = c(1),
  indVarying = FALSE,
  randomScaleTime = c(1, 1),
  activateRPB = FALSE,
  checkSingleStudyResults = FALSE,
  customPar = FALSE,
  T0means = 0,
  manifestMeans = 0,
  manifestVars = NULL,
  CoTiMAStanctArgs = NULL,
  scaleTime = NULL
)
}
\arguments{
\item{primaryStudies}{list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}}

\item{activeDirectory}{activeDirectory}

\item{problemStudy}{number (position in list) where the problem study in primaryStudies is found}

\item{reFits}{how many reFits should be done}

\item{finishsamples}{number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).}

\item{n.latent}{number of latent variables of the model (hast to be specified)!}

\item{coresToUse}{if neg., the value is subtracted from available cores, else value = cores to use}

\item{indVarying}{control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means}

\item{randomScaleTime}{lower and upper limit of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)}

\item{activateRPB}{set to TRUE to receive push messages with 'CoTiMA' notifications on your phone}

\item{checkSingleStudyResults}{displays estimates from single study 'ctsem' models and waits for user input to continue.
Useful to check estimates before they are saved.}

\item{customPar}{logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)}

\item{T0means}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{manifestMeans}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{manifestVars}{define the error variances of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.}

\item{CoTiMAStanctArgs}{parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function}

\item{scaleTime}{scale time (interval) - sometimes desirable to improve fitting}
}
\value{
returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
best fit.
}
\description{
Initial fitting (i.e., applies \code{\link{ctmaInit}}) to a primary study reFit times to capitalize on chance for obtaining
a hard-to-find optimal fit. This could be very helpful if a primary yields out-of-range estimates, which could happen if the fitting
algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
ctmaOptimizeInit is like gambling, hoping that at leas one set of starting values (the number is tries is specified in the reFits argument)
eneables finding the global optimal fit. On unix-like machines (e.g. MacOS), this could be done in parallel mode if coresToUse > 1.
}
\note{
All but one of multiple cores are used on unix-type machines for parallel fitting

During fitting, not output is generated. Be patient.
}
\examples{
\dontrun{
optimFit313 <- ctmaOptimizeInit(primaryStudies=CoTiMAstudyList_3,
                                activeDirectory="/Users/tmp/",  # adapt!
                                problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
                                reFits=10,
                                n.latent=2)
summary(optimFit313)
}

}
