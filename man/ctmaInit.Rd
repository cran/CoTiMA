% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmaInit.R
\name{ctmaInit}
\alias{ctmaInit}
\title{ctmaInit}
\usage{
ctmaInit(
  activateRPB = FALSE,
  activeDirectory = NULL,
  binaries = NULL,
  chains = NULL,
  checkSingleStudyResults = FALSE,
  cint = 0,
  coresToUse = c(2),
  CoTiMAStanctArgs = NULL,
  customPar = FALSE,
  diff = NULL,
  digits = 4,
  doPar = 1,
  drift = NULL,
  experimental = FALSE,
  finishsamples = NULL,
  indVarying = FALSE,
  indVaryingT0 = NULL,
  iter = NULL,
  lambda = NULL,
  loadSingleStudyModelFit = c(),
  manifestMeans = 0,
  manifestVars = NULL,
  n.latent = NULL,
  n.manifest = 0,
  nopriors = FALSE,
  optimize = TRUE,
  posLL = TRUE,
  primaryStudies = NULL,
  priors = FALSE,
  sameInitialTimes = FALSE,
  saveRawData = list(),
  saveSingleStudyModelFit = c(),
  scaleTI = NULL,
  scaleTime = NULL,
  silentOverwrite = FALSE,
  T0means = 0,
  T0var = "auto",
  useSV = FALSE,
  verbose = NULL
)
}
\arguments{
\item{activateRPB}{set to TRUE to receive push messages with 'CoTiMA' notifications on your phone}

\item{activeDirectory}{defines another active directory than the one used in \code{\link{ctmaPrep}}}

\item{binaries}{which manifest is a binary. Still experimental}

\item{chains}{number of chains to sample, during HMC or post-optimization importance sampling.}

\item{checkSingleStudyResults}{Displays estimates from single study ctsem models and waits for user input to continue. Useful to check estimates before they are saved.}

\item{cint}{default 'auto' (= 0). Are set free if random intercepts model with varying cints is requested (by indvarying='cint')}

\item{coresToUse}{if neg., the value is subtracted from available cores, else value = cores to use}

\item{CoTiMAStanctArgs}{parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function}

\item{customPar}{logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)}

\item{diff}{labels for diffusion effects. Have to be either of the character strings of the type "diff_eta1" or "diff_eta2_eta1" (= freely estimated) or values (e.g., 0 for effects to be excluded, which is usually not recommended)}

\item{digits}{number of digits used for rounding (in outputs)}

\item{doPar}{parallel and multiple fitting if single studies. A value > 1 will fit each study doPar times in parallel mode during which no output is generated (screen remains silent). Useful to obtain best fit.}

\item{drift}{labels for drift effects. Have to be either of the character strings of the type V1toV2 (= freely estimated) or values (e.g., 0 for effects to be excluded, which is usually not recommended)}

\item{experimental}{set TRUE to try new pairwise N function}

\item{finishsamples}{number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).}

\item{indVarying}{control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means}

\item{indVaryingT0}{(default = NULL). Automatically set to TRUE if not set to FALSE if indVarying ist set TRUE. indVaryingT0=TRUE fits the regular random intercept models.}

\item{iter}{number of interation (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation}

\item{lambda}{R-type matrix with pattern of fixed (=1) or free (any string) loadings.}

\item{loadSingleStudyModelFit}{load the fit of single study ctsem models}

\item{manifestMeans}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{manifestVars}{define the error variances of the manifests within a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.}

\item{n.latent}{number of latent variables of the model (hast to be specified)!}

\item{n.manifest}{number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).}

\item{nopriors}{Deprecated, but still working. If TRUE, any priors are disabled – sometimes desirable for optimization}

\item{optimize}{if set to FALSE, Stan's Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .}

\item{posLL}{logical. Allows (default = TRUE) of positive loglik (neg -2ll) values}

\item{primaryStudies}{list of primary study information created with \code{\link{ctmaPrep}}}

\item{priors}{if FALSE, any priors are disabled – sometimes desirable for optimization}

\item{sameInitialTimes}{Only important for raw data. If TRUE (default=FALSE), T0MEANS occurs for every subject at the same time, rather than just at the earliest observation.}

\item{saveRawData}{save (created pseudo) raw date. List: saveRawData$studyNumbers, $fileName, $row.names, col.names, $sep, $dec}

\item{saveSingleStudyModelFit}{save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)}

\item{scaleTI}{scale TI predictors}

\item{scaleTime}{scale time (interval) - sometimes desirable to improve fitting}

\item{silentOverwrite}{overwrite old files without asking}

\item{T0means}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{T0var}{(default = 'auto')}

\item{useSV}{if TRUE (default=FALSE) start values will be used if provided in the list of primary studies}

\item{verbose}{integer from 0 to 2. Higher values print more information during model fit - for debugging}
}
\value{
ctmaFit returns a list containing some arguments supplied, the fitted models, different elements summarizing the main results,
model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
coresToUse, n.latent, n.manifest, and primaryStudyList. The study count is returned as n.studies, the created matrix of loadings of
manifest on latent factors is returned as lambda, and a re-organized list of primary studies with some information ommited is returned as
studyList. The fitted models for each primary study are found in studyFitList, which is a large list with many elements (e.g., the ctsem
model specified by CoTiMA, the rstan model created by ctsem, the fitted rstan model etc.). Further results returned are emprawList
(containing the pseudo raw data created), statisticsList (comprising baisc stats such as average sample size, no. of measurement points,
etc.), a list with modelResults (i.e., DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef,
CINT=model_Cint_Coef), and the paramter names internally used. The summary list,  which is printed if the summary function is applied to the
returned object, comprises "estimates" (the aggregated effects), possible randomEffects (not yet fully working),  confidenceIntervals, the
minus2ll value and its n.parameters, and possible warning messages (message). Plot type is plot.type=c("drift") and model.type="stanct"
("omx" was deprecated).
}
\description{
Fits ctsem models to each primary study in the supplied list of primary studies prepared by \code{\link{ctmaPrep}}.
}
\examples{
# Fit a ctsem model to all three primary studies summarized in
# CoTiMAstudyList_3 and save the three fitted models
\dontrun{
CoTiMAInitFit_3 <- ctmaInit(primaryStudies=CoTiMAstudyList_3,
                            n.latent=2,
                            checkSingleStudyResults=FALSE,
                            activeDirectory="/Users/tmp/") # adapt!
summary(CoTiMAInitFit_3)
}

}
