% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.gllvm.R
\name{predict.gllvm}
\alias{predict.gllvm}
\alias{predict}
\title{Predict Method for gllvm Fits}
\usage{
\method{predict}{gllvm}(
  object,
  newX = NULL,
  newTR = NULL,
  newLV = NULL,
  type = "link",
  level = 1,
  offset = TRUE,
  ...
)
}
\arguments{
\item{object}{an object of class 'gllvm'.}

\item{newX}{A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.}

\item{newTR}{A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.}

\item{newLV}{A new matrix of latent variables.  If omitted, the original matrix of latent variables is used. Note that number of rows/sites must be the same for \code{newX} (if X covariates are included in the model).}

\item{type}{the type of prediction required. The default (\code{"link"}) is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable. that is, the predictions for the binomial model are predicted probabilities. In case of ordinal data, \code{type = "response"} gives predicted probabilities for each level of ordinal variable.}

\item{level}{specification for how to predict. Level one (\code{level = 1}) attempts to use the predicted site scores from variational approximations or laplace approximation or given site scores in \code{newLV}. Level 0 sets the latent variable to zero. Defaults to 1.}

\item{offset}{specification whether of not offset values are included to the predictions in case they are in the model, defaults to \code{TRUE} when offset values that are used to fit the model are included to the predictions. Alternatives are matrix/vector (number of rows must match with the \code{newX}) of new offset values or \code{FALSE}, when offsets are ignored.}

\item{...}{not used.}
}
\value{
A matrix containing requested predictor types.
}
\description{
Obtains predictions from a fitted generalized linear latent variable model object.
}
\details{
If \code{newX}, \code{newTR} and \code{newLV} are omitted the predictions are based on the data used for fitting the model. Notice that \code{newTR} need to match with the number of species in the original data.
Instead, new sites can be specified in \code{newX}. If predictors \code{newX} (and \code{newTR}) are given, and \code{newLV} is not, latent variables are not used in the predictions.
}
\examples{
\donttest{
# Load a dataset from the mvabund package
data(antTraits, package = "mvabund")
y <- as.matrix(antTraits$abund)
X <- scale(antTraits$env[, 1:3])
# Fit gllvm model
fit <- gllvm(y = y, X, family = poisson())
# fitted values
predfit <- predict(fit, type = "response")

# linear predictors
predlin <- predict(fit)
# Predict new sites:
# Generate matrix of environmental variables for 10 new sites
xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
colnames(xnew) <- colnames(X)
predfit <- predict(fit, newX = xnew, type = "response", level = 0)

TR <- (antTraits$tr[, 1:3])
fitt <- gllvm(y = y, X, TR, family = poisson())
# linear predictors
predlin <- predict(fitt)
# Predict new sites:
# Generate matrix of environmental variables for 10 new sites
xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
colnames(xnew) <- colnames(X)
# Generate matrix of traits for species
trnew <- data.frame(Femur.length = rnorm(41), No.spines = rnorm(41),
 Pilosity = factor(sample(0:3, 41, replace = TRUE)))
predfit <- predict(fitt, newX = xnew, newTR = trnew, type = "response", level = 0)
}
}
\author{
Jenni Niku <jenni.m.e.niku@jyu.fi>,  David Warton
}
