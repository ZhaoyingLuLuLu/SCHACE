\name{main.SCHACE}
\alias{main.SCHACE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Curve Fitting
}
\description{
Fits the curve with lambda and degrees of freedom of B-splines (user-specified or tuned by CV procedure or information criteria) and provides the estimated locations of change points and the jump sizes.
}
\usage{
main.SCHACE(
  y,
  folds = 3,
  Bdf_set = seq(3, 15, 1),
  lambda = NULL,
  Bdf = NULL,
  tuning = c("IC", "crossvalidation"),
  IC = c("BIC", "BIC.Chen"),
  methods = c("MSE", "Abs", "Trim"),
  clambda = "lambdamin",
  percent = 0.1
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Reponse variable
}
  \item{folds}{
Number of folds
}
  \item{Bdf_set}{
Degrees of freedom of B-splines
}
  \item{tuning}{
Tuning procedures
}
  \item{IC}{
Information Criteria
}
  \item{methods}{
Loss function used in CV procedure
}
  \item{clambda}{
lambda.min or lambda.1se for CV procedure, where lambda sequence is chosen by R package glmnet.
}
  \item{percent}{
The proportion of the largest squared differences of the observed and predicted values of observations that are trimmed when methods = "Trim".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
  \item{df}{
Degrees of freedom of B-splines chosen by tuning procedure or user-supplied
}
  \item{lambda}{
Lambda chosen by tuning procedure or user-supplied
}
  \item{predicted y}{
Predicted y
}
  \item{beta}{
Predicted beta
}
  \item{number of CPs}{
Estimated number of change points
}
  \item{locations of detected CPs}{
Estimated locations of change points
}
  \item{jump size}{
Estimated jump sizes
}
}
\references{
Zhaoying Lu, Ning Hao, and Hao Helen Zhang, Simultaneous Change-point Detection and Curve Estimation.
}
\author{
Zhaoying Lu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(splines)
library(glmnet)
data_blip <- function(n = 256, range = 1, rsnr = 4){

  sample_x <- seq(0, range, length.out = n)
  true_y <- (0.32 + 0.6*sample_x + 0.3*exp(-100*(sample_x - 0.3)^2))*as.numeric(identity(sample_x <= 0.8)) + (-0.28 + 0.6*sample_x + 0.3*exp(-100*(sample_x - 1.3)^2))*as.numeric(identity(sample_x > 0.8))
  sdf <- sd(true_y)
  sigma <- sdf/rsnr

  epsilon <- rnorm(n, 0, sigma)

  sample_y <- true_y + epsilon

  position = which(as.numeric(identity(sample_x <= 0.8) == 0) > 0)[1] - 1

  true_CP <- sample_x[position]

  data = list("true_y" = true_y, "sample_x" = sample_x, "sample_y" = sample_y, "true CPs" = true_CP, "noise" = sigma)
  return(data)

}

data_sample <- data_blip(rsnr = 8)
y <- data_sample$"sample_y"

main.SCHACE(
  y,
  folds = 3,
  Bdf_set = seq(3, 15, 1),
  lambda = NULL,
  Bdf = NULL,
  tuning = "IC",
  IC = "BIC.Chen",
  methods = c("MSE", "Abs", "Trim"),
  clambda = "lambdamin",
  percent = 0.1
  )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
