\name{CV.SCHACE}
\alias{CV.SCHACE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation procedure
}
\description{
Runs k-fold cross-validation for tuning paramters lambda and degrees of freedom of B-splines.
}
\usage{
CV.SCHACE(A, N, sample_x, sample_y, folds, Bdf_set, methods, clambda, percent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Design matrix X
}
  \item{N}{
Number of observations in sample_y
}
  \item{sample_x}{
Ordered design points, automatically generated in main.SCHACE()
}
  \item{sample_y}{
Reponse variable
}
  \item{folds}{
Number of folds
}
  \item{Bdf_set}{
Degrees of freedom of B-splines
}
  \item{methods}{
Loss function used in CV procedure
}
  \item{clambda}{
lambda.min or lambda.1se
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
%%  If it is a LIST, use
\item{df}{
  Degrees of freedom of B-splines chosen by CV procedure
}
\item{lambda}{
  Lambda chosen by CV procedure
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
