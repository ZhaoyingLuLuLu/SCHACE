\name{trimmedmse}
\alias{trimmedmse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trimmed Mean Squared Errors
}
\description{
Calculates the trimmed Mean Squared Error Loss in cross-validation procedure.
}
\usage{
trimmedmse(test.y, pred.y, percent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.y}{
Actual values of observations in test data set.
}
\item{pred.y}{
Predicted values of observations in test data set.
}
\item{percent}{
The proportion of the largest squared differences of the observed and predicted values of observations that are trimmed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Trimmed Mean Squared Errors
}
\references{
Andreas Alfons (2011). cvTools: Cross-validation tools for regression models. R package version 0.3.2, https://cran.r-project.org/web/packages/cvTools
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[cvTools]{tmspe}
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
