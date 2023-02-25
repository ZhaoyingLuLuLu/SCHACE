\name{IC.SCHACE}
\alias{IC.SCHACE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Information Criteria
}
\description{
Calculates tuning parameters, lambda and degrees of freedom of B-splines, based on information criteria.
}
\usage{
IC.SCHACE(A, N, sample_x, sample_y, Bdf_set, IC)
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
  \item{Bdf_set}{
Degrees of freedom of B-splines
}
  \item{IC}{
Information Criteria
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{df}{
  Degrees of freedom of B-splines chosen by Information Criteria
}
\item{lambda}{
  Lambda chosen by Information Criteria
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
