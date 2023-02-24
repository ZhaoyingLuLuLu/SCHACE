\name{foldK}
\alias{foldK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
foldK
}
\description{
Given the sample size (n) and the number of groups (K), divides the index set \{1, ..., n\} into K subsets based on the remainder of the division of n by K.
}
\usage{
foldK(n,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations
}
\item{K}{
number of groups
}
}
\details{

}
\value{
It returns order-preserved data splitting for K-fold cross-validation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zhaoying Lu and Ning Hao
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
foldK(50,3)
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
