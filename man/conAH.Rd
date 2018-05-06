\name{conAH}
\alias{conAH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conAH
}
\description{
Test Allelic Heterogeneity with summary statistics and intersection-union tests.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
conAH(Z,ld,n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
    a q*1 matrix containing the Z-scores for q SNPs.
  }
  \item{ld}{
    a q*q correlation matrix of the q SNPs.
  }
  \item{n}{
    a q*1 matrix containing the sample sizes used to get the Z-scores for the q SNPs. Can be a number if the q sample sizes are the same.
  }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
conAH returns the p-value of testing AH.
}
\references{
%Schaid, D. J., Tong, X., Larrabee, B., Kennedy, R. B., Poland, G. A., & Sinnwell, J.P. (2016). Statistical Methods for Testing Genetic Pleiotropy. Genetics 204(2): 483-497.
%% ~put references to the literature/web site here ~
}
\author{
Yangqing Deng and Wei Pan.
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
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
