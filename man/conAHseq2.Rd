\name{conAHseq2}
\alias{conAHseq2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conAHseq2
}
\description{
A faster version of conAHseq, which considers less cases but gives the causal locations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
conAHseq2(Z,ld,n,k=6,alpha=0.05)
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
  \item{k}{
    a maximum number of causal SNPs.
  }
  \item{alpha}{
    the significance threshold.
  }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
conAHseq2 returns the significant causal SNPs, the p-value of each step and the predicted number of causal SNPs.
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
