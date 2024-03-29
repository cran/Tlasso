% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal.R
\name{signal}
\alias{signal}
\title{Regression Parameter of Conditional Linear Model}
\usage{
signal(Omega.list, i = 1, k = 1)
}
\arguments{
\item{Omega.list}{list of precision matrices of tensor, i.e., \code{Omega.list[[k]]} is the kth precision matrix.
Omega.list can be either true precision matrices or output of \code{Tlasso.fit}.   
for the kth tensor mode, \eqn{k \in \{1 , \ldots, K\}}{ 1 <= k <= K }.}

\item{i}{index of interested regression parameter, default is 1. See details in \href{https://arxiv.org/abs/1609.04522}{Lyu et al. (2019)}.}

\item{k}{index of interested mode, default is 1.}
}
\value{
A vector of regression paramter.
}
\description{
Compute regression parameter of conditional linear model of separable tensor normal distribution described in \href{https://arxiv.org/abs/1609.04522}{Lyu et al. (2019)}.
}
\details{
This function computes regression parameter and is fundamental for sample covariance of residuals and 
bias correction. See details in \href{https://arxiv.org/abs/1609.04522}{Lyu et al. (2019)}.
}
\examples{

m.vec = c(5,5,5)  # dimensionality of a tensor 
n = 5   # sample size 
k=1 # index of interested mode
lambda.thm = 20*c( sqrt(log(m.vec[1])/(n*prod(m.vec))), 
                   sqrt(log(m.vec[2])/(n*prod(m.vec))), 
                   sqrt(log(m.vec[3])/(n*prod(m.vec))))
DATA=Trnorm(n,m.vec,type='Chain') 
# obersavations from tensor normal distribution
out.tlasso = Tlasso.fit(DATA,T=1,lambda.vec = lambda.thm)   
# output is a list of estimation of precision matrices
signal(out.tlasso, i=2 , k=k )
# the regression parameter for conditional linear model of 2rd row in 1st mode

}
\seealso{
\code{\link{covres}}, \code{\link{biascor}}
}
\author{
Xiang Lyu, Will Wei Sun, Zhaoran Wang, Han Liu, Jian Yang, Guang Cheng.
}
