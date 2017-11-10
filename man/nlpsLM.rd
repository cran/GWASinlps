\name{nlpsLM}

\alias{nlpsLM}
 
\title{Nonlocal prior based single-step SNP selection for Genome-wide association study data with continuous phenotype}

\description{Performs variable selection in a single iteration with data from Genome-wide association study (GWAS) with continuous phenotype using non-local prior based Bayesian general linear model as described in Sanyal et al. (2017)(submitted).}

\usage{
nlpsLM(x, y, prior, tau, priorDelta = modelbbprior(1,1), k0, rxx, 
       niter = 2000, verbose = FALSE)
}

\arguments{
  \item{x}{The SNP genotype matrix with subjects represented by rows and SNPs represented by columns. For continuous y (see below), the elements of x should be of numeric type. Missing values are currently not accepted.}
  \item{y}{The vector of phenotype values. Can be continuous or binary (0 and 1) phenotype.}
  \item{prior}{"mom" for pMOM prior, or "imom" for piMOM prior.}
  \item{tau}{the value of the scale parameter tau of the nonlocal prior. }
  \item{priorDelta}{Prior for model space. Defaults to modelbbprior(1,1).}
  \item{k0}{GWASinlps tuning parameter, denoting the number of leading SNPs (see References).}
  \item{rxx}{GWASinlps tuning parameter, denoting the correlation threshold to determine leading sets (see References).}
  \item{niter}{Number of MCMC iterations for nonlocal prior based Bayesian variable selection. Defaults to 2000.}
  \item{verbose}{If TRUE, prints some details. FALSE by default.}
}

\details{The nlpsLM function performs SNP selection in one iteraion for continuous phenotypes. The GWASinlps function repeatedly calls the nlpsLM function. The nlpsLM procedure starts by determining the k0 \emph{leading SNPs} having the highest association with y. The measure of association is Pearson's correlation coefficient. These k0 leading SNPs, in turn, determine the k0 \emph{leading sets}, where each leading set consists of all SNPs with absolute correlation coefficient more than or equal to rxx with the correspondng leading SNP. Then non-local prior based Bayesian variable selection is run within each leading set, and the predictors appearing in the HPPM are considered selected. For more details, see the References.
}

\value{
A list with elements
  \item{hppm }{The set of predictors appearing in the HPPM of at least one leading set. }
  \item{not.selected }{The set of predictors appearing in at least one leading set, but in none of the HPPMs.}
}

\references{
Sanyal et al. (2017), "GWASinlps: Nonlocal prior based SNP selection tool for genome-wide association studies" (submitted to Bioinformatics).
}

\author{Nilotpal Sanyal <nisanyal@ucsd.edu>}
% \note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GWASinlps}}
}

\examples{
%\dontrun{
p = 1000; n = 100
m = 10
# Generate SNP genotype matrix
set.seed(1) 
f = runif( p, .1, .2 ) # simulate minor allele frequency
x = matrix( nrow = n, ncol = p )
colnames(x) = 1:p
for(j in 1:p)
  x[,j] = rbinom( n, 2, f[j] )
# Generate data
causal_snps = sample( 1:p, m )
beta = rep( 0, p )
set.seed(1)
beta[causal_snps] = rnorm(m, mean = 0, sd = 2 )
y = x \%*\% beta + rnorm(n, 0, 1) 
# Fix scale parameter tau 
tau = 0.022
# Perform GWASinlps
out = nlpsLM(x, y, prior = "mom", tau = tau, k0 = 10, rxx = .5, niter = 10000, verbose = TRUE) 
out
%}
}