\name{MixMAP-class}
\alias{MixMAP-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MixMAP-class
}
\description{
Objects of class \code{MixMAP} are returned as output from the function \code{MixMAP}.  
}
\arguments{
%- maybe also 'usage' for other objects documented here.

  \item{output}{
A \code{data.frame} containing posterior estimates for all genes, as well as other associated gene level information.  
}
  \item{num.genes.detected}{
A numeric vector containing the number of genes detected and the total number of genes.    
}
  \item{detected.genes}{
A \code{data.frame} containing information on the detected genes.  This includes the posterior estimate, variance, and upper bound of the prediction interval fo all genes.  Other gene level information that is included is the gene location including chromosome and base pair along with the total number of SNPs within each gene (snpCount), and some summary information regarding the SNP p-values within each gene.  
}

  \item{num.genes.detected}{
A numeric vector containing the number of genes detected and the total number of genes.    
}

   \item{lmer.out}{A \code{mer} object containing all of the model output information, including parameter estimates, from the \code{lmer} function.}

	}
	
	
	

	

\author{
Gregory J. Matthews
}

\examples{
showClass("MixMAP")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GWAS }
\keyword{ Mixed Models }% __ONLY ONE__ keyword per line
\keyword{ Genetics }% __ONLY ONE__ keyword per line