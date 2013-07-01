\name{MixMAP-package}
\alias{MixMAP-package}
\alias{lme4-package}
\docType{package}
\title{
MixMAP: Mixed Modeling of Meta-Analysis P-values
}
\description{
Mixed Modeling approach based on Meta-analysis p-values for detecting gene level associations from genome wide assocation studies (GWAS) or candidate gene studies. 

This package uses raw p-values from previous GWAS, and information about gene-level groupings to search for gene-level assocations between complex disease phenotypes and genetic loci. 
}
\details{
\tabular{ll}{
Package: \tab MixMAP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-08-13\cr
License: \tab GPL-3\cr
}

This goal of this package is to implement the MixMAP algorithm.  The aim of the algorithm is to search for assocations between genes and complex diseases by using individual SNP-level p-values.  The user must provides a file with SNP name and SNP p-value as well as gene name, chromosome, and basepair location.  The output of the function MixMAP is an object of class 'MixMAP', which contains information on genes that have been detected as being associated with the phenotype of interest.  An object of class 'MixMAP' has a plot method associated with it to visually display the result of the MixMAP algorithm in a Manhattan style plot as well as a summary method.  
}
\author{
Gregory J. Matthews and Andrea Foulkes

Maintainer: Gregory J. Matthews <gjm112@gmail.com>

}
\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  (2013) ``Mixed Modeling of Meta-Analysis P-Values (MixMAP) Suggests Multiple Novel Gene Loci for Low Density Lipoprotein Cholesterol". PLoS ONE 8(2): e54812.
}

\keyword{ package }
\seealso{
\code{\link{lme4-package}}
}
\examples{
library(MixMAP)  
#Load data
#This data has been prepared to be used as input to the MixMAP function directly
data(MixMAP_example)
str(MixMAP_example)
#Run MixMAP
MixOut<-MixMAP(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",
  chr="Chr",coord="Coordinate",gene="Gene")
#Display first ten detected genes
summary(MixOut)
#MixManhattan Plot
plot(MixOut)
}
