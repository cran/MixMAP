\name{MixMAP_example}
\alias{MixMAP_example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A working example of the MixMAP algorithm.}
\description{
Data containing the results of a meta-analysis consisting of many GWAS whose phenotype of interest is LDL cholesterol.  This file contains all of the information necessary to run the MixMAP function directly without any pre-processing.  
}
\usage{
data(MixMAP_example)
}
%- maybe also 'usage' for other objects documented here.
\format{
A data frame with 31825 observations and the following five variables:

MarkerName: a character with SNP name

GC.Pvalue: a numeric with p-values for each SNP

Coordinate: integer containing the base pair location of the SNP within the chromosome

Chr: integer containing the chromsome number of the SNP

Gene: A character with gene symbol


}
\details{
This file contains all of the information needed to run the MixMAP function directly.   
}
\source{Nature 466, 707--713 (05 August 2010) 
The p-values in this file are from the paper Teslovich Et Al. (2009) ``Biological, clinical and population relevance of 95 loci for blood lipids" Nature 466: 707--713.  The phenotype of interest was low-density lipoprotein (LDL) cholesterol.  The SNPs chosen are a subset that are found in the IBC array and have a unique gene name associated with them.  
}

\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  (2013) ``Mixed Modeling of Meta-Analysis P-Values (MixMAP) Suggests Multiple Novel Gene Loci for Low Density Lipoprotein Cholesterol". PLoS ONE 8(2): e54812.

Matthews, G.J. and Foulkes, A.S. (2015) ``{MixMAP}: An {R} Package for Mixed Modeling of Meta-Analysis $p$ Values in Genetic Association Studies." Journal of Statistical Software.  66(3): 1-11.

}
%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(MixMAP)	
#Load data
#This data has been prepared to be used as input to the MixMAP function directly
data(MixMAP_example)
str(MixMAP_example)

#Run mixmapTest
MixOut<-mixmapTest(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",
  chr="Chr",coord="Coordinate",gene="Gene")
#Display first ten detected genes
summary(MixOut)
#MixManhattan Plot
plot(MixOut)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GWAS }
\keyword{ Mixed Models }% __ONLY ONE__ keyword per line
\keyword{ Genetics }% __ONLY ONE__ keyword per line


