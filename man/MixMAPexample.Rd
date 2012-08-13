\name{MixMAP_example}
\alias{MixMAP_example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MixMAP_example}
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
\source{
The p-values in this file are from the 2010 paper Teslovich Et Al.  ``Biological, clinical and population relevance of 95 loci for blood lipids" in Nature where the 
 phenotype of interest was low-density lipoprotein (LDL) cholesterol.  The SNPs chosen are a subset that are found in the IBC array and have a unique gene name associated with them.  
}

\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  ``Mixed modeling of Meta-Analysis P-values (MixMAP) suggests multiple novel gene loci for low density lipoprotein cholesterol." Under Review; Available upon request. %%PLOS one.  

}
%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(MixMAP)	
#Load data
#This data has been prepared to be used as input to the MixMAP function directly
data(MixMAP_example)
str(MixMAP_example)
#Run MixMAP
MixOut<-MixMAP(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",chr="Chr",bp="Coordinate",gene="Gene",use.post.var=TRUE)
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


