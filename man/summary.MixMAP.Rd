\name{Summary}
\alias{summary,MixMAP-method}
\docType{methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary
}
\description{
Displays a summary of the MixMAP object.
}
\usage{
\S4method{summary}{MixMAP}(mixmap.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mixmap.object}{
An object of class MixMAP, which is output from the function MixMAP.  
}
 }
\details{
This method returns up to the top ten detected genes as well as information on the total number of genes and the number of genes detected.
}
\value{
Returns a summary of the MixMAP object.
}
\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  (2013) ``Mixed Modeling of Meta-Analysis P-Values (MixMAP) Suggests Multiple Novel Gene Loci for Low Density Lipoprotein Cholesterol". PLoS ONE 8(2): e54812. 
}
\author{
Gregory J. Matthews
}

\examples{
data(MixMAP_example)
#Run mixmapTest
MixOut<-mixmapTest(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",
  chr="Chr",coord="Coordinate",gene="Gene")
summary(MixOut)

#Run mixmapPI
MixOutPI<-mixmapPI(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",
  chr="Chr",coord="Coordinate",gene="Gene")
summary(MixOutPI)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GWAS }
\keyword{ Genetics }% __ONLY ONE__ keyword per line
\keyword{ Mixed Models}% __ONLY ONE__ keyword per line
