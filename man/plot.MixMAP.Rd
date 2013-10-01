\name{Plot}
\alias{plot,MixMAP-method}
\docType{methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot
}
\description{
The purpose of plot is to graphically display the results of MixMAP in a similar way as a Manhattan plot that is used for Genome Wide Association Studies.    
}
\usage{
\S4method{plot}{MixMAP}(x,col.genes=c("black","gray"),col.detected=c("blue","violet"),
  col.text="black",title="MixMAP Manhattan Plot",display.text=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class mixmap, which is output from the function MixMAP.  
}
  \item{col.genes}{
A character vector containing colors for alternate chromosomes.  Default is black and gray.
}
  \item{col.detected}{
A character vector containing with the names of color to display genes detected by MixMAP.  Genes detected by MixMAP but not by single SNP analysis will be dispalyed using the first color.  Genes deteced in both MixMAP and single SNP analysis are displayed using the second colors.  Default is blue and violet.
}
  \item{col.text}{
A character string with the name of a color to display text of the names of detected genes.  Default is black.
}
  \item{title}{
A character string used a title for the MixManhattan plot.  Default is "MixMAP Manhattan Plot".
}
  \item{display.text}{
Either TRUE or FALSE.  FALSE will not label detected genes; TRUE will add text.  Default is TRUE.  
}
}
\details{
Since MixMAP searches for genes with small empirical Bayes estimates, the resulting graph plots all genes with a positive empirical Bayes estimate as zero, and genes with a negative Bayes estimate are displayed as their absolute values.  
}
\value{
Returns a MixManhattan plot.
}
\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  (2013) ``Mixed Modeling of Meta-Analysis P-Values (MixMAP) Suggests Multiple Novel Gene Loci for Low Density Lipoprotein Cholesterol". PLoS ONE 8(2): e54812.   
}
\author{
Gregory J. Matthews
}

\examples{
data(MixMAP_example)
#Run MixMAP
MixOut<-MixMAP(MixMAP_example,pval="GC.Pvalue",snp="MarkerName",
  chr="Chr",coord="Coordinate",gene="Gene")
plot(MixOut)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GWAS }
\keyword{ Genetics }% __ONLY ONE__ keyword per line
\keyword{ Mixed Models}% __ONLY ONE__ keyword per line
