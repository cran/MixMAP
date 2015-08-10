\name{mixmapTest}
\alias{mixmapTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implements the MixMAP algorithm using hypothesis testing framework.  
}
\description{
This function implements the MixMAP algorithm, which performs gene-level tests of association using data from a previous GWAS or data from a meta-analysis of several GWAS.  Conceptually, genes are detected as significant if the collection of p-values within a gene are determined to be collectively smaller than would be observed by chance.  
}
\usage{
mixmapTest(data.set, pval="pval", snp="snp", gene="gene",
  coord="coord", chr="chr", alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	  \item{data.set}{
A data.frame containing the input data.  Each observation in this data set is a SNP.  This file must contain, at least, the SNP name, p-value for a SNP, group name (usually a gene name) where SNP is located, SNP location coordinate, and chromosome number of the SNP.  
}

  \item{pval}{
A character string with the name of the variable containing the p-values in the data.set.  Default is "pval". The data contained in this variable must be numeric.  
}
  \item{snp}{
A character string with the name of the variable containing the SNP name in the data.set.  Default is "snp".
}
  \item{gene}{
A character string with the name of the variable containing the gene name in the data.set.  Default is "gene".
}
  \item{coord}{
A character string with the name of the variable containing the location coordinate of the SNP in the data.set.  Default is "coord".  The data contained in this variable must be numeric.   
}
  \item{chr}{
A character string with the name of the variable containing the chromosome number of the SNP location in the data.set.  Default is "chr".  The data contained in this variable must be numeric.  
}

  \item{alpha}{
A numeric scalar indicating the level of significance the user chooses to use for detection.  The default is 0.05.
}

}
\details{
The user must provide a file that includes SNP name, SNP p-value, and a group name (assumed to be gene, but the user can use any grouping they choose), among other inputs.  

Note about SNPs in genes: 
It is possible for individual SNPs to be located in two overlapping genes.  In this case, the user can choose, when creating the input file, to list that SNP in both genes, only one of the genes, or simply ignore SNPS that are not in a unique gene.  If the user chooses to list a SNP in more than one gene, the SNP must have multiple rows in the input file with a different gene in each row.  %%It is not recommended that a SNP be listed in multiple genes.  

If the user only has SNP names and p-values, gene name must be appended to the file.  The user can either use their own file to append gene name, base pair, and chromosome.  The R package biomaRt located in bioconductor is a good source for linking SNPs to genes.  

}
\value{
An object of class \code{MixMAP}.
%%  If it is a LIST, use
  \item{output}{A \code{data.frame} with a row for each gene containing gene symbol, posterior estimates for all gene level effects, variance used in intervals, upper bound of one sided interval, the number of SNPs in each gene, the chromsome of the SNP, the location coordinate of the SNP, the gene-level p-value, the Bonferroni adjusted gene-level p-value, and the q-value based on Banjamini-Hochberg fasle discovery rate.}
 \item{num.genes.detected }{A vector containing the number of SNPs detected and the total number of genes}
  \item{detected.genes }{A \code{data.frame} with a row for each gene containing gene symbol, posterior estimates for all gene level effects, variance used in intervals, upper bound of one sided interval, the number of SNPs in each gene, the chromsome of the SNP, the location coordinate of the SNP, the gene-level p-value, the adjusted gene-level p-value, the name of the SNP with the smallest p-value, the minimum p-value in the gene, and a 5 number summary of the p-values within each gene}
   \item{lmer.out }{A \code{mer} object containing all of the model output information, including parameter estimates, from the \code{lmer} function.}
%% ...
}
\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  (2013) ``Mixed Modeling of Meta-Analysis P-Values (MixMAP) Suggests Multiple Novel Gene Loci for Low Density Lipoprotein Cholesterol". PLoS ONE 8(2): e54812.

Matthews, G.J. and Foulkes, A.S. (2015) ``{MixMAP}: An {R} Package for Mixed Modeling of Meta-Analysis $p$ Values in Genetic Association Studies." Journal of Statistical Software.  66(3): 1-11.
}
\author{
Gregory J. Matthews
}

\examples{
library(MixMAP)	
#Load data
#This data has been prepared to be used as input to the MixMAP function
data(MixMAP_example)
#Run MixMAP
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
