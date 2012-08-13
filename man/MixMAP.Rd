\name{MixMAP}
\alias{MixMAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MixMAP
MixMAP
}
\description{
This function implements the MixMAP algorithm, which performs gene-level tests of association using data from a previous GWAS or data from a meta-analysis of several GWAS.  Conceptually, genes are detected as significant if the collection of p-values within a gene are determined to be collectively smaller than would be observed by chance.  
}
\usage{
MixMAP(data.set, pval="PVAL", snp="SNP", gene="GENE", bp="BP", chr="CHR", alpha = 0.05, use.post.var=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	  \item{data.set}{
A data.frame containing the input data.  Each observation in this data set is a SNP.  This file must contain, at least, the SNP name, p-value for a SNP, group name (usually a gene name) where SNP is located, SNP base pair location, and chromosome number of the SNP.  
}

  \item{pval}{
A character string with the name of the variable containing the p-values in the data.set.  Default is "PVAL".
}
  \item{snp}{
A character string with the name of the variable containing the SNP name in the data.set.  Default is "SNP".
}
  \item{gene}{
A character string with the name of the variable containing the gene name in the data.set.  Default is "GENE".
}
  \item{bp}{
A character string with the name of the variable containing the base pair location of the SNP in the data.set.  Default is "BP".
}
  \item{chr}{
A character string with the name of the variable containing the chromosome number of the SNP location in the data.set.  Default is "CHR".
}

  \item{alpha}{
A numeric scalar indicating the level of significance the user chooses to use for detection.  The default is 0.05.
}
  \item{use.post.var}{
	A logical value.  TRUE will use the posterior estimate of variance for prediction intervals (this is the default in the package lme4).  FALSE will use the variance of the estimated BLUP minus the random effect.  The default is FALSE.
}
}
\details{
The user must provide a file that includes SNP name, SNP p-value, and a group name (assumed to be gene, but the user can use any grouping they choose), among other inputs.  

Note about SNPs in genes: 
It is possible for individual SNPs to be located in two overlapping genes.  In this case, the user can choose, when creating the input file, to list that SNP in both genes, only one of the genes, or simply ignore SNPS that are not in a unique gene.  If the user chooses to list a SNP in more than one gene, the SNP must have multiple rows in the input file with a different gene in each row.  %%It is not recommended that a SNP be listed in multiple genes.  

If the user only has SNP names and p-values, gene name must be appended to the file.  The user can either use their own file to append gene name, base pair, and chromosome.  The R package biomaRt located in bioconductor is a good source for linking SNPs to genes.  

%%A link between SNPs and genes can be found by using ensembl in the R package biomaRt located at bioconductor.  

Note about variance for prediction intervals:
By default in lme4, the lmer function returns the posterior variance estimate for the best linear unbiased estimators of the random effects.  The user can choose to use this variance by chooisng the option use.post.var=TRUE.  By default in MixMAP, use.post.var=FALSE, which will a slightly different variance calculation in the prediction interval that is always larger than the posterior variance.  Practically though, when the number of groups (genes, in this setting) is large, these two variances are very close to each other.  For more details see: Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  ``Mixed modeling of Meta-Analysis P-values (MixMAP) suggests multiple novel gene loci for low density lipoprotein cholesterol." %%PLOS one.  

}
\value{
An object of class 'MixMAP'.
%%  If it is a LIST, use
  \item{output}{data.frame with a row for each gene containing gene symbol, posterior estimates for all gene level effects, variance used in prediction intervals, upper bound of one sided prediction interval, and number of SNPs in each gene.}
 \item{num.genes.detected }{A vector containing the number of SNPs detected and the total number of genes}
  \item{detected.genes }{A data.frame with a row for each detected gene containing gene symbol, posterior estimates, posterior variances of estimates, upper bound of one-sided prediction interval, number of SNPs in gene, min p-value in gene, and a 5 number summary of the p-values within each gene}
   \item{lmer.out }{A mer object containing all of the model output information, including parameter estimates, from the lmer function.}
%% ...
}
\references{
Foulkes, A.S., Matthews, G.J., Das, U., Ferguson, J., Reilly, M.  ``Mixed modeling of Meta-Analysis P-values (MixMAP) suggests multiple novel gene loci for low density lipoprotein cholesterol." Under Review; Available upon request. %%PLOS one.  

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
