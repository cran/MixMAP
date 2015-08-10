MixMAP
======

This package contains functions to implement the MixMAP algorithm, which performs gene level tests of association using data from a previous GWAS or data from a meta-analysis of several GWAS.  Conceptually, genes are detected as significant if the collection of p-values within a gene are determined to be collectively smaller than would be observed by chance.

The function mixmapPI implements the version of MixMAP described in Foulkes et al. 2013.  However, the function mixmapTest implements an updated version of MixMAP that is less sensitive to the number of SNPs within a gene.  Therefore, we suggest using mixmapTest rather than mixmapPI.  
