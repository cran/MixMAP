MixMAP <-
function(data.set,pval="PVAL",snp="SNP",gene="GENE",bp="BP",chr="CHR",alpha=0.05,use.post.var=FALSE){
############################
#defining errors
############################
#names must be specified
#p-values
if (!pval%in%names(data.set)) stop(gettextf(paste0('Variable "',pval,'" not found in input data.frame.  Please specify variable name for p-values')))
#snp
if (!snp%in%names(data.set)) stop(gettextf(paste0('Variable "',snp,'" not found in input data.frame.  Please specify variable name for SNPs.')))
#basepair
if (!bp%in%names(data.set)) stop(gettextf(paste0('Variable "',bp,'" not found in input data.frame.  Please specify variable name for basepair location.')))
#chromosome
if (!chr%in%names(data.set)) stop(gettextf(paste0('Variable "',chr,'" not found in input data.frame.  Please specify variable name for chromosome')))
#gene
if (!gene%in%names(data.set)) stop(gettextf(paste0('Variable "',gene,'" not found in input data.frame.  Please specify variable name for genes.')))

#Lengths of input must be the same
if (pval%in%names(data.set) & length(data.set[[pval]])!=length(data.set[[gene]])) stop(gettextf(paste('Lengths differ: Length of pval is ',length(data.set[[pval]]),'; Length of gene is ',length(data.set[[gene]]),sep="")))

############################
#Warnings
############################
#are pvalues numeric?
if (!is.numeric(data.set[[pval]])) stop(gettextf('p-values must be numeric'))
if (sum(is.na(data.set[[pval]]))>0) stop(gettextf('Some p-values are missing'))
if (sum(is.na(data.set[[gene]]))>0) stop(gettextf('Some gene names are missing'))

############################
#Pull out the subset of data that will be used
############################
#Pull out the data that we need from the bigger data file
dat.temp<-data.frame(pval.temp=data.set[[pval]],gene.temp=as.character(data.set[[gene]]),SNP.temp=as.character(data.set[[snp]]))
fret<-dat.temp[dat.temp$gene.temp!="",]

#How many SNPs per gene?
tab<-data.frame(table(dat.temp$gene.temp))
names(tab)<-c(gene,"SNP.COUNT")

#Inverse normal transformation of the p-values after ranking
dat.temp$probit.rank.transform<-qnorm((rank(dat.temp$pval.temp)-0.5)/length(dat.temp$pval.temp))

#Run lmer function
fm.rawg=lmer(probit.rank.transform ~ 1+(1|gene.temp),data=dat.temp)
aa=ranef(fm.rawg,postVar=TRUE)
beta<-fixef(fm.rawg)
post.est=aa$gene.temp[,1,]
post.var=attr(aa$gene.temp,"postVar")[1,1,]
n.i<-as.vector(table(dat.temp$gene.temp))
sigma.sq.b<-VarCorr(fm.rawg)$gene.temp[1,1]
sigma.sq<-attr(VarCorr(fm.rawg),"sc")

JSJ<-(n.i/(sigma.sq*sigma.sq.b))/(n.i/sigma.sq+1/sigma.sq.b)
A<-(sum(JSJ)-JSJ)/sum(JSJ)
var.bi<-sigma.sq.b^2*JSJ*A
var.bihat.minus.bi<-sigma.sq.b-var.bi

############################################################
#Calculating upper limit of prediction interval depending on which variance to use.
if (use.post.var==FALSE)
	{
		pred.upper<-post.est+sqrt(var.bihat.minus.bi)*qnorm(1-alpha/(length(post.est)))##Bonferroni correction
		var.out<-var.bihat.minus.bi
	}  
if (use.post.var==TRUE)
	{
		pred.upper<-post.est+sqrt(post.var)*qnorm(1-alpha/(length(post.est)))  ##Bonferroni correction
		var.out<-post.var
	}
############################################################

############################
#Defining Output
############################
out<-data.frame(GENE=as.character(rownames(aa$gene)),post.est=post.est,var=var.out,pred.upper=as.numeric(as.character(pred.upper)))
names(out)[1]<-gene
out<-merge(out,tab,by.x=gene,by.y=gene,all.x=TRUE)

data.set.g<-data.set[!duplicated(data.set[[gene]]),c(gene,chr,bp)]
out<-merge(out,data.set.g,by.x=gene,by.y=gene,all.x=TRUE)
out[[bp]]<-as.numeric(as.character(out[[bp]]))
out[[chr]]<-as.numeric(as.character(out[[chr]]))

cutoff<-0
num<-c("number detected"=sum(out$pred.upper<cutoff),"total number of genes"=dim(out)[1])
detected<-out[out$pred.upper<cutoff,]

############################
#If any genes detected
############################
if (num[1]>0){
genes.detect<-as.character(detected[[gene]])
SNP.temp<-dat.temp[dat.temp$gene.temp%in%genes.detect,]
SNP.temp$gene.temp<-as.character(SNP.temp$gene.temp)
SNP.temp$SNP.temp<-as.character(SNP.temp$SNP.temp)

#Pull out the min SNP 
snp.tmp.list<-list()
for (g in unique(SNP.temp$gene.temp)){
tmp<-SNP.temp[SNP.temp$gene.temp==g,]
snp.tmp.list[[g]]<-unlist(c(tmp[min(tmp$pval.temp)==tmp$pval.temp,][1,-1],summary(tmp$pval.temp)))
}

snp.min<-data.frame(do.call(rbind,snp.tmp.list))
names(snp.min)<-c(gene,"min.SNP","probit.rank.transform","pval.min","pval.Q1","pval.median","pval.mean","pval.Q3","pval.max")
snp.min[[gene]]<-as.character(snp.min[[gene]])
snp.min$min.SNP<-as.character(snp.min$min.SNP)
snp.min$probit.rank.transform<-as.numeric(as.character(snp.min$probit.rank.transform))
snp.min$pval.min<-as.numeric(as.character(snp.min$pval.min))
snp.min$pval.Q1<-as.numeric(as.character(snp.min$pval.Q1))
snp.min$pval.median<-as.numeric(as.character(snp.min$pval.median))
snp.min$pval.mean<-as.numeric(as.character(snp.min$pval.mean))
snp.min$pval.Q3<-as.numeric(as.character(snp.min$pval.Q3))
snp.min$pval.max<-as.numeric(as.character(snp.min$pval.max))


#merge on the min pvalue and SNP name
detected<-merge(detected,snp.min,by.x=gene,by.y=gene,all.x=TRUE)

#merge in the CHR and BP
#CHR.BP.temp<-data.set[as.character(data.set[[snp]])%in%snp.min$min.SNP,c("GENE","CHR","BP")]
#detected<-merge(detected,CHR.BP.temp,by.x="GENE",by.y="GENE",all.x=TRUE)

names(out)[c(1,6,7)]<-names(detected)[c(1,6,7)]<-c("GENE","CHR","BP")
out[["GENE"]]<-as.character(out[["GENE"]])
detected[["GENE"]]<-as.character(detected[["GENE"]])

#detected.merg<-merge(detected,gene.location.file.num,by.x="gene",by.y="external_gene_id",all.x=TRUE)
MixMAP.out<-new("MixMAP",output=out,num.genes.detected=num,detected.genes=detected,lmer.out=fm.rawg)

}
############################
#If no genes are detected
############################
if (num[1]==0) {MixMAP.out<-new("MixMAP",output=out,num.genes.detected=num,detected.genes=detected,lmer.out=fm.rawg)}

#return MixMAP object
MixMAP.out
}
