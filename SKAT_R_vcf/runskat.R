args<-commandArgs(TRUE)
myskat <- function(pheno_file,gene_matrix) {
  library(SKAT)
  set.seed(1)
  d<-read.table(pheno_file,header=TRUE,sep=" ")
  pheno = d$PLq
  Ztmp<-read.table(gene_matrix,header=TRUE,sep="\t")
  Z<-data.matrix(Ztmp,rownames.force = NA)
  y.c<-d$PLq
  
  #obj<-SKAT_Null_Model(y.c ~ 1, out_type="D")
  #print (SKAT(Z, obj,kernel = "linear.weighted")$p.value)
  
  X<-matrix(c(d$age,d$sex),nrow=123,ncol=2,byrow=FALSE)
  obj<-SKAT_Null_Model(y.c ~ X, out_type="D",n.Resampling=0 , type.Resampling="bootstrap", Adjustment=TRUE)
  skt = SKAT(Z, obj,kernel = "linear.weighted")
  #print (skt$p.value)
  return(skt$p.value)
}
#pheno_file<-"/data/Rstudio/data/pheno_all1.ped"
#gene_matrix<-"/data/Rstudio/data/VDR.csv"
pheno_file<-args[1]
gene_matrix<-args[2]
p<-myskat(pheno_file,gene_matrix)
gene_name<-args[3]
gene<-data.frame(g=gene_name, pp = p)
write.table(gene, file = args[4], sep= "\t",row.names = FALSE,col.names =FALSE, append = TRUE, quote=FALSE)
print(p)




