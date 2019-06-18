args<-commandArgs(TRUE)
# input
# args[1]: pheno_file
# args[2]: gene_matrix
# args[3]: gene_name
# args[4]: file output
myskat <- function(pheno_file, gene_matrix) {
  library(SKAT)
  set.seed(1)
  d<-read.table(pheno_file, header=TRUE, sep=" ")
  # phenotype
  pheno<-d$PLq
  # import covariates to X
  X<-matrix(c(d$age, d$sex), nrow=123, ncol=2, byrow=FALSE)
  
  # genotype matrix
  Ztmp<-read.table(gene_matrix, header=TRUE, sep="\t")
  Z<-data.matrix(Ztmp, rownames.force=NA)
  
  # run SKAT
  obj<-SKAT_Null_Model(pheno ~ X, out_type="D", n.Resampling=0, type.Resampling="bootstrap", Adjustment=TRUE)
  skt = SKAT(Z, obj, kernel="linear.weighted")
  return(skt$p.value)
}

pheno_file<-args[1]
gene_matrix<-args[2]
p<-myskat(pheno_file, gene_matrix)
gene_name<-args[3]
gene<-data.frame(g=gene_name, pp=p)
write.table(gene, file=args[4], sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
print(p)




