methyAge<-function(beta,fastImputation=FALSE,normalize=TRUE,nCores=2)
{
#source("mage_norm_function.R")
#require(sqldf)
#require(impute)
#require(flashClust)
#require(dynamicTreeCut)
#require(RPMM)

if(min(beta,na.rm=TRUE)<0 | max(beta,na.rm=TRUE)>1){
stop("Warning: Methylation beta value input should be within [0,1]")}

horvath=data.frame()
hannum=data.frame()
phenoage=data.frame()

load(system.file("mage_ref.RData",package="ENmix"))
names(refmeth)[which(names(refmeth)=="meth_mean")]="goldstandard"
#load(system.file("mPOA_Models.RData",package="ENmix"))
PACE=DunedinPACE(betas=beta, proportionOfProbesRequired = 0.8)

missing_horvath=as.vector(horvath$cg[-1])[!(as.vector(horvath$cg[-1]) %in% rownames(beta))]
missing_hannum=as.vector(hannum$cg)[!(as.vector(hannum$cg) %in% rownames(beta))]
missing_phenoage=as.vector(phenoage$cg[-1])[!(as.vector(phenoage$cg[-1]) %in% rownames(beta))]
missing=FALSE
if(length(missing_horvath)>0){missing=TRUE
cat(length(missing_horvath), "CpG missed methy age calculation using Horvath method: ",missing_horvath)}
if(length(missing_hannum)>0){missing=TRUE
cat(length(missing_hannum), "CpG missed methy age calculation using Hannum method: ",missing_hannum)}
if(length(missing_phenoage)>0){missing=TRUE
cat(length(missing_phenoage), "CpG missed methy age calculation using phenoage method: ",missing_phenoage)}

tmp=NULL
if(missing){
cat("Warning: Methylation values for missing CpGs will be imputed with reference values, but the results will be less acurate")
clockcg=unique(c(as.vector(horvath$cg)[-1],as.vector(hannum$cg),as.vector(phenoage$cg)[-1]))
missedcg=clockcg[!(clockcg %in%  rownames(beta))]
if(length(missedcg)>0){
  reference=refmeth[as.vector(refmeth$cg) %in% missedcg,]
  tmp=matrix(rep(reference$goldstandard,ncol(beta)),nrow=length(missedcg))
  rownames(tmp)=as.vector(reference$cg)
  colnames(tmp)=colnames(beta)
# beta=rbind(beta,tmp)
tmp=t(tmp)
}
}


refmeth=refmeth[as.vector(refmeth$cg) %in% rownames(beta),]
beta=beta[as.vector(refmeth$cg),]

#STEP 2: Imputing
set.seed(1)
nSamples=ncol(beta)
beta= t(beta)
noMissingPerSample=apply(as.matrix(is.na(beta)),1,sum)
if(!fastImputation & nSamples>1){
if(max(noMissingPerSample,na.rm=TRUE)>0 ){
dimnames1=dimnames(beta)
beta= data.frame(t(impute::impute.knn(t(beta))$data))
dimnames(beta)=dimnames1
}}
if(fastImputation | nSamples==1){
if(max(noMissingPerSample,na.rm=TRUE)>0){
dimnames1=dimnames(beta)
for (i in which(noMissingPerSample>0)){
selectMissing1=is.na(beta[i,])
beta[i,selectMissing1] = as.numeric(refmeth$goldstandard[selectMissing1])
}
dimnames(beta)=dimnames1
}}

#normalization using Horvath modified BMIQ method
if(normalize){
nCores=min(floor(nrow(beta)/2),nCores)
if(nCores>detectCores()){nCores <- detectCores()}
N=ceiling(nrow(beta)/(nCores))
parts=rep(1:nCores,each = N)[1:nrow(beta)]
bbb<-function(s,parts,dat){
  id=which(parts==s)
  dat1=dat[id,]
  BMIQcalibration(datM=dat1,goldstandard.beta= refmeth$goldstandard,plots=FALSE)
}
resu <- mclapply(1:nCores,bbb,parts=parts,dat=beta, mc.cores=nCores)
beta=as.matrix(do.call(rbind, resu))
}
beta=cbind(beta,tmp)

methyAge=data.frame(SampleID=rownames(beta))
datMethClock=beta[,as.character(horvath$cg[-1])]
#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
mAge=anti.trafo(horvath$coef[1]+as.matrix(datMethClock)%*% as.numeric(horvath$coef[-1]))
methyAge$mAge_Horvath=mAge

datMethClock=beta[,as.character(hannum$cg)]
mAge=as.matrix(datMethClock) %*% as.numeric(hannum$coef)
methyAge$mAge_Hannum=mAge

datMethClock=beta[,as.character(phenoage$cg[-1])]
mAge=phenoage$coef[1]+as.matrix(datMethClock)%*% phenoage$coef[-1]
methyAge$PhenoAge=mAge

methyAge$PACE=PACE[as.vector(methyAge$SampleID)]

return(methyAge)
}


