
methscore<-function(datMeth,datPheno=NULL,fastImputation=FALSE,normalize=TRUE,
		    GrimAgeComponent=3,AgeResid=FALSE)
{
#input datMeth value matrix must be in between 0 to 1, percentage of cg methylated
datMeth<-as.matrix(datMeth)
if(min(datMeth,na.rm=TRUE)<0 | max(datMeth,na.rm=TRUE)>1){
    stop("Warning: Methylation datMeth value should be within [0,1]")}

#if GrimAge Components flag
if (!GrimAgeComponent %in% c(1, 2, 3)) {
  stop("Error: GrimAgeComponent must be 1, 2, or 3")
}

#remove CpGs with missing values for all samples
datMeth=datMeth[rowSums(is.na(datMeth)) != ncol(datMeth),]

#check pheno 
if (is.null(datPheno)) {stop("Error: datPheno is missing.")}
datPheno <- as.data.frame(datPheno)
required_cols <- c("Age", "Female", "SampleID")
missing_cols <- setdiff(required_cols, colnames(datPheno))
if (length(missing_cols) > 0) {
  stop(paste0("Error: datPheno must have the following columns: ",
    paste(missing_cols, collapse = ", ")))
}
if (!all(colnames(datMeth) %in% datPheno$SampleID)) {
  stop("Error: some samples in datMeth do not have phenotype in datPheno")
}
rownames(datPheno) <- datPheno$SampleID
datPheno <- datPheno[colnames(datMeth), , drop = FALSE]
if (any(is.na(datPheno))) {
  stop("Error: missing value in phenotype data is not allowed")
}

#load in all model and reference data
load(system.file("mage_ref.RData",package="ENmix"))
load(system.file("mPOA_Models.RData",package="ENmix"))
#load(system.file("pcclock_model.RData",package="ENmix"))
#pcc1=pcc2=pcc3=pcc4=NULL
#load in PC related clocks
load(system.file("pcclock_model1.RData",package="ENmix"))
load(system.file("pcclock_model2.RData",package="ENmix"))
load(system.file("pcclock_model3.RData",package="ENmix"))
load(system.file("pcclock_model4.RData",package="ENmix"))
pcc=c(pcc1,pcc2,pcc3,pcc4)

#check whether datMeth probe names include suffix, such as in EPICv2 array
if(sum(as.character(refmeth$cg) %in% rownames(datMeth))<50){datMeth=rm.cgsuffix(datMeth)}

#combine reference CpG
#normalize to Horvath refmeth before combine
refcg <- unique(c(
  as.character(refmeth$cg),
  as.character(mPOA_Models$gold_standard_probes),
  as.character(episcore_model$CpG_Site),
  as.character(frailty_model$cg[-1]),
  names(pcc$pcCpGs),
  as.character(refmeth2$cg),
  as.character(cAge_CpG$means$cpg),
  as.character(bAge_CpG$cpgs$CpG_Site),
  as.character(DNAmFitnessModels$AllCpGs)
))

#remove unnecessary probes for estimation
datMeth=datMeth[rownames(datMeth) %in% refcg,,drop=FALSE]

#imputation to fill missing values
noMissingPerCpG <- rowSums(is.na(datMeth))
if(max(noMissingPerCpG)>0 & ncol(as.matrix(datMeth))>1){
message("Imputation to fill missing values...")
if(!fastImputation){
    set.seed(1)
    datMeth <- t(impute::impute.knn(t(datMeth))$data)
}else{
    for (i in which(noMissingPerCpG>0)){
        idx <- is.na(datMeth[i,])
        datMeth[i,idx] <- mean(datMeth[i,],na.rm=TRUE)
}}}

#Count of model required CpGs that are missing in user provided data
## Count of model-required CpGs that are missing in user-provided data
cgcheck <- data.frame(
  predictor = character(),
  nCpG_required = integer(),
  nCpG_present = integer(),
  nCpG_missing = integer(),
  stringsAsFactors = FALSE
)
add_check <- function(name, required, present) {
  cgcheck <<- rbind(
    cgcheck,
    data.frame(
      predictor = name,
      nCpG_required = required,
      nCpG_present = present,
      stringsAsFactors = FALSE
    )
  )
}
rownames_meth <- rownames(datMeth)
add_check("HorvathAge", nrow(horvath) - 1,
          sum(as.character(horvath$cg) %in% rownames_meth))
add_check("PhenoAge", nrow(phenoage) - 1,
          sum(as.character(phenoage$cg) %in% rownames_meth))
add_check("HannumAge",
          nrow(hannum),
          sum(as.character(hannum$cg) %in% rownames_meth))
add_check("PACE",length(mPOA_Models$model_probes),
          sum(as.character(mPOA_Models$model_probes) %in% rownames_meth))
add_check("Frailty",nrow(frailty_model) - 1,
          sum(as.character(frailty_model$cg[-1]) %in% rownames_meth))
add_check("DNAmTL",nrow(DNAmTL_CpGs),
          sum(as.character(DNAmTL_CpGs$ID) %in% rownames_meth))
add_check("pcgtAge",length(EpiToc_CpGs),
          sum(EpiToc_CpGs %in% rownames_meth))
add_check("TNSC",nrow(EpiToc2_CpGs),
          sum(rownames(EpiToc2_CpGs) %in% rownames_meth))
add_check("TNSC2",nrow(EpiToc2_CpGs),
          sum(rownames(EpiToc2_CpGs) %in% rownames_meth))
add_check("Zhang10CpG",nrow(Zhang_10_CpG),
          sum(as.character(Zhang_10_CpG$Marker) %in% rownames_meth))
add_check("Horvath2",nrow(Horvath2_CpGs),
          sum(as.character(Horvath2_CpGs$ID) %in% rownames_meth))
add_check("MiAge", nrow(MiAgedat$MiAge_CpGs),
          sum(as.character(MiAgedat$MiAge_CpGs$CpGs) %in% rownames_meth))
add_check("PedBE",nrow(PEDBE_CpGs), 
	  sum(as.character(PEDBE_CpGs$ID) %in% rownames_meth))
add_check("GACPC",nrow(GACPC_CpGs),
          sum(as.character(GACPC_CpGs$CpG) %in% rownames_meth))
add_check("GARPC",nrow(GARPC_CpGs),
          sum(as.character(GARPC_CpGs$CpG) %in% rownames_meth))
add_check("GARRPC",nrow(GARRPC_CpGs),
          sum(as.character(GARRPC_CpGs$CpG) %in% rownames_meth))
add_check("BohlinGAge", nrow(Bohlin_CpGs),
          sum(as.character(Bohlin_CpGs$CpG) %in% rownames_meth))
add_check("KnightGAge", nrow(Knight_CpGs),
          sum(as.character(Knight_CpGs$CpG) %in% rownames_meth))
#pcclock
pcvar <- c(
  "PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge",
  "PCDNAmTL", "PCPACKYRS", "PCADM", "PCB2M", "PCCystatinC",
  "PCGDF15", "PCLeptin", "PCPAI1", "PCTIMP1", "PCGrimAge"
)
pc_cpgs <- names(pcc$pcCpGs)
for (vv in pcvar) {
  add_check(vv,length(pc_cpgs), sum(pc_cpgs %in% rownames_meth))
}

## cAge
add_check("cAge",nrow(cAge_CpG$means),
  sum(as.character(cAge_CpG$means$cpg) %in% rownames_meth))
## bAge
add_check("bAge",nrow(bAge_CpG$cpgs),
  sum(as.character(bAge_CpG$cpgs$CpG_Site) %in% rownames_meth))
add_check("bAge_Years", nrow(bAge_CpG$cpgs),
  sum(as.character(bAge_CpG$cpgs$CpG_Site) %in% rownames_meth))
## episcore
for (ms in unique(episcore_model$Predictor)) {
  tmp <- episcore_model[episcore_model$Predictor == ms, , drop = FALSE]
  add_check(ms, nrow(tmp), sum(as.character(tmp$CpG_Site) %in% rownames_meth))
}
## DNAm Fit Age
fa <- c(
  "DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max",
  "DNAmGait_wAge", "DNAmGrip_wAge", "DNAmFEV1_wAge",
  "DNAmFitAge"
)
for (i in seq_along(fa)) {
  add_check(fa[i], length(DNAmFitnessModels$AllCpGs),
    sum(as.character(DNAmFitnessModels$AllCpGs) %in% rownames_meth)
  )
}

#GrimAge
add_check("DNAmPACKYRS", length(GrimAgeV1_model$PACKYRS.model[-1]),
          sum(names(GrimAgeV1_model$PACKYRS.model) %in% rownames_meth))
add_check("DNAmADM", length(GrimAgeV1_model$ADM.model[-1]),
          sum(names(GrimAgeV1_model$ADM.model) %in% rownames_meth))
add_check("DNAmB2M", length(GrimAgeV1_model$B2M.model[-1]),
          sum(names(GrimAgeV1_model$B2M.model) %in% rownames_meth))
add_check("DNAmCystatinC", length(GrimAgeV1_model$CystatinC.model[-1]),
          sum(names(GrimAgeV1_model$CystatinC.model) %in% rownames_meth))
add_check("DNAmGDF15", length(GrimAgeV1_model$GDF15.model[-1]),
          sum(names(GrimAgeV1_model$GDF15.model) %in% rownames_meth))
add_check("DNAmLeptin", length(GrimAgeV1_model$Leptin.model),
          sum(names(GrimAgeV1_model$Leptin.model) %in% rownames_meth))
add_check("DNAmPAI1", length(GrimAgeV1_model$PAI1.model),
          sum(names(GrimAgeV1_model$PAI1.model) %in% rownames_meth))
add_check("DNAmTIMP1", length(GrimAgeV1_model$TIMP1.model[-1]),
          sum(names(GrimAgeV1_model$TIMP1.model) %in% rownames_meth))
add_check("DNAmlogCRP", length(GrimAgeV2_model$logCRP.model),
          sum(names(GrimAgeV2_model$logCRP.model) %in% rownames_meth))
add_check("DNAmlogA1C", length(GrimAgeV2_model$logA1C.model[-1]),
          sum(names(GrimAgeV2_model$logA1C.model) %in% rownames_meth))
add_check("GrimAgeV1", length(GrimAge_CpG),
	  sum(GrimAge_CpG %in% rownames_meth))
add_check("GrimAgeV2", length(GrimAge_CpG),
          sum(GrimAge_CpG %in% rownames_meth))

#AdaptAge_CpGs
add_check("AdaptAge", length(AdaptAge_model$CpG),
          sum(AdaptAge_model$CpG %in% rownames_meth))

#AdaptAge_CpGs
add_check("CausAge", length(CausAge_model$CpG),
          sum(CausAge_model$CpG %in% rownames_meth))

#DamAge_CpGs
add_check("DamAge", length(DamAge_model$CpG),
          sum(DamAge_model$CpG %in% rownames_meth))


message("Number of CpGs missing for estimates\n")
cgcheck$nCpG_missing=cgcheck$nCpG_required-cgcheck$nCpG_present
print(cgcheck)
if(max(cgcheck$nCpG_missing)>0){
	message("Missing probes will be imputed with reference values\n")}
message("See file summary_methscore_CpG.csv for the CpG count summary 
	and citation for each predictor")
cgcheck=merge(methscore_dict,cgcheck,by="predictor")
rownames(cgcheck)=cgcheck$predictor
cgcheck=cgcheck[as.character(methscore_dict$predictor),]
write.csv(cgcheck,file="summary_methscore_CpG.csv",row.names=FALSE)

#function for normalization using a modified RCP method, and imputing missing CpGs
norm_impute <-function(dat,refdat,normalize = TRUE, missedcg = character(0)){
   if(normalize){dat=rcp2(dat,refdat)}
    #replace missing probes with reference values
    missedcg <- missedcg[!(missedcg %in% rownames(dat))]
    if(length(missedcg)>0){
        refsub=refdat[as.character(refdat$cg) %in% missedcg, , drop = FALSE]
        missedcg.dat=matrix(rep(refsub$meth_mean,ncol(dat)),
			    nrow=length(missedcg),ncol = ncol(dat))
        rownames(missedcg.dat)=as.character(refsub$cg)
        colnames(missedcg.dat)=colnames(dat)
        dat=rbind(dat,missedcg.dat)
     }
     return(dat)    
}

#Age transformation and probe annotation functions for Hovath age
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); 
    y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) {
       	ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }


###Methylation predictors
#methylation ages
message("Calculating methylation scores ...")
mScore=data.frame(SampleID=colnames(datMeth))

modelcg=unique(c(as.character(horvath$cg[-1]),as.character(phenoage$cg[-1]),
		 as.character(hannum$cg)))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=refmeth
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

#HorvathAge
intercept=horvath$coef[horvath$cg=="(Intercept)"]
horvath=horvath[as.character(horvath$cg) %in% rownames(datMeth2), , drop = FALSE]
mAge=anti.trafo(
    colSums(as.numeric(horvath$coef) * datMeth2[as.character(horvath$cg), , drop = FALSE],
		na.rm=T)+intercept)
mScore$HorvathAge=mAge[as.character(mScore$SampleID)]
#Hannum Age
hannum=hannum[hannum$cg %in% rownames(datMeth2), , drop = FALSE]
mAge=colSums(hannum$coef * datMeth2[as.character(hannum$cg), , drop = FALSE],na.rm=TRUE)
mScore$HannumAge=mAge[as.character(mScore$SampleID)]
#PhenoAge
intercept=phenoage$coef[phenoage$cg=="intercept"]
phenoage=phenoage[phenoage$cg %in% rownames(datMeth2), , drop = FALSE]
mAge=colSums(phenoage$coef * datMeth2[as.character(phenoage$cg), , drop = FALSE],na.rm=T)+intercept
mScore$PhenoAge=mAge[as.character(mScore$SampleID)]

#PACE estimation
modelcg=unique(as.character(mPOA_Models$model_probes))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=mPOA_Models$gold_standard_probes,meth_mean=mPOA_Models$gold_standard_means)
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

PACE=DunedinPACE(betas=datMeth2, proportionOfProbesRequired = 0.8)
mScore$PACE=PACE[as.character(mScore$SampleID)]

#Frailty
modelcg=unique(as.character(frailty_model$cg[-1]))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=frailty_model[-1,c("cg","meth_mean")]
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

intercept=frailty_model$coef[frailty_model$cg=="intercept"]
frailty_model=frailty_model[frailty_model$cg %in% rownames(datMeth2), , drop = FALSE]
eFRS=colSums(frailty_model$coef * datMeth2[as.character(frailty_model$cg), , drop = FALSE],na.rm=T)+intercept
mScore$Frailty=eFRS[as.character(mScore$SampleID)]


#PEDBE,EpiToc,EpiToc2,Zhang10CpG,Horvath2,MiAge,DNAmTL,PEDBE,GACPC,GARPC,GARRPC,Bohlin,Knight
modelcg <- unique(c(
  as.character(PEDBE_CpGs$ID),
  as.character(EpiToc_CpGs),
  rownames(EpiToc2_CpGs),
  as.character(Zhang_10_CpG$Marker),
  as.character(Horvath2_CpGs$ID),
  as.character(MiAgedat$MiAge_CpGs$CpGs),
  as.character(DNAmTL_CpGs$ID),
  as.character(GACPC_CpGs$CpG),
  as.character(GARPC_CpGs$CpG),
  as.character(GARRPC_CpGs$CpG),
  as.character(Bohlin_CpGs$CpG),
  as.character(Knight_CpGs$CpG)
))
missedcg <- setdiff(modelcg, rownames(datMeth))
refdat=refmeth2; names(refdat)=c("cg","meth_mean")
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

#DNAmTL
DNAmTL_CpGs=DNAmTL_CpGs[DNAmTL_CpGs$ID %in% rownames(datMeth2), , drop = FALSE]
DNAmTL=colSums(DNAmTL_CpGs$Coef * datMeth2[as.character(DNAmTL_CpGs$ID), , drop = FALSE],na.rm=TRUE)+7.924780053
mScore$DNAmTL=DNAmTL[as.character(mScore$SampleID)]

#EpiTOC
pcgtAge <- colMeans(datMeth2[rownames(datMeth2) %in% EpiToc_CpGs, , drop = FALSE],na.rm=TRUE)
mScore$pcgtAge=pcgtAge[as.character(mScore$SampleID)]

#EpiTOC2
EpiToc2_CpGs=EpiToc2_CpGs[rownames(EpiToc2_CpGs) %in% rownames(datMeth2), , drop = FALSE]
tmp.m=datMeth2[rownames(EpiToc2_CpGs), , drop = FALSE]
TNSC <- 2*colMeans((1/(EpiToc2_CpGs[,1]*(1-EpiToc2_CpGs[,2]))) * (tmp.m - EpiToc2_CpGs[,2]),na.rm=TRUE)
#approximated = T
TNSC2 <- 2*colMeans((1/EpiToc2_CpGs[,1]) * tmp.m,na.rm=TRUE)
mScore$TNSC=TNSC[as.character(mScore$SampleID)]
mScore$TNSC2=TNSC2[as.character(mScore$SampleID)]

#Zhang 10 CpG clock
Zhang_10_CpG=Zhang_10_CpG[as.character(Zhang_10_CpG$Marker) %in% rownames(datMeth2), , drop = FALSE]
Zhang10CpG=colSums(Zhang_10_CpG$coef * datMeth2[as.character(Zhang_10_CpG$Marker), , drop = FALSE],na.rm=T)
mScore$Zhang10CpG=Zhang10CpG[as.character(mScore$SampleID)]

#Horvath2
Horvath2_CpGs=Horvath2_CpGs[as.character(Horvath2_CpGs$ID) %in% rownames(datMeth2), , drop = FALSE]
Horvath2=anti.trafo(colSums(Horvath2_CpGs$Coef * datMeth2[as.character(Horvath2_CpGs$ID), , drop = FALSE],na.rm=T)-0.447119319)
mScore$Horvath2=Horvath2[as.character(mScore$SampleID)]

#MiAge mitotic age
#function used by the original MiAge code http://www.columbia.edu/~sw2206/softwares.htm
MiAge_fr2 <- function (x, b, c, d, betaj)
{
    nj = x
    return(sum((c + b^(nj - 1) * d - betaj)^2, na.rm = T))
}
MiAge_grr2 <-function (x, b, c, d, betaj)
{
    nj = x
    return(2 * sum((c + b^(nj - 1) * d - betaj) * b^(nj - 1) *
        log(b) * d, na.rm = T))
}
mitotic.age <- function (beta, b, c, d)
{
#    library(methods)
    upperage = 10000
    lowerage = 10
    n = rep(500, ncol(beta));names(n)=colnames(beta)
    no.initial.n = 5
    for (j in 1:ncol(beta)) {
        current.value = MiAge_fr2(n[j], b, c, d, beta[, j])
        columnoptim = vector("list", no.initial.n)
        val = rep(NA, no.initial.n)
        for (jj in 1:(no.initial.n - 1)) {
            temp = try(optim(par = lowerage + jj * (upperage -
                lowerage)/no.initial.n, fn = MiAge_fr2, gr = MiAge_grr2,
                b = b, c = c, d = d, betaj = beta[, j], method = "L-BFGS-B",
                lower = lowerage, upper = upperage, control = list(factr = 1)),
                silent = T)
            if (!is(temp, "try-error")) {
                columnoptim[[jj]] = temp
                val[jj] = temp$value
            }
        }
        temp = try(optim(par = n[j], fn = MiAge_fr2, gr = MiAge_grr2,
            b = b, c = c, d = d, betaj = beta[, j], method = "L-BFGS-B",
            lower = lowerage, upper = upperage, control = list(factr = 1)),
            silent = T)
        if (!is(temp, "try-error")) {
            columnoptim[[no.initial.n]] = temp
            val[no.initial.n] = temp$value
        }
        temp = columnoptim[[which(val == min(val, na.rm = T))[1]]]
        if (!is(temp, "try-error")) {
            n[j] = temp$par
        }
        else {
            print(2)
            print(temp)
        }
    }
    return(n)
}

MiAge=mitotic.age(
		  datMeth2[na.omit(match(MiAgedat$MiAge_CpGs$CpGs,
		  rownames(datMeth2))),,drop=FALSE],MiAgedat$MiAge_parameters[[1]],
		  MiAgedat$MiAge_parameters[[2]],
		  MiAgedat$MiAge_parameters[[3]])
mScore$MiAge=MiAge[as.character(mScore$SampleID)]

#PedBE The Pediatric-Buccal-Epigenetic (PedBE) clock
PEDBE_CpGs=PEDBE_CpGs[PEDBE_CpGs$ID %in% rownames(datMeth2),]
PedBE=anti.trafo(colSums(PEDBE_CpGs$Coef * datMeth2[as.character(PEDBE_CpGs$ID),,drop=FALSE],na.rm=T)-2.10)
mScore$PedBE=PedBE[as.character(mScore$SampleID)]

#Placental epigenetic clocks
#Control placental clock (CPC)
GACPC_CpGs=GACPC_CpGs[GACPC_CpGs$CpG %in% rownames(datMeth2),]
GACPC=colSums(GACPC_CpGs$coef * datMeth2[as.character(GACPC_CpGs$CpG),,drop=FALSE],na.rm = T)+13.06182
mScore$GACPC=GACPC[as.character(mScore$SampleID)]

#Robust placental clock (RPC)
GARPC_CpGs=GARPC_CpGs[GARPC_CpGs$CpG %in% rownames(datMeth2),]
GARPC=colSums(GARPC_CpGs$coef * datMeth2[as.character(GARPC_CpGs$CpG),,drop=FALSE], na.rm = T)+24.99772
mScore$GARPC=GARPC[as.character(mScore$SampleID)]

#Refined robust placental clock for uncomplicated term pregnancies
GARRPC_CpGs=GARRPC_CpGs[GARRPC_CpGs$CpG %in% rownames(datMeth2),]
GARRPC=colSums(GARRPC_CpGs$coef * datMeth2[as.character(GARRPC_CpGs$CpG),,drop=FALSE], na.rm = T)+ 30.74966
mScore$GARRPC=GARRPC[as.character(mScore$SampleID)]

#Bohlin Gestational age
Bohlin_CpGs=Bohlin_CpGs[Bohlin_CpGs$CpG %in% rownames(datMeth2),]
BohlinGAge=colSums(Bohlin_CpGs$coef * datMeth2[as.character(Bohlin_CpGs$CpG),,drop=FALSE],na.rm=T)+ 277.2421
mScore$BohlinGAge=BohlinGAge[as.character(mScore$SampleID)]

#Knight Gestational age
Knight_CpGs=Knight_CpGs[Knight_CpGs$CpG %in% rownames(datMeth2),]
KnightGAge=colSums(Knight_CpGs$coef * datMeth2[as.character(Knight_CpGs$CpG),,drop=FALSE],na.rm=T) + 41.7
mScore$KnightGAge=KnightGAge[as.character(mScore$SampleID)]

#PC clocks
modelcg=names(pcc$pcCpGs)
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=names(pcc$pcCpGs),meth_mean=pcc$pcCpGs)
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

if(!is.null(datPheno)){
    pcclock=calcPCClocks(methdat=datMeth2,phenodat=datPheno,pcc)
    pcclock=pcclock[,!(names(pcclock) %in% c("Age","Female"))]
    mScore <- cbind(mScore, pcclock[match(mScore$SampleID, pcclock$SampleID), setdiff(colnames(pcclock), "SampleID"),drop=FALSE])
}

# Calculate GrimAge Clocks
# adopted code from methylCIPHER
modelcg <- GrimAge_CpG
missedcg <- setdiff(modelcg, rownames(datMeth))
datMeth2=norm_impute(datMeth,refdat=refmeth2,normalize,missedcg)
if(!is.null(datPheno)){
DNAm <- datMeth2[GrimAge_CpG, , drop = FALSE]
pp <- datPheno
rownames(pp)=pp$SampleID
DNAm <- DNAm [,rownames(pp),drop=FALSE]
DNAm <- rbind(DNAm, Age= pp[colnames(DNAm), "Age"])
pp$DNAmPACKYRS <- drop(GrimAgeV1_model$PACKYRS.model %*% DNAm[names(GrimAgeV1_model$PACKYRS.model),,drop = FALSE]) + GrimAgeV1_model$PACKYRS.intercept
pp$DNAmADM <- drop(GrimAgeV1_model$ADM.model %*% DNAm[names(GrimAgeV1_model$ADM.model), , drop = FALSE]) +  GrimAgeV1_model$ADM.intercept
pp$DNAmB2M <- drop(GrimAgeV1_model$B2M.model %*% DNAm[names(GrimAgeV1_model$B2M.model), , drop = FALSE]) +  GrimAgeV1_model$B2M.intercept
pp$DNAmCystatinC <- drop(GrimAgeV1_model$CystatinC.model %*% DNAm[names(GrimAgeV1_model$CystatinC.model), , drop = FALSE]) +  GrimAgeV1_model$CystatinC.intercept
pp$DNAmGDF15 <- drop(GrimAgeV1_model$GDF15.model %*% DNAm[names(GrimAgeV1_model$GDF15.model), , drop = FALSE]) +  GrimAgeV1_model$GDF15.intercept
pp$DNAmLeptin <- drop(GrimAgeV1_model$Leptin.model %*% DNAm[names(GrimAgeV1_model$Leptin.model), , drop = FALSE]) +  GrimAgeV1_model$Leptin.intercept
pp$DNAmPAI1 <- drop(GrimAgeV1_model$PAI1.model %*% DNAm[names(GrimAgeV1_model$PAI1.model), , drop = FALSE]) +  GrimAgeV1_model$PAI1.intercept
pp$DNAmTIMP1 <- drop(GrimAgeV1_model$TIMP1.model %*% DNAm[names(GrimAgeV1_model$TIMP1.model), , drop = FALSE]) +  GrimAgeV1_model$TIMP1.intercept
pp$GrimAgeV1 <- drop(as.matrix(pp[, GrimAgeV1_model$components, drop = FALSE]) %*%    GrimAgeV1_model$GrimAge.model)
y <- pp$GrimAgeV1
pp$GrimAgeV1 <- (((y - GrimAgeV1_model$GrimAge.transform[3]) /  GrimAgeV1_model$GrimAge.transform[4]) *  GrimAgeV1_model$GrimAge.transform[2]) +  GrimAgeV1_model$GrimAge.transform[1]
pp$DNAmlogA1C <- drop(GrimAgeV2_model$logA1C.model %*% DNAm[names(GrimAgeV2_model$logA1C.model), , drop = FALSE]) +  GrimAgeV2_model$logA1C.intercept
pp$DNAmlogCRP <- drop(GrimAgeV2_model$logCRP.model %*% DNAm[names(GrimAgeV2_model$logCRP.model), , drop = FALSE]) +  GrimAgeV2_model$logCRP.intercept
pp$GrimAgeV2 <- drop(as.matrix(pp[, GrimAgeV2_model$components, drop = FALSE]) %*%    GrimAgeV2_model$GrimAge.model)
y <- pp$GrimAgeV2
pp$GrimAgeV2 <- (((y - GrimAgeV2_model$GrimAge.transform[3]) /  GrimAgeV2_model$GrimAge.transform[4]) *  GrimAgeV2_model$GrimAge.transform[2]) +  GrimAgeV2_model$GrimAge.transform[1]
pp <- pp[match(mScore$SampleID, rownames(pp)), , drop = FALSE]
mScore <- cbind( mScore, pp[, !(colnames(pp) %in% c("Age", "Female", "SampleID")), drop = FALSE])
}


#reference value for AdaptAge CausAge DamAge are in refmeth2
modelcg <- c(AdaptAge_model$CpG,CausAge_model$CpG,DamAge_model$CpG)
missedcg <- setdiff(modelcg, rownames(datMeth))
datMeth2=norm_impute(datMeth,refdat=refmeth2,normalize,missedcg)
#AdaptAge
AdaptAge <- drop(AdaptAge_model$Beta %*% datMeth2[AdaptAge_model$CpG,,drop=FALSE]) - 511.9742762
mScore$AdaptAge <- AdaptAge[match(mScore$SampleID,names(AdaptAge))]

#CausAge
CausAge <- drop(CausAge_model$Beta %*% datMeth2[CausAge_model$CpG,,drop=FALSE]) + 86.80816381
mScore$CausAge <- CausAge[match(mScore$SampleID,names(CausAge))]

#DamAge
DamAge <- drop(DamAge_model$Beta %*% datMeth2[DamAge_model$CpG,,drop=FALSE]) + 543.4315887
mScore$DamAge <- DamAge[match(mScore$SampleID,names(DamAge))]


#cAge and bAge
modelcg=unique(as.character(cAge_CpG$means$cpg))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=as.character(cAge_CpG$means$cpg),meth_mean=cAge_CpG$means$mean)
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)
mScore$cAge=calc_cAge(datMeth2,cAge_CpG)[as.character(mScore$SampleID)]

#use episcore reference data for bAge
modelcg=unique(as.character(bAge_CpG$cpgs$CpG_Site))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
tmp0=episcore_model[!duplicated(episcore_model$CpG_Site),]
refdat=data.frame(cg=tmp0$CpG_Site,meth_mean=tmp0$Mean_Beta_Value)
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)
tmp=calc_bAge(datMeth2,datPheno,mScore,bAge_CpG,GrimAgeComponent)[as.character(mScore$SampleID),]
mScore$bAge=tmp$bAge
mScore$bAge_Years=tmp$bAge_Years

#episcore
modelcg=unique(as.character(episcore_model$CpG_Site))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
tmp0=episcore_model[!duplicated(episcore_model$CpG_Site),]
refdat=data.frame(cg=tmp0$CpG_Site,meth_mean=tmp0$Mean_Beta_Value)
datMeth2=norm_impute(datMeth,refdat=refdat,normalize,missedcg)

out <- data.frame(SampleID=colnames(datMeth2))
for(i in unique(episcore_model$Predictor)){
  tmp_coef = episcore_model[episcore_model$Predictor %in% i, ]
  tmp_coef=tmp_coef[tmp_coef$CpG_Site %in% rownames(datMeth2),]
  if(nrow(tmp_coef) > 1) {
    out[,i]=colSums(tmp_coef$Coefficient * datMeth2[as.character(tmp_coef$CpG_Site),],na.rm=T)
  } else {
    out[,i] = tmp_coef$Coefficient * datMeth2[as.character(tmp_coef$CpG_Site),]
  }
}
out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
mScore <- cbind(mScore, out[match(mScore$SampleID, out$SampleID), setdiff(colnames(out), "SampleID"), drop = FALSE])

#DNAmFitAge
FitAge=calc_DNAmFitAge(datMeth,datPheno,mScore,GrimAgeComponent,DNAmFitnessModels)
mScore <- cbind(mScore, FitAge[match(mScore$SampleID, FitAge$SampleID), setdiff(colnames(FitAge), "SampleID"), drop = FALSE])



#column order
mScore=mScore[,c("SampleID",as.character(methscore_dict$predictor))]
###################Age residuals
if(AgeResid){
if(!is.null(datPheno)){
    rownames(datPheno)=datPheno$SampleID
    datPheno=datPheno[as.character(mScore$SampleID),]
    clockColumns=colnames(mScore)[!(colnames(mScore)=="SampleID")]
    mScore[,"Age"] = datPheno$Age
    for (i in clockColumns){
        mScore[,paste0(i,"Resid")] = resid(lm(mScore[,i] ~ datPheno$Age))
    }
}}
return(mScore)
}

#calculation of PC Clocks
#Higgins-Chen et al, PMID 36277076,https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9586209/
calcPCClocks <- function(methdat, phenodat,pcc=NULL){

  if(!("Age" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named Age")
  }
  if(!("Female" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named Female")
  }
  if(!("SampleID" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named SampleID")
  }else if(sum(colnames(methdat) %in% phenodat$SampleID)<ncol(methdat)){
     stop("Error: some samples in methdat do not have phenotype in phenodat")
  }else{
     rownames(phenodat)=phenodat$SampleID
     phenodat=phenodat[colnames(methdat),]
  }

anti.trafo <-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

  #check number of missing CpGs
  flag=apply(is.na(methdat),1,sum)!=ncol(methdat)
  methdat=methdat[flag,]
#  load(file = paste(path_to_PCClocks_directory,"pcclock_model.RData", sep = ""))
#  if(is.null(pcc)){load(system.file("pcclock_model.RData",package="ENmix"))}
  CpGs=names(pcc$pcCpGs)
  methdat=methdat[rownames(methdat) %in% CpGs,]
  missingCpGs=CpGs[!(CpGs %in% rownames(methdat))]
#  message(paste0("PC clocks: Missing ",length(missingCpGs)," out of ",length(CpGs)," CpGs required for PC clocks calculation"))
  #impute missing probes and values
  tmp=matrix(rep(pcc$pcCpGs[missingCpGs],ncol(methdat)),ncol=ncol(methdat))
  rownames(tmp)=missingCpGs
  methdat=rbind(methdat,tmp)
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  methdat <- apply(methdat,1,meanimpute)
  methdat=methdat[,CpGs]
  phenodat=phenodat[rownames(methdat),]
  #Calculate PC Clocks
  #Initialize a data frame for PC clocks
  DNAmAge <- data.frame(phenodat)
#  message("Calculating PC Clocks now")
  DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(methdat),2,pcc$PCHorvath1$center) %*% pcc$PCHorvath1$model + pcc$PCHorvath1$intercept))
  DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(methdat),2,pcc$PCHorvath2$center) %*% pcc$PCHorvath2$model + pcc$PCHorvath2$intercept))
  DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCHannum$center) %*% pcc$PCHannum$model + pcc$PCHannum$intercept)
  DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCPhenoAge$center) %*% pcc$PCPhenoAge$model + pcc$PCPhenoAge$intercept)
  DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCDNAmTL$center) %*% pcc$PCDNAmTL$model + pcc$PCDNAmTL$intercept)
  pp=cbind(Female = DNAmAge$Female,Age = DNAmAge$Age)
  DNAmAge$PCPACKYRS <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCPACKYRS$model + pp %*% pcc$PCPACKYRS$pfa + pcc$PCPACKYRS$intercept)
  DNAmAge$PCADM <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCADM$model + pp %*% pcc$PCADM$pfa + pcc$PCADM$intercept)
  DNAmAge$PCB2M <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCB2M$model + as.matrix(pp[,2]) %*% pcc$PCB2M$pfa + pcc$PCB2M$intercept)
  DNAmAge$PCCystatinC <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCCystatinC$model + pp %*% pcc$PCCystatinC$pfa + pcc$PCCystatinC$intercept)
  DNAmAge$PCGDF15 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCGDF15$model + as.matrix(pp[,2]) %*% pcc$PCGDF15$pfa + pcc$PCGDF15$intercept)
  DNAmAge$PCLeptin <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCLeptin$model + pp %*% pcc$PCLeptin$pfa + pcc$PCLeptin$intercept)
  DNAmAge$PCPAI1 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCPAI1$model + as.matrix(pp[,1]) %*% pcc$PCPAI1$pfa + pcc$PCPAI1$intercept)
  DNAmAge$PCTIMP1 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCTIMP1$model + pp %*% pcc$PCTIMP1$pfa + pcc$PCTIMP1$intercept)
  DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,pcc$PCGrimAge$components]) %*% pcc$PCGrimAge$model + pcc$PCGrimAge$intercept)

 return(DNAmAge)
}

#cAge
calc_cAge <- function(data,cAge_CpG){
#load(system.file("cAge.RData",package="ENmix"))
coef_linear=cAge_CpG$coef_linear
coef_log=cAge_CpG$coef_log
intercept=cAge_CpG$intercept
intercept_log=cAge_CpG$intercept_log
means=cAge_CpG$means
rownames(coef_linear)=coef_linear$CpG_Site
rownames(coef_log)=coef_log$CpG_Site

data <- data[rownames(data) %in% rownames(means),]
## impute missing CpGs
if (nrow(data)!=nrow(means)) {
    missing_cpgs <- means[!rownames(means) %in% rownames(data),]
    mat=matrix(rep(missing_cpgs$mean,ncol(data)),ncol=ncol(data))
    rownames(mat)=rownames(missing_cpgs)
    colnames(mat)=colnames(data)
    data <- rbind(data,mat)
}
## impute missing values
na_to_mean <-function(x) {x[is.na(x)] <- mean(x, na.rm=T);x}
data <- t(apply(data,1,na_to_mean))

## Prep for linear predictor
cg1=as.character(coef_linear$CpG_Site[-grep('_2', coef_linear$CpG_Site)])
cg2=as.character(coef_linear$CpG_Site[grep('_2', coef_linear$CpG_Site)])
cg2=gsub("_2","",cg2)
tmp <- data[cg2,]^2;rownames(tmp)=paste0(rownames(tmp),"_2")
scores_linear <- rbind(data[cg1,],tmp)

## Prep for log predictor
cg1=as.character(coef_log$CpG_Site[-grep('_2', coef_log$CpG_Site)])
cg2=as.character(coef_log$CpG_Site[grep('_2', coef_log$CpG_Site)])
cg2=gsub("_2","",cg2)
tmp <- data[cg2,]^2;rownames(tmp)=paste0(rownames(tmp),"_2")
scores_log <- rbind(data[cg1,],tmp)

## Calculate cAge with linear model
coef_linear <- coef_linear[rownames(scores_linear),]
pred_linear_pp <- colSums(scores_linear * coef_linear$Coefficient) + intercept

## Identify any individuals predicted as under 20s, and re-run with model trained on log(age)
over20s <- names(pred_linear_pp[pred_linear_pp > 20])
pred_linear_pp1 <- pred_linear_pp[over20s]
under20s <- names(pred_linear_pp[pred_linear_pp < 20])

## Now re-run model for those
coef_log <- coef_log[rownames(scores_log),]
pred_log_pp <- colSums(scores_log * coef_log$Coefficient) + intercept_log
pred_log_pp <- exp(pred_log_pp[under20s])

c(pred_log_pp, pred_linear_pp1)
}

#bAge
calc_bAge<-function(data,phenodat,grim,bAge_CpG,GrimAgeComponent){
grim$DNAmGrimAge=grim$PCGrimAge
rownames(grim)=grim$SampleID
rownames(phenodat)=phenodat$SampleID
cpgs=bAge_CpG$cpgs
coefficients=bAge_CpG$coefficients
##standardize CpG by CpG
coef <- data[rownames(data) %in% as.character(cpgs$CpG_Site),]
ids <- colnames(coef)
coef <- t(apply(coef, 1, scale))
colnames(coef) <- ids

##impute missing CpGs
if (nrow(coef) != length(unique(cpgs$CpG_Site))) {
    missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)), c("CpG_Site", "Mean_Beta_Value")]
    missing_cpgs=missing_cpgs[!duplicated(as.character(missing_cpgs$CpG_Site)),]
    mat=matrix(rep(missing_cpgs$Mean_Beta_Value,ncol(coef)),ncol=ncol(coef))
    rownames(mat)=as.character(missing_cpgs$CpG_Site)
    colnames(mat)=colnames(coef)
    coef = rbind(coef,mat)
}
## impute NA with CpG means
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}
coef <- t(apply(coef,1,function(x) na_to_mean(x)))

#Calculate Episcores
loop <- unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  tmp=coef[as.character(tmp_coef$CpG_Site),]
  if(nrow(tmp_coef) > 1) {
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp,na.rm=T)
  } else {
    out[colnames(coef),i] = tmp*tmp_coef$Coefficient
  }
}

###### Standadize GrimAge components
samples <- rownames(out)
grim_pred <- grim[samples, c("DNAmGrimAge"), drop = FALSE]
if(GrimAgeComponent %in% c(1,2)){
  grim <- grim[samples, c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")]
}else{
    grim <- grim[samples, c("PCADM", "PCB2M", "PCCystatinC", "PCGDF15", "PCLeptin", "PCPACKYRS", "PCPAI1", "PCTIMP1")]
}
colnames(grim)=c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")
grim <- scale(grim)
## Calculate bAge
scores <- cbind(Age=phenodat[samples, c("Age")], grim, out)
scores <- scores[, coefficients$Variable]
pred_pp <- colSums(coefficients[,"Coefficient"] * t(scores))

## Scale to same scale as age
scale_pred <- function(x, mean_pred, sd_pred, mean_test, sd_test) {
  scaled <- mean_test + (x - mean_pred)*(sd_test/sd_pred)
  return(scaled)
}
# Scale to same Z scale
scale_Z <- function(x, mean_pred, sd_pred) {
  scaled <- (x - mean_pred)/sd_pred
  return(scaled)
}

mean_pred <- mean(pred_pp)
mean_test <- mean(phenodat$Age) # Mean age in testing data
sd_pred <- sd(pred_pp)
sd_test <- sd(phenodat$Age) # SD age in testing data

bAge <- scale_Z(pred_pp, mean_pred, sd_pred)
bAge_Years <- scale_pred(pred_pp, mean_pred, sd_pred, mean_test, sd_test)

pred=data.frame(SampleID=names(bAge),bAge=bAge,bAge_Years=bAge_Years)
rownames(pred)=pred$SampleID
return(pred)
}


#######DNAmFitAge begin
calc_DNAmFitAge<-function(datMeth,datPheno,mScore,GrimAgeComponent,DNAmFitnessModels){

	if (GrimAgeComponent == 1) {
           grim1 <- mScore[, c("SampleID", "GrimAgeV1")]
        } else if (GrimAgeComponent == 2) {
           grim1 <- mScore[, c("SampleID", "GrimAgeV2")]
        } else {
           grim1 <- mScore[, c("SampleID", "PCGrimAge")]
	}
        names(grim1)[2] <- "DNAmGrimAge"

	pheno=merge(datPheno,grim1,by="SampleID")
	pheno=pheno[,c("SampleID","Age","Female","DNAmGrimAge")]
	rownames(pheno)=pheno$SampleID
	pheno=pheno[colnames(datMeth),]

data_prep <- function(dataset,pheno){
  dataset <- dataset[rownames(dataset) %in% DNAmFitnessModels$AllCpGs,]
  if(nrow(dataset) != length(DNAmFitnessModels$AllCpGs)){
      cpgs_toadd <- DNAmFitnessModels$AllCpGs[!DNAmFitnessModels$AllCpGs %in% rownames(dataset)]

      #separate by sex to impute missing CpGs
      output=NULL
      rownames(pheno)=pheno$SampleID
      cid=intersect(colnames(dataset), rownames(pheno))
      dataset=dataset[,cid];pheno=pheno[cid,]
      if(sum(pheno$Female==1)>0){
         dat <- dataset[,pheno$Female == 1]
         dat=rbind(as.matrix(dat),matrix(rep(DNAmFitnessModels$Female_Medians_All[cpgs_toadd],ncol(dat)),ncol=ncol(dat),dimnames=list(cpgs_toadd,colnames(dat))))
         output=cbind(output,dat)
      }
      if(sum(pheno$Female==0)>0){
         dat <- dataset[,pheno$Female == 0]
         dat=rbind(as.matrix(dat),matrix(rep(DNAmFitnessModels$Male_Medians_All[cpgs_toadd],ncol(dat)),ncol=ncol(dat),dimnames=list(cpgs_toadd,colnames(dat))))
         output=cbind(output,dat)
      }

#      print(paste0("Total ", length(cpgs_toadd)," Missing CpGs that are assigned median values from training data"))
      dataset=output
  }
  return(dataset)
}
# Function to provide estimates for any 1 DNAm fitness models
# TidyModel is a specific model in DNAmFitnessModels list
DNAmEstimatorAnyModel <- function(dataset, TidyModel){
  intercept=matrix(rep(1.0, ncol(dataset)),ncol=ncol(dataset),dimnames=list("(Intercept)",colnames(dataset)))
  dataset=rbind(intercept,dataset)
  dataset <- dataset[as.character(TidyModel$term),]
  dm=dimnames(dataset)
  dataset=matrix(as.numeric(dataset),nrow=nrow(dataset))
  dimnames(dataset)=dm
  estimate <- colSums(TidyModel$estimate  * dataset)
  return(estimate)
}

# Function to calculate all DNAm fitness estimates #
DNAmFitnessEstimators <- function(data, pheno){
  rownames(pheno)=pheno$SampleID
  pheno=pheno[colnames(data),]
  data=rbind(t(pheno),data)

  if(sum(pheno$Female ==1)>0){
  data_fem <- data[,pheno$Female ==1]
  fem_est1 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_noAge_Females) # gait without age
  fem_est2 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_noAge_Females) # grip
  fem_est3 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$VO2maxModel) # vo2max
  fem_est4 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_wAge_Females) # gait w age
  fem_est5 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_wAge_Females) # grip w age
  fem_est6 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$FEV1_wAge_Females) # fev1 w age
  }
  if(sum(pheno$Female ==0)>0){
  data_male <- data[,pheno$Female ==0]
  male_est1 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_noAge_Males) # gait
  male_est2 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_noAge_Males) # grip
  male_est3 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$VO2maxModel) # vo2max
  male_est4 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_wAge_Males) # gait
  male_est5 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_wAge_Males) # grip
  male_est6 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$FEV1_wAge_Males) # fev1
  }
  if(sum(pheno$Female ==1)==0){est1=male_est1;est2=male_est2;est3=male_est3;est4=male_est4;est5=male_est5;est6=male_est6
    }else if(sum(pheno$Female ==0)==0){est1=fem_est1;est2=fem_est2;est3=fem_est3;est4=fem_est4;est5=fem_est5;est6=fem_est6
    }else{
          est1 <- c(fem_est1, male_est1)
          est2 <- c(fem_est2, male_est2)
          est3 <- c(fem_est3, male_est3)
          est4 <- c(fem_est4, male_est4)
          est5 <- c(fem_est5, male_est5)
          est6 <- c(fem_est6, male_est6)
    }
  data_and_est=data.frame(DNAmGait_noAge=est1,DNAmGrip_noAge=est2,DNAmVO2max=est3,DNAmGait_wAge=est4,DNAmGrip_wAge=est5,DNAmFEV1_wAge=est6)
  data_and_est$SampleID=rownames(data_and_est)
  return(data_and_est)
}


# Function to calculate DNAmFitAge
FitAgeEstimator <- function(data){
  # prep dataset for fitage- dataset by sex
  fem <- data[data$Female == 1,]
  male <- data[data$Female == 0,]

  # can only estimate FitAge if all variables are present, remove those without them #
  fem_comcase <- fem[complete.cases(fem), ]
  male_comcase <- male[complete.cases(male), ]

  # Female FitAge
  female_fitest <- 0.1044232 * ((fem_comcase$DNAmVO2max - 46.825091) / (-0.13620215)) +
                0.1742083 * ((fem_comcase$DNAmGrip_noAge - 39.857718) / (-0.22074456)) +
                0.2278776 * ((fem_comcase$DNAmGait_noAge - 2.508547) / (-0.01245682))  +
                0.4934908 * ((fem_comcase$DNAmGrimAge - 7.978487) / (0.80928530))

  # Male FitAge
  male_fitest <- 0.1390346 * ((male_comcase$DNAmVO2max - 49.836389) / (-0.141862925)) +
                0.1787371 * ((male_comcase$DNAmGrip_noAge - 57.514016) / (-0.253179827)) +
                0.1593873 * ((male_comcase$DNAmGait_noAge - 2.349080) / (-0.009380061))  +
                0.5228411 * ((male_comcase$DNAmGrimAge - 9.549733) / (0.835120557))

  fem <- data.frame(fem_comcase, DNAmFitAge = female_fitest)
  male <- data.frame(male_comcase, DNAmFitAge = male_fitest)

  resu <- rbind(fem, male)
  resu=resu[,!names(resu) %in% c("Age","Female","DNAmGrimAge")]
  return(resu)
}

sample_data_prep <- data_prep(dataset = datMeth,pheno=pheno)
sample_data_FitnessEst <- DNAmFitnessEstimators(sample_data_prep, pheno)
sample_data_FitAge_prep <- merge(pheno, sample_data_FitnessEst, by = "SampleID")
FitAge <- FitAgeEstimator(sample_data_FitAge_prep)
return(FitAge)
}
##### end DNAmFitAge

