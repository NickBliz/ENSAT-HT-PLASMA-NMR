#remove peaks not present in at least 80% of samples in each group
setwd("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/Publication 3")
library(readr)
UF_PILOT4_3_withexcl <- read_csv("ENP_final_3_new_FINAL_peaks_excluded_new2_FM.csv")
UF_PILOT4_3_withexcl$`2.022404917`=NULL
UF_PILOT4_3_withexcl$`1.903997967`=NULL
UF_PILOT4_3_withexcl$`2.214349772`=NULL
UF_PILOT4_3_withexcl$`3.185692061`=NULL
UF_PILOT4_3_withexcl$`2.262080201`=NULL
UF_PILOT4_3_withexcl$`2.275711999`=NULL
UF_PILOT4_3_withexcl$`2.313307518`=NULL
UF_PILOT4_3_withexcl$`2.37023004`=NULL
UF_PILOT4_3_withexcl$`2.388646009`=NULL
UF_PILOT4_3_withexcl$`4.132772125`=NULL
#PEAKS EXCLUDED FOR APPROACH B:
#UF_PILOT4_3_withexcl$`5.219707859`=NULL
#UF_PILOT4_3_withexcl$`5.227224754`=NULL
#UF_PILOT4_3_withexcl$`4.080120359`=NULL
#UF_PILOT4_3_withexcl$`4.093821199`=NULL
#UF_PILOT4_3_withexcl$`4.107654959`=NULL
#UF_PILOT4_3_withexcl$`4.120725314`=NULL
#UF_PILOT4_3_withexcl$`1.320775629`=NULL
#UF_PILOT4_3_withexcl$`1.307020131`=NULL
#UF_PILOT4_3_withexcl$`2.324851518`=NULL
#UF_PILOT4_3_withexcl$`2.332132451`=NULL
#UF_PILOT4_3_withexcl$`2.341139339`=NULL
#UF_PILOT4_3_withexcl$`2.046730765`=NULL
#UF_PILOT4_3_withexcl$`2.059835912`=NULL
#UF_PILOT4_3_withexcl$`2.074880492`=NULL
#UF_PILOT4_3_withexcl$`2.095209445`=NULL
#UF_PILOT4_3_withexcl$`2.102557302`=NULL
#UF_PILOT4_3_withexcl$`2.107835865`=NULL
#UF_PILOT4_3_withexcl$`2.112807515`=NULL
#UF_PILOT4_3_withexcl$`2.12163218`=NULL
#UF_PILOT4_3_withexcl$`2.13235951`=NULL
#UF_PILOT4_3_withexcl$`2.139543856`=NULL
#UF_PILOT4_3_withexcl$`2.144890749`=NULL
#UF_PILOT4_3_withexcl$`2.418146656`=NULL
#UF_PILOT4_3_withexcl$`2.42836363`=NULL
#UF_PILOT4_3_withexcl$`2.433390098`=NULL
#UF_PILOT4_3_withexcl$`2.443586856`=NULL
#UF_PILOT4_3_withexcl$`2.448745467`=NULL
#UF_PILOT4_3_withexcl$`2.460335404`=NULL
#UF_PILOT4_3_withexcl$`3.041061161`=NULL
#UF_PILOT4_3_withexcl$`3.057273107`=NULL
#UF_PILOT4_3_withexcl$`3.554511486`=NULL
#UF_PILOT4_3_withexcl$`3.567271011`=NULL
#UF_PILOT4_3_withexcl$`2.355572113`=NULL
#UF_PILOT4_3_withexcl$`3.283913914`=NULL
#UF_PILOT4_3_withexcl$`3.548012672`=NULL
#UF_PILOT4_3_withexcl$`3.021469175`=NULL
#UF_PILOT4_3_withexcl$`3.916752776`=NULL
#UF_PILOT4_3_withexcl$`3.346299278`=NULL
#UF_PILOT4_3_withexcl$`3.136864883`=NULL
#UF_PILOT4_3_withexcl$`3.177317058`=NULL
jack=UF_PILOT4_3_withexcl
require("plyr")
plyr::count(jack, c("GROUP"))
jack=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="QC"),]))>=26  | colSums(is.na(jack[which(jack$GROUP=="HV"),]))>=19)]
pearl=jack#[,-which(colSums(is.na(jack[which(jack$GROUP=="POST"),]))>=7 & colSums(is.na(jack[which(jack$GROUP=="PPGL"),]))>=7)]
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$`ENSAT-HT ID`, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.996102403`)
sparrow$`5.996102403`<- NULL
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA.csv", sparrow)

#PQN
james=sparrow[,3:ncol(sparrow)]
james[james==0]=NA
black=sparrow[sparrow$`captain$GROUP`=="HV",3:ncol(sparrow)]
dauntless=apply(black, 2, median, na.rm=T)
check=t(t(james)/dauntless)
check2=apply(check, 1, median, na.rm=T)
checkmate=james/check2
checkmate[is.na(checkmate)]=0
beckett=cbind(sparrow[,1:2], checkmate)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN.csv", beckett)

#remove peaks with a RSD>0.3 in replicate samples (in this case both PA, PPGL)
cutler=beckett
library(sjstats)
hay=function(x) {
  cv(cutler[which(cutler$`captain$GROUP`=="QC"),3:(as.numeric(ncol(cutler)))][,x])
}
cvQC=sapply((1:(as.numeric(ncol(cutler))-2)), hay)
norrington=checkmate[,which(cvQC<0.3)]
commodore=cbind(cutler[,1:2], norrington)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3.csv", commodore)
#before exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"
#after exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=cbind(commodore[,1:2], ncheckO.imputed$data)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3_knn10.csv", commodore)

#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(ENP_vanilla_forGLOG$`captain$GROUP`=="QC"),]
#PLASMA_QC_PQN=PLASMA_QC_PQN[-c(92:105, 148:166),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCled_PQN=QCled_PQN[-2,]
colnames(QCled_PQN)=QCled_PQN[1,]
QCled_PQN=QCled_PQN[-1,]
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2)),1))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
led_PQN=led_PQN[-2,]
colnames(led_PQN)=led_PQN[1,]
led_PQN=led_PQN[-1,]
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(cutler[,1:2], final[,1:as.numeric(ncol(final))])
final=final[,-1]
#write.csv(final, file="ENP_final_3_bothQCexp_16072019_excl_final_initial_sample_exclusions_20_MA_PQN_0.3_knn10_GLOG.csv")

#robust PCA
library(rospca)
X=final[which(final$`captain$GROUP`=='QC'),-1]
res=robpca(X, k=0, skew = F)
diagPlot((res))
which(res$flag.all==F)
EXCL_QC=which(res$flag.all==F)

library(rospca)
X=final[which(final$`captain$GROUP`=='HV'),-1]
res=robpca(X, k=0, skew = F)
diagPlot((res))
which(res$flag.all==F)
EXCL_HV=which(res$flag.all==F)

EXCL=c(EXCL_QC, EXCL_HV)


ix <- which(UF_PILOT4_3_withexcl$`ENSAT-HT ID` %in% c(names(EXCL)))
jack=UF_PILOT4_3_withexcl[-ix,]
require("plyr")
plyr::count(jack, c("GROUP"))
jack=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="QC"),]))>=21  | colSums(is.na(jack[which(jack$GROUP=="HV"),]))>=15)]
pearl=jack#[,-which(colSums(is.na(jack[which(jack$GROUP=="POST"),]))>=7 & colSums(is.na(jack[which(jack$GROUP=="PPGL"),]))>=7)]
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$`ENSAT-HT ID`, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.996102403`)
sparrow$`5.996102403`<- NULL

#PQN
james=sparrow[,3:ncol(sparrow)]
james[james==0]=NA
black=sparrow[sparrow$`captain$GROUP`=="HV",3:ncol(sparrow)]
dauntless=apply(black, 2, median, na.rm=T)
check=t(t(james)/dauntless)
check2=apply(check, 1, median, na.rm=T)
checkmate=james/check2
checkmate[is.na(checkmate)]=0
beckett=cbind(sparrow[,1:2], checkmate)

#remove peaks with a RSD>0.3 in replicate samples (in this case both PA, PPGL)
cutler=beckett
library(sjstats)
hay=function(x) {
  cv(cutler[which(cutler$`captain$GROUP`=="QC"),3:(as.numeric(ncol(cutler)))][,x])
}
cvQC=sapply((1:(as.numeric(ncol(cutler))-2)), hay)
norrington=checkmate[,which(cvQC<0.3)]
commodore=cbind(cutler[,1:2], norrington)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3.csv", commodore)
#before exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"
#after exclusion: "1.860408043" "1.884843241" "1.9553774"   "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "3.995651938" "4.102385996" "4.115225666" "5.171586437"

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=cbind(commodore[,1:2], ncheckO.imputed$data)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3_knn10.csv", commodore)

#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(ENP_vanilla_forGLOG$`captain$GROUP`=="QC"),]
#PLASMA_QC_PQN=PLASMA_QC_PQN[-c(92:105, 148:166),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCled_PQN=QCled_PQN[-2,]
colnames(QCled_PQN)=QCled_PQN[1,]
QCled_PQN=QCled_PQN[-1,]
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2))))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
led_PQN=led_PQN[-2,]
colnames(led_PQN)=led_PQN[1,]
led_PQN=led_PQN[-1,]
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(cutler[,1:2], final[,1:as.numeric(ncol(final))])
final=final[,-1]
row.names(final)=commodore[,1]
#write.csv(final, file="PVP_final_no2.02_noFRPAHVs.csv")

ENP2=final[-which(final$`captain$GROUP`=="HV" | final$`captain$GROUP` =="QC"),-1]
ENP1=read_csv("ENP_final_3_new_FINAL_peaks_excluded_new2_new.csv")
ENP1=ENP1[which(ENP1$`ENSAT-HT ID` %in% row.names(ENP2)),]

ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = F)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]

PAVPHT1=list(mitch, all.BER, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 7, type.measure = "class", standardize = F)
    #fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F, lambda = cvfit1$lambda.min)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=FALSE)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich

PAVPHT3=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)


ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]

sex=as.factor(WHYF$GENDER)
levels(sex)=c("0", "1")
levels(sex)=as.numeric(levels(sex))
sex=as.numeric(sex)
sex[sex==1]=0
sex[sex==2]=1
age=as.numeric(WHYF$`PATIENT AGE`)
colnames(ex)=round(as.numeric(colnames(ex)),3)
ex=cbind(ex, sex, age)

X=ex
Y=WHYF$GROUP

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = T)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
plsda.select <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = T)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]

PAVPHT5=list(mitch, all.BER, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 7, type.measure = "class")
    #fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F, lambda = cvfit1$lambda.min)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich

PAVPHT7=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA" & clus==0)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA" & clus==0),])))
    folds2=split(sample(which(WHYF$GROUP=="PA" & clus==1)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA" & clus==1),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT" & clus==0))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT" & clus==0)),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT" & clus==1))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT" & clus==1)),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)

WHYF=cbind(WHYF, clus)

hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    clus=as.factor(jack$clus)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    enp.fac=as.factor(clus)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    cl=as.matrix(enp.fac)
    cl=as.numeric(cl)
    cl=as.matrix(cl)
    colnames(cl)=c("cl")
    WHYFA=cbind(group, cl)
    X=as.matrix(X)
    ASCA_UB = AnovaUnbalancedDataUpdated(X, WHYFA, Interact = c("12"), PREINT =T)
    DI<- ASCA_UB$BalDesMat[,grepl("&", colnames(ASCA_UB$BalDesMat))]
    BI<- ASCA_UB$BetaMAT[grepl("&", colnames(ASCA_UB$BalDesMat)),]
    XI<-DI%*%BI
    ex3=X-XI
    X=ex3
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    clus=as.factor(jack$clus)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    enp.fac=as.factor(clus)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    cl=as.matrix(enp.fac)
    cl=as.numeric(cl)
    cl=as.matrix(cl)
    colnames(cl)=c("cl")
    FA=cbind(group, cl)
    Interact= c("12")
    FAC<-c()
    for (i in 1:ncol(FA)) 
    {
      nrn<-paste(rep(colnames(FA)[i],nrow(FA)),FA[,i],sep =":")
      FAC<-cbind(FAC,nrn)
    }
    colnames(FAC)<-colnames(FA)
    
    Int.FA<-strsplit(Interact,split = "")
    INT<-c()
    for (i in 1:length(Int.FA)) 
    {
      nrn<-c()
      for(j in 1:nrow(FAC))
      {
        nrn<-rbind(nrn,paste(FAC[j,as.numeric(Int.FA[[i]])],collapse ="&"))
      }
      colnames(nrn)<-paste(colnames(FAC)[as.numeric(Int.FA[[i]])],collapse ="&")
      INT<-cbind(INT,nrn)
    }
    
    ALL.FAC<-cbind(FAC,INT)
    
    ###recode the factors
    D<-c()
    for(j in 1:ncol(ALL.FAC))
    {
      D<-cbind(D,decodeClassLabels(ALL.FAC[,j]))
    }
    DES<-cbind(rep(1,nrow(X)),D)
    colnames(DES)[1]<-"OVERALL"
    
    #######colnames for intraction and factors
    ##factors
    COL.names.FA<-vector("list",length =ncol(FAC) )
    COL.names1<-c()
    for(i in 1:ncol(FAC)) 
    {
      COL.names.FA[[i]]<-unique(FAC[,i])
      COL.names1<-c(COL.names1,unique(FAC[,i]))
    }
    ##interactions
    COL.names.INT<-vector("list",length(Int.FA))
    COL.names2<-c()
    for (i in 1:length(Int.FA)) 
    {
      L<- vector("list",length(Int.FA[[i]]))
      for (j in 1:length(Int.FA[[i]])) 
      {
        L[[j]]<-unique(FAC[,as.numeric(Int.FA[[i]][j])])
      }
      L<-as.matrix(expand.grid(L))
      
      COLN<-c()
      for (j in 1:nrow(L)) 
      {
        COLN<-c(COLN,paste(c(L[j,]),collapse ="&"))
      }
      
      COL.names.INT[[i]]<-COLN
      COL.names2<-c(COL.names2,COLN)
    }
    ##colnames
    COL.names<-c("OVERALL",COL.names1,COL.names2)
    
    ##design matrix
    DES.MAT<-array(0,c(nrow(DES),length(COL.names)))
    colnames(DES.MAT)<-COL.names
    for (j in 1:ncol(DES.MAT)) 
    {
      if(colnames(DES.MAT)[j]%in%colnames(DES))
      {
        DES.MAT[,j]<-DES[,which(colnames(DES)%in%colnames(DES.MAT)[j])]
      }else{
        DES.MAT[,j]<-DES.MAT[,j]
      }
    }
    
    ################################################
    ###Balanced Design matrix#######################
    BD.MAT<-array(1,dim = dim(DES.MAT))
    for (j in 2:ncol(DES.MAT)) 
    {
      n1<-length(which(DES.MAT[,j]==1))
      n0<-length(which(DES.MAT[,j]==0))
      BD.MAT[which(DES.MAT[,j]==1),j]<- n0
      BD.MAT[which(DES.MAT[,j]==0),j]<- -n1
    }
    colnames(BD.MAT)<-colnames(DES.MAT)
    DI2<- BD.MAT[,grepl("&", colnames(BD.MAT))]
    XI<-DI2%*%BI
    ex3=X-XI
    X=ex3
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

jack=WHYF
ENP=ex
dataset_train=ENP
group=as.factor(jack$GROUP)
clus=as.factor(jack$clus)
ex1=as.matrix(ENP)
X=ex
Y=group
enp.fac=as.factor(Y)
levels(enp.fac)=c("0", "1")
levels(enp.fac)=as.numeric(levels(enp.fac))
group=as.matrix(enp.fac)
group=as.numeric(group)
group=as.matrix(group)
colnames(group)=c("GROUP")
Y=as.factor(group)
clus=rep(1,as.numeric(nrow(WHYF)))
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="GBGL")]=0
clus=as.factor(clus)
levels(clus)=as.numeric(levels(clus))
clus=as.matrix(clus)
clus=as.numeric(clus)
clus=as.matrix(clus)
colnames(clus)=c("clus")
enp.fac=as.factor(clus)
levels(enp.fac)=as.numeric(levels(enp.fac))
cl=as.matrix(enp.fac)
cl=as.numeric(cl)
cl=as.matrix(cl)
colnames(cl)=c("cl")
WHYFA=cbind(group, cl)
X=as.matrix(X)
set.seed(1)
ASCA_UB = AnovaUnbalancedDataUpdated(X, WHYFA, Interact = c("12"), PREINT =T)
DI<- ASCA_UB$BalDesMat[,grepl("&", colnames(ASCA_UB$BalDesMat))]
BI<- ASCA_UB$BetaMAT[grepl("&", colnames(ASCA_UB$BalDesMat)),]
XI<-DI%*%BI
ex3=X-XI
ex=ex3
colnames(ex)=round(as.numeric(colnames(ex)),3)
Y=WHYF$GROUP

plsda.select <- plsda(ex, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(ex, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = F)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(ex, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]


PAVPHT9=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

ex=X
f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA" & clus==0)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA" & clus==0),])))
    folds2=split(sample(which(WHYF$GROUP=="PA" & clus==1)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA" & clus==1),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT" & clus==0))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT" & clus==0)),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT" & clus==1))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT" & clus==1)),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)

WHYF=cbind(WHYF, clus)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    clus=as.factor(jack$clus)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    enp.fac=as.factor(clus)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    cl=as.matrix(enp.fac)
    cl=as.numeric(cl)
    cl=as.matrix(cl)
    colnames(cl)=c("cl")
    WHYFA=cbind(group, cl)
    X=as.matrix(X)
    ASCA_UB = AnovaUnbalancedDataUpdated(X, WHYFA, Interact = c("12"), PREINT =T)
    DI<- ASCA_UB$BalDesMat[,grepl("&", colnames(ASCA_UB$BalDesMat))]
    BI<- ASCA_UB$BetaMAT[grepl("&", colnames(ASCA_UB$BalDesMat)),]
    XI<-DI%*%BI
    ex3=X-XI
    X=ex3
    Y=as.factor(group)
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 7, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    clus=as.factor(jack$clus)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    enp.fac=as.factor(clus)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    cl=as.matrix(enp.fac)
    cl=as.numeric(cl)
    cl=as.matrix(cl)
    colnames(cl)=c("cl")
    FA=cbind(group, cl)
    Interact= c("12")
    FAC<-c()
    for (i in 1:ncol(FA)) 
    {
      nrn<-paste(rep(colnames(FA)[i],nrow(FA)),FA[,i],sep =":")
      FAC<-cbind(FAC,nrn)
    }
    colnames(FAC)<-colnames(FA)
    
    Int.FA<-strsplit(Interact,split = "")
    INT<-c()
    for (i in 1:length(Int.FA)) 
    {
      nrn<-c()
      for(j in 1:nrow(FAC))
      {
        nrn<-rbind(nrn,paste(FAC[j,as.numeric(Int.FA[[i]])],collapse ="&"))
      }
      colnames(nrn)<-paste(colnames(FAC)[as.numeric(Int.FA[[i]])],collapse ="&")
      INT<-cbind(INT,nrn)
    }
    
    ALL.FAC<-cbind(FAC,INT)
    
    ###recode the factors
    D<-c()
    for(j in 1:ncol(ALL.FAC))
    {
      D<-cbind(D,decodeClassLabels(ALL.FAC[,j]))
    }
    DES<-cbind(rep(1,nrow(X)),D)
    colnames(DES)[1]<-"OVERALL"
    
    #######colnames for intraction and factors
    ##factors
    COL.names.FA<-vector("list",length =ncol(FAC) )
    COL.names1<-c()
    for(i in 1:ncol(FAC)) 
    {
      COL.names.FA[[i]]<-unique(FAC[,i])
      COL.names1<-c(COL.names1,unique(FAC[,i]))
    }
    ##interactions
    COL.names.INT<-vector("list",length(Int.FA))
    COL.names2<-c()
    for (i in 1:length(Int.FA)) 
    {
      L<- vector("list",length(Int.FA[[i]]))
      for (j in 1:length(Int.FA[[i]])) 
      {
        L[[j]]<-unique(FAC[,as.numeric(Int.FA[[i]][j])])
      }
      L<-as.matrix(expand.grid(L))
      
      COLN<-c()
      for (j in 1:nrow(L)) 
      {
        COLN<-c(COLN,paste(c(L[j,]),collapse ="&"))
      }
      
      COL.names.INT[[i]]<-COLN
      COL.names2<-c(COL.names2,COLN)
    }
    ##colnames
    COL.names<-c("OVERALL",COL.names1,COL.names2)
    
    ##design matrix
    DES.MAT<-array(0,c(nrow(DES),length(COL.names)))
    colnames(DES.MAT)<-COL.names
    for (j in 1:ncol(DES.MAT)) 
    {
      if(colnames(DES.MAT)[j]%in%colnames(DES))
      {
        DES.MAT[,j]<-DES[,which(colnames(DES)%in%colnames(DES.MAT)[j])]
      }else{
        DES.MAT[,j]<-DES.MAT[,j]
      }
    }
    
    ################################################
    ###Balanced Design matrix#######################
    BD.MAT<-array(1,dim = dim(DES.MAT))
    for (j in 2:ncol(DES.MAT)) 
    {
      n1<-length(which(DES.MAT[,j]==1))
      n0<-length(which(DES.MAT[,j]==0))
      BD.MAT[which(DES.MAT[,j]==1),j]<- n0
      BD.MAT[which(DES.MAT[,j]==0),j]<- -n1
    }
    colnames(BD.MAT)<-colnames(DES.MAT)
    DI2<- BD.MAT[,grepl("&", colnames(BD.MAT))]
    XI<-DI2%*%BI
    ex3=X-XI
    X=ex3
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)


jack=WHYF
ENP=ex
dataset_train=ENP
group=as.factor(jack$GROUP)
clus=as.factor(jack$clus)
ex1=as.matrix(ENP)
X=ex
Y=group
enp.fac=as.factor(Y)
levels(enp.fac)=c("0", "1")
levels(enp.fac)=as.numeric(levels(enp.fac))
group=as.matrix(enp.fac)
group=as.numeric(group)
group=as.matrix(group)
colnames(group)=c("GROUP")
Y=as.factor(group)
clus=rep(1,as.numeric(nrow(WHYF)))
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="GBGL")]=0
clus=as.factor(clus)
levels(clus)=as.numeric(levels(clus))
clus=as.matrix(clus)
clus=as.numeric(clus)
clus=as.matrix(clus)
colnames(clus)=c("clus")
enp.fac=as.factor(clus)
levels(enp.fac)=as.numeric(levels(enp.fac))
cl=as.matrix(enp.fac)
cl=as.numeric(cl)
cl=as.matrix(cl)
colnames(cl)=c("cl")
WHYFA=cbind(group, cl)
X=as.matrix(X)
set.seed(1)
ASCA_UB = AnovaUnbalancedDataUpdated(X, WHYFA, Interact = c("12"), PREINT =T)
DI<- ASCA_UB$BalDesMat[,grepl("&", colnames(ASCA_UB$BalDesMat))]
BI<- ASCA_UB$BetaMAT[grepl("&", colnames(ASCA_UB$BalDesMat)),]
XI<-DI%*%BI
ex3=X-XI
ex=ex3
colnames(ex)=round(as.numeric(colnames(ex)),3)
Y=WHYF$GROUP

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=FALSE)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich


PAVPHT11=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)



ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
ex=ex[-which(WHYF$`CENTER ID`=="FRPA" & WHYF$GROUP=="PHT"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" & WHYF$GROUP=="PHT"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)
Original_PCA=mixOmics::pca(X, scale = F, center = T, ncomp = 10)
plotIndiv(Original_PCA, group = Y, legend = T, style = "3d", pch="sphere")
plotVar(Original_PCA)
clus=rep(1,as.numeric(nrow(WHYF)))
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP=="PA"&WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="GBGL")]=0
clus=as.factor(clus)
levels(clus)=as.numeric(levels(clus))
clus=as.matrix(clus)
clus=as.numeric(clus)
clus=as.matrix(clus)
colnames(clus)=c("clus")
length(which(clus=="0"))
length(which(clus=="1"))
plotIndiv(Original_PCA, group = clus, legend = T, style = "3d", pch="sphere")

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(clus=="0")), rep(1:f, length=nrow(X[which(clus=="0"),])))
    folds2=split(sample(which(clus=="1")), rep(1:f, length=nrow(X[which(clus=="1"),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=clus[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    plsda.train <- plsda(X, Y, ncomp = ncomp, scale = F)
    jack=clus[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


allconf=Reduce('+', allconf)
((allconf[2,2]/(allconf[2,2]+allconf[2,1]))+(allconf[1,1]/(allconf[1,1]+allconf[1,2])))/2


vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]][,as.numeric(ncol(mitch[[x]][[y]][[3]]))]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

allv=function(x){
  do.call(cbind, allvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checker=checke2[order(checke2)]
checker[which(checker>1)]

C1VC21=list(all.BER, checker, checke2)

ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)
Original_PCA=mixOmics::pca(X, scale = F, center = T, ncomp = 10)
plotIndiv(Original_PCA, group = Y, legend = T, style = "3d", pch="sphere")
plotVar(Original_PCA)
SA=as.numeric(WHYF$`SAMPLE AGE`)
SA[SA<median(SA, na.rm=T)]=0
SA[!SA==0]=1
enp.fac=as.factor(SA)
SA=as.matrix(enp.fac)
SA=as.numeric(SA)
SA=as.matrix(SA)
colnames(SA)=c("SA")
length(which(SA=="0"))
length(which(SA=="1"))
plotIndiv(Original_PCA, group = SA, legend = T, style = "3d", pch="sphere")

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(SA=="0")), rep(1:f, length=nrow(X[which(SA=="0"),])))
    folds2=split(sample(which(SA=="1")), rep(1:f, length=nrow(X[which(SA=="1"),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=SA[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    plsda.train <- plsda(X, Y, ncomp = ncomp, scale = F)
    jack=SA[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


allconf=Reduce('+', allconf)
((allconf[2,2]/(allconf[2,2]+allconf[2,1]))+(allconf[1,1]/(allconf[1,1]+allconf[1,2])))/2


vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]][,as.numeric(ncol(mitch[[x]][[y]][[3]]))]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

allv=function(x){
  do.call(cbind, allvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checker=checke2[order(checke2)]
checker[which(checker>1)]

C1VC23=list(all.BER, checker, checke2)

ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
ex=ex[-which((WHYF$`CENTER ID`=="FRPA" & WHYF$GROUP=="PA")|WHYF$`CENTER ID`=="GBGL"|WHYF$`CENTER ID`=="ITPD"),]
WHYF=WHYF[-which((WHYF$`CENTER ID`=="FRPA" & WHYF$GROUP=="PA")|WHYF$`CENTER ID`=="GBGL"|WHYF$`CENTER ID`=="ITPD"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)
Original_PCA=mixOmics::pca(X, scale = F, center = T, ncomp = 10)
plotIndiv(Original_PCA, group = WHYF$`CENTER ID`, legend = T, style = "3d", pch="sphere")
plotVar(Original_PCA)


ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
ex=ex[which(WHYF$`CENTER ID`=="GYDR" | WHYF$`CENTER ID`=="ITTU"),]
WHYF=WHYF[which(WHYF$`CENTER ID`=="GYDR" | WHYF$`CENTER ID`=="ITTU"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 6, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 6, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = F)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]

PAVPHT_GYDRITTU1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 5, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 6, type.measure = "class", standardize=FALSE)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich

PAVPHT_GYDRITTU3=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)



ex=ENP2[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
WHYF=ENP1[which(ENP1$GROUP=="PA" | ENP1$GROUP=="PHT"),]
ex=ex[which(WHYF$`CENTER ID`=="GYDR" | WHYF$`CENTER ID`=="ITTU"),]
WHYF=WHYF[which(WHYF$`CENTER ID`=="GYDR" | WHYF$`CENTER ID`=="ITTU"),]

sex=as.factor(WHYF$GENDER)
levels(sex)=c("0", "1")
levels(sex)=as.numeric(levels(sex))
sex=as.numeric(sex)
sex[sex==1]=0
sex[sex==2]=1
age=as.numeric(WHYF$`PATIENT AGE`)
colnames(ex)=round(as.numeric(colnames(ex)),3)
ex=cbind(ex, sex, age)

X=ex
Y=WHYF$GROUP

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = T)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
plsda.select <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 6, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 6, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = T)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]

PAVPHT_GYDRITTU5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 5, type.measure = "class", standardize=T)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 6, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich


PAVPHT_GYDRITTU7=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)

