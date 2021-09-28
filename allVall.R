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

ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

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

allVall1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 7, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 8, type.measure = "class", standardize=FALSE)
richCS=coef(cvfit1, s = "lambda.min")$CS
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richCS, richPA, richPHT, richPPGL)
colnames(rich)=c("CS", "PA", "PHT", "PPGL")
rich=rich[-which(rowSums(rich)==0),]
rich

allVall3=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, richCS, richPA, richPHT, richPPGL)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]

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
colnames(X)=round(as.numeric(colnames(X)),3)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

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

allVall5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds3=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds4=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 7, type.measure = "class", standardize = T)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.CS=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.CS=lapply((1:n), each.CS)
all.CS=unlist(all.CS)

each.PA=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

each.PPGL=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
richCS=coef(cvfit1, s = "lambda.min")$CS
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richCS, richPA, richPHT, richPPGL)
colnames(rich)=c("CS", "PA", "PHT", "PPGL")
rich=rich[-which(rowSums(rich)==0),]
rich

allVall7=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, richCS, richPA, richPHT, richPPGL)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
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

each.PPGL=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

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
plotIndiv(splsda.select, legend=T, pch=19, ellipse = T, cex=5)
plotIndiv(Original_PCA,group=Y, legend=T, pch=19, ellipse = T, cex=5)

ALLVALL_GYDRITTU1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)



f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 5, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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

each.PPGL=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 6, type.measure = "class", standardize=F)
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richPA, richPHT, richPPGL)
colnames(rich)=c("PA", "PHT", "PPGL")
rich=rich[-which(rowSums(rich)==0),]
rich

ALLVALL_GYDRITTU3=list(mitch, all.BER, checker, checke2, WHYF, allMC, allCC, cvfit1, rich)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]

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
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
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

each.PPGL=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

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

ALLVALL_GYDRITTU5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
    folds3=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2, folds3)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "multinomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 5, type.measure = "class", standardize = T)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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

each.PPGL=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "multinomial", alpha=1, nfolds = 6, type.measure = "class", standardize=T)
richPA=coef(cvfit1, s = "lambda.min")$PA
richPHT=coef(cvfit1, s = "lambda.min")$PHT
richPPGL=coef(cvfit1, s = "lambda.min")$PPGL
rich=cbind(richPA, richPHT, richPPGL)
colnames(rich)=c("PA", "PHT", "PPGL")
rich=rich[-which(rowSums(rich)==0),]
rich


ALLVALL_GYDRITTU7=list(mitch, all.BER, checker, checke2, WHYF, allMC, allCC, cvfit1, rich)



ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD"),]

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)


f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
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

each.PPGL=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)


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



PAVPPGL1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 5, type.measure = "class", standardize = F)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.PPGL=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 6, type.measure = "class", standardize=F)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-which(rowSums(rich)==0),]
rich=as.matrix(rich)
rich=rich[-1,]
rich=rich[order(rich)]
rich=as.matrix(rich)
rich


PAVPPGL3=list(mitch, all.BER, checker, checke2, WHYF, allMC, allCC, cvfit1, rich)



ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD"),]

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
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
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

each.PPGL=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

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



PAVPPGL5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)



f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
    folds2=split(sample(which(WHYF$GROUP=="PPGL")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PPGL"),])))
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
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 5, type.measure = "class", standardize = T)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack$GROUP
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
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


each.PPGL=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PPGL=lapply((1:n), each.PPGL)
all.PPGL=unlist(all.PPGL)

X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 6, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-which(rowSums(rich)==0),]
rich=as.matrix(rich)
rich=rich[-1,]
rich=rich[order(rich)]
rich=as.matrix(rich)
rich


PAVPPGL7=list(mitch, all.BER, checker, checke2, WHYF, allMC, allCC, cvfit1, rich)

