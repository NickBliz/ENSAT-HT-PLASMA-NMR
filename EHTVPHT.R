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


ENP2=final[-which(final$`captain$GROUP`=="HV" | final$`captain$GROUP` =="QC"),-1]
ENP1=read_csv("ENP_final_3_new_FINAL_peaks_excluded_new2_new.csv")             #load metadata
ENP1=ENP1[which(ENP1$`ENSAT-HT ID` %in% row.names(ENP2)),]                    #match metadata

ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF$GROUP2=WHYF$GROUP
WHYF$GROUP[which(WHYF$GROUP=="CS" | WHYF$GROUP=="PA" |WHYF$GROUP=="PPGL")]="AHT"
X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)
Original_PCA=mixOmics::pca(X, scale = F, center = T, ncomp = 10)
group2=as.factor(WHYF$GROUP2)
group1=as.factor(WHYF$GROUP)
levels(group1)=c("EHT", "PHT")
plotIndiv(Original_PCA, group = group1, legend = T, pch=c(19), cex = 5, title="(a): PCA DISEASE GROUP")
plotIndiv(Original_PCA, group = group2, legend = T, pch=c(19), cex = 5)
cluster_plot=plotIndiv(Original_PCA, group = cluster, ellipse = T, legend = T, pch=c(19), cex=5, title="(b): PCA SAMPLE CENTER OF ORIGIN")
center_plot=plotIndiv(Original_PCA, group = WHYF$`CENTER ID`, legend = T, pch=c(19), cex=5, title="(b): PCA SAMPLE CENTER OF ORIGIN")
center_plot=plotIndiv(Original_PCA, group = WHYF$`CENTER ID`, col.per.group = c("darkblue", "darkviolet", "blue", "darkgrey", "cyan", "azure4", "burlywood3", "gold", "darkgoldenrod1", "gold4", "darkseagreen", "palegreen", "black"), legend = T, pch=c(19), cex=5, title="(b): PCA SAMPLE CENTER OF ORIGIN")
center_plot$graph$layers[[14]]=cluster_plot$graph$layers[[3]]
center_plot$graph$layers[[15]]=cluster_plot$graph$layers[[4]]
center_plot
pV=plotVar(Original_PCA, cutoff = 0.5)

sa=as.numeric(WHYF$`SAMPLE AGE`)
kruskal.test(sa~WHYF$`CENTER ID`)
#p-value < 2.2e-16
#PPGL GYDRITTU: p-value = 7.072e-10
sa[which(sa<median(sa))]="<median"
sa[-which(sa=="<median")]=">=median"
pI=plotIndiv(Original_PCA, group = sa, legend = T, pch=c(19), cex = 5, title="(c): PCA SAMPLE AGE")
cor.test(pI$df$x, sa, method = "kendall")
#p-value = 0.0001889
#gydrittu: p-value = 0.04014
#PPGL: p-value = 8.415e-13
cor.test(pI$df$x, sa, method = "spearman")
#p-value = 0.0001845
#gydrittu: p-value = 0.0288
#PPGL: p-value = 5.245e-14
cor.test(pI$df$y, sa, method = "kendall")
#p-value = 1.517e-13
#excl: p-value = 0.06625
#gydrittu: p-value = 4.896e-06
#PPGL: p-value = 0.0004957
#PPGL EXCL: p-value = 0.002942
#PPGL GYDRITTU: p-value = 0.01303
cor.test(pI$df$y, sa, method = "spearman")
#p-value = 7.629e-14
#excl: p-value = 0.07204
#gydrittu: p-value = 2.406e-06
#PPGL: p-value = 9.628e-05
#PPGL EXCL: p-value = 0.003378
#PPGL GYDRITTU: p-value = 0.01186
SEQ=as.numeric(WHYF$`SEQUENCE NUMBER`)
kruskal.test(SEQ~WHYF$`CENTER ID`)
SEQ[which(SEQ<median(SEQ))]="<median"
SEQ[-which(SEQ=="<median")]=">=median"
pI=plotIndiv(Original_PCA, group = SEQ, legend = T, pch=c(19), cex = 5, title="(c): PCA RUN ORDER")
cor.test(pI$df$x, SEQ, method = "kendall")
cor.test(pI$df$x, SEQ, method = "spearman")
cor.test(pI$df$y, SEQ, method = "kendall")
cor.test(pI$df$y, SEQ, method = "spearman")
RUN=as.numeric(WHYF$RUN)
kruskal.test(RUN~WHYF$`CENTER ID`)
#p-value = 7.513e-10
#excl: p-value = 7.513e-10
#gydrittu: p-value = 7.112e-08
#PPGL: p-value = 4.662e-11
#PPGL GYDR-ITTU: p-value = 1.891e-06
RUN[which(RUN<median(RUN))]="<median"
RUN[-which(RUN=="<median")]=">=median"
pI=plotIndiv(Original_PCA, group = RUN, legend = T, pch=c(19), cex = 5, title="(c): PCA BATCH")
cor.test(pI$df$x, RUN, method = "kendall")
#gydrittu: p-value = 0.005102
#IN PPGL-PHT: p-value = 0.07233
#PPGL GYDRITTU: p-value = 0.002841
cor.test(pI$df$x, RUN, method = "spearman")
#gydrittu: p-value = 0.006237
#IN PPGL-PHT: p-value = 0.05966
#PPGL GYDRITTU: p-value = 0.003859
cor.test(pI$df$y, RUN, method = "kendall")
#p-value = 0.05031
#gydrittu: p-value = 0.01047
#PPGL EXCL: p-value = 0.09404
cor.test(pI$df$y, RUN, method = "spearman")
#p-value = 0.06207
#gydrittu: p-value = 0.01146
#PPGL EXCL: p-value = 0.08891

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_AHTVPHT1=test.predict$B.hat[,,ncomp]

AHTVPHT1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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

each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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

AHTVPHT3=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA" | ENP1$GROUP=="CS"),]
WHYF$GROUP2=WHYF$GROUP
WHYF$GROUP[which(WHYF$GROUP=="CS" | WHYF$GROUP=="PA" |WHYF$GROUP=="PPGL")]="AHT"

X=ex
Y=WHYF$GROUP
Original_PCA=mixOmics::pca(X, scale = F, center = T, ncomp = 10)
plotIndiv(Original_PCA, group = Y, legend = T, style = "3d", pch="sphere")

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
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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

AHTVPHT5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)

X=as.matrix(ex)
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich

AHTVPHT7=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)

clus=rep(1,as.numeric(nrow(WHYF)))
clus[which(WHYF$GROUP2=="PPGL"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="CS"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="IRGA")]=0
clus[which(WHYF$`CENTER ID`=="GBGL")]=0
clus=as.factor(clus)
levels(clus)=as.numeric(levels(clus))
clus=as.matrix(clus)
clus=as.numeric(clus)
clus=as.matrix(clus)
colnames(clus)=c("clus")


f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT" & clus==0)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT" & clus==0),])))
    folds2=split(sample(which(WHYF$GROUP=="AHT" & clus==1)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT" & clus==1),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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
clus[which(WHYF$GROUP2=="PPGL"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="CS"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="IRGA")]=0
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

AHTVPHT9=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


ex=X
f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT" & clus==0)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT" & clus==0),])))
    folds2=split(sample(which(WHYF$GROUP=="AHT" & clus==1)), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT" & clus==1),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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
clus[which(WHYF$GROUP2=="PPGL"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="PA"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$GROUP2=="CS"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$`CENTER ID`=="ITPD")]=0
clus[which(WHYF$`CENTER ID`=="IRGA")]=0
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
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich

AHTVPHT11=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF$GROUP2=WHYF$GROUP
WHYF$GROUP[which(WHYF$GROUP=="PA" |WHYF$GROUP=="PPGL")]="AHT"

X=ex
Y=WHYF$GROUP
colnames(X)=round(as.numeric(colnames(X)),3)

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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

AHTVPHT_GYDRITTU1=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)



MC.mitch=function(x) {
  fold.MC=function(y){
    mitch[[x]][[y]][[8]]
  }
  fold.MCs=lapply((1:f), fold.MC)
}

allMC=lapply((1:50), MC.mitch)


allMCC=function(x){
  do.call(rbind, allMC[[x]])
}

allMCCC=lapply(1:50, allMCC)
allMCCC=do.call(rbind, allMCCC)
nMC=plyr::count(allMCCC, "`ENSAT-HT ID`")
ch=WHYF
WHYF=WHYF[order(WHYF$`ENSAT-HT ID`),]
WHYF$MC=0
WHYF$MC[which(WHYF$`ENSAT-HT ID` %in% nMC$ENSAT.HT.ID)]=nMC$freq
wilcox.test(WHYF$MC[which(WHYF$GENDER=="MALE")], WHYF$MC[which(WHYF$GENDER=="FEMALE")])
t.test(WHYF$MC[which(WHYF$GENDER=="MALE")], WHYF$MC[which(WHYF$GENDER=="FEMALE")])
wilcox.test(WHYF$MC[which(SA=="0")], WHYF$MC[which(SA=="1")])
wilcox.test(WHYF$MC[which(RUN=="0")], WHYF$MC[which(RUN=="1")])
wilcox.test(WHYF$MC[which(SEQ=="0")], WHYF$MC[which(SEQ=="1")])
wilcox.test(WHYF$MC[which(AGE=="0")], WHYF$MC[which(AGE=="1")])
wilcox.test(WHYF$MC[which(clus=="0")], WHYF$MC[which(clus=="1")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="GBGL")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="FRPA")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")])
t.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")], WHYF$MC[which(WHYF$`CENTER ID`=="GBGL")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")], WHYF$MC[which(WHYF$`CENTER ID`=="FRPA")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="FRPA")], WHYF$MC[which(WHYF$`CENTER ID`=="GBGL")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")])
t.test(WHYF$MC[which(WHYF$`CENTER ID`=="ITTU")], WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")], WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")])
t.test(WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")], WHYF$MC[which(WHYF$`CENTER ID`=="GYDR")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")], WHYF$MC[which(WHYF$`CENTER ID`=="GBGL")])
wilcox.test(WHYF$MC[which(WHYF$`CENTER ID`=="PLWW")], WHYF$MC[which(WHYF$`CENTER ID`=="FRPA")])
wilcox.test(WHYF$MC[which(WHYF$GROUP2=="PA")], WHYF$MC[which(WHYF$GROUP2=="PPGL")])
wilcox.test(WHYF$MC[which(WHYF$GROUP2=="PA")], WHYF$MC[which(WHYF$GROUP2=="CS")])
wilcox.test(WHYF$MC[which(WHYF$GROUP2=="CS")], WHYF$MC[which(WHYF$GROUP2=="PPGL")])

mc=as.numeric(WHYF$MC)
age=as.numeric(WHYF$`PATIENT AGE`)
cor.test(age, mc, method = "pearson")
cor.test(age, mc, method = "kendall")
cor.test(age, mc, method = "spearman")
wilcox.test(age[which(mc==0)], age[which(mc>0)])
t.test(age[which(mc==0)], age[which(mc>0)])
sa=as.numeric(WHYF$`SAMPLE AGE`)
cor.test(sa, mc, method = "pearson")
cor.test(sa, mc, method = "kendall")
cor.test(sa, mc, method = "spearman")
wilcox.test(sa[which(mc==0)], sa[which(mc>0)])
t.test(sa[which(mc==0)], sa[which(mc>0)])
run=as.numeric(WHYF$RUN)
cor.test(run, mc, method = "pearson")
cor.test(run, mc, method = "kendall")
cor.test(run, mc, method = "spearman")
wilcox.test(run[which(mc==0)], run[which(mc>0)])
t.test(run[which(mc==0)], run[which(mc>0)])
seq=as.numeric(WHYF$`SEQUENCE NUMBER`)
cor.test(seq, mc, method = "pearson")
cor.test(seq, mc, method = "kendall")
cor.test(seq, mc, method = "spearman")
wilcox.test(seq[which(mc==0)], seq[which(mc>0)])
t.test(seq[which(mc==0)], seq[which(mc>0)])
WHYF=ch

CC.mitch=function(x) {
  fold.CC=function(y){
    mitch[[x]][[y]][[9]]
  }
  fold.CCs=lapply((1:f), fold.CC)
}

allCC=lapply((1:50), CC.mitch)


X=ex
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 6, type.measure = "class", standardize=F)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich[order(rich)])
rich


AHTVPHT_GYDRITTU3=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)


ex=ENP2[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]
WHYF=ENP1[which(ENP1$GROUP=="PPGL" | ENP1$GROUP=="PHT" | ENP1$GROUP=="PA"),]

ex=ex[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF=WHYF[-which(WHYF$`CENTER ID`=="FRPA" | WHYF$`CENTER ID`=="ITPD" | WHYF$`CENTER ID`=="GBGL"),]
WHYF$GROUP2=WHYF$GROUP
WHYF$GROUP[which(WHYF$GROUP=="PA" |WHYF$GROUP=="PPGL")]="AHT"

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
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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

AHTVPHT_GYDRITTU5=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)


f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="AHT")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="AHT"),])))
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


each.AHT=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.AHT=lapply((1:n), each.AHT)
all.AHT=unlist(all.AHT)

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


AHTVPHT_GYDRITTU7=list(mitch, all.BER, checker, checke2, perm1, WHYF, allMC, allCC, cvfit1, rich)
