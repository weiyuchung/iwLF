#已設定提前停止規則,通常不會迭代到10次

library(caret)
library("logicFS")
library(pacman)
pacman::p_load(tidyverse, doParallel, foreach)
library(LogicForest)
library(parallel)
library(PRROC)
library(e1071)

X <- make.snp.dummy(as.matrix(X))

###

weighted_LF3 <-
function(resp, Xs, nBSXVars, anneal.params, nBS=100, h=0.5,probs)
{
 
  tree_pre_X_in_count=rep(0,ncol(Xs))
  pred<-ncol(Xs)
  if (missing(anneal.params)) {anneal.params<-logreg.anneal.control(start=2, end=-1, iter=50000)}
  if (missing(nBSXVars)) {nBSXVars<-floor((pred)^0.7)   }


  n<-nrow(Xs)
  if (is.null(colnames(Xs))) {x.nms<-paste("X", 1:pred, sep="")}else {x.nms<-colnames(Xs)}
  fitlist<-vector("list", nBS)
  IBdata<-vector("list", nBS)
  OOBdata<-vector("list", nBS)
  OOBpred<-matrix(nrow=n, ncol=nBS)
  single.predimport<-vector("list",nBS)
  vimp.import<-vector("list", nBS)
  treepreds.list<-vector("list", nBS)
  pre_X_in.list<-vector("list", nBS)

##

##
##


getDoParWorkers()
resu=foreach(b = 1:nBS, .combine = rbind) %dopar% { 

     BSindices<-sample(1:n, n, replace=TRUE) 
     OOBindices<-(1:n)[!is.element(1:n, BSindices)]
     BS<-Xs[BSindices, ]
     OOB<-Xs[OOBindices, ]
     BSY<-resp[BSindices] 
     OOBY<-resp[OOBindices]

     XVarIndices<-sort(sample(1:pred, nBSXVars, replace=FALSE,prob=probs)) 
     tree_pre_X_in_count[XVarIndices]=tree_pre_X_in_count[XVarIndices]+1
     rsBSX<-BS[ ,XVarIndices]
     rsOOBX<-OOB[,XVarIndices]
     FinalBS<-cbind(rsBSX, BSY) 
     FinalOOB<-cbind(rsOOBX, OOBY)  
     colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
     fit <- logreg(resp = BSY, bin = FinalBS[,1:nBSXVars], 
            type = 1, select = 1, ntrees =1 , nleaves = 10,   anneal.control = anneal.params)
     oobpred<-predict.logreg(fit, newbin=as.matrix(FinalOOB[,1:nBSXVars]))
misclass=sum(abs(oobpred-OOBY))/length(OOBY)

     if (sum(fit$model$tree[[1]]$trees[,3])!=0) {
        pred.import<-pimp.import(fit=fit, data=as.data.frame(FinalBS), 
        testdata=as.data.frame(FinalOOB), BSpred=length(XVarIndices),  pred=pred, Xs=XVarIndices ) 

	list(pred.import$pimp.vimp,
	pred.import$single.vimp,
	pred.import$Xids,
	fit,BSindices,
	OOBindices,
	XVarIndices,
	oobpred,
	OOBindices,misclass)
        }else {
	list(0,
	0,
	0,
	fit,BSindices,
	OOBindices,
	XVarIndices,
	oobpred,
	OOBindices,misclass)

}

}

vimp.import=resu[,1]
names(vimp.import)=NULL
single.predimport=resu[,2]
names(single.predimport)=NULL
treepreds.list=resu[,3]
names(treepreds.list)=NULL
fitlist=resu[,4]
names(fitlist)=NULL
IBdata=resu[,5]
names(IBdata)=NULL
OOBdata=resu[,6]
names(OOBdata)=NULL
pre_X_in.list=resu[,7]
names(pre_X_in.list)=NULL
for (b in 1:nBS){
OOBpred[as.vector(resu[,9][[b]]),b]<-resu[,8][[b]]
}
OOBmisclass.list=resu[,10]
names(OOBmisclass.list)=NULL
weight_metrix=1-as.numeric(unlist(OOBmisclass.list))

  OOB.pred<-matrix(0, nrow=n, ncol=2)
  for (i in 1:n)
     {
     pred.ids<-which(OOBpred[i,]==1|OOBpred[i,]==0)
     pred.vec<-OOBpred[i,c(pred.ids)]
     OOBprop<-sum(pred.vec*weight_metrix[pred.ids])/sum(weight_metrix[pred.ids])
     OOBi.pred<-ifelse(OOBprop>h, 1, 0)
     OOB.pred[i,]<-c(OOBi.pred, OOBprop)
     }


  colnames(OOB.pred)<-c("predicted_class","proportion")
  OOBmisclass<-sum(abs(OOB.pred[,1]-resp))/n
  pred.importmat<-matrix(0, nrow=nBS, ncol=pred)
  colnames(pred.importmat)<-x.nms
  for (i in 1:nBS)
    {
    pred.ind<-treepreds.list[[i]]
    m<-length(pred.ind)
    for (j in 1:m)
      {
      col<-pred.ind[j]
      pred.importmat[i,col]<-single.predimport[[i]][j]
      }
    }
weight_metrix=1-as.matrix(unlist(OOBmisclass.list))
for (i in 1:(ncol(pred.importmat)-1)){
weight_metrix=cbind(weight_metrix,1-as.matrix(unlist(OOBmisclass.list)))
}
  pred.imp<-colSums(pred.importmat*weight_metrix)
  names(pred.imp)<-x.nms
  freq.table<-table(names(unlist(single.predimport)))
  all.pimps<-unique(names(unlist(vimp.import)))
  npimps<-length(unique(names(unlist(vimp.import))))
  pimp.freqtable<-table(names(unlist(vimp.import)))
  pimptable<-matrix(0, nrow=nBS, ncol=npimps)
  colnames(pimptable)<-all.pimps
  for (i in 1:nBS)
    {
    npimps.tree<-length(vimp.import[[i]])
    col.ids<-which(colnames(pimptable)%in%names(vimp.import[[i]]))
    n.ids<-length(col.ids)
    for (j in 1:n.ids)
      {
      pimptable[i,col.ids[j]]<- vimp.import[[i]][j]
      }
    }

MMM=1-as.matrix(unlist(OOBmisclass.list))

MM=matrix(rep(MMM,ncol(pimptable)),nrow=nrow(pimptable),ncol=ncol(pimptable))
  pimpsum<-colSums(pimptable*MM)
  t5PIs<-names(sort(pimpsum, decreasing=TRUE)[1:5]) 
  which_X_not_in_model=which(tree_pre_X_in_count==0)
  tree_pre_X_in_count[which_X_not_in_model]=1
 # pred.imp2=pred.imp/tree_pre_X_in_count
 # pred.imp3=  pred.imp2
 # pred.imp3[which_X_not_in_model]=probs[which_X_not_in_model]
 # pred.imp=  pred.imp/sum(1-as.numeric(unlist(OOBmisclass.list)))
 # pimpsum=  pimpsum/sum(1-as.numeric(unlist(OOBmisclass.list)))

  ans<-list(AllFits=fitlist, Top5.PI=t5PIs, Predictor.importance=pred.imp, Predictor.frequency=freq.table,
            PI.frequency=pimp.freqtable, PI.importance=pimpsum,  ModelPI.import=vimp.import,
            OOBmisclass=OOBmisclass, OOBprediction=OOB.pred, IBdata=IBdata, OOBdata=OOBdata, predictors=pred, Xs=Xs,
            pre_X_in.list=pre_X_in.list,OOBmisclass.list=OOBmisclass.list)
 class(ans)<-"logforest"
 ans
}
####
iLF3 <-
function(resp, Xs, nBSXVars, 
	      anneal.params, nBS=100, h=0.5, maxK=10,err_conv,dect_order=1)
{
err_list=list()
weight_list=list()
dect_list=list()
rank_list=list()
AUPRC_list=list()
f1_score_list=list()
WFP_list=list()
TN_list=list()
FN_list=list()
TP_list=list()
FP_list=list()
time_list=list()

true_inter="x1"
if(dect_order>1){
for(i in 2:dect_order){
true_inter=paste(true_inter,"&",gsub(" ","",paste("x",i)))
}
}

pre_err=1
X_weight=rep(1,ncol(Xs))
X_weight=X_weight/sum(X_weight)

for (i in 1:maxK){
weight_list[[i]]=X_weight
ptm <- proc.time()
wLF=weighted_LF3(resp=resp, Xs=Xs, nBSXVars=nBSXVars,  
	      anneal.params=anneal.params, nBS=nBS, h=h, probs=X_weight)
time_list[[i]]=(proc.time() - ptm)[3]
pred_prob <- wLF$OOBprediction [,2]
true_label <- resp
pr <- pr.curve(scores.class0 =pred_prob, weights.class0 = true_label)
AUPRC_list[[i]]=pr$auc.integral 

actual_labels <- resp
predicted_labels <- wLF$OOBprediction [,1]
conf_matrix <- confusionMatrix(factor(predicted_labels), factor(actual_labels))
TP <- conf_matrix$table[2,2]
FP <- conf_matrix$table[1,2]
FN <- conf_matrix$table[2,1]
TN <- conf_matrix$table[1,1]
FP_list[[i]]=FP
TP_list[[i]]=TP
FN_list[[i]]=FN
TN_list[[i]]=TN
precision <- TP / (TP + FP)
recall <- TP / (TP + FN)
f1_score <- 2 * ((precision * recall) / (precision + recall))
f1_score_list[[i]]=f1_score

WFP_list[[i]]=FP*sum(resp)/(length(resp)-sum(resp))/length(resp)

now_err=wLF$OOBmisclass
err_list[[i]]=now_err
if(sum(names(wLF$PI.frequency)==true_inter)>0){dect_list[[i]]=wLF$PI.frequency[true_inter]}else{dect_list[[i]]=0}
if(sum(names(wLF$PI.importance)==true_inter)>0){rank_list[[i]]=which(names(sort(wLF$PI.importance,decreasing=T))==true_inter)}else{rank_list[[i]]=0}

if((now_err-pre_err)/pre_err>0.05){
break
}
if(now_err==0){
best_model=wLF
best_k=i
break
}
if(now_err<=pre_err){
best_model=wLF
best_k=i
pre_err=now_err
}


X_weight=ifelse(wLF$Predictor.importance>=(10^(-50)),wLF$Predictor.importance,(10^(-50)))
X_weight=X_weight/sum(X_weight)
}

list(fit=best_model,err_list=err_list,weight_list=weight_list,dect_list=dect_list,rank_list=rank_list,nBS=nBS,AUPRC_list=AUPRC_list,f1_score_list=f1_score_list,FP_list=FP_list,TN_list=TN_list,TP_list=TP_list,FN_list=FN_list,
time_list=time_list,best_k=best_k)
}
#####

predict.iLF <-
function(object, newdata, cutoff=0.5,...)
 {
 if (class(object$fit)!= "logforest")
    stop("object not of class logforest")  
  nBS<-length(object$fit$AllFits)
  trees<- object$fit$AllFits
X_list=object$fit$pre_X_in.list

weight_matrix=as.matrix(1-unlist(object$fit$OOBmisclass.list ))
  h<-length(trees)
  if (missing(newdata)) {LFprediction<-object$fit$OOBprediction[,1]
      proportion_one<-object$fit$OOBprediction[,2]
      ans<-list(LFprediction=LFprediction, proportion_one=proportion_one)
      }
  if (!missing(newdata))
     {
     pred<-ncol(newdata)
     if (pred!=object$fit$predictors)
       stop("the predictors in newdata do not match the original predictors")
     size<-nrow(newdata)
     predict.new<-vector("list", nBS)
     for (i in 1:nBS)
       {   
       newX<-newdata[,1:pred]
       newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX[,object$fit$pre_X_in.list[[i]]]))
       predict.new[[i]]<- newpredict
       }
     predictlist<-unlist(predict.new)
     predictvector<-as.vector(predictlist)
     predictmatrix<-matrix(c(predictvector), nrow=size, ncol=h, byrow=FALSE)
predict.pos=predictmatrix%*%weight_matrix/sum(weight_matrix)
status=predict.pos>cutoff
predictions=cbind(predict.pos,status)


     predmatrix<-cbind(predictmatrix, predictions)
     predframe<-as.data.frame(predmatrix)
     names(predframe)[1:h]<-paste("tree", 1:h, sep="")
     names(predframe)[h+1]<-paste("proportion_one")
     names(predframe)[h+2]<-paste("prediction")
     ans<-list(LFprediction=predframe$prediction, proportion_one=predframe$proportion_one,
             AllTrees=predframe)
     } 
  class(ans)<-"LFprediction"
  ans
}

(core <- makeCluster(detectCores()))
clusterExport(core, c( "logreg", "predict.logreg", "pimp.import","prime.imp"))
registerDoParallel(core)
##
newanneal<-logreg.anneal.control(start=2, end=-1, iter=50000)

vv=iLF3(resp=Y, Xs=as.data.frame(X),  
	     nBS=100, h=0.5, maxK=10,nBSXVars=floor(ncol(X)^0.7),
	      anneal.params=newanneal,err_conv=0.05,dect_order=8)

vv$err_list[[vv$best_k]]


pred_Y=predict.iLF(vv,test_X)$LFprediction
1-sum(abs(pred_Y-test_Y))/length(test_Y)







