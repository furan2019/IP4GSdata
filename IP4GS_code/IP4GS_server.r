#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(DT)
library(shinybusy)
library(umap)
library(Rcpp)
library(tsne)
# library(RColorBrewer)
library(GGally)
library(rrBLUP)
library(BGLR)
library(pbapply)
library(reshape)
library(reshape2)

## functions ##################################
cvSampleIndex <- function( sampleNum, cross = 5, seed = 1,randomSeed = FALSE ) {
  if(randomSeed == TRUE){
    seed <- randomSeed()
  }
  cv <- cross
  resList <- list()
  
  # leave-one-out
  if( cv == sampleNum ){
    vec <- 1:sampleNum
    for( i in 1:sampleNum ){
      resList[[i]] <- list( trainIdx = vec[-i], testIdx = i, cvIdx = i)
    }
  }else {
    #random samples 
    set.seed(seed)
    index <- sample(1:sampleNum, sampleNum, replace = FALSE )
    step = floor( sampleNum/cv )
    
    start <- NULL
    end <- NULL
    train_sampleNums <- rep(0, cv)
    for( i in c(1:cv) ) {
      start <- step*(i-1) + 1
      end <- start + step - 1
      if( i == cv ) 
        end <- sampleNum
      
      testIdx <- index[start:end]
      trainIdx <- index[-c(start:end)]
      resList[[i]] <- list( trainIdx = trainIdx, testIdx = testIdx, cvIdx = i)
    }
  }
  names(resList) <- paste0("cv",1:cross)
  resList
}


fit.BGLR <- function( trainMarkerMat, trainPheVal, predictMarkerMat, modelMethods,trainX = NULL, predX = NULL,  
                      outputModel = FALSE, nIter = 1500, burnIn = 500, thin = 5,
                      saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                      verbose = FALSE, rmExistingFiles = TRUE, groups=NULL){
  require("BGLR")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  numT <- nrow(predictMarkerMat)
  MarkerMat <- rbind(predictMarkerMat,trainMarkerMat)
  pheVal <- c(rep(NA,numT),trainPheVal)
  fix <- rbind(predX,trainX)
  
  checkModel <- modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")
  if( ! checkModel ) {
    stop("Error: not defined model!")
  }
  
  if (!is.null(trainX) & !is.null(predX)) {
    ETA=list(list(X = MarkerMat, model= modelMethods),
             list(X = fix, model = "FIXED"))
  }else{
    ETA <- list(list(X=MarkerMat, model= modelMethods))
  }
  
  BGLRModel.fit <- BGLR( y = pheVal, ETA=ETA,  nIter = nIter, burnIn = burnIn, thin = thin, 
                         saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                         verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups )
  
  Res <- BGLRModel.fit$yHat[1:numT]
  if(outputModel){
    Res <- list(model = BGLRModel.fit,predictRes = Res)
  }
  Res
}
##
trainModel_RRBLUP <- function( markerMat, phenVec,X = NULL){
  require("rrBLUP")
  phen_answer<-mixed.solve(phenVec, Z=markerMat, K=NULL, SE = FALSE, return.Hinv=FALSE,X = X)
  beta <- phen_answer$beta
  phD <- phen_answer$u
  e <- as.matrix(phD)
  return( list(beta = beta, e = e, phD = phD) )
}


GSReModel <- function(markers, pheVal, modelMethods,NAImpute = F,
                      K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,
                      alpha = 1,
                      X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8, weights = NULL,
                      ...){  
  
  if (sum(is.na(pheVal)) != 0) {
    NAIdx <- which(is.na(pheVal) == T)
    if (NAImpute) {
      pheVal[NAIdx] <- mean(pheVal,na.rm =T)
    }else{
      pheVal <- pheVal[-NAIdx]
      markers <- markers[-NAIdx,]
      if (!is.null(X)) {
        X <- X[-NAIdx,]
      }
    }  
  }
  checkModel <- modelMethods %in% c("RRBLUP","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GSReModel!")
  }
  #   if (modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")){
  #     BGLRmethods <- modelMethods
  #     modelMethods <- "BGLRModel"
  #   }
  switch(modelMethods,
         #        BGLRModel  = trainedPredictModel_BGLR(trainMarkerMat = markers, trainedPhenVec = pheVal, modelMethods = BGLRmethods ,nIter = nIter, burnIn = burnIn, thin = thin, 
         #                                                saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
         RRBLUP  = trainModel_RRBLUP(markerMat = markers, phenVec = pheVal,X = X),
         LASSO = trainModel_LASSO(markers,pheVal,alpha = alpha, ...),
         SPLS = trainModel_spls(markers,pheVal,K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,maxstep = maxstep, ...),
         bigRR = trainModel_bigRR(markers = markers, X = X, pheVal = pheVal, weight = weights,family = family, lambda = lambda, tol.err = tol.err,tol.conv = tol.conv, ...)
  )
  
}

################################  LASSO ##############################
trainModel_LASSO <- function(markers,pheVal,alpha = 1, ...){
  require("glmnet")
  #glmnet fits a lasso model when we specify that alpha=1
  LASSO.fit <- glmnet(y=pheVal,x=markers,alpha=1, ...)
  #cv.glmnet finds the lambda value such that the the cvm value is the minimum
  cv <- cv.glmnet(y = pheVal, x=markers)
  LASSO_Res <- list(LASSO.fit = LASSO.fit,cv = cv)
  LASSO_Res
}
####################### LASSO pred ######################### 
pred_LASSO <- function(trainModel,testMat){
  predict(object = trainModel$LASSO.fit,testMat,s = trainModel$cv$lambda.min)
}

GSmachine <- function(markers, pheVal, modelMethods ="SVC", posPercentage = 0.4, BestIndividuals = c("top"),
                      ntree = 500,NAImpute = T,
                      nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9), ...){
  
  if (sum(is.na(pheVal)) != 0) {
    NAIdx <- which(is.na(pheVal) == T)
    if (NAImpute) {
      pheVal[NAIdx] <- mean(pheVal,na.rm =T)
    }else{
      pheVal <- pheVal[-NAIdx]
      markers <- markers[-NAIdx,]
    }  
  }
  
  if(is.null(nodesize)){
    if(modelMethods == "RFC"){
      nodesize <- 1
    }else if(modelMethods == "RFR"){
      nodesize <- 5
    }
  }else{
    nodesize <- nodesize
  }
  require("e1071")
  require("randomForest") 
  if( !modelMethods%in% c("SVR","SVC","RFR","RFC") ) {
    stop("Error: not defined category")
  }
  if (modelMethods %in%  c("SVC","RFC")){
    posNegSampleList <- sampleClassify(phenotype = pheVal ,posPercentage = posPercentage ,BestIndividuals = BestIndividuals )
    markers <- markers[c(posNegSampleList$posSampleIndex,posNegSampleList$negSampleIndex),]
    pheVal <-  as.factor( c( rep("1", length(posNegSampleList$posSampleIndex)), rep("0", length(posNegSampleList$negSampleIndex)) ) )
    
  }
  
  if(modelMethods %in%  c("SVR","SVC")){
    modelMethods <- "svm"
  }
  
  if(modelMethods %in% c("RFR","RFC")){
    modelMethods <- "randomforest"
  }
  switch(modelMethods,
         svm = svm(x= markers, y = pheVal,kernel = kernel,cost=cost,gamma = gamma,probability = TRUE, ...),
         randomforest = randomForest( x = markers, y = pheVal, ntree = ntree, importance = F,nodesize = nodesize,...))
}

predictGS <- function(testMat, trainModel,predX = NULL,modelMethods = "SVC",type = "fit"){
  ########## check the methods of GS
  checkModel <- modelMethods %in% c("RRBLUP","SVR","SVC","RFR","RFC","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GS Model")
  }
  
  Methods <- modelMethods
  
  ####### check testset 
  if(!is.matrix(testMat)) {
    testMat <- rbind(testMat,testMat)
    testnum <- 1 
  }else{
    testnum <- nrow(testMat)
  }
  
  #############
  if (!is.null(predX)) {
    if (modelMethods %in% c("RRBLUP","SVR","RFR","LASSO","SPLS","bigRR")){
      predresult <- switch(Methods,
                           #BGLRModel = {mkr.effs <- as.numeric(trainModel$ETA[[1]]$b); testMat %*% mkr.effs},
                           RRBLUP = {pred_phenVec <-testMat %*% trainModel$e;beta <- matrix(trainModel$beta,nrow = ncol(predX));beta <- predX %*% beta;as.numeric(pred_phenVec[,1]) + as.numeric(beta)},
                           SVR = {testMat <- cbind(testMat,predX);predict( trainModel, testMat)},
                           RFR = {testMat <- cbind(testMat,predX);predict( trainModel,testMat )},
                           LASSO = {testMat <- cbind(testMat,predX);pred_LASSO(trainModel,testMat)},
                           SPLS = {testMat <- cbind(testMat,predX);pred_SPLS(trainModel,testMat,type= type)},
                           bigRR = {beta <- matrix(trainModel$beta,nrow = ncol(predX));beta <- predX %*% beta;as.numeric(trainModel$beta) + as.matrix(testMat %*% trainModel$u)}
      )
    }
    else if (modelMethods %in% c("SVC","RFC")){
      predresult <- switch(modelMethods,
                           SVC = {testMat <- cbind(testMat,predX);obj_pred <- predict(trainModel,testMat, probability = TRUE); as.matrix(attr(obj_pred, "probabilities")[,"1"])},
                           RFC = {testMat <- cbind(testMat,predX);predict(trainModel, testMat, type= "vote" )[,"1"]})
    }
  }else{
    if (modelMethods %in% c("RRBLUP","SVR","RFR","LASSO","SPLS","bigRR")){
      predresult <- switch(Methods,
                           #BGLRModel = {mkr.effs <- as.numeric(trainModel$ETA[[1]]$b); testMat %*% mkr.effs},
                           RRBLUP = {pred_phenVec <-  testMat %*% trainModel$e; as.numeric(pred_phenVec[,1]) + as.numeric(trainModel$beta)},
                           SVR = predict( trainModel, testMat),
                           RFR = predict( object = trainModel,testMat ),
                           LASSO = pred_LASSO(trainModel,testMat),
                           SPLS = pred_SPLS(trainModel,testMat,type= type),
                           bigRR = as.numeric(trainModel$beta + testMat %*% trainModel$u)
      )
    }
    else if (modelMethods %in% c("SVC","RFC")){
      predresult <- switch(modelMethods,
                           SVC = {obj_pred <- predict(trainModel,testMat, probability = TRUE); as.matrix(attr(obj_pred, "probabilities")[,"1"])},
                           RFC = predict(trainModel, testMat, type= "vote" )[,"1"])
    }
  }
  predresult[1:testnum]
}

G2P <- function(trainMarker,trainPheno,testMarker,testPheno = NULL,modelMethods ="BayesA",outputModel = FALSE,  # main parameters
                nIter = 1500, burnIn = 500, thin = 5, 
                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                posPercentage = 0.4,BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9),  # machine learing parameters
                K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,trace = FALSE,maxstep = 100, # SPLS parameters
                alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 20,
                ...){   # LASSO parameters
  if(is.null(testPheno)){
    ntestSample <- nrow(testMarker)
    if (is.null(ntestSample)) {
      testPheno <- rep(0,1)
    }else{
      testPheno <- rep(0,ntestSample)
    }
  }
  allModelMethod <- c('SPLS','LASSO','SVC','SVR','RFC','RFR','RRBLUP','bigRR','BayesA','BayesB','BayesC','BL','BRR','RKHS','RR','BRNN')
  modelMethods <- modelMethods
  ## judge 
  if (!all(modelMethods %in% allModelMethod)) {
    stop("The provided list of modelMethods is out of G2P model list, please check!")
  }
  
  require("pbapply")
  pboptions(type = "timer")
  result <- pblapply(modelMethods, function(ii){cat(ii,"is modeling ...","\n");
    switch(ii,
           RR = fit.RR(trainMarker,trainPheno,testMarker,outputModel = TRUE),
           BRNN = fit.BRNN(trainMarker,trainPheno,testMarker,outputModel = TRUE,verbose = verbose, neurons = neurons, epochs = epochs, ...),
           RKHS = fit.RKHS(trainMarker,trainPheno,testMarker,outputModel = TRUE, nIter = nIter, burnIn = burnIn,thin = thin,
                           saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesA = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesA",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesB = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesB",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesC = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesC",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BL = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BL",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                         saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BRR = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BRR",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                          saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           RRBLUP = {resModel <- GSReModel(modelMethods = "RRBLUP",markers = trainMarker,pheVal = trainPheno,
                                           K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
                                           alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RRBLUP");
           list(model = resModel,predictRes = predictRes)},
           LASSO = {resModel <- GSReModel(modelMethods = "LASSO",markers = trainMarker,pheVal = trainPheno,
                                          K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
                                          alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "LASSO");
           list(model = resModel,predictRes = predictRes)},
           SPLS = {resModel <- GSReModel(modelMethods = "SPLS",markers = trainMarker,pheVal = trainPheno,
                                         K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
                                         alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SPLS");
           list(model = resModel,predictRes = predictRes)},
           bigRR = {resModel <- GSReModel(modelMethods = "bigRR",markers = trainMarker,pheVal = trainPheno,
                                          K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
                                          alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "bigRR");
           list(model = resModel,predictRes = predictRes)},
           SVC = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "SVC" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVC");
           list(model = resModel,predictRes = predictRes)},
           
           RFC = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "RFC" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFC");
           list(model = resModel,predictRes = predictRes)},
           
           SVR = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "SVR" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVR");
           list(model = resModel,predictRes = predictRes)},
           
           RFR = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "RFR" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFR");
           list(model = resModel,predictRes = predictRes)}
    )
  })
  
  ModelList <- lapply(result, function(x){a <- x[[1]]})
  PredResList <- lapply(result, function(x){b <- x[[2]]})
  lengtList <- length(result)
  realPheno <- as.matrix(testPheno)
  resMat <- as.matrix(realPheno)
  for(i in 1:lengtList){
    resMat <- cbind(resMat,as.matrix(PredResList[[i]]))
  }
  colnames(resMat) <- c("realPhenScore",modelMethods)
  rownames(resMat) <- rownames(testMarker)
  names(ModelList) <- modelMethods
  finalRes <- list(model = ModelList, predScore = resMat)
  # predscores <- predscores[,c("realPhenScore",modelMethods)]
  ####### output result
  if(outputModel){
    return(finalRes)
  }
  else{
    return(resMat)
  }
}

##################### G2P CV #######################
G2PCrossValidation <-function(cvSampleList = NULL,cross = 10,times = 1,seed = 1,cpus = 1, markers, pheVal, modelMethods ="SVC", outputModel = FALSE,
                              nIter = 1500, burnIn = 500, thin = 5, 
                              saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                              verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
                              posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500, nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9),
                              K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,  # SPLS parameters
                              alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                              epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 20,
                              ...){
  require("parallel")
  require("pbapply")
  
  sampleNum <- nrow(markers)
  phenotypeNum <- length(pheVal)
  cvSampleList = cvSampleList;
  cross = cross; seed = seed; cpus = cpus; markers = markers;pheVal = pheVal;
  modelMethods = modelMethods; BestIndividuals = BestIndividuals;
  nIter = nIter; burnIn = burnIn; thin = thin; saveAt = saveAt; S0 = S0; df0 =df0; R2 = R2; weights = weights;posPercentage = posPercentage;
  verbose = verbose; rmExistingFiles = rmExistingFiles; groups=groups;ntree = ntree  ;nodesize = nodesize ;importance = importance;
  kernel = kernel;gamma = gamma; cost = cost;outputModel = outputModel;
  K = K;eta = eta;select = select;fit = fit;scale.x = scale.x;scale.y = scale.y;eps = eps;trace = trace;
  alpha = alpha;X = X;family = family; lambda = lambda; tol.err = tol.err; tol.conv = tol.conv;maxstep = maxstep;
  epochs = epochs; neurons = neurons; Z = Z;  effectConcider = effectConcider; mmerIters = mmerIters;
  
  if(sampleNum != phenotypeNum) {
    stop("Marker count is not equal to phenotype count!")
  }
  
  if(!is.numeric(markers) | !is.numeric(pheVal)){
    stop("Marker or phenotype is not numeric, please check it!")  
  }
  
  
  cl <- makeForkCluster(cpus)
  cat(cpus," cores were used for cross validation ... \n")
  cat("Start cross validation ... \n")
  
  if(!is.null(cvSampleList)){
    cvSampleList <- cvSampleList
  }else{
    if (times == 1) {
      cvSampleList <- cvSampleIndex(sampleNum = sampleNum,cross = cross,seed = seed,randomSeed = F)
    }else if (times > 1){
      cvSampleListToatal <- c()
      for (i in 1:times) {
        cvSampleList <- cvSampleIndex(sampleNum = sampleNum,cross = cross,seed = seed,randomSeed = T)
        names(cvSampleList) <- paste0("Times_",i,"_",names(cvSampleList))
        cvSampleListToatal <- c(cvSampleListToatal,cvSampleList)
        cvSampleListToatal
      }
      cvSampleList <- cvSampleListToatal
    }
  }
  
  pboptions(type = "timer")
  results <- pblapply(1:length(cvSampleList),function(x){
    # library("BGLR")
    # library("G2P")
    library("brnn")
    library("glmnet")
    library("spls")
    library("pls")
    library("e1071")
    library("BGLR")
    library("rrBLUP")
    library("randomForest")
    library("sommer")
    library("hglm")
    # source("./code/G2P.r")
    
    cat("All needed package have loaded in all cores! \n")
    
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <-  cvSampleList[[x]]$testIdx
    trainMarker <- markers[trainIdx,]
    trainPheno <- pheVal[trainIdx]
    testMarker <- markers[testIdx,]
    testPheno <- pheVal[testIdx]
    G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker,testPheno =testPheno ,modelMethods = modelMethods ,BestIndividuals = BestIndividuals,
        nIter = nIter, burnIn = burnIn, thin = thin,  saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,posPercentage = posPercentage,
        verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,ntree = ntree  ,nodesize = nodesize ,importance = importance,
        kernel = kernel,gamma = gamma, cost = cost,outputModel = outputModel,
        K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
        alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep,
        epochs = epochs, neurons = neurons, Z = Z,  effectConcider = effectConcider, mmerIters = mmerIters,
        ...)
  },cl=cl) # lapply
  stopCluster(cl)
  cat(times," times",cross," fold cross validation is done! \n")
  if (times == 1) {
    names(results) <- paste0("CV",1:cross)
  }else{
    final_res <- list()
    length(final_res) <- times
    names(final_res) <- paste0("Rep",1:times)
    for (i in 1:times) {
      final_res[[i]] <- results[((i-1)*cross):(i*cross)]
      names(final_res[[i]]) <- paste0("CV",1:cross)
    }
    results <- final_res
  }
  results
}

################################################ eval ###############################################
meanNDCG <- function( realScores, predScores, topK = 10 ){
  resVec <- rep(0, topK )
  for( idx in 1:topK ){
    resVec[idx] <- NDCG( realScores, predScores, topK = idx )
  }
  meanNDCG <- mean(resVec)
  names(meanNDCG ) <- paste0("top",topK)
  return (meanNDCG)
}
##from plos one, 2015, A Ranking Approach to Genomic Selection
NDCG <- function( realScores, predScores, topK = 10){
  
  if( length(realScores) != length(predScores)){
    stop("Error: different length between realScores and predScores")
  }
  if( length(realScores) < topK ) {
    stop("Error: too large topK")
  }
  
  scoreMat <- cbind(realScores,predScores)
  scoreMatSortbyPred <- scoreMat[order(scoreMat[,2],decreasing = TRUE),]
  scoreMatSortByReal <- scoreMat[order(scoreMat[,1],decreasing = TRUE),]
  
  DCG <- rep(0, topK)
  IDCG <- rep(0, topK)
  for(idx in 1:topK){
    DCG[idx] <-  scoreMatSortbyPred[idx,1]/log(idx+1,2)
    IDCG[idx] <- scoreMatSortByReal[idx,1]/log(idx+1,2)
  }
  
  NDCG <- sum(DCG)/sum(IDCG) 
  names(NDCG) <- paste0("top",topK)
  return(NDCG)
}
##evaluation method: pearson, spearman, kendall, MSE
corEvaluation <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),BestIndividuals = "top",probIndex = NULL){
  # Probability handle
  if (!is.null(probIndex)) {
    if(BestIndividuals == "top"){
      realScores <- realScores
      predScores[,probIndex] <- predScores[,probIndex]
    }else if(BestIndividuals == "buttom"){
      realScores <- realScores
      predScores[,probIndex] <- 1 - predScores[,probIndex]
    }else if(BestIndividuals == "middle"){
      realScores <- abs(realScores)
      predScores[,probIndex] <- 1 - predScores[,probIndex] 
    }
  }else{
    realScores <- realScores
    predScores <- predScores
  }
  
  if( length(method) > 1 ){
    method <- method[1]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  res <- ""
  if( (method == "pearson") | (method == "kendall") | (method == "spearman") ){
    res <- cor( realScores, predScores,  use="complete", method = method  )
  }else if(method == "MSE") {
    res <- apply(predScores,2,function(ii){
      deltaVec <- abs( realScores - ii)
      deltaVec <- deltaVec^2
      mean(deltaVec)})
  }else if(method == "R2"){
    res <-apply(predScores,2,function(ii){
      R2 <- summary(lm(realScores ~ ii))$r.squared })
  }
  
  res <- matrix(res,nrow = 1,dimnames = list(method,colnames(predScores)))
  res
}
## evaluation method alpha : pearson, spearman, kendall, MSE
corEvaluationAlpha <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha = 15,BestIndividuals = "top",probIndex = NULL){
  topNum <- 1: round(length(realScores)*(topAlpha/100))
  if( length(method) > 1 ){
    method <- method[1]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  mat <- cbind(realScores,predScores)
  
  if(BestIndividuals == "top"){
    mat <- mat[order(realScores,decreasing = T),]
  }else if(BestIndividuals == "buttom"){
    mat <- mat[order(realScores,decreasing = F),]
  }else if(BestIndividuals == "middle"){
    mat <- mat[order(abs(realScores),decreasing = F),]
  }
  
  mat <- mat[topNum,]
  
  realScores <- as.matrix(mat[,1])
  predScores <- as.matrix(mat[,-1])
  
  res <- corEvaluation(realScores, predScores, method = method,BestIndividuals = BestIndividuals,probIndex = probIndex)
  res
}

multiAlphaCor <- function(realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha = 1:15,BestIndividuals = "top",probIndex = NULL){
  res <- lapply(method,function(meth){
    # cormat <- apply(predScores,2,function(x){
    cormat <- sapply(topAlpha,function(ii){
      corEvaluationAlpha(realScores = realScores,predScores = predScores,method = meth,topAlpha = ii,BestIndividuals = BestIndividuals, probIndex = probIndex)})
    cormat <- matrix(cormat,ncol = length(topAlpha))
    cormat <- t(cormat)
    rownames(cormat) <- paste0("top",topAlpha)
    colnames(cormat) <- colnames(as.matrix(predScores))
    cormat
  })
  names(res) <- method
  res
}

classEvaluation <- function(realScores, predScores, topAlpha = 15, Beta = 1,Probability = FALSE,evalMethod = "AUC",BestIndividuals = c("top", "middle", "buttom") ) {
  if( length(realScores) != length(predScores) ){
    stop("Error: different length between realScores and predScores")
  }
  
  if( length(BestIndividuals) > 1 ) {
    BestIndividuals <- BestIndividuals[1]
  }
  
  realScores <- as.numeric(realScores)
  predScores <- as.numeric(predScores)
  total <- length(realScores)
  topK <- round( total*topAlpha/100 )
  classVec <- c( rep(1, topK), rep(0, total-topK))
  ##
  if(BestIndividuals == "top"){
    decreaseReal <- TRUE
    decreasePred <- TRUE
  }else if(BestIndividuals == "buttom"){
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }else if(BestIndividuals == "middle"){
    realScores <- abs(realScores)
    predScores <- abs(predScores)
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }
  
  scoreMat <- cbind( realScores, predScores )
  newScoreMat <- scoreMat[order(scoreMat[,1], decreasing = decreaseReal),] 
  newScoreMat <- cbind( newScoreMat, classVec )
  topRealMean <- mean( newScoreMat[1:topK,1] )
  #
  threshold <- newScoreMat[topK,1]
  #
  
  ##### RE ,kappa
  newScoreMat <- newScoreMat[order(newScoreMat[,2], decreasing = decreasePred),]
  # classVec <- c(rep(1,length(which(newScoreMat[,2] <= threshold))),rep(0, total-(length(which(newScoreMat[,2] <= threshold)))))
  newScoreMat <- cbind( newScoreMat, classVec )
  colnames(newScoreMat) <- c("real", "predicted", "realClass", "predClass")
  PredRealMean <- mean( newScoreMat[1:topK,1] ) 
  
  TP <- sum(newScoreMat[1:topK, 3])
  FP <- topK - TP
  FN <- topK - TP
  TN <- total - topK - FN
  Po <- (TP+TN)/total
  Pe <- ((FP+TN)/total)*((FN+TN)/total) + ((TP+FN)/total)*((TP+FP)/total)
  allRealMean <- mean( scoreMat[,1] )
  precision = TP/(TP + FP)
  recall = TP/(TP +FN)
  ###  the area under the receiver operating characteristics curve(AUC),and  the area under the precision-recall  curve (AUCpr)
  result <- sapply(evalMethod,function(one_method){
    switch(one_method,
           Kappa =  (Po-Pe)/(1-Pe),
           RE = ( PredRealMean - allRealMean)/(topRealMean - allRealMean),
           #AUC = roc.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE)$AUC ,
           AUC = roc(newScoreMat[,3],newScoreMat[,4],quiet = T)$auc[[1]],
           AUCpr = pr.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE,sorted = TRUE)$auc.integral,
           accuracy = (TP + TN)/total,
           F1 = (1 + Beta^2)*precision *recall/(precision + recall )
    )
  })
  names(result) <- paste(evalMethod,"_top", topAlpha, sep = "" )
  return(result) 
}
##################################### to evaluate for multiple parameters set
multiParameters <- function(realScores,predScores, Probability = TRUE,evalMethod = c("RE"),topAlpha ,Beta = 1,BestIndividuals = c("top"),probIndex = NULL){
  
  ############ calculate one or multiple topAlpha 
  multiTopAlpha <- function(realScores,predScores, Probability ,evalMethod ,topAlpha,Beta,BestIndividuals){
    ######################### one or multiple topAlpha for classifiction evaluation
    if (length(intersect(evalMethod,c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" ))) != 0){
      result <-  sapply(topAlpha,function(ii)
        classEvaluation(realScores = realScores,predScores = predScores,Probability = Probability,
                        evalMethod = evalMethod,topAlpha = ii,Beta = Beta,BestIndividuals = BestIndividuals)
      ) }
    
    ######################### one or multiple topAlpha for NDCG evaluation
    ################ set format of output
    if(length(evalMethod) > 1){
      result <- t(result)
    }else{
      result <- as.matrix(result)
    }
    dimnames(result) <- list(paste0("top",topAlpha),evalMethod )
    return(result)
  }
  
  
  predScores <- as.matrix(predScores)
  evalNum <- 1:ncol(predScores)
  
  ## class ##
  if(!is.null(probIndex)){
    classMethSite <- probIndex
    evalNum1 <- evalNum[-classMethSite]
    multiPredresultS <- lapply(evalNum1,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    classPredresult <- lapply(classMethSite,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = TRUE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    multiPredresult <- list()
    length(multiPredresult) <- ncol(predScores)
    multiPredresult[evalNum1] <- multiPredresultS
    multiPredresult[classMethSite] <- classPredresult
  }else{
    ############ evaluate one or multiple prediction result
    multiPredresult <- lapply(evalNum,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
  }
  ######### set format of output
  result <- lapply(evalMethod,function(ii){
    one_evalMethod <- sapply(evalNum,function(jj) multiPredresult[[jj]][,ii])
    if(length(topAlpha) > 1){
      colnames(one_evalMethod) <- colnames(predScores)
    }else{
      one_evalMethod <- matrix(one_evalMethod,nrow = 1,dimnames = list(paste0("top",topAlpha),colnames(predScores)))
    }
    one_evalMethod
  })
  names(result) <- evalMethod
  return(result)
}

evaluateGS <- function(realScores, predScores, evalMethod = "RE", Beta = 1, BestIndividuals = "top", topAlpha = 1:90, allNDCG = F, globalAlpha = F, probIndex = NULL){
  require(pROC)
  require(PRROC)
  ## data process
  ## remove NAs
  evalMat <- cbind(as.matrix(realScores),as.matrix(predScores))
  evalMat <- na.omit(evalMat)
  ## warning prediction results not converging
  check <- apply(as.matrix(evalMat[,2:ncol(evalMat)]),2,function(x){
    if (length(unique(x)) == 1) {
      1
    }else{
      0
    }
  })
  
  if (1 %in% check) {
    warning(paste0("The prediction resluts of model ",paste0(names(check)[check == 1]),collapse = ",")," is abnormal. Please check.")
    abn.idx <- which(check == 1)
    if (!is.null(probIndex)){
      for(i in 1:length(probIndex)){
        diff.idx <- probIndex[i] - abn.idx
        ## minus
        probIndex[i] <- probIndex[i] - length(diff.idx[diff.idx > 0])
      }
    }
    evalMat <- evalMat[,-(which(check == 1)+1)]
    if (!is.matrix(evalMat)) {
      stop("There is no predictions value and cannot be assessed!")
    }
  }
  
  # probMat <- transValMat2ProbMat(evalMat = evalMat,BestIndividuals = BestIndividuals)
  realScores <- evalMat[,1]
  predScores <- evalMat[,-1]
  ## Evaluation
  selectFun <- NULL
  predScores <- as.matrix(predScores)
  globalMethods <-  c("pearson", "kendall", "spearman", "MSE","R2")
  thresholdMethods <- c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" )
  NDCGMethods <- c("meanNDCG", "NDCG")
  ###################### correlation methods
  if (length(intersect(evalMethod,globalMethods)) != 0){
    if(globalAlpha){
      selectFun <- c(selectFun,"globalEvaluation","globalEvaluationAlpha")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }else{
      selectFun <- c(selectFun,"globalEvaluation")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }
  }
  
  ################# classification evaluation
  if (length(intersect(evalMethod,thresholdMethods)) != 0){
    selectFun <- c(selectFun,"thresholdEvaluation")
    thresholdMethods <- intersect(evalMethod,thresholdMethods)
  }
  
  result <- lapply(selectFun,function(one_fun){
    switch(one_fun,
           globalEvaluation = {corResult <- sapply(globalMethods,function(one_methods){corEvaluation( realScores, predScores, method = one_methods,BestIndividuals = BestIndividuals,probIndex = probIndex)});
           matrix(t(corResult),ncol= ncol(predScores),dimnames = list(globalMethods,colnames(predScores)))},
           globalEvaluationAlpha = multiAlphaCor(realScores,predScores,topAlpha = topAlpha,method = globalMethods,BestIndividuals = BestIndividuals,probIndex = probIndex),
           thresholdEvaluation = multiParameters(realScores, predScores, topAlpha = topAlpha, Probability = Probability,evalMethod = thresholdMethods,Beta =Beta , BestIndividuals = BestIndividuals, probIndex = probIndex)
    )})
  ############ the format of output
  finalresult <- list()
  nList <- length(result)
  nList <- c("one","two","three")[nList]
  if(length(result) != 0){
    #     id <- 1
    #     if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];id <- 2};finalresult <- c(finalresult,result[[id]])
    switch(nList,
           one = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]]}else{finalresult <- c(finalresult,result[[1]])}},
           two = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]])}else{finalresult <- c(finalresult,result[[1]],result[[2]])}},
           three = {finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]],result[[3]])})
  }else{finalresult <- list()}
  #################### NDCG evaluation
  if (length(intersect(evalMethod,NDCGMethods)) != 0){
    selectFun <- intersect(evalMethod,NDCGMethods)
    result2 <- lapply(selectFun,function(one_fun){
      switch (one_fun,
              NDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "NDCG",topAlpha = topAlpha,realScores = realScores,predScores = x,allNDCG = allNDCG)}),
              meanNDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "meanNDCG",topAlpha = topAlpha,realScores = realScores,predScores = x,allNDCG = allNDCG)})
      )})
    names(result2) <- intersect(evalMethod,NDCGMethods)
    finalresult <- c(finalresult,result2)
  }
  finalresult
}
NDCGEvaluation <- function(method = NULL,topAlpha,realScores,predScores, allNDCG = T){
  
  if(allNDCG == T){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
  }else{
    ndcg <- round(length(realScores)*(topAlpha/100))
  }
  
  if(method == "NDCG"){
    result <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
      
    })
    if(allNDCG == F){
      names(result) <- paste0("top",topAlpha)
    }
  }else if(method == "meanNDCG"){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
    NDCGEval <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
    })
    result <- sapply(topAlpha,function(ii){
      iii <- round(ii*length(realScores)/100)
      sum(NDCGEval[1:iii])/iii
    })
    names(result) <- paste0("top",topAlpha)
  }
  result
}


## server overall #################
server <- shinyServer(function(input, output) {
    ## set the maximum size of upload 
    color <- function(){
      color = colorRampPalette(rev(palette()))(8000)
      color
    }
    options(shiny.maxRequestSize=100*1024^2)
    
    ##################################### QC module and pre-processing module #####################################
    apply_pb <- function(X, MARGIN, FUN, ...) {
        env <- environment()
        pb_Total <- sum(dim(X)[MARGIN])
        counter <- 0
        pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
        wrapper <- function(...) {
            curVal <- get("counter", envir = env)
            assign("counter", curVal + 1, envir = env)
            setTxtProgressBar(get("pb", envir = env), curVal + 1)
            FUN(...)
        }
        res <- apply(X, MARGIN, wrapper, ...)
        close(pb)
        res
    }
    
    apply_pb_shiny <- function(X, MARGIN, message, FUN, ...){
      env <- environment()
      pb_Total <- dim(X)[MARGIN]
      counter <- 0
      # percentage <- 0
      wrapper <- function(...) {
        # Sys.sleep(0.1);
        # percentage <<- percentage + 1/pb_Total*100
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        details <- paste0("Progress: \n",round((curVal+1)/pb_Total,4)*100,"%")
        incProgress(1/pb_Total, detail = details)
        FUN(...)
      }
      res <- withProgress(message = message, value=0, {
        apply(X, MARGIN, wrapper, ...)
      })
      res
    }
    
    apply_pb_busy <- function(X, MARGIN, FUN, ...) {
      env <- environment()
      pb_Total <- dim(X)[MARGIN]
      counter <- 0
      # pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
      show_modal_progress_line()
      wrapper <- function(...) {
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        # setTxtProgressBar(get("pb", envir = env), curVal + 1)
        # percentage <<- percentage + 1/pb_Total*100
        update_modal_progress(
          value = (curVal + 1)/ pb_Total,
          text = paste("Process", trunc( (curVal + 1)/pb_Total), sprintf("(%02d%%)",(curVal + 1) ))
        )
        FUN(...)
      }
      res <- apply(X, MARGIN, wrapper, ...)
      remove_modal_progress()
      res
    }
    
    GSImputation <- function(G,imputeMethod = "median",silent = FALSE){
        if(is.numeric(G)){
            if(imputeMethod == "mean"){
                if(silent){
                    res <- apply(G,2,function(x){
                        x[which(is.na(x))] <- mean(x,na.rm=TRUE)
                        x
                    })
                }else{
                    res <- apply_pb_shiny(G,2,'Imputation is doing, method is "mean" ...',function(x){
                        x[which(is.na(x))] <- mean(x,na.rm=TRUE)
                        x
                    })
                }
            }else if(imputeMethod=="median"){
                if(silent){
                    res <- apply(G,2,function(x){
                        tt <- table(x)
                        x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
                        x
                    })
                }else{
                    res <- apply_pb_shiny(G,2,'Imputation is doing, method is "median" ...',function(x){
                        tt <- table(x)
                        x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
                        x
                    })
                }
            }
            class(res) <- "numeric"
        }else{
            if(imputeMethod =="mean"){
                stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
            }else {
                if(silent){
                    res <- apply(G,2,function(x){
                        tt <- table(x)
                        x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
                        x
                    })
                }else{
                    res <- apply_pb_shiny(G,2,'Imputation is doing, method is "median" ...',function(x){
                        tt <- table(x)
                        x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
                        x
                    })
                }
            }
        }
        return(res)
    }
    
    MAFSummary <- function(G,format = "letter",hete = 1,silent = FALSE){
        siSNP2biSNP <- function(x) {
            y <- rep(NA, length(x))
            y[which(x == "A")] <- "AA"
            y[which(x == "T")] <- "TT"
            y[which(x == "C")] <- "CC"
            y[which(x == "G")] <- "GG"
            y[which(x == "R")] <- "AG"
            y[which(x == "Y")] <- "CT"
            y[which(x == "S")] <- "CG"
            y[which(x == "W")] <- "AT"
            y[which(x == "K")] <- "GT"
            y[which(x == "M")] <- "AC"
            y[which(x == "+")] <- NA
            y[which(x == "0")] <- NA
            y[which(x == "-")] <- NA
            y[which(x == "N")] <- NA
            return(y)
        }
        if (format == "letter") {
            judge <- apply(G[,sample(1:ncol(G),round(ncol(G)*0.02))],2,function(x){
                x <- x[!is.na(x)]
                if (length(x) == 0) {
                    y <- 0.5
                }else{
                    allele <- strsplit(sample(x,1),split = "")[[1]]
                    allele.len <- length(allele)
                    if (allele.len == 1) {
                        y <- 1
                    }else{
                        y <- 0
                    }
                }
            })
            if (sum(judge)/(round(ncol(G)*0.02)) > 0.9) {
                message("Input genotype is coded by single-letter, if not, please check your data!")
                G.rownames <- rownames(G)
                G.colnames <- colnames(G)
                message("Input genotype is coded by single-letter, converting single-letter to bi-letter ...")
                if (silent) {
                    G <- apply(G,2,siSNP2biSNP)
                } else {
                    G <- apply_pb_shiny(G,2,'Converting single-letter to bi-letter ...',siSNP2biSNP)
                }
                rownames(G) <- G.rownames
                colnames(G) <- G.colnames
            }
            
            ## maf calculation
            message("Calculating minor allele frequency (MAF) ...")
            if (silent) {
                maf.info <- apply(G,2,function(x){
                    x <- na.omit(x)
                    fasta <- paste0(x,collapse = "")
                    allele.array <- strsplit(fasta,split = "")
                    allele.info <- sort(table(na.omit(allele.array)),decreasing = T)
                    if(length(allele.info) < 2){
                        af <- allele.info[1]/sum(allele.info)
                        maf <- 1-af
                        af.allele <- names(allele.info)[1]
                        if (is.null(af.allele)) {
                            af.allele <- NA
                        }
                        maf.allele <- NA
                    }else if(length(allele.info) > 2){
                        warning("Data have multi-allele loci, these loci will be removed!")
                        af <- maf <- maf.allele <- af.allele <- "remove"
                    }else{
                        af <- allele.info[1]/sum(allele.info)
                        maf <- allele.info[2]/sum(allele.info)
                        maf.allele <- names(allele.info)[2]
                        af.allele <- names(allele.info)[1]
                    }
                    res <- c(af.allele,maf.allele,af,maf)
                    res
                })
            }else{
                maf.info <- apply_pb_shiny(G,2,'Statistic for genotype information ...',function(x){
                    x <- na.omit(x)
                    fasta <- paste0(x,collapse = "")
                    allele.array <- strsplit(fasta,split = "")
                    allele.info <- sort(table(na.omit(allele.array)),decreasing = T)
                    if(length(allele.info) < 2){
                        af <- allele.info[1]/sum(allele.info)
                        maf <- 1-af
                        af.allele <- names(allele.info)[1]
                        if (is.null(af.allele)) {
                            af.allele <- NA
                        }
                        maf.allele <- NA
                    }else if(length(allele.info) > 2){
                        warning("Data have multi-allele loci, these loci will be removed!")
                        af <- maf <- maf.allele <- af.allele <- "remove"
                    }else{
                        af <- allele.info[1]/sum(allele.info)
                        maf <- allele.info[2]/sum(allele.info)
                        maf.allele <- names(allele.info)[2]
                        af.allele <- names(allele.info)[1]
                    }
                    res <- c(af.allele,maf.allele,af,maf)
                    res
                })
            }
            ## remove multi-allele loci 
            if(length(which(maf.info == "remove",arr.ind = T)) > 0){
                message("Checking for markers with more than 2 alleles. If found will be removed.")
                remove.idx <- unique(which(maf.info == "remove",arr.ind = T)[,1])
                maf.info <- maf.info[-remove.idx,]
                G <- G[,-remove.idx]
            }
            maf.info <- t(maf.info)
            maf.info <- as.data.frame(maf.info)
            
            class(maf.info$V3) <- class(maf.info$V4) <- "numeric"
            names(maf.info) <- c("Ref","Alt","AF","MAF")
            
        } else if (format == "number"){
            if(!is.numeric(G)){
                stop("The input format is not equal with defined! please check your data!")
            }else{
                message("Calculating minor allele frequency (MAF) ...")
                if (silent) {
                    maf.info <- apply(G,2,function(x){
                        x <- na.omit(x)
                        if (hete %in% x) {
                            if (length(table(x)) == 2) {
                                maf <- length(which(x == hete))/(2*length(x))
                            }else{
                                a <- length(which(x == hete))
                                b <- sort(table(x[-which(x == hete)]))
                                maf <- (b[[1]]*2+a)/(length(x)*2)
                                maf
                            }
                        }else{
                            a <- table(x)
                            b <- sort(a)
                            if (length(a) == 0) {
                                maf <- NA
                            }else{
                                maf <- b[1]/length(x)
                                if (maf > 0.5) {
                                    maf <- 1- maf
                                }
                            }
                        }
                        maf.info <- maf
                    })
                }else{
                    maf.info <- apply_pb_shiny(G,2,'Calculating minor allele frequency (MAF) ...',function(x){
                        x <- na.omit(x)
                        if (hete %in% x) {
                            if (length(table(x)) == 2) {
                                maf <- length(which(x == hete))/(2*length(x))
                            }else{
                                a <- length(which(x == hete))
                                b <- sort(table(x[-which(x == hete)]))
                                maf <- (b[[1]]*2+a)/(length(x)*2)
                                maf
                            }
                        }else{
                            a <- table(x)
                            b <- sort(a)
                            if (length(b) == 0) {
                                maf <- NA
                            }else{
                                maf <- b[1]/length(x)
                                if (maf > 0.5) {
                                    maf <- 1- maf
                                }
                            }
                        }
                        maf.info <- maf
                    })
                }
            }
        }
        maf.info
    }
    
    transLetter2number <- function(G,maf=0,mr=0,ref.alt = NULL,impute = TRUE,imputeMethod = "median"){
      siSNP2biSNP <- function(x) {
        y <- rep(NA, length(x))
        y[which(x == "A")] <- "AA"
        y[which(x == "T")] <- "TT"
        y[which(x == "C")] <- "CC"
        y[which(x == "G")] <- "GG"
        y[which(x == "R")] <- "AG"
        y[which(x == "Y")] <- "CT"
        y[which(x == "S")] <- "CG"
        y[which(x == "W")] <- "AT"
        y[which(x == "K")] <- "GT"
        y[which(x == "M")] <- "AC"
        y[which(x == "+")] <- NA
        y[which(x == "0")] <- NA
        y[which(x == "-")] <- NA
        y[which(x == "N")] <- NA
        return(y)
      }
      
      G <- as.matrix(G)
      judge <- apply(G[,sample(1:ncol(G),round(ncol(G)*0.02))],2,function(x){
        x <- x[!is.na(x)]
        if (length(x) == 0) {
          y <- 0.5
        }else{
          allele <- strsplit(sample(x,1),split = "")[[1]]
          allele.len <- length(allele)
          if (allele.len == 1) {
            y <- 1
          }else{
            y <- 0
          }
        }
      })
      if (sum(judge)/(round(ncol(G)*0.02)) > 0.9) {
        message("Input genotype is coded by singe-letter, if not, please check your data!")
        G.rownames <- rownames(G)
        G.colnames <- colnames(G)
        message("Input genotype is coded by singe-letter, converting single-letter to bi-letter ...")
        G <- apply(G,2,siSNP2biSNP)
        rownames(G) <- G.rownames
        colnames(G) <- G.colnames
      }
      
      ## missing rate
      if (mr > 0) {
        message("Calculating missing rate (MR) ...")
        missing.rate <- apply_pb_shiny(G,2,"Calculating missing rate (MR) ...",function(x){
          NA.idx <- which(is.na(x))
          NA.rate <- length(NA.idx)/length(x)
          NA.rate
        })
        if (sum(missing.rate > mr) > 0) {
          G <- G[,-which(missing.rate > mr)]
          message(paste0(sum(missing.rate > mr)," marker(s) was(were) filtered by missing rate > ",mr,"."))
        }else{
          message("None marker was removed by missing rate with ",mr,".")
        }
      }
      
      ## MAF
      message("Calculating minor allele frequency (MAF) ...")
      maf.info <- apply_pb_shiny(G,2,"Calculating minor allele frequency (MAF) ...",function(x){
        x <- na.omit(x)
        if (length(x) == 0) {
          warning("Data have column(s) with NA completely, these loci will be removed!")
          af <- maf <- maf.allele <- af.allele <- "remove"
        }else{
          fasta <- paste0(x,collapse = "")
          allele.array <- strsplit(fasta,split = "")
          allele.info <- sort(table(allele.array),decreasing = T)
          if(length(allele.info) < 2){
            af <- allele.info[1]/sum(allele.info)
            maf <- 1-af
            af.allele <- names(allele.info)[1]
            if (is.null(af.allele)) {
              af.allele <- NA
            }
            maf.allele <- NA
          }else if(length(allele.info) > 2){
            warning("Data have multi-allele loci, these loci will be removed!")
            af <- maf <- maf.allele <- af.allele <- "remove"
          }else{
            af <- allele.info[1]/sum(allele.info)
            maf <- allele.info[2]/sum(allele.info)
            maf.allele <- names(allele.info)[2]
            af.allele <- names(allele.info)[1]
          }
        }
        res <- c(af.allele,maf.allele,af,maf)
        res
      })
      
      maf.info <- t(maf.info)
      
      ## remove multi-allele loci 
      if(length(which(maf.info == "remove",arr.ind = T)) > 0){
        message("Checking for markers with more than 2 alleles. If found will be removed.")
        remove.idx <- unique(which(maf.info == "remove",arr.ind = T)[,1])
        maf.info <- maf.info[-remove.idx,]
        G <- G[,-remove.idx]
      }
      
      maf.info <- as.data.frame(maf.info)
      class(maf.info$V3) <- class(maf.info$V4) <- "numeric"
      names(maf.info) <- c("Ref","Alt","AF","MAF")
      
      ## filter G by maf
      if(maf > 0){
        maf.filter.idx <- which(maf.info$MAF < maf)
        if (length(maf.filter.idx) > 0) {
          G <- G[,-c(which(maf.info$MAF < maf))]
          maf.info <- maf.info[-c(which(maf.info$MAF < maf)),]
          message(paste0(length(maf.filter.idx)," marker(s) was(were) filtered by maf < ",maf,"."))
        }else{
          message("None marker was removed by maf with ",maf,".")
        }
      }
      
      ## Transformation of letter to number
      message("Tranformation is begainning ...")
      Ref <- maf.info$Ref
      lines.ID <- rownames(G)
      markers.ID <- colnames(G)
      
      G <- t(G)
      G.number <- apply_pb_shiny(cbind(Ref, G), 1,"Tranformation is begainning ...",function(x) {
        tmp <- gregexpr(pattern = x[1], text = x[-1], 
                        fixed = T)
        res <- as.integer(lapply(tmp, function(z) {
          ifelse(z[1] < 0,2,2 - length(z))
        }))
        return(res)
      })
      
      
      ## impute 
      if (impute) {
        missing <- which(is.na(G.number))
        if (length(missing) > 0) {
          message(paste0("Imputing missing data with mode ",imputeMethod,"."))
          G.number <- GSImputation(G.number,imputeMethod = imputeMethod,silent = silent)
        }
      }else {
        message("Imputation not required. Be careful using non-imputed matrices in mixed model solvers.")
      }
      
      rownames(G.number) <- lines.ID
      colnames(G.number) <- markers.ID
      final.res <- list(G.n = G.number,maf.info = maf.info)
      final.res
    } 
    
    ###############################################################################
    ## ++1.1 read data in ################
    ## ++++1.1.1 marker data input options ###########
    ## ++++2.1.1 marker data reactive options ###########
    data.exp0 <- reactive({
        req(input$markersData)
        infile.data <- input$markersData
        data.exp0 <- switch(input$markerFormat,
                            "Hapmap" = read.table(file = infile.data$datapath,header = input$genome.headerCK,sep = "\t",check.names = F,comment.char = ""),
                            "Matrix" = read.table(file = infile.data$datapath,header = input$genome.headerCK,sep = "\t",check.names = F,row.names = 1),
                            "Custom" = read.table(file = infile.data$datapath,header = input$genome.headerCK,sep = "\t",check.names = F,row.names = 1)
        )
        return(data.exp0)
    })
    
    ## markers data
    D.markers.raw <- reactive({
        req(data.exp0())
      
        D.markers <- data.exp0()
        
        D.markers <- as.matrix(D.markers)
        
        # rownames(D.markers) <- D.markers[,1]
        # D.markers <- D.markers[,-1]
        # class(D.markers) <- "numeric"
        
        return(D.markers)
        
    })
    
    D.markers.raw.show <- reactive({
      req(D.markers.raw())
      D.markers <- D.markers.raw()
      
      n.nrow <- dim(D.markers)[1]
      n.col <- dim(D.markers)[2]
      
      n.row <- ifelse(n.nrow >= 100,100,nrow)
      n.col <- ifelse(n.col >= 100,100,ncol)
      
      D.markers <- D.markers[1:n.row,1:n.col]
      
      return(D.markers)
      
    })
    ## marker view before processing
    output$markerView <- renderDT(
        D.markers.raw.show(),
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = FALSE),
        rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    
    ##++++2.1.2 data QC processing #####################################
    D.markers.pp <- reactive({
            req(input$dataPreprocessing)
            req(D.markers.raw())
            if(input$markerFormat != "Custom"){
                if(input$markerFormat == "Hapmap"){
                    G <- D.markers.raw()
                    rownames(G) <- G[,1]
                    G <- G[,-c(1:11)]
                    G <- t(G)
                } else if (input$markerFormat == "matrix"){
                    G <- D.markers.raw()
                }

                impute <- input$filter_imputation
                imputeMethod <- input$filter_imputation_method
                maf <- input$filter_maf
                mr <- input$filte_mr
                GSDataQC.res <- transLetter2number(G,maf = maf,mr = mr,impute = impute,imputeMethod = imputeMethod)
            }else{
              GSDataQC.res <- list(G.n = D.markers.raw(),maf.info = data.frame(Ref = rep(NA,ncol( D.markers.raw())),
                                                                               Alt = rep(NA,ncol( D.markers.raw())),
                                                                               AF = rep(NA,ncol( D.markers.raw())),
                                                                               MAF = rep(NA,ncol( D.markers.raw()))))
            }
              GSDataQC.res
          
        })
        
        ## observation for statistic and processed marker  #######
    observeEvent(input$dataPreprocessing,{output$sta_maf_Info <- renderDT(
      D.markers.pp()[[2]],
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    })
    
    D.markers.pro.show <- reactive({
      req(D.markers.pp())
      D.markers <- D.markers.pp()[[1]]
      
      n.nrow <- dim(D.markers)[1]
      n.col <- dim(D.markers)[2]
      
      n.row <- ifelse(n.nrow >= 100,100,nrow)
      n.col <- ifelse(n.col >= 100,100,ncol)
      
      D.markers <- D.markers[1:n.row,1:n.col]
      
      return(D.markers)
      
    })
    
    output$pro_marker <- renderDT(
      D.markers.pro.show(),
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    ##++++2.1.2 marker data download ######################
    output$download.MAF.Info <- downloadHandler(
      filename=function(){
        paste0("MAF-data-", Sys.Date(), ".csv", sep="")
      },
      content=function(filename){
        write.table(D.markers.pp()[[2]],filename,row.names = T,col.names = T,sep = ",",quote = F)
      },"csv"
    )
    
    output$download.numeric.data <- downloadHandler(
      filename=function(){
        paste0("Marker-data-", Sys.Date(), ".txt", sep="")
      },
      content=function(filename){
        write.table(D.markers.pp()[[1]],filename,row.names = T,col.names = T,sep = "\t",quote = F)
      },"text"
    )
  
## pheno ----------------------------------------------------------------------------
    ##++phenotype data input ###############################
    D.phenotype <-reactive({
      infile.target<-input$phenotypeData
      phenoFormat <- input$phenoFormat
      if (is.null(infile.target))
        return(NULL)
      else{
        if (input$phe.rownames) {
          D.phenotype <- switch(phenoFormat,
                                "CSV" = read.table(file=infile.target$datapath,header=input$phe.headerCK,sep = ",",check.names = F,stringsAsFactors = F,row.names = 1),
                                "table text" = read.table(file=infile.target$datapath,header=input$phe.headerCK,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
          )
        }else{
          D.phenotype <- switch(phenoFormat,
                                "CSV" = read.table(file=infile.target$datapath,header=input$phe.headerCK,sep = ",",check.names = F,stringsAsFactors = F),
                                "table text" = read.table(file=infile.target$datapath,header=input$phe.headerCK,sep = "\t",check.names = F,stringsAsFactors = F)
          )
        }
      }
      return(D.phenotype)
    })
    
    ##++phenotype datatable #######################
    output$pheRaw <- renderDT(
      D.phenotype(),
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    ##++phenotype interact selection UI ################
    observe({
      # infile.target<-input$phenotypeData
      if (is.null(input$phenotypeData)){
        output$traitSelectionPanel <- renderUI(textInput("traitSelection", label = h4("Select trait"), value = "No trait input ..."))
      }else{
      a <- colnames(D.phenotype())
      names(a) <- a
      output$traitSelectionPanel <- renderUI(selectInput("traitSelection", h4("Selection trait"), 
                                                                choices = as.list(a)))
      }
    })
    
    observe({
      # infile.target<-input$phenotypeData
      if (is.null(input$phenotypeData)){
        output$traitSelectionPanel1 <- renderUI(textInput("traitSelection1", label = "Select trait:", value = "No trait input ..."))
      }else{
        a <- colnames(D.phenotype())
        names(a) <- a
        output$traitSelectionPanel1 <- renderUI(selectInput("traitSelection1", "Selection trait:",
                                                           choices = as.list(a)))
      }
    })
    
    observe({
      # infile.target<-input$phenotypeData
      if (is.null(input$phenotypeData)){
        output$traitSelectionPanel2 <- renderUI(textInput("traitSelection2", label = h4("Select trait"), value = "No trait input ..."))
      }else{
        a <- colnames(D.phenotype())
        names(a) <- a
        output$traitSelectionPanel2 <- renderUI(selectInput("traitSelection2", h4("Selection trait"),
                                                           choices = as.list(a)))
      }
    })
    
    #++++1. phenotype visulization ##########################
    observe({
      # infile.target<-input$phenotypeData
        output$phenotypePlot <- renderUI(plotlyOutput("phenotypePlotR",width = paste0(input$plotPanelWidth,"%"), height = paste0(input$plotPanelHeight,"px")) %>% withSpinner(color="#0dc5c1"))
        output$genotypePlot <- renderUI(plotlyOutput("genotypePlotR",width = paste0(input$plotPanelWidth,"%"), height = paste0(input$plotPanelHeight,"px")) %>% withSpinner(color="#0dc5c1"))
    })
    #++++1.1  seperate plot #################################
    output$phenotypePlotR <- renderPlotly({
      req(D.phenotype())
      req(input$traitSelection)
      y <- D.phenotype()[,input$traitSelection]
      req(class(y) == "numeric")
      names.y <- input$traitSelection
      data2plot <- data.frame(phenotype = y,a = rep("a",length(y)))

      plot.type <- input$phenotypePlotType
    if (plot.type == "density") {
      g <- ggplot(data2plot, aes(y)) +
        geom_density(size =input$densityLineSizeD,color = input$densityLineColorD)+
        labs(title="Phenotype density plot", x=input$traitSelection, y="Density") +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",color = "black"),
              axis.text = element_text(size = 12),axis.title=element_text(size=12,face = "bold"),
              legend.position = "bottom",legend.text = element_text(size = 10),
              legend.title = element_text(size = 12,face = "bold"),legend.background = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
              # axis.text.x = element_blank(),axis.ticks.x = element_blank(),
              strip.background = element_blank(),strip.text = element_text(size =12),
              panel.spacing = unit(0.5, "lines"))
      plotly::ggplotly(g)
    } else if (plot.type == "histogram"){
      g <- ggplot(data2plot, aes(y)) +
        geom_histogram(aes(y = ..count..),bins = input$binsN, closed = "left",fill = input$binsColor,color = input$binsColor) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title="Phenotype bar density plot", x=input$traitSelection, y="Count") +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",color = "black"),
              axis.text = element_text(size = 12),axis.title=element_text(size=12,face = "bold"),
              legend.position = "bottom",legend.text = element_text(size = 10),
              legend.title = element_text(size = 12,face = "bold"),legend.background = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
              # axis.text.x = element_blank(),axis.ticks.x = element_blank(),
              strip.background = element_blank(),strip.text = element_text(size =12),
              panel.spacing = unit(0.5, "lines"))
      if (input$showDensityLine)
        g <- g + geom_freqpoly(binwidth = (range(y,na.rm = T)[2]-range(y,na.rm = T)[1])/input$binsN, color = input$densityLineColor,
                               size = input$densityLineSize)
      plotly::ggplotly(g)
    } else {
      g <- ggpairs(D.phenotype())+theme_bw() 
      plotly::ggplotly(g)
    }
    })
    
    #++++2. genotype visulization ##########################
    #++++++ structure func ############
    geneStructureCompute <- function(G.n,which,k=3,d=3){
      ##1. PCA
      res <- lapply(which,function(x){ switch (x,
                                               PCA = {
                                                 PCA_comp <- prcomp(G.n)
                                                 # get_eigenvalue(PCA_comp)
                                                 PCA_res <- as.data.frame(PCA_comp$x[,1:d])
                                                 PCA_res
                                               },
                                               UMAP = {
                                                 require(umap) 
                                                 umap <- umap(G.n, n_components = d) 
                                                 umap_res <- umap[["layout"]] 
                                                 umap_res <- as.data.frame(umap_res)
                                                 names(umap_res) <- paste0("Dimension",1:d)
                                                 umap_res
                                               },
                                               TSNE = {
                                                 require(Rcpp)
                                                 require(tsne)
                                                 tsne <- tsne(G.n, initial_dims = 10,k = d)
                                                 tsne_res <- data.frame(tsne)
                                                 names(tsne_res) <- paste0("Dimension",1:d)
                                                 tsne_res
                                               }
      )
      })
      names(res) <- which
      
      # if (is.null(label)) {
        if (k > 1) {
          hc = hclust(dist(G.n))
          groups <- as.factor(cutree(hc,k=k))
        }else{
          groups <- as.factor(rep(1,nrow(G.n)))
        }

        res <- lapply(res, function(x){
          res <- data.frame(x,groups = groups)
          res
        })
# 
#       }else{
#         res <- lapply(res, function(x){
#           res <- data.frame(x,groups = label[,1],colors = label[,2])
#           res
#         })
#       }
      return(res)
    }
    #++++++ read in label data ##############
    labelData <- reactive({
      infile.target <- input$groupInformation
      if(is.null(infile.target))
        return(NULL)
      else{
        labelData <- read.table(file = infile.data$datapath,header = F,sep = "\t",check.names = F,comment.char = "",sep = "\t")
      }
      return(labelData)
    })
  
    structureRes <- reactive({
      req(input$genotypeComputStart)
      req(D.markers.pp())
      req(input$genotypeComputType)
      # req(labelData())
    
      G.n <- D.markers.pp()[[1]]
      
      res <- geneStructureCompute(G.n,input$genotypeComputType, k = input$cutTreeCount,
                                  d = input$dimensionSelect)
      # names(res) <- input$genotypeComputType
      
      if(!is.null(labelData())){
        res <- lapply(res,function(x){
          # req(labelData())
          x <- data.frame(x[1:input$cutTreeCount],labelData())
          names(x)[(input$cutTreeCount+1):(input$cutTreeCount+2)] <- c("groups","color")
        })
      }else{
        res <- res
      }
      cat("PCA ect ... done!")
      return(res)
    })
    
    observe({
      a <- input$genotypeComputType
      names(a) <- a
      output$genotypePlotType <- renderUI(radioButtons("genotypePlotTypeSelection", h4("Select plot data:"), 
                                                       choices = as.list(a),inline = T))
    })
    
    output$PCARes <- renderDT(
      structureRes()[["PCA"]],
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    ## UMAP data ###################
    output$UMAPRes <- renderDT(
      structureRes()[["UMAP"]],
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    ## TSNE data #####################
    output$TSNERes <- renderDT(
      structureRes()[["TSNE"]],
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    #++++2.1. plot ##########################
    output$genotypePlotR <- renderPlotly({
      req(structureRes())
      req(input$genotypePlotTypeSelection)
      plotData <- structureRes()[[input$genotypePlotTypeSelection]]

      if (!is.null(labelData())) {
        colors_list <- unique(plotData$color)
        names(colors_list) <- plotData$groups[match(colors_list,plotData$color)]
      }else if (input$genotypePointColor != "") {
        color_list <- strsplit(input$genotypePointColor,",| ")[[1]]
      }else{
        color_list <- rainbow(n=3)
      }
      
      ## readin label data ###########
      
      color_list <- color_list
      
      fig1 <- plot_ly(plotData, x = plotData[,1], y = plotData[,2],
                     color = plotData$groups,
                     colors = color_list,
                     type = 'scatter', mode = 'markers',size = 2)
      fig1 <- fig1 %>% layout(
        xaxis = list( 
          title = names(plotData)[1]),  
        yaxis = list( 
          title =  names(plotData)[2])) 
      
      fig2 <- plot_ly(plotData, x = plotData[,1], y = plotData[,2], z = plotData[,3], color =  plotData$groups, 
                      colors = color_list,
                      size = 2,type = "scatter3d") 
      fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = names(plotData)[1]), 
                                           yaxis = list(title = names(plotData)[2]), 
                                           zaxis = list(title = names(plotData)[3]))) 
      
      subplot(fig1, fig2,margin = 0.1) %>% layout(scene = list(domain = list(x = c(0.5, 1), y = c(0,1))))
    })
    # extra matrix
    D.Z<-reactive({
        infile.file.Z<-input$file.Z
        
        if (is.null(infile.target))
            return(NULL)
        if (input$fix.rownamesCK) {
            D.Z<-as.matrix(read.table(file=infile.target$datapath,header=input$fix.headerCK,row.names = 1,stringsAsFactors = F))
        } else {
            D.Z<-as.matrix(read.table(file=infile.target$datapath,header=input$fix.headerCK,stringsAsFactors = F))
        }
        return(D.Z)
    })

    ## panel2 G2P and CV parameter setting #######################
    ##++++ G2P parameters manual ################
    modelsParameterModal <- function(){
      modalDialog(easyClose = TRUE,
                  wellPanel(
                            HTML('The details of parameters:
                                                     <ul><li><b>nIter, burnIn, thin:</b>(integer)The number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.</li>
                                                     <li><b>S0, df0:</b>(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or <a href="http://genomics.cimmyt.org/BGLR-extdoc.pdf">BGLR</a>.Default S0 NULL,df0 5.</li>
                                                     <li><b>R2:</b>(numeric, (0,1)) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or <a href="http://genomics.cimmyt.org/BGLR-extdoc.pdf">BGLR</a>.Defult 0.5.</li>
                                                     <li><b>ntree:</b>RandomForest parameter:Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.Defualt 500.</li>
                                                     <li><b>nodesize:</b>Randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).</li>
                                                     <li><b>kernel:</b>svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial)Default "linear".</li>                                                     
                                                     <li><b>gamma:</b>svm parameter parameter needed for all kernels except linear (default: 1/(data dimension)).</li>
                                                     <li><b>cost:</b>svm cost,default 2^(-9).</li>
                                                     </ul>')
                  ),
                  footer = tagList(
                    actionButton("okShowModelModal", "OK")
                  )
                  
      )
    }
    
    observeEvent(input$showModelModal, {
      showModal(modelsParameterModal())
    })
    
    observeEvent(input$okShowModelModal, {
      removeModal()
    })
    
    ##++++ eval parameters manual ################
    evalParameterModal <- function(){
      modalDialog(easyClose = TRUE,
                  wellPanel(
                           HTML('The details of parameters:
                                                     <ul><li><b>Threshold: top alpha(%):</b>(numeric,(0,100])A vector is the proportion of excellent individuals,default 1:90.</li>
                                                     <li><b>Bata:</b>(numeric)The parameter of "F1".</li>
                                                     <li><b>BestIndividuals:</b>(character)The position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".</li>
                                                     <li><b>Probability:</b>(logical)Whether the predScores is probability? Default FALSE.</li>
                                                     <li><b>AllNDCG:</b>(logical)Compute all NDCG results one by one?default FALSE.</li>
                                                     <li><b>probIndex:</b>(integer)indicates the column index which prediction result is probability.</li>
                                                     </ul>')
                  ),
                  footer = tagList(
                    actionButton("okShowEvalModal", "OK")
                  )
                  
      )
    }
    
    observeEvent(input$showEvalModal, {
      showModal(evalParameterModal())
    })
    
    observeEvent(input$okShowEvalModal, {
      removeModal()
    })
    
    
    observe({
      # req(input$model.group.R)
      # req(input$model.group.C)
      a <- input$modelsSelection
      req(a)
      names(a) <- a
      output$showSelectedModel <- renderUI(radioButtons("selectedModels", h4("The following model(s) have been selected:"), 
                                                       choices = as.list(a),inline = T))
      
    })
    
    ## evaluation render UI ###################
    observe({
      # req(input$model.group.R)
      # req(input$model.group.C)
      a <- input$measuresSelection
      req(a)
      names(a) <- a
      output$showSelectedMeasures <- renderUI(radioButtons("selectedMeasures", h4("The following evaluation metric(s) have been selected:"), 
                                                        choices = as.list(a),inline = T))
    })
    
    ## render selected Parameters ###############
    modelingDataInfoShow <- reactive({
      req(D.markers.pp())
      req(D.phenotype())
      
      dim.G <- dim(D.markers.pp()[[1]])
      length.p <- nrow(D.phenotype())
    
      res <- paste0("<h4><b>The marker matrix with ",dim.G[1]," rows and ",dim.G[2], " columns. <br />", "The phenotypic data have ",length.p," lines. <br /></h4>")
    })
    
    
    output$selectedParametersGPorCV <- renderPrint({
      HTML('<h4><b>Genotype to phenotype prediction (G2P) or CrossValidation (CV):</b>',input$GPorCV,'</h4>')
    })
    
    output$dataInformations <- renderPrint({
      HTML(modelingDataInfoShow())
    })
    
    observeEvent(input$GPorCV,{
      if (input$GPorCV == "CV"){
        output$selectedParametersCV <- renderPrint({
          HTML('<ul><li><b>Selected trait:</b>',input$traitSelection2,'</li>',
               '<li><b>Cross validation Methods:</b>',input$CVMethods,'</li>',
               '<li><b>Proportion of testing-sets:</b>',input$holdoutFold,'</li>',
               '<li><b>Repeats of CV (Holdout):</b>',input$holdoutRepeat,'</li>',
               '<li><b>Fold of K-fold CV:</b>',input$KfoldFold,'</li>',
               '<li><b>Repeats of CV (K-fold):</b>',input$KfoldRepeat,'</li>',
               '</ul>')})
      }else if (input$GPorCV == "G2P"){
        output$selectedParametersCV <- renderPrint({
          HTML('<ul><li><b>Selected trait:</b>',input$traitSelection1,'</li>',
               '<li><b>Training set index:</b>',input$trainIdx,'</li>',
               '<li><b>Testing set index:</b>',input$testIdx,'</li>',
               '</ul>')
        })
      }
    })
    
    output$selectedParametersModel <- renderPrint({
      HTML('<ul><li><b>The selected methods:</b>',input$modelsSelection,'</li>',
           '<li><b>nIter:</b>',input$G2P.nIter,'</li>',
           '<li><b>burnIn:</b>',input$G2P.burnIn,'</li>',
           '<li><b>thin:</b>',input$G2P.thin,'</li>',
           '<li><b>S0:</b>',input$G2P.S0,'</li>',
           '<li><b>df0:</b>',input$G2P.df0,'</li>',
           '<li><b>R2:</b>',input$G2P.R2,'</li>',
           '<li><b>alpha:</b>',input$G2P.alpha,'</li>',
           '<li><b>gamma:</b>',input$G2P.gamma,'</li>',
           '<li><b>cost:</b>',input$G2P.cost,'</li>',
           '<li><b>kernel:</b>',input$G2P.kernel,'</li>',
           '<li><b>ntree:</b>',input$G2P.ntree,'</li>',
           '<li><b>nodesize:</b>',input$G2P.nodesize,'</li>',
           '</ul>')
      
    })
    output$selectedParametersEval <- renderPrint({
      HTML('<ul><li><b>The selected evaluation metrics:</b>',input$measuresSelection,'</li>',
           '<li><b>Top alpha(%):</b>',input$eval.topAlpha,'</li>',
           '<li><b>Best individuals:</b>',input$eval.BestIndividuals,'</li>',
           '<li><b>Beta:</b>',input$eval.Beta,'</li>',
           '<li><b>All NDCG:</b>',input$eval.allNDCG,'</li>',
           '</ul>')
      
    })
    
    ## G2P prediciton module ########################
    G2P.res <- eventReactive(input$G2P.run,{
      req(D.markers.pp())
      req(D.phenotype())
      req(input$traitSelection1)
      req(input$trainIdx)
      req(input$testIdx)
      
      G <- D.markers.pp()[[1]]
      p <- D.phenotype()[,input$traitSelection1]
      # G <- G.n
      # p <- phe.test
      
      trainIdx <- eval(parse(text = input$trainIdx))
      testIdx <- eval(parse(text = input$testIdx))
      trainMarker <- G[trainIdx,]
      trainPheno <- p[trainIdx]
      testMarker <- G[testIdx,]
      testPheno <- p[testIdx]
      modelMethods <- input$modelsSelection
      
      # modelMethods <- gsub("RRBLUP","rrBLUP",modelMethods)
      
      # if(is.na(input$G2P.lambda)){lambda <- NULL}else{lambda <- input$G2P.lambda}
      if(is.na(input$G2P.S0)){S0 <- NULL}else{S0 <- input$G2P.S0}
      
      G2P.res <- G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker, testPheno = testPheno,
                     modelMethods = modelMethods,
                     # outputModel = input$G2P.outputModel,  # main parameters
                     nIter = input$G2P.nIter, burnIn = input$G2P.burnIn, thin = input$G2P.thin,
                     saveAt = "", S0 = S0, df0 = input$G2P.df0, R2 = input$G2P.R2, weights = NULL,
                     verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                     ntree = input$G2P.ntree,nodesize = input$G2P.nodesize,
                     kernel = input$G2P.kernel,gamma = input$G2P.gamma, cost = input$G2P.cost  # machine learing parameters
                     )
      if(is.null(rownames(G2P.res))){
        rownames(G2P.res) <- 1:nrow(testMarker)
      }else{
        rownames(G2P.res) <- rownames(testMarker)
      }
      # output$G2PCVprocess <- renderPrint(HTML('<h5><b>G2P is Done !!!</h5></b>'))
      G2P.res
    })
    
    output$G2PPredRes <- renderDT(
      G2P.res(),
      # predRes,
      extensions = 'Buttons',
      filter = 'top',
      options = list(dom = "Blfrtip",
                     buttons = list("copy",list(extend = "collection",
                                                buttons = c("csv", "excel"),
                                                text = "Download")),
                     lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                     pageLength = 10, autoWidth = FALSE),
      rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    ## G2P evalRes ######################

    G2PEval.res <- eventReactive(input$eval.run.G2P,{
      req(input$measuresSelection)
      req(G2P.res())
      predRes <- G2P.res()

      evalMethod <- input$measuresSelection
      evalMethod <- gsub("Pearson","pearson",evalMethod)
      evalMethod <- gsub("Kendall","kendall",evalMethod)
      evalMethod <- gsub("Spearman","spearman",evalMethod)
      evalMethod <- gsub("F-score","F1",evalMethod)
      evalMethod <- gsub("Accuracy","accuracy",evalMethod)

      evaluareTest <- evaluateGS(realScores = predRes[,1], predScores = predRes[,2:ncol(predRes)],
                               # evalMethod = c("pearson", "kendall","spearman","RE","Kappa","NDCG","meanNDCG",
                               #                "MSE","R2","F1","accuracy"),
                               evalMethod = evalMethod,
                               topAlpha = input$eval.topAlpha[1]:input$eval.topAlpha[2],
                               # allNDCG = input$eval.allNDCG,
                               BestIndividuals = input$eval.BestIndividuals,
                               Beta = input$eval.Beta,
                               globalAlpha = input$eval.globalAlpha)

      # evaluareTest
      namesEval <- names(evaluareTest)
      namesEval <- gsub("Pearson","pearson",namesEval)
      namesEval <- gsub("Kendall","kendall",namesEval)
      namesEval <- gsub("Spearman","spearman",namesEval)
      namesEval <- gsub("F-score","F1",namesEval)
      namesEval <- gsub("Accuracy","accuracy",namesEval)

      names(evaluareTest) <- namesEval
      evaluareTest
    })
    
    # overall arguments 
    observe({
      req(input$eval.topAlpha)
      output$showSelectAlpha <- renderUI(
      numericInput("visualOverallAlpha", label = h5("which threshold in overview ? (%)"), 
                   value = input$eval.topAlpha[2],min = input$eval.topAlpha[1],max = input$eval.topAlpha[2])
      )
    })
    
    # G2PEval.res <- eventReactive()
    # show eval data 
    # observeEvent(req(G2PEval.res()),{
    observe({
      req(G2PEval.res())
      req(input$visualOverallAlpha)
      a <- input$visualOverallAlpha
      evalres <- G2PEval.res()
      namesEval <- names(evalres)
      
      if(length(evalres) > 1 & "corMethods" %in% namesEval){
        overall <- do.call(rbind,lapply(evalres[-1],function(x){x[paste0("top",a),]}))
        rownames(overall) <- paste0(rownames(overall),"_top",a)
        overall <- rbind(evalres[["corMethods"]],overall)
      }else if(length(evalres) == 1 & "corMethods" %in% namesEval){
        overall <- evalres[["corMethods"]]
      }else{
        overall <- do.call(rbind,lapply(evalres,function(x){x[paste0("top",a),]}))
        rownames(overall) <- paste0(rownames(overall),"_top",a)
      }
      
      ## render DT 
      output$evalResOverall <- renderDT(
        overall,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 16, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("pearson" %in% namesEval & input$eval.globalAlpha){a1 <- evalres[["pearson"]]}else{a1 <- NULL}
      output$evalResPearson <- renderDT(
        # ifelse(("pearson" %in% namesEval) & input$eval.globalAlpha, evalres[["pearson"]],NULL),
        a1,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("kendall" %in% namesEval & input$eval.globalAlpha){a2 <- evalres[["kendall"]]}else{a2 <- NULL}
      output$evalResKendall <- renderDT(
        a2,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("spearman" %in% namesEval& input$eval.globalAlpha){a3 <- evalres[["spearman"]]}else{a3 <- NULL}
      output$evalResSpearman <- renderDT(
        a3,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("MSE" %in% namesEval& input$eval.globalAlpha){a4 <- evalres[["MSE"]]}else{a4 <- NULL}
      output$evalResMSE <- renderDT(
        a4,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("R2" %in% namesEval& input$eval.globalAlpha){a5 <- evalres[["R2"]]}else{a5 <- NULL}
      output$evalResR2 <- renderDT(
        a5,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("RE" %in% namesEval){a6 <- evalres[["RE"]]}else{a6 <- NULL}
      output$evalResRE <- renderDT(
        a6,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')

      if ("Kappa" %in% namesEval){a7 <- evalres[["Kappa"]]}else{a7 <- NULL}
      output$evalResKappa <- renderDT(
        a7,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("NDCG" %in% namesEval){a8 <- evalres[["NDCG"]]}else{a8 <- NULL}
      output$evalResNDCG <- renderDT(
        a8,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("meanNDCG" %in% namesEval){a9 <- evalres[["meanNDCG"]]}else{a9 <- NULL}
      output$evalResMeanNDCG <- renderDT(
        a9,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("F1" %in% namesEval){a10 <- evalres[["F1"]]}else{a10 <- NULL}
      output$evalResF1 <- renderDT(
        a10,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
      
      if ("accuracy" %in% namesEval){a11 <- evalres[["accuracy"]]}else{a11 <- NULL}
      output$evalResAccuracy <- renderDT(
        a11,
        extensions = 'Buttons',
        filter = 'top',
        options = list(dom = "Blfrtip",
                       buttons = list("copy",list(extend = "collection",
                                                  buttons = c("csv", "excel"),
                                                  text = "Download")),
                       lengthMenu = list(c(10, 20, -1), c(10, 20, "All")),
                       pageLength = 10, autoWidth = TRUE),
        rownames = TRUE, escape = FALSE, class = 'compact nowrap stripe hover')
    
    })
    
    ## G2P visualization #######################
    ## ++1.scatter plot ########
    ## function
    scatterPlot <- function (predmat,x1,x2,col.low = "blue",col.high = "red",col.mid = NULL,col = "blue",show.line = FALSE,
                             color_szie = FALSE,alpha = 0.8,make.plotly = TRUE,sizeRange) 
    {
      plot.dataframe <- as.data.frame(predmat[,c(x1,x2)])
      namesData <- colnames(plot.dataframe)
      colnames(plot.dataframe) <- c("x1","x2")
      if (color_szie) {
        g <- ggplot(plot.dataframe, aes(x = x1, y = x2,color = x1,size = x2)) + 
          geom_point(alpha = alpha) + 
          scale_colour_gradientn(colours = c(col.low,col.mid,col.high))+
          labs(x = namesData[1],y = namesData[2],color = namesData[1], size = namesData[2])+
          scale_size(name = namesData[2],range = sizeRange)
      }else{
        g <- ggplot(plot.dataframe, aes(x = x1, y = x2)) + 
          geom_point(alpha = alpha,color = col,size = sizeRange[1])+
          labs(x = namesData[1],y = namesData[2],color = namesData[1], size = namesData[2])
      }
      if(show.line){
        g <- g + geom_smooth(method = "lm",se = F)
      }
      g
      if (make.plotly){ 
        plotly::ggplotly(g)}else{g}
    }
    
    ## port
    output$visualization.scatter <- renderPlotly({
      if(is.null(G2P.res()) | input$method1 == " " | input$method2 == " "){
        return(NULL)
      }else{
        if(input$scatter.col.mid == 0){
          col.mid = NULL}else{
            col.mid = color()[input$scatter.col.mid]
          }
        if(input$scatter.col.low == 0){col.low = NULL}else{col.low = color()[input$scatter.col.low]}
        if(input$scatter.col.high == 0){col.high = NULL}else{col.high = color()[input$scatter.col.high]}
        predRes <- G2P.res()
      scatterPlot(predmat = predRes,x1 = input$method1,x2= input$method2,col.low = col.low,
                  col.high = col.high,col.mid = col.mid,col = input$scatter.color,show.line = input$scatter.showline,color_szie = input$scatterColor_size,
                  alpha = input$scatter.alpha,sizeRange = input$scatter.sizeRange,make.plotly = T
      )
      }
    })
    
    # UI
    observe({
      # infile.target<-input$phenotypeData      # input$modelsSelection
      if (is.null(input$modelsSelection)){
        output$scatter.selectM1 <- renderUI(selectInput("method1", label = h4("Method 1:"), value = " "))
      }else{
        output$scatter.selectM1 <- renderUI(selectInput("method1", h4("Method 1:"),
                                                           choices = as.list(c("realPhenScore",input$modelsSelection))))
      }
    })

    observe({
      # infile.target<-input$phenotypeData
      if (is.null(input$modelsSelection)){
        output$scatter.selectM2 <- renderUI(selectInput("method2", label = h4("Method 2:"), value = " "))
      }else{
        output$scatter.selectM2 <- renderUI(selectInput("method2", h4("Method 2:"),
                                                        choices = as.list(c("realPhenScore",input$modelsSelection))))
      }
    })
    
    ## ++2.bar plot ####
    ## function
    barPlot <- function(data,
                        xlab = "GS Methods",ylab = "Value",legendTittle = "Measures",
                        make.plotly = FALSE,other = "normal"){
      color <- rainbow(nrow(data))
      plot.dataframe <- melt(data)
      g <- ggplot(plot.dataframe, aes(x = X2, y = value, fill=X1)) + 
        geom_bar(stat="identity", position="dodge")+
        labs(x = xlab,y = ylab,fill = legendTittle)
      scale_colour_manual(values= color)
      switch (other,
              normal = {g <- g},
              parallel = {g <- g + coord_flip()},
              sector = {g <- g + coord_polar(theta = "x")},
              pie = {g <- g + coord_polar(theta = "y")}
      )
      
      if (make.plotly){ 
        plotly::ggplotly(g)}else{g}
    }
    ## port
    
    overall.plotdata <- reactive({
      req(G2PEval.res())
      req(input$visualOverallAlpha)
      a <- input$visualOverallAlpha
      evalres <- G2PEval.res()
      namesEval <- names(evalres)
      
      if(length(evalres) > 1 & "corMethods" %in% namesEval){
        overall <- do.call(rbind,lapply(evalres[-1],function(x){x[paste0("top",a),]}))
        rownames(overall) <- paste0(rownames(overall),"_top",a)
        overall <- rbind(evalres[["corMethods"]],overall)
      }else if(length(evalres) == 1 & "corMethods" %in% namesEval){
        overall <- evalres[["corMethods"]]
      }else{
        overall <- do.call(rbind,lapply(evalres,function(x){x[paste0("top",a),]}))
        rownames(overall) <- paste0(rownames(overall),"_top",a)
      }
    })
    
    output$visualization.bar <- renderPlotly({
      if(is.null(overall.plotdata())){
        return(NULL)
      }else{
        data <- overall.plotdata()
      }
      barPlot(data = data,xlab = input$globalBar.xlab,ylab = input$globalBar.ylab,
              legendTittle = input$globalBar.legend.title,
              other = input$globalBar.type,
              make.plotly = T
      )
    })
    
    ## ++3.line plot ######
    ## function 
    linePlot <- function(evalMat,xlab = "Top alpha (%)",ylab = "Value",legendTittle = "GS Method",
                         size = 2,alpha = 0.8,make.plotly = FALSE){
      color <- rainbow(ncol(evalMat))
      times <- nrow(evalMat)
      plot.dataframe <- melt(evalMat)
      g <- ggplot(plot.dataframe, aes(x = rep(1:times,ncol(evalMat)), y = value, colour=X2,group=X2)) + 
        geom_line(size=size,alpha = alpha)+
        labs(x = xlab,y = ylab,colour = legendTittle)+
        scale_colour_manual(values = color)
      if (make.plotly){ 
        plotly::ggplotly(g)}else{g}
    }
    ## port
    output$visualization.line <- renderPlotly({
      if(is.null(G2PEval.res()[[input$thresholdLine.measure]])){
        return(NULL)
      }else{
        evalMat <- G2PEval.res()[[input$thresholdLine.measure]]
      }
      linePlot(evalMat = evalMat,xlab = input$thresholdLine.xlab,ylab = input$thresholdLine.ylab,
               legendTittle = input$thresholdLine.legend.title,size = input$thresholdLine.width,
               alpha = input$thresholdLine.alpha,
               make.plotly = T
      )
    })
    
    ## UI
    observe({
      if (is.null(input$modelsSelection)){
        output$lineplotEval <- renderUI(selectInput("thresholdLine.measure",h4("Which metric is used to plot?"), value = " "))
      }else{
        output$lineplotEval <- renderUI(selectInput("thresholdLine.measure",h4("Which metric is used to plot?"), value = as.list(input$measuresSelection)))
      }
    })
    
    ## G2P CV module ###########################
    G2PCV.res <- eventReactive(input$G2PCV.run,{
      req(D.markers.pp())
      req(D.phenotype())
      req(input$traitSelection2)
      req(input$trainIdx)
      req(input$testIdx)
      
      G <- D.markers.pp()[[1]]
      p <- D.phenotype()[,input$traitSelection1]
      # G <- G.n
      # p <- phe.test
      
      trainIdx <- eval(parse(text = input$trainIdx))
      testIdx <- eval(parse(text = input$testIdx))
      trainMarker <- G[trainIdx,]
      trainPheno <- p[trainIdx]
      testMarker <- G[testIdx,]
      testPheno <- p[testIdx]
      modelMethods <- input$modelsSelection
      
      # modelMethods <- gsub("RRBLUP","rrBLUP",modelMethods)
      
      # if(is.na(input$G2P.lambda)){lambda <- NULL}else{lambda <- input$G2P.lambda}
      if(is.na(input$G2P.S0)){S0 <- NULL}else{S0 <- input$G2P.S0}
      
      G2P.res <- G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker, testPheno = testPheno,
                     modelMethods = modelMethods,
                     # outputModel = input$G2P.outputModel,  # main parameters
                     nIter = input$G2P.nIter, burnIn = input$G2P.burnIn, thin = input$G2P.thin,
                     saveAt = "", S0 = S0, df0 = input$G2P.df0, R2 = input$G2P.R2, weights = NULL,
                     verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                     ntree = input$G2P.ntree,nodesize = input$G2P.nodesize,
                     kernel = input$G2P.kernel,gamma = input$G2P.gamma, cost = input$G2P.cost  # machine learing parameters
      )
      if(is.null(rownames(G2P.res))){
        rownames(G2P.res) <- 1:nrow(testMarker)
      }else{
        rownames(G2P.res) <- rownames(testMarker)
      }
      # output$G2PCVprocess <- renderPrint(HTML('<h5><b>G2P is Done !!!</h5></b>'))
      G2P.res
    })
})
