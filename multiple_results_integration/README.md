## Code of integration of multi-model results

### 1. Download all R scripts in the path
### 2. Prepare a multi-model results from IP4GS or your own results like test file
### 3. Run integration like following sample code 
```R
source("Merge.r")
source("Ensemble.r")
source("Evaluate.r")

predRes <- as.matrix(read.table("predRes_test.csv",sep = ",",header = T,row.names = 1))

## merge 2 model prediction results ## 
MergeRes <- Merge(predResMat = predRes[,c("realPhenScore","RRBLUP","RFR")], autoOptimize = F,ratio = c(3,7))

## merge 2 model prediction results automatically ## 
MergeRes <- Merge(predResMat = predRes, autoOptimize = T)

## multiple predResults integraiton
EnsembleRes <- Ensemble(predMat = predRes, nrandom = 10, evalMethods = "pearson",
                            by = 0.1, evaluation = T)
```

### 4. Results of integration
* MergeRes：a matrix involve observed value of trait, merge models and the merge result.
* EnsembleRes
  * $BestWeight：the best weight of methods in all repeat}
  * $finalMat：the final matrix cbind predMat with final ensemble score}
  * $evalRes：the evaluation results of finalMat with evalMethods}
  * $weightMat：a weight matrix including all repeats}
  * $evalMat：a evaluation results matrix including all repeats}
