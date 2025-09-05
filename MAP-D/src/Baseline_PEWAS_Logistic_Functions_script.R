## Sivateja Tangirala
## 06/17/2020
## Baseline Functions (logistic regression)

# Load required packages (individual tidyverse components)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(gdata)
  library(RNOmni)
  library(biglm)
  library(broom)
  library(pROC)
})


constructFormulaLogistic <- function (exposure, data, depvar, adjustments) {
  cat(sprintf("constructing formula for %s regressed on %s and adjusted for %s \n",
              depvar, exposure, paste(adjustments, collapse = ' + ')))
  formula <- NULL
  if (!is.null(adjustments)) {
    formula <-   {`if` (is.numeric(data[,exposure]),
                        sprintf("%s ~ I(scale(%s)) + %s", depvar, exposure, adjustments = paste(adjustments, collapse = ' + ')),
                        sprintf("%s ~ %s + %s", depvar, exposure, adjustments = paste(adjustments, collapse = ' + ')))}
  } else {
    formula <-   {`if` (is.numeric(data[,exposure]),
                        sprintf("%s ~ I(scale(%s))", depvar, exposure),
                        sprintf("%s ~ %s", depvar, exposure))}
  }
  cat(sprintf("Formula: %s \n", formula))
  return (list(formula,exposure,adjustments))
}

executeModelLogistic <- function (formula,data,exposure,depvar,adjustments) {
  mod <- NULL
  results <- NULL
  
  df_full <- data[,c(exposure,adjustments,depvar)]
  
  df_filter_NA_full <- df_full %>% drop_na()
  
  df_filter_NA_2_full <- df_filter_NA_full[, !duplicated(colnames(df_filter_NA_full), fromLast = TRUE)]
  
  # df_filter_NA_3_full <- drop.levels(df_filter_NA_2_full)
  
  data1 <- df_filter_NA_2_full
  
  cat(sprintf("Executing formula: %s \n", formula))
  cat(exposure)
  print(sprintf("\n The dimension of data1: %i",dim(data1)))
  tryCatch(
    mod1 <- glm(as.formula(formula),data=data1,family = "binomial"),
    error = function(e) {
      cat(sprintf("Failed on: %s \n", formula))
      rm(formula)
    }) 
  tryCatch(
    mod2 <- glm(as.formula(sprintf("%s ~ %s", depvar,adjustments = paste(adjustments, collapse = ' + '))),data=data1,family="binomial") ,
    error = function(e) {
      cat(sprintf("mod2 Failed on: %s \n", formula))
      rm(formula)
    })
  if(exists("mod1") & exists("mod2")){
    if (!is.null(mod1) & !is.null(mod2)){
      tryCatch({
        results1 <- mod1 %>% tidy()
      },
      error = function(e) {
        cat(sprintf("Failed on: %s  \n", formula))
        
      })
      if (!is.null(results1)) {
        tryCatch({results2 <- mod2 %>% tidy()
        probs_1 <- predict(mod1,type=c("response"))
        data1$probs_1 <- probs_1
        xwas_roc_1 <- pROC::roc(as.formula(sprintf("%s ~ probs_1", depvar)), data = data1) 
        auc_1 <- gsub('.*: ',"",xwas_roc_1$auc)
        
        probs_2 <- predict(mod2,type=c("response"))
        data1$probs_2 <- probs_2
        xwas_roc_2 <- pROC::roc(as.formula(sprintf("%s ~ probs_2", depvar)), data = data1) 
        auc_2 <- gsub('.*: ',"",xwas_roc_2$auc)
        
        results_bind_1 <- cbind(results1,AUC = rep(auc_1,times=nrow(results1)),AUCadjVariables = rep(auc_2,times=nrow(results1)),Phenotype= rep(depvar,times=nrow(results1)),SampleSize = rep(nrow(data1),times=nrow(results1)),Exposure=rep(exposure,times=nrow(results1)))
        rm(results1,results2)},error = function(e){
          cat(sprintf("Failed on: roc --  %s or %s  \n", sprintf("%s ~ probs_1", depvar),sprintf("%s ~ probs_2", depvar)))
          results_bind_1 <- cbind(results1,AUC = rep(NA,times=nrow(results1)),AUCadjVariables = rep(NA,times=nrow(results1)),Phenotype= rep(depvar,times=nrow(results1)),SampleSize = rep(nrow(data1),times=nrow(results1)),Exposure=rep(exposure,times=nrow(results1)))
          
        })
      }
      
    }
    if(exists("results_bind_1")){
      return(as.data.frame(results_bind_1))
    }
  }
  
}

executeEWASLogistic.map <- function  (data, depvar,adjustments,
                                      exposures) {
  
  
  return_df  <- exposures %>%
    map(constructFormulaLogistic, data=data, depvar, adjustments) %>%
    map_dfr(.f= function(x){executeModelLogistic(x[[1]],data,exposure =x[[2]],depvar,adjustments=x[[3]])})
  
  return(return_df)
}



EWASLogistic <- function  (data, depvar,adjustments,exposures,
                           outFileName) {
  
  cat(sprintf("executing EWAS (logistic) for %s", depvar))
  
  cols <- ncol(data)
  
  executeEWASLogistic (data, depvar, adjustments,
                       exposures, outFileName)
  
  
  
}


executeEWASLogistic <- function (data, depvar, adjustments,
                                 exposures, outFileName) {
  ## call executeEWASLogistic.map and save results
  
  results <- executeEWASLogistic.map(data,depvar=depvar,
                                     adjustments=adjustments,
                                     exposures=exposures)
  
  saveRDS(results,sprintf("%s.RDS",outFileName))
  
  
}





