## Sivateja Tangirala

## Baseline Functions 

library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(quantreg)


constructFormula <- function (exposure, data, depvar, adjustments) {
  cat(sprintf("constructing formula for %s regressed on %s and adjusted for %s \n",
              depvar, exposure, paste(adjustments, collapse = ' + ')))
  colIndices <- match(exposure,names(data))
  formula <- NULL
  if (is.null(adjustments)) {
    formula <-   {`if` (is.numeric(select(data,all_of(colIndices))[[1]]),
                        sprintf("I(scale(%s)) ~ I(scale(%s))", depvar, exposure),
                        sprintf("I(scale(%s)) ~ %s", depvar, exposure))}
  } else {
    formula <-   {`if` (is.numeric(select(data,all_of(colIndices))[[1]]),
                        sprintf("I(scale(%s)) ~ I(scale(%s)) + %s", depvar, exposure, paste(adjustments, collapse = ' + ')),
                        sprintf("I(scale(%s)) ~ %s + %s", depvar, exposure, paste(adjustments, collapse = ' + ')))}
  }
  cat(sprintf("Formula: %s \n", formula))
  return (list(formula,exposure,adjustments))
}

executeModel <- function (formula,data,exposure,depvar,adjustments) {
  mod <- NULL
  results <- NULL
  
  df_full <- data[,c(exposure,adjustments,depvar)]
  
  df_filter_NA_full <- df_full %>% drop_na()
  
  df_filter_NA_2_full <- df_filter_NA_full[, !duplicated(colnames(df_filter_NA_full), fromLast = TRUE)]
  
  data1 <- df_filter_NA_2_full
  
  cat(sprintf("Executing formula: %s \n", formula))
  cat(exposure)
  print(sprintf("\n The dimension of data1: %i",dim(data1)))
  tryCatch(
    mod1 <- lm(as.formula(formula),data = data1),
    error = function(e) {
      cat(sprintf("Failed on: %s \n", formula))
      rm(formula)
    }) 

  tryCatch(
      results1 <- mod1 %>% tidy(),
      error = function(e) {
        cat(sprintf("Failed on: %s  \n", formula))
        
      })
  tryCatch(
    {if (!is.null(results1)) {
      results_bind_1 <- cbind(results1,Phenotype= rep(depvar,times=nrow(results1)),SampleSize = rep(nrow(data1),times=nrow(results1)),Exposure=rep(exposure,times=nrow(results1)))
      rm(results1)
    }},error = function(e) {
      cat(sprintf("Failed on: %s  \n", formula))
      
    })
    
  tryCatch(
    return (as.data.frame(results_bind_1)),error = function(e) {
      cat(sprintf("Failed on: %s  \n", formula))
      
    }
    
  )
  
  
}

executeEWAS.map <- function  (data, depvar,adjustments,
                          exposures) {
  
  
  return_df  <- exposures %>%
    map(constructFormula, data=data, depvar, adjustments) %>%
    map_dfr(.f= function(x){executeModel(x[[1]],data,exposure =x[[2]],depvar,adjustments=x[[3]])})
  
  return(return_df)
}



EWAS <- function  (data, depvar,adjustments,exposures,
                   outFileName) {
  
  cat(sprintf("executing EWAS for %s", depvar))
  
  cols <- ncol(data)
 
  executeEWAS (data, depvar, adjustments,
           exposures, outFileName)
  
  
  
}


executeEWAS <- function (data, depvar, adjustments,
                     exposures, outFileName) {
  ## call executeEWAS.map and save results
  
  results <- executeEWAS.map(data,depvar=depvar,
                         adjustments=adjustments,
                         exposures=exposures)
  
  saveRDS(results,sprintf("%s.RDS",outFileName))
  
  
}


## code inspired from chiragjp/dhs_india/xwas



