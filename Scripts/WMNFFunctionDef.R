#### Define Function to Create List of Models (linear and GAM):

## Function to write out all possible model terms, to be combined in later funcitons:
build_model_terms <- function(x, by, op) {
  s_terms <- character()
  int_terms <- character()
  ## for loops to write all possible smooth terms:
  ## l is dummy variable specifying inclusion of id = 1:
  for (l in 1:2) {
    if (l == 1) {
      ## k is a dummy variable indicating whether or not theres an interaction term:
      for (k in 1:2) {
        if (k == 1) {
          for (i in 1:length(x)) {
            for (j in 1:length(by)) {
              s_termtext <- paste0("s(", x[i], ", k = 4", ", by =", by[j], " ,id = 1", ")")
              s_terms <- c(s_terms, s_termtext)
            }
          } 
        }
        else if (k == 2) {
          for (j in 1:length(by)) {
            s_termtext <- paste0("s(", x[1], ", ", x[2], ", k = 4", ", by =", by[j], " ,id = 1", ")")
            s_terms <- c(s_terms, s_termtext)
          }
        }
      }
    }
    else if (l == 2) {
      ## k is a dummy variable indicating whether or not theres an interaction term:
      for (k in 1:2) {
        if (k == 1) {
          for (i in 1:length(x)) {
            for (j in 1:length(by)) {
              s_termtext <- paste0("s(", x[i], ", k = 4", ", by =", by[j], ")")
              s_terms <- c(s_terms, s_termtext)
            }
          } 
        }
        else if (k == 2) {
          for (j in 1:length(by)) {
            s_termtext <- paste0("s(", x[1], ", ", x[2], ", k = 4", ", by =", by[j], ")")
            s_terms <- c(s_terms, s_termtext)
          }
        }
      }
    }
  }
  ## k is a dummy variable indicating whether there are 1 or 2 intercept terms
  for (k in 1:2) {
    if (k == 1) {
      for (i in 1:length(by)) {
        by1 <- by[! by %in% " "]
        if (!is.na(by1[i])) {
          int_termtext <- paste0(by[i])
          int_terms <- c(int_terms, int_termtext)
        }
        else {}
      }  
    }
    else if (k == 2) {
      for (l in 1:2) {
        int_termtext <- paste0(by[1], op[l], by[2])
        int_terms <- c(int_terms, int_termtext)
      }
    }
  }
  model_terms <- list(int_terms, s_terms)
  return(model_terms)
}

## Function builds models with one s() term - including those with two explanatory variables and an 
## interaction:
build_models_1 <- function(y, x, by, int_terms, s_terms) {
  models <- character()
  for (i in 1:length(s_terms)) {
    for (j in 1:length(int_terms)) {
      ## k is a dummy variable that indicates whether there are 1 or two intercept terms
      for (k in 1:2) {
        if (k == 1) {
          ## When we have only 1 intercept term:
          if (!grepl("\\+", int_terms[j]) && !grepl("\\*", int_terms[j])) {
            if (grepl(int_terms[j], s_terms[i]) && !grepl("SPPSTAGE", s_terms[i])) {
              newmodeltext <- paste0(y, "~ ", s_terms[i], " + ", int_terms[j])
              newmodel <- as.formula(parse(text = newmodeltext)[[1]])
              models <- c(models, newmodel)
            }
            else if (grepl("SPPSTAGE", s_terms[i])) {
              newmodeltext <- paste0(y, " ~ ", s_terms[i], " + ", "SPPSTAGE")
              newmodel <- as.formula(parse(text = newmodeltext)[[1]])
              models <- c(models, newmodel)
            }
            else if (grepl("by = )", s_terms[i]) && !grepl("SPPSTAGE", int_terms[j])) {
              newmodeltext <- paste0(y, " ~ ", s_terms[i], " + ", int_terms[j])
              newmodel <- as.formula(parse(text = newmodeltext)[[1]])
              models <- c(models, newmodel)
            }
          }  
        }
        else if (k == 2) {
          if (grepl("SPPSTAGE", s_terms[i])) {}
          else if (grepl("\\+", int_terms[j]) | grepl("\\*", int_terms[j])) {
            newmodeltext <- paste0(y, " ~ ", s_terms[i], " + ", int_terms[j])
            newmodel <- as.formula(parse(text = newmodeltext)[[1]])
            models <- c(models, newmodel)
          }
        }
      }
    }  
  }
  models <- unique(models)
  return(models)
}

## Function to build models for case with two s() terms, no interactions:
## I got a bit carried away with this one I think.
build_models_2 <- function(y, x, by, int_terms, s_terms) {
  models <- character()
  ## create a vector of DBH s() terms and one of HT s() terms
  dbh_terms <- character()
  ht_terms <- character()
  for (i in 1:length(s_terms)) {
    if (grepl("logDBHcm", s_terms[i]) && !grepl("logHTcm", s_terms[i])) {
      dbh_terms <- c(dbh_terms, s_terms[i])
    }
    else if (grepl("logHTcm", s_terms[i]) && !grepl("logDBHcm", s_terms[i])) {
      ht_terms <- c(ht_terms, s_terms[i])
    }
  }
  ## Make sure there are no duplicates
  dbh_terms <- unique(dbh_terms)
  ht_terms <- unique(ht_terms)
  
  for (i in 1:length(dbh_terms)) {
    for (l in 1:length(ht_terms)) {
      for (j in 1:length(int_terms)) {
        for (k in 1:2) {
          if (k == 1) {
            ## When we have only 1 intercept term (except with special case of SPP4STANDAGE):
            if (!grepl("\\+", int_terms[j]) && !grepl("\\*", int_terms[j])) {
              ## for the case without SPP4STANDAGE, require that both dbh_term and ht_term are represented in int_term:
              if ((grepl(int_terms[j], dbh_terms[i]) && grepl(int_terms[j], ht_terms[l])) | (grepl(int_terms[j], paste0(dbh_terms[i], ht_terms[l])) && grepl("by = )", paste0(dbh_terms[i], ht_terms[l]))) && !grepl("SPPSTAGE", paste0(dbh_terms[i], ht_terms[l]))) {
                newmodeltext <- paste0(y, " ~ ", dbh_terms[i], " + ", ht_terms[l], " + ", int_terms[j])
                newmodel <- as.formula(parse(text = newmodeltext)[[1]])
                models <- c(models, newmodel)
              }
              ## for the case with SPPSTAGE:
              else if (grepl("SPPSTAGE", paste0(dbh_terms[i], ht_terms[i]))) {
                ## for the case with SPP4STANDAGE and different by term, require that diff by term is represented in int_term
                if (grepl(int_terms[j], paste0(dbh_terms[i], ht_terms[l])) && int_terms[j] != "SPPSTAGE") {
                  newmodeltext <- paste0(y, " ~ ", dbh_terms[i], " + ", ht_terms[l], " + ", "SPPSTAGE", " + ", int_terms[j])
                  newmodel <- as.formula(parse(text = newmodeltext)[[1]])
                  models <- c(models, newmodel)
                }
                ## for the case with only SPPSTAGE as a by term, only SPPSTAGE as an int_term
                else if (grepl(int_terms[j], dbh_terms[i]) && grepl(int_terms[j], ht_terms[l])) {
                  newmodeltext <- paste0(y, " ~ ", dbh_terms[i], " + ", ht_terms[l], " + ", "SPPSTAGE")
                  newmodel <- as.formula(parse(text = newmodeltext)[[1]])
                  models <- c(models, newmodel)
                }
              }
              ## for the case with empty by terms
              else if (grepl("by = )", dbh_terms[i]) && grepl("by = )", ht_terms[l]) && !grepl("SPPSTAGE", int_terms[j])) {
                newmodeltext <- paste0(y, " ~ ", dbh_terms[i], " + ", ht_terms[l], " + ", int_terms[j])
                newmodel <- as.formula(parse(text = newmodeltext)[[1]])
                models <- c(models, newmodel)
              }
            }  
          }
          else if (k == 2) {
            if (grepl("SPPSTAGE", paste0(int_terms[j], dbh_terms[i], ht_terms[l]))) {}
            else if (grepl("\\+", int_terms[j]) | grepl("\\*", int_terms[j])) {
              newmodeltext <- paste0(y, " ~ ", dbh_terms[i], " + ", ht_terms[l], " + ", int_terms[j])
              newmodel <- as.formula(parse(text = newmodeltext)[[1]])
              models <- c(models, newmodel)
            }
          }
        }
      }  
    }
  }
  models <- unique(models)
  return(models)
}

## Function to put it all together:
allom_gams <- function(data1, y, x, by, rank = "AIC") {
  ## y, x and by must be input as concatenated string variables
  x <- c("logDBHcm", "logHTcm", "logPVcm3")
  by <- c("SPP4", "STANDAGE", "SPPSTAGE")
  by <- append(by, " ")
  op <- c(" + ", " * ")
  
  ## call function to build models for case with one explanatory variable
  model_terms <- build_model_terms(x = x, by = by, op = op)
  int_terms <- model_terms[[1]]
  s_terms <- model_terms[[2]]
  models_1 <- build_models_1(y = y, x = x, by = by, int_terms = int_terms, s_terms = s_terms)
  models_2 <- build_models_2(y = y, x = x, by = by, int_terms = int_terms, s_terms = s_terms)
  models <- c(models_1, models_2)
  
  ## fit gam models:
  for (i in 1:length(models)) {
    assign(paste0("model", i), gam(models[[i]], data = data1))
  }
  
  model_list <- character()
  for (i in 1:length(models)) {
    if (i < max(length(models))) {
      model_list <- paste0(model_list, "model", i, ", ")
    }
    else {
      model_list <- paste0(model_list, "model", i)
    }
  }
  results <- eval(parse(text = paste("model.sel(", model_list, ", rank =", rank, ")")))
  return(results)   
}

## Function to plot predictions:
plot_pred <- function(data, model, pred.length = 50, plotcov = c("logDBHcm, logHTcm"), tolerance = 3) {
  ## Drop unused levels from covariate vectors:
  data$SPP4 <- droplevels(data$SPP4)
  data$STANDAGE <- droplevels(data$STANDAGE)
  
  ## Build predicted values into dataframe:
  numspp <- length(levels(data$SPP4))
  numstage <- length(levels(data$STANDAGE))
  
  ## write in STANDAGE levels for each SPP
  stagevec <- rep(levels(data$STANDAGE), numspp)
  sppvec <- rep(levels(data$SPP4), each = numstage)
  
  ## combine into dataframe
  fac_data <- data.frame(sppvec, stagevec)
  
  ## get x scale for each STANDAGE and covariate
  STAGE <- levels(data$STANDAGE)
  minDBH <- numeric(length = length(STAGE))
  maxDBH <- numeric(length = length(STAGE))
  minHT <- numeric(length = length(STAGE))
  maxHT <- numeric(length = length(STAGE))
  ranges <- data.frame(STAGE, minDBH, maxDBH, minHT, maxHT)
  
  for (i in 1:length(levels(data$STANDAGE))) {
    newsub <- subset(data, STANDAGE == levels(data$STANDAGE)[i])
    ranges$minDBH[i] <- min(newsub$logDBHcm)
    ranges$maxDBH[i] <- max(newsub$logDBHcm)
    ranges$minHT[i] <- min(newsub$logHTcm)
    ranges$maxHT[i] <- max(newsub$logHTcm)
  }
  
  ## set up ylims 
  minlogABOVEg <- numeric()
  maxlogABOVEg <- numeric()
  for (j in 1:length(levels(data$SPP4))) {
    newsub <- subset(data, SPP4 == levels(data$SPP4)[j])
    for (i in 1:length(levels(data$STANDAGE))) {
      newsub1 <- subset(newsub, STANDAGE == levels(data$STANDAGE)[i])
      minlogABOVEg <- c(minlogABOVEg, min(newsub1$logABOVEg))
      maxlogABOVEg <- c(maxlogABOVEg, max(newsub1$logABOVEg))
    }
  }
  minlogABOVEg <- rep(minlogABOVEg, each = pred.length)
  maxlogABOVEg <- rep(maxlogABOVEg, each = pred.length)
  
  ## write in regularly spaces covariate values to predict upon
  logDBHcm <- numeric()
  logHTcm <- numeric()
  for (i in 1:nrow(fac_data)) {
    range <- subset(ranges, STAGE == stagevec[i])
    logDBHcm <- c(logDBHcm, seq(range$minDBH, range$maxDBH, length.out = pred.length))
    logHTcm <- c(logHTcm, seq(range$minHT, range$maxHT, length.out = pred.length))
  }
  
  SPP4 <- rep(fac_data$sppvec, each = pred.length)
  STANDAGE <- rep(fac_data$stagevec, each = pred.length)
  
  new_data <- data.frame(SPP4, STANDAGE, logDBHcm, logHTcm, minlogABOVEg, maxlogABOVEg)
  
  ## Create SPPSTAGE variable in new data for prediction
  for (i in 1:length(new_data$SPP4)) {
    new_data$SPPSTAGE[i] <- paste(new_data$SPP4[i], new_data$STANDAGE[i])
  }
  new_data$SPPSTAGE <- as.factor(new_data$SPPSTAGE)
  
  ## Eliminate covariate combinations in predictive data for which there was no observed data
  del_vec <- character()
  for (i in 1:nrow(new_data)) {
    if (!any(grepl(new_data$SPPSTAGE[i], levels(data$SPPSTAGE)))) {
      new <- rownames(new_data)[i]
      del_vec <- c(del_vec, new)
    }
  }
  
  del_vec <- as.numeric(del_vec[!is.na(del_vec)])
  new_data <- new_data[-del_vec,]
  
  ## write predictions to object
  predicted_data <- predict(model, newdata = new_data, se.fit = TRUE)
  
  ## write in predictions and standard errors to dataframe
  new_data$logABOVEg <- predicted_data$fit
  new_data$SE <- predicted_data$se.fit
  new_data$CI.low <- new_data$logABOVEg - predicted_data$se.fit
  new_data$CI.high <- new_data$logABOVEg + predicted_data$se.fit
  
  for (i in 1:nrow(new_data)) {
    if (new_data$CI.high[i] > new_data$maxlogABOVEg[i]+tolerance) {
      new_data$CI.high[i] <- NA
    }
    if (new_data$CI.low[i] < new_data$minlogABOVEg[i]-tolerance | new_data$CI.low[i] > new_data$maxlogABOVEg[i]+tolerance) {
      new_data$CI.low[i] <- NA
    }
    if (new_data$logABOVEg[i] > new_data$maxlogABOVEg[i]+tolerance) {
      new_data$logABOVEg[i] <- NA
    }
    else if (new_data$logABOVEg[i] < new_data$minlogABOVEg[i]-tolerance) {
      new_data$logABOVEg[i] <- NA
    }
  }
  
  plots <- list()
  if (grepl("logDBHcm", plotcov) && grepl("logHTcm", plotcov)) {
    DBHplot <- ggplot(new_data, aes(logDBHcm, logABOVEg)) +
      geom_line(size = 1) +
      geom_line(aes(logDBHcm, CI.high), color = "palevioletred1") +
      geom_line(aes(logDBHcm, CI.low), color = "palevioletred1") +
      geom_jitter(data = data, aes(logDBHcm, logABOVEg), size = 0.75) +
      facet_grid(SPP4 ~ STANDAGE, scales = "free")
    
    HTplot  <- ggplot(new_data, aes(logHTcm, logABOVEg)) +
      geom_line(size = 1) +
      geom_line(aes(logHTcm, CI.high), color = "palevioletred1") +
      geom_line(aes(logHTcm, CI.low), color = "palevioletred1") +
      geom_jitter(data = data, aes(logHTcm, logABOVEg), size = 0.75) +
      facet_grid(SPP4 ~ STANDAGE, scales = "free")
    print(DBHplot)
    print(HTplot)
  }
  else if (grepl("logHTcm", plotcov) && !grepl("logDBHcm", plotcov)) {
    HTplot  <- ggplot(new_data, aes(logHTcm, logABOVEg)) +
      geom_line(size = 1) +
      geom_line(aes(logHTcm, CI.high), color = "palevioletred1") +
      geom_line(aes(logHTcm, CI.low), color = "palevioletred1") +
      geom_jitter(data = allom_data_trunc, aes(logHTcm, logABOVEg), size = 0.75) +
      facet_grid(SPP4 ~ STANDAGE, scales = "free")
    print(HTplot)
  }
  else if (grepl("logDBHcm", plotcov) && !grepl("logHTcm", plotcov)) {
    DBHplot <- ggplot(new_data, aes(logDBHcm, logABOVEg)) +
      geom_line(size = 1) +
      geom_line(aes(logDBHcm, CI.high), color = "palevioletred1") +
      geom_line(aes(logDBHcm, CI.low), color = "palevioletred1") +
      geom_jitter(data = data, aes(logDBHcm, logABOVEg), size = 0.75) +
      facet_grid(SPP4 ~ STANDAGE, scales = "free")
    print(DBHplot)
  }
  else {
    print("error: plot cov must contain either logHTcm, logDBHcm or both")
  }
}
