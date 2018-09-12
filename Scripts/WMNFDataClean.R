
## load required packages:

library(here)
library(mgcv)
library(plyr)
library(visreg)
library(ggplot2)
library(MuMIn)
library(reticulate)

## read in the main data
allom_data <- read.csv("WMNF_allom5.csv", header = TRUE)
str(allom_data)

allom_data[, "STANDAGE"] <- as.factor(allom_data[, "STANDAGE"])

## Remove the data from the 98 year old stand so that all of the data are comparable:
allom_data_Y <- subset(allom_data, STANDAGE != "98years")
str(allom_data_Y)

max(allom_data_Y$DBHcm)

## get subset of data which includes trees in 98year stand which are smaller than the 
## largest trees in the younger stands.

allom_data_trunc <- subset(allom_data, DBHcm <= 66)

## log transform aboveground biomass:
allom_data_trunc$logABOVEg <- log(allom_data_trunc$ABOVEg)
ggplot(allom_data_trunc, aes(x = DBHcm, y = logABOVEg)) +
  geom_jitter() +
  geom_smooth(size = 0.5)

## lets try taking the log transformation of both variables:
allom_data_trunc$logDBHcm <- log(allom_data_trunc$DBHcm)
ggplot(allom_data_trunc, aes(x = logDBHcm, y = logABOVEg)) +
  geom_jitter() +
  geom_smooth(size = 0.5)

## Create SPPSTAGE in order to fully cross in GAM
for (i in 1:length(allom_data_trunc$SPP4)) {
  allom_data_trunc$SPPSTAGE[i] <- paste(allom_data_trunc$SPP4[i], allom_data_trunc$STANDAGE[i])
}

allom_data_trunc$SPPSTAGE <- as.factor(allom_data_trunc$SPPSTAGE)

## log transform HTcm:
allom_data_trunc$logHTcm <- log(allom_data_trunc$HTcm)

## log transform PVcm3
allom_data_trunc$logPVcm3 <- log(allom_data_trunc$PVcm3)

## remove NA values of biomass:
allom_data_trunc <- subset(allom_data_trunc, logABOVEg != "NA")

## look at how many data points we have for each species.  ACRU and PIRU seem low. 
## Let's create a dataset with them removed to compare to one with them included.

count(allom_data_trunc$SPP4)

allom_data_truncSPP4 <- subset(allom_data_trunc, SPP4 != "ACRU")
allom_data_truncSPP4 <- subset(allom_data_truncSPP4, SPP4 != "PIRU")
count(allom_data_truncSPP4$SPP4)

## make another subset of the data with only SPPSTAGE levels that have more than 10 data points.

SPPSTAGE_counts <- data.frame(count(allom_data_trunc$SPPSTAGE))
SPPSTAGE_counts <- subset(SPPSTAGE_counts, freq >= 10)

SPPSTAGE_list <- as.factor(SPPSTAGE_counts$x)

allom_data_truncSPPSTAGE <- subset(allom_data_trunc, (allom_data_trunc$SPPSTAGE %in% SPPSTAGE_list) == TRUE)

allom_data_truncSPPSTAGE$SPPSTAGE

## Try STANDAGE as numerical

allom_data_truncSTANDAGENUM <- allom_data_trunc
allom_data_truncSTANDAGENUM$STANDAGE <- as.numeric(allom_data_trunc$STANDAGE)

