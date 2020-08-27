## Set up data and libraries needed
set.seed(1)

## Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(fields)
library(fdasrvf) 
library(here)
library(foreach)
library(doParallel)
library(wrassp) 
library(RColorBrewer)
library(ggplot2)
library(gganimate)
library(scales)
library(kernlab)
library(rworldxtra)
library(INLA) #cite
library(mgcv) #cite
library(signal, warn.conflicts = F, quietly = T) # signal processing functions
library(oce, warn.conflicts = F, quietly = T) # image plotting functions and nice color maps
library(tuneR, warn.conflicts = F, quietly = T) # nice functions for reading and manipulating .wav files
library(glmnet)

# 9 colours for categorical variables
cols <- brewer.pal(9, 'Set1') %>% scales::alpha(alpha=0.5)
trcols <- brewer.pal(9, 'Set1') %>% scales::alpha(alpha=0.1)
opcols <- brewer.pal(9, 'Set1')

## Data
## NSCV data
load('data/nscv_smooth_formant.RData') # Smoothed NSCV formants
load('data/nscv_smooth_mfcc.RData') # Smoothed NSCV MFCCs
load('data/nscv_log.RData') # NSCV metadata

# BNC data 
load('data/locations.RData') # BNC location coordinates
load('data/bnc_aligned.RData') # Aligned BNC MFCCs and formants
load('data/bnc_log.RData') # BNC metadata

# Maps
data(countriesHigh)
uk.id <- grep("Kingdom", countriesHigh$NAME)
uk.nodes.full <- countriesHigh@polygons[[uk.id]]@Polygons[[1]]@coords
# remove last node (duplicated)
uk.nodes <- uk.nodes.full <- uk.nodes.full[-nrow(uk.nodes.full), ] 

### values == coarser approximation. 
### 0.01 creates a mesh that contains all the data
idx <- inla.simplify.curve(loc=uk.nodes.full, idx=1:nrow(uk.nodes.full), 
                           eps=0.01) 
border.nodes <- uk.nodes.full[idx,]
bn <- data.frame(lon=border.nodes[,1], lat=border.nodes[,2])

# Knots used in the soap film smoother for spatial smoothing
knotsfew <- read.csv('data/knotsfew.csv')

rm(list = c('border.nodes',  ls(pattern='uk'), 'countriesHigh', 'idx'))
