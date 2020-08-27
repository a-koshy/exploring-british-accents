# Clean, extract, smooth and align NSCV and BNC MFCCs and formants
# This code takes around 3.5 hours to run through on a laptop with 4 cores.
# The sounds/NSCV/ folder should contain the NSCV files which are available from http://wrap.warwick.ac.uk/138368/ 
source('analysis/functions.R')

library(magrittr)
library(dplyr)
library(tuneR, warn.conflicts = F, quietly = T)
library(wrassp)
library(fdasrvf)
library(here)

set.seed(1)

# METADATA ---------- ####
## Clean BNC log ####
load('data/bnc_class_sounds.rda')
load('data/aa-index-to-remove.RData') # indices with bad sound quality
current.bnc$soundnum <- 1:nrow(current.bnc)
bnc_log <- current.bnc %>% 
    dplyr::filter(!grepl('India|European|NA|Zealand|Turkish|Chinese|
                  Canad|United States',
                         dialect),
                  age > 10,
                  duration >= 0.2 & duration <= 1,
                  !grepl('radio|television', activity),
                  !grepl('telvision|presenter|commentator|newsreader|sports reporter', 
                         occupation),  # telvision is not a typo
                  !is.na(placenamecleaned)) %>% 
    anti_join(data.frame(index.to.remove), by=c('index'='index.to.remove')) %>% 
    mutate(soc=factor(soc),
           persname=factor(persname),
           placenamecleaned=factor(placenamecleaned))
save(bnc_log, file='data/bnc_log.RData')


# First write wav files from audio object into a BNC directory
start <- Sys.time()
if(!dir.exists('sounds/BNC')) {dir.create('sounds/BNC')}
for (i in 1:nrow(bnc_log)){
    obs <- bnc_log$soundnum[i]
    wavfile <- paste0('sounds/BNC/bnc', i, '.wav')
    if (!file.exists(wavfile)){
        writeWave(sounds[[obs]], wavfile, extensible = FALSE)
    }
}
stop <- Sys.time()
stop-start


## BNC formants ####
bnc_formant <- list()
start <- Sys.time()
print(paste('Started extracting', nrow(bnc_log), 'BNC formants at:', start))
pb <- txtProgressBar(min = 0, max = nrow(bnc_log), style = 3)
for (i in 1:nrow(bnc_log)) {
    setTxtProgressBar(pb, i)
    obs <- bnc_log$soundnum[i]
    filename <- paste0('sounds/BNC/bnc', i, '.wav')
    f <- forest(filename, toFile = F)
    bnc_formant[[i]] <- f$fm
}
stop <- Sys.time()
stop-start
save(bnc_formant, file=here('data/bnc_raw_formant.RData'))


## BNC MFCCs ####
# takes around 40 min
bnc_mfcc <- list()
start <- Sys.time()
print(paste('Started extracting', nrow(bnc_log), 'BNC MFCCs at:', start))
pb <- txtProgressBar(min = 0, max = nrow(bnc_log), style = 3)
for (i in 1:nrow(bnc_log)) {
    setTxtProgressBar(pb, i)
    obs <- bnc_log$soundnum[i]
    wavfile <- paste0('sounds/BNC/bnc', i, '.wav')
    m <- sound_to_mfcc(wavfile)
    bnc_mfcc[[i]] <- m$mfcc
}
stop <- Sys.time()
stop-start
save(bnc_mfcc, file=here('data/bnc_raw_mfcc.RData'))


# Smooth BNC formants ####
# load('data/bnc_raw_formants.RData')
bnc_smooth_formant <- smooth_robust(bnc_formant, r=0.4, degree=1)
save(bnc_smooth_formant, file='data/bnc_smooth_formant.RData')

# Smooth BNC MFCCs ####
# load('data/bnc_raw_mfcc.RData')
bnc_smooth_mfcc <- smooth_mfccs(bnc_mfcc)
save(bnc_smooth_mfcc, file='data/bnc_smooth_mfcc.RData')

# Aligning all BNC words within speakers using MFCC 1
alignedbnc_mfcc <- bnc_warpings <- array(NA, dim=dim(bnc_smooth_mfcc$mfcc_norm))
alignedbnc_formant <- array(NA, dim=c(40, nrow(bnc_log), 4))
counter = 0
nspeaker <- bnc_log %>% 
    group_by(id) %>% 
    summarise(n=n(),
              soundnums=list(soundnum)) %>% 
    arrange(desc(n))

# Takes 45 min
start <- Sys.time()
for(i in nspeaker$id[lengths(nspeaker$soundnums)>1]){
    counter = counter+1
    print(paste0('Speaker ', i, ': ', counter,'/', 408))
    
    # align MFCC 1 within all words from speaker i
    indices <- which(bnc_log$id==i)
    alignspeaker <- align_fPCA(bnc_smooth_mfcc$mfcc_norm[,indices,1],
                               time=seq(0, 1, length.out = 40),
                               parallel = T, cores=4)
    # warp whole MFCC and formant matrix using the warping functions
    bnc_warpings[,indices,] <- alignspeaker$gam
    alignedbnc_mfcc[,indices,] <- align_matrix(alignspeaker$gam, 
                                               bnc_smooth_mfcc$mfcc_norm[,indices,])
    alignedbnc_formant[,indices,] <- align_matrix(alignspeaker$gam, 
                                                  bnc_smooth_formant$formant_norm[,indices,])
}
stop <- Sys.time()
stop-start

# single observations remain the same
for(i in nspeaker$id[lengths(nspeaker$soundnums)==1]){
    indices <- which(bnc_log$id==i)
    alignedbnc_mfcc[,indices,] <- bnc_smooth_mfcc$mfcc_norm[,indices,]
    alignedbnc_formant[,indices,] <- bnc_smooth_formant$formant_norm[,indices,]
}
save(alignedbnc_formant, alignedbnc_mfcc, bnc_warpings, file='data/bnc_aligned.RData')


# Preprocessing NSCV -------- ####

# create metadata
ids <- c('p1', 'p2', 'p3', 'p4') %>% rep(each=100)
mics <- c('dict', 'mic') %>% rep(each=50)
word <- c('class', 'grass', 'last', 'fast', 'pass') %>% rep(each=10)
accent <- c('S', 'N') %>% rep(each=5)
rep <- 1:5
nscv_log <- data.frame(id = ids, mic = mics, word, accent, rep) %>% 
    mutate(filename = paste0(here('sounds/NSCV/'), paste(id, word, tolower(accent), rep, mic, sep='-'),
                             '.wav'))

# Get durations
for(i in 1:nrow(nscv_log)){
    wavobj <- readWave(nscv_log$filename[i], header=T)
    nscv_log$duration[i] <- round(wavobj$samples/wavobj$sample.rate, digits=2)
}

# extract formants
print('Extracting formants')
nscv_formant <- list()
start <- Sys.time()
for (i in 1:nrow(nscv_log)){
    s <- forest(nscv_log$filename[i], toFile = F)
    nscv_formant[[i]] <- s$fm
}

# extract MFCCs
nscv_mfcc <- list()
start <- Sys.time()
print(paste('Started extracting', nrow(nscv_log), 'sound MFCCs at:', start))
pb <- txtProgressBar(min = 0, max = nrow(nscv_log), style = 3)
for (i in 1:nrow(nscv_log)) {
    setTxtProgressBar(pb, i)
    m <- sound_to_mfcc(nscv_log$filename[i])
    nscv_mfcc[[i]] <- m$mfcc
}
stop <- Sys.time()
print(stop-start)

# Smooth formants
print('Smoothing formants')
nscv_smooth_formant <- smooth_robust(nscv_formant, r=0.4, degree=1)

# Smooth mfccs
print('Smoothing MFCCs')
nscv_smooth_mfcc <- smooth_mfccs(nscv_mfcc)

save(nscv_formant, file=here('data/nscv_raw_formant.RData'))
save(nscv_mfcc, file=here('data/nscv_raw_mfcc.RData'))
save(nscv_smooth_formant, file=here('data/nscv_smooth_formant.RData'))
save(nscv_smooth_mfcc, file=here('data/nscv_smooth_mfcc.RData'))
save(nscv_log, file=here('data/nscv_log.RData'))


