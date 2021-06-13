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


# Smooth BNC formants ####
# load('data/bnc_raw_formants.RData')
bnc_smooth_formant <- smooth_robust(bnc_formant, r=0.4, degree=1)

# Smooth BNC MFCCs ####
# load('data/bnc_raw_mfcc.RData')
bnc_smooth_mfcc <- smooth_mfccs(bnc_mfcc)

save(bnc_log, file='data/bnc_log.RData')
save(bnc_formant, file=here('data/bnc_raw_formant.RData'))
save(bnc_mfcc, file=here('data/bnc_raw_mfcc.RData'))
save(bnc_smooth_formant, file='data/bnc_smooth_formant.RData')
save(bnc_smooth_mfcc, file='data/bnc_smooth_mfcc.RData')


# Preprocessing NSCV -------- ####

# create metadata
ids <- c('p1', 'p2', 'p3', 'p4') %>% rep(each=100)
mics <- c('dict', 'mic') %>% rep(each=50)
word <- c('class', 'grass', 'last', 'fast', 'pass') %>% rep(each=10)
accent <- c('S', 'N') %>% rep(each=5)
rep <- 1:5
nscv_log <- data.frame(id = ids, mic = mics, word, accent, rep) %>% 
    mutate(filename = paste0('sounds/NSCV/', paste(id, word, tolower(accent), rep, mic, sep='-'),
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

# Align NSCV within accent and speaker
# align MFCC 1 - only one microphone together
# the other microphone for the same sound will use the same warping function
mic_obs = which(nscv_log$mic == 'mic')
micalign = time_warping(f=nscv_smooth_mfcc$mfcc_norm[,mic_obs,1], time=tseq, MaxItr = 100)
nscv_warpings = micalign$gam

nscv_aligned_mfcc <- array(NA, dim = c(40,400,40))
nscv_aligned_mfcc[,mic_obs,] = align_matrix(micalign$gam,
                                           nscv_smooth_mfcc$mfcc_norm[,mic_obs,])
nscv_aligned_mfcc[,-mic_obs,] = align_matrix(micalign$gam,
                                            nscv_smooth_mfcc$mfcc_norm[,-mic_obs,])
nscv_aligned_formant <- array(NA, dim=c(40,400,4))
nscv_aligned_formant[,mic_obs,] = align_matrix(micalign$gam,
                                              nscv_smooth_formant$formant_norm[,mic_obs,])
nscv_aligned_formant[,-mic_obs,] = align_matrix(micalign$gam,
                                               nscv_smooth_formant$formant_norm[,-mic_obs,])

# align BNC MFCC 1 curves to mean aligned MFCC1 curve from NSCV 
nscvmean = micalign$fmean
load('data/bnc_smooth_formant.RData')
load('data/bnc_smooth_mfcc.RData')

bnc_warpings = matrix(NA, nrow=40,
                      ncol=nrow(bnc_log))
# takes around 4 min
for (i in 1:nrow(bnc_log)) {
    bnc_warpings[,i] = pair_align_functions(f1=nscvmean, 
                                            f2=bnc_smooth_mfcc$mfcc_norm[,i,1], 
                                            time=tseq)$gam
}

# align whole MFCC matrix and formant matrix using BNC warpings
bnc_aligned_mfcc = align_matrix(bnc_warpings, bnc_smooth_mfcc$mfcc_norm)
bnc_aligned_formant = align_matrix(bnc_warpings, bnc_smooth_formant$formant_norm)


save(nscv_formant, file=here('data/nscv_raw_formant.RData'))
save(nscv_mfcc, file=here('data/nscv_raw_mfcc.RData'))
save(nscv_smooth_formant, file=here('data/nscv_smooth_formant.RData'))
save(nscv_smooth_mfcc, file=here('data/nscv_smooth_mfcc.RData'))
save(nscv_log, file=here('data/nscv_log.RData'))
save(nscv_aligned_formant, nscv_aligned_mfcc, nscv_warpings, nscvmean,
     file='data/nscv_aligned.RData')
save(bnc_aligned_formant, bnc_aligned_mfcc, bnc_warpings,
     file='data/bnc_aligned_to_nscv.RData')

# Prepare alignment within CV folds for cross validating models
s <- split(1:400, nscv_log$id)
tseq = seq(0, 1, length.out=40)
nscv_cv_aligned = list()
for (i in 1:4) {
    testind <- s[[i]]
    
    train = list(formant = nscv_smooth_formant$formant_norm[, -testind,],
                 mfcc = nscv_smooth_mfcc$mfcc_norm[, -testind,],
                 id = nscv_log$id[-testind],
                 mic = nscv_log$mic[-testind],
                 accent=nscv_log$accent[-testind])
    test = list(formant = nscv_smooth_formant$formant_norm[,testind,],
                mfcc = nscv_smooth_mfcc$mfcc_norm[,testind,],
                id = nscv_log$id[testind],
                mic = nscv_log$mic[testind],
                accent=nscv_log$accent[testind])
    
    # align curves within training set
    # align MFCC 1 - only one microphone together
    mic_obs = which(train$mic == 'mic')
    micalign = time_warping(f=train$mfcc[,mic_obs, 1], time=tseq, MaxItr = 20)
    train_warpings = micalign$gam
    
    train_aligned_formant <- array(NA, dim=c(40,300,4))
    train_aligned_formant[,mic_obs,] = align_matrix(micalign$gam,
                                                    train$formant[,mic_obs,])
    train_aligned_formant[,-mic_obs,] = align_matrix(micalign$gam,
                                                     train$formant[,-mic_obs,])
    train_aligned_mfcc <- array(NA, dim = c(40,300,40))
    train_aligned_mfcc[,mic_obs,] = align_matrix(micalign$gam,
                                                 train$mfcc[,mic_obs,])
    train_aligned_mfcc[,-mic_obs,] = align_matrix(micalign$gam,
                                                  train$mfcc[,-mic_obs,])
    
    # align test curves
    train_mean = micalign$fmean
    test_warpings = matrix(NA, nrow=40,
                           ncol=50)
    testmic = which(nscv_log$mic[testind] == 'mic')
    for (j in 1:length(testmic)) {
        test_warpings[,j] = pair_align_functions(f1=train_mean, 
                                                 f2=test$mfcc[,testmic[j],1], 
                                                 time=tseq)$gam
    }
    test_aligned_formant <- array(NA, dim=c(40,100,4))
    test_aligned_mfcc <- array(NA, dim=c(40,100,40))
    test_aligned_formant[,testmic,] = align_matrix(test_warpings, test$formant[,testmic,])
    test_aligned_formant[,-testmic,] = align_matrix(test_warpings, test$formant[,-testmic,])
    
    test_aligned_mfcc[,testmic,] = align_matrix(test_warpings, test$mfcc[,testmic,])
    test_aligned_mfcc[,-testmic,] = align_matrix(test_warpings, test$mfcc[,-testmic,])
    
    # save aligned test and train formants
    nscv_cv_aligned[[i]] <- list()
    nscv_cv_aligned[[i]]$train_aligned_formant = train_aligned_formant
    nscv_cv_aligned[[i]]$test_aligned_formant = test_aligned_formant
    nscv_cv_aligned[[i]]$train_aligned_mfcc = train_aligned_mfcc
    nscv_cv_aligned[[i]]$test_aligned_mfcc = test_aligned_mfcc
    
}

save(nscv_cv_aligned, file='data/nscv_aligned_cv.RData')
