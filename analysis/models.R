# Models
source('analysis/setup.R')
source('analysis/functions.R')

# Functional linear model with formants ####
classgap <- list(gap=t(nscv_smooth_formant$formant[,,2]-nscv_smooth_formant$formant[,,1]),
                 time=matrix(seq(0,1, length.out = 40), ncol=40, nrow=400, byrow=T),
                 accent=nscv_log$accent,
                 f1=apply(nscv_smooth_formant$formant[,,1], 2, mean))

gapmodel <- gam((accent=='S') ~ s(time, by=gap, k=15, bs='cr')+f1,
                data=classgap, family=binomial())
summary(gapmodel)
plot(gapmodel, shade=TRUE)
abline(h=0)
AIC(gapmodel) # 57 

# Get confidence intervals
betas <- coef(gapmodel)
vb <- vcov(gapmodel, unconditional = T)
se <- sqrt(diag(vb))
betas[1] + (c(-1,1)*qnorm(0.975)*se[1]) # CI for intercept
betas[2] + (c(-1,1)*qnorm(0.975)*se[2]) # CI for F1

# Cross validate FLM using NSCV data
s <- split(1:400, nscv_log$id)
gacc <- c()
predrespflm <- predlinkflm <- gpreds <- matrix(NA, 100, 4)
testfit <- list(time=seq(0, 1, length.out = 40) %>% matrix(data=., nrow=40, ncol=40, byrow=F),
                gap= rep(1, times=40) %>% matrix(data=., nrow=40, ncol=40),
                f1=rep(0, times=40))
testbeta <- matrix(NA, nrow=40, ncol=4)
print('Cross validating FLM model:')
for (i in 1:length(s)) {
    print(paste('Fold', i))
    testind <- s[[i]]
    train <- list(gap=classgap$gap[-testind,], 
                  time=classgap$time[-testind,], 
                  accent=classgap$accent[-testind],
                  f1=classgap$f1[-testind])
    test <- list(gap=classgap$gap[testind,], 
                 time=classgap$time[testind,], 
                 accent=classgap$accent[testind],
                 f1=classgap$f1[testind])
    m <- gam((accent=='S') ~ s(time, by=gap, bs='cr', k=15)+f1, 
             data=train, family=binomial())
    
    testbeta[,i] <- predict(m, testfit, type='terms')[,2]/40 #discretised beta curve
    plot(m)
    abline(h=0)
    predrespflm[,i] <- predict(m, newdata=test, type='response') #response scale
    predlinkflm[,i] <- predict(m, newdata=test, type='link') #link scale
    
    gpreds[,i] <- ifelse(predrespflm[,i] > 0.5, 'S', 'N')
    gacc[i] <- sum(gpreds[,i]==test$accent)/length(testind)
}
print('Done')
gacc
mean(gacc) # 90.25%
range(predlinkflm) #-55, 161


#confusion matrix - truth on columns
predicted <- as.vector(gpreds)
# Confusion matrices for FLM predictions
table(predicted, classgap$accent)

# Penalised logistic regression with MFCCs ####
# PCA on MFCC curves together
nscv_centred_mfcc <- nscv_smooth_mfcc$mfcc_norm
nscv_centred_mfcc[,,1] <- apply(nscv_centred_mfcc[,,1], 2, function(x) x-mean(x))

# Stack MFCC curves into a single row
mfccstacked <- matrix(NA, nrow=400, ncol=40*40)
for(i in 1:40){
  mfccstacked[,(40*(i-1)+1):(40*i)] <- nscv_centred_mfcc[,,i] %>% t()
}

mfccpm <- prcomp(mfccstacked) # FPCA

# Full model with PC scores
pccvfull <- cv.glmnet(mfccpm$x, nscv_log$accent=='S', alpha=1, 
                      family=binomial())
pcamodfull <- glmnet(mfccpm$x, nscv_log$accent=='S', alpha=1, family=binomial(),
                     lambda=pccvfull$lambda.min)
pcamodfull$beta[1:20] %>% abs() %>% plot() # first 20 coefficients
(mfccpm$sdev^2/sum(mfccpm$sdev^2)) %>% plot() # eigenvalues for the eigenfunctions

# Cross validate PLR using NSCV data
pcacc <- c()
predrespplr <- predlinkplr <- pcpreds <- matrix(NA, nrow=100, ncol=4)
print('Cross validating PLR model:')
for (i in 1:4){
  print(paste('Fold', i))
  testind <- s[[i]]
  train <- mfccstacked[-testind,]
  test <- mfccstacked[testind,]
  pcatrain <- prcomp(train)
  
  pccv <- cv.glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                    alpha=1, family=binomial())
  pcamod <- glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                   alpha=1, family=binomial(), lambda=pccv$lambda.min)
  
  testpcscores <- predict(pcatrain, test)
  predrespplr[,i] <- predict(pcamod, testpcscores, type='response')
  predlinkplr[,i] <- predict(pcamod, testpcscores, type='link')
  pcpreds[,i] <- ifelse(predrespplr[,i] > 0.5, 1, 0)
  pcacc[i] <- sum(pcpreds[,i]==(nscv_log$accent[testind]=='S'))/
    length(testind)
}
print('Done')
pcacc
mean(pcacc) # 98.5 
range(predlinkplr) # -12, 15
hist(predlinkplr)

# confusion matrix
predicted <- as.vector(pcpreds)
# Confusion matrices for PLR predictions
table(predicted, nscv_log$accent)

# Combined effect of included PCs: using betas
nonzero <- which(pcamodfull$beta != 0)
betas <- matrix(data=pcamodfull$beta[nonzero], nrow=length(nonzero), ncol=400,
                byrow = F)

contrib <- array(NA, dim=c(40,40,length(nonzero)))
for (i in 1:length(nonzero)){
  contrib[,,i] <- pcamodfull$beta[nonzero[i]] * 
    mfccpm$rotation[,nonzero[i]] %>% 
    matrix(nrow=40, ncol=40, byrow=F)
}
sumcontrib <- apply(contrib, 1:2, sum)
#save(sumcontrib, file='output/sumcontrib.RData')


# Geographical variation model ####
# set up BNC data
nspeaker <- bnc_log %>% 
  group_by(id) %>% 
  summarise(n=n(),
            soundnums=list(soundnum)) %>% 
  arrange(desc(n))

# Formants
bncgap <- list(gap=t(alignedbnc_formant[,,2]-alignedbnc_formant[,,1]),
               time=matrix(seq(0, 1, length.out=40), ncol=40, nrow=nrow(bnc_log), byrow=T),
               f1=apply(alignedbnc_formant[,,1], 2, mean))

# Predictions for each sound
bncfpreds <- predict(gapmodel, bncgap, type='response') # FLM predictions
bncflink <- predict(gapmodel, bncgap, type='link')
bncfterms <- predict(gapmodel, bncgap, type='terms')
range(as.numeric(bncflink)) # -52, 35

# MFCCs
bncmfccstacked <- matrix(NA, nrow=nrow(bnc_log), ncol=40*40)
# Centre MFCC 1
bncmfccstacked[,1:40] <- apply(alignedbnc_mfcc[,,1], 2, function(x) x-mean(x)) %>% 
  t()
for(i in 2:40){
    bncmfccstacked[,(40*(i-1)+1):(40*i)] <- alignedbnc_mfcc[,,i] %>% t()
}

# Project on FPCs to get the scores
bncpcs <- predict(mfccpm, bncmfccstacked)

# Predictions for each sound
bncppreds <- predict(pcamodfull, bncpcs, type='response') # PLR predicted probs
bncplink <- predict(pcamodfull, bncpcs, type='link')

hist(bncplink)
hist(predlinkplr) #cross validated
hist(bncppreds)
range(bncplink) #-12,14 

table(bncppreds>0.5, bncfpreds>0.5) # agreement between FLM and PLR model predictions


# Compute predictions for each speaker using average sound
nspeakerpreds <- bnc_log %>% 
  mutate(ppred = as.numeric(bncppreds),
         fpred = as.numeric(bncfpreds),
         plink = as.numeric(bncplink),
         flink = as.numeric(bncflink)) %>% 
  group_by(id) %>% 
  summarise(n=n(),
            soundnums=list(soundnum),
            avgppred = plogis(mean(plink)),
            avgfpred = plogis(mean(flink)),
            ) %>% 
  arrange(desc(n))

# Constructing dataset of probs for the soap film
c <- bnc_log %>% 
  group_by(id, placenamecleaned) %>% 
  summarise(n_eachloc = n()) %>% 
  group_by(id) %>% 
  summarise(majloc = max(n_eachloc))

d <- bnc_log %>% 
  group_by(id, placenamecleaned) %>% 
  summarise(n_eachloc = n()) %>% 
  inner_join(c, by=c('id', 'n_eachloc' = 'majloc')) %>% 
  distinct(id, .keep_all = TRUE)

bnclocpred <- nspeakerpreds %>% 
  inner_join(d, by='id') %>% 
  inner_join(locations, by=c('placenamecleaned'='location.names'))

lf.inside <- point.in.polygon(bnclocpred$lon, bnclocpred$lat,
                              bn$lon, bn$lat)

# Beta regressions with soap film smoother
# FLM predictions
locationflm <- gam(avgfpred ~ s(lon, lat, k=40, bs='so', xt=list(bnd=list(bn))),
                   knots=knotsfew, data=bnclocpred[lf.inside==1,], 
                   family = betar())
summary(locationflm)

# PLR predictions
locationplr <- gam(avgppred ~ s(lon, lat, k=40, bs='so', 
                                        xt=list(bnd=list(bn))),
                   knots=knotsfew, data=bnclocpred[lf.inside==1,], 
                   family = betar())
summary(locationplr)
