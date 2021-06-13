# Models
source('analysis/setup.R')
source('analysis/functions.R')

# Functional linear model with formants ####

# Set up data
shorter = c(10, 20, 30)
classgap <- list(gap=t(nscv_aligned_formant[,,2]-nscv_aligned_formant[,,1]),
                 time=matrix(seq(0,1, length.out = 40), ncol=40, nrow=400, byrow=T),
                 accent=as.factor(nscv_log$accent),
                 id=as.factor(nscv_log$id),
                 mic=as.factor(nscv_log$mic),
                 word=as.factor(nscv_log$word),
                 f2=apply(nscv_aligned_formant[,,2], 2, mean),
                 f1=apply(nscv_aligned_formant[,,1], 2, mean),
                 f1f = t(nscv_aligned_formant[,,1]),
                 f2f = t(nscv_aligned_formant[,,2]),
                 f1mid = (nscv_aligned_formant[20,,1]),
                 f2mid = (nscv_aligned_formant[20,,2]),
                 f1shorter = t(nscv_aligned_formant[shorter, , 1]),
                 f2shorter = t(nscv_aligned_formant[shorter, , 2]),
                 timeshorter = matrix(seq(0,1, length.out = length(shorter)), ncol=length(shorter), nrow=400, byrow=T),
                 gapshorter = t(nscv_aligned_formant[shorter,,2] - nscv_aligned_formant[shorter,,1]))

classgap_unaligned <- list(gap=t(nscv_smooth_formant$formant_norm[,,2]-nscv_smooth_formant$formant_norm[,,1]),
                           time=matrix(seq(0,1, length.out = 40), ncol=40, nrow=400, byrow=T),
                           accent=as.factor(nscv_log$accent),
                           id=as.factor(nscv_log$id),
                           mic=as.factor(nscv_log$mic),
                           word=as.factor(nscv_log$word),
                           idaccent= as.factor(paste0(nscv_log$id, nscv_log$accent)),
                           f2=apply(nscv_smooth_formant$formant_norm[,,2], 2, mean),
                           f1=apply(nscv_smooth_formant$formant_norm[,,1], 2, mean),
                           f1f = t(nscv_smooth_formant$formant_norm[,,1]),
                           f2f = t(nscv_smooth_formant$formant_norm[,,2]),
                           f1mid = (nscv_smooth_formant$formant_norm[20,,1]),
                           f2mid = (nscv_smooth_formant$formant_norm[20,,2]),
                           f1shorter = t(nscv_smooth_formant$formant_norm[shorter, , 1]),
                           f2shorter = t(nscv_smooth_formant$formant_norm[shorter, , 2]),
                           timeshorter = matrix(seq(0,1, length.out = length(shorter)), ncol=length(shorter), nrow=400, byrow=T),
                           gapshorter = t(nscv_smooth_formant$formant_norm[shorter,,2] - nscv_smooth_formant$formant_norm[shorter,,1]))

formant_model_formula = (accent=='S') ~ s(time, by=f2f, k=15, bs='cr')+s(id, bs='re')

formant_model = gam(formant_model_formula, data=classgap, family = binomial())
summary(formant_model)
plot(formant_model, shade=T)
AIC(formant_model)

betas = coef(formant_model)
vb = vcov(formant_model, unconditional = T)
se <- sqrt(diag(vb))
betas[1] + (c(-1,1)*qnorm(0.975)*se[1]) # CI for intercept
gam.vcomp(formant_model)^2 # variance component from random effect is 0.006, 
# 2.126 from functional term

# Cross validate FLM using NSCV data
formant_model_cv = check_model(formant_model_formula, fulldata=classgap, align=T,
                               optimizer=c('outer', 'bfgs'))
mean(formant_model_cv$cv_accuracy) # accuracy 98
range(formant_model_cv$cv_link) # -900 to 866 on linear predictor scale

# Confusion matrix for FLM predictions
predicted_flm = as.vector(formant_model_cv$cv_preds)
table(predicted_flm, nscv_log$accent)

# Penalised logistic regression with MFCCs ####
# PCA on MFCC curves together
nscv_centred_mfcc <- nscv_aligned_mfcc
nscv_centred_mfcc[,,1] <- apply(nscv_centred_mfcc[,,1], 2, function(x) x-mean(x))

# Stack MFCC curves into a single row
mfccstacked <- matrix(NA, nrow=nrow(nscv_log), ncol=40*40)
for(i in 1:40){
  mfccstacked[,(40*(i-1)+1):(40*i)] <- nscv_centred_mfcc[,,i] %>% t()
}

mfccpm <- prcomp(mfccstacked) # FPCA

# Model with PC scores with full dataset
pccvfull <- cv.glmnet(mfccpm$x, nscv_log$accent=='S', alpha=1, 
                      family=binomial())
pcamodfull <- glmnet(mfccpm$x, nscv_log$accent=='S', alpha=1, family=binomial(),
                     lambda=pccvfull$lambda.min)
pcamodfull$beta[1:20] %>% abs() %>% plot() # first 20 coefficients
(mfccpm$sdev^2/sum(mfccpm$sdev^2)) %>% plot() # eigenvalues for the eigenfunctions

# Cross validate PLR using NSCV data
pcacc <- c()
predrespplr <- predlinkplr <- pcpreds <- matrix(NA, nrow=100, ncol=4)
nonzero_cv <- betas_cv <- pcs_cv <- list()
s <- split(1:400, nscv_log$id)
print('Cross validating PLR model:')
for (i in 1:4){
  print(paste('Fold', i))
  testind <- s[[i]]
  
  train_centred_mfcc <- nscv_cv_aligned[[i]]$train_aligned_mfcc
  train_centred_mfcc[,,1] <- apply(train_centred_mfcc[,,1], 2, function(x) x-mean(x))
  
  # Stack training MFCC curves into a single row
  trainstacked <- matrix(NA, nrow=300, ncol=40*40)
  for(j in 1:40){
    trainstacked[,(40*(j-1)+1):(40*j)] <- train_centred_mfcc[,,j] %>% t()
  }
  
  test_centred_mfcc <- nscv_cv_aligned[[i]]$test_aligned_mfcc
  test_centred_mfcc[,,1] <- apply(test_centred_mfcc[,,1], 2, function(x) x-mean(x))
  
  # Stack testing MFCC curves into a single row
  teststacked <- matrix(NA, nrow=100, ncol=40*40)
  for(j in 1:40){
    teststacked[,(40*(j-1)+1):(40*j)] <- test_centred_mfcc[,,j] %>% t()
  }
  
  train <- trainstacked
  test <- teststacked
  pcatrain <- prcomp(train)
  
  pccv <- cv.glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                    alpha=1, family=binomial())
  pcamod <- glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                   alpha=1, family=binomial(), lambda=pccv$lambda.min)
  
  testpcscores <- predict(pcatrain, test)
  predrespplr[,i] <- predict(pcamod, testpcscores, type='response')
  predlinkplr[,i] <- predict(pcamod, testpcscores, type='link')
  pcpreds[,i] <- ifelse(predrespplr[,i] > 0.5, 1, 0)
  pcacc[i] <- sum(pcpreds[,i]==(nscv_log$accent[testind]=='S'))/length(testind)
  
  # Save chosen FPCs, and weighted loadings from each fold
  nonzero_cv[[i]] <- which(pcamod$beta != 0)
  pcs_cv[[i]] <- array(NA, dim=c(40,40,length(nonzero_cv[[i]])))
  for (j in 1:length(nonzero_cv[[i]])){
    pcs_cv[[i]][,,j] <- pcamod$beta[nonzero_cv[[i]][j]] * pcatrain$rotation[,nonzero_cv[[i]][j]] %>% 
  matrix(nrow=40, ncol=40, byrow=F)
  }
  
}
print('Done')
pcacc
mean(pcacc) # 95.25
range(predlinkplr) # -13, 15 after alignment
hist(predlinkplr)

# Confusion matrices for PLR predictions
predicted_plr <- as.vector(pcpreds)
table(predicted_plr, nscv_log$accent)

# Confusion matrix comparing predictions from the two models
table(predicted_flm, predicted_plr) # 93.75% of predictions agree

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

save(sumcontrib, file='output/sumcontrib.RData')

# Combined model with MFCCs and formants ####
all_fm <- array(NA, dim=c(40, 400, 44))
all_fm[,,1:4] <- nscv_aligned_formant
all_fm[,,5:44] <- nscv_aligned_mfcc
# centre first MFCC
all_fm[,,5] <- apply(all_fm[,,5], 2, function(x) x-mean(x))

# Stack MFCC and formant curves from each observation into a single row
allstacked <- matrix(NA, nrow=400, ncol=40*44)
for(i in 1:44){
  allstacked[,(40*(i-1)+1):(40*i)] <- all_fm[,,i] %>% t()
}
all_pca <- prcomp(allstacked) # FPCA

# Model with PC scores with full dataset of formants and MFCCs
pccv_all <- cv.glmnet(all_pca$x, nscv_log$accent=='S', alpha=1, 
                      family=binomial())
pcamod_all <- glmnet(all_pca$x, nscv_log$accent=='S', alpha=1, family=binomial(),
                     lambda=pccv_all$lambda.min)
pcamod_all$beta[1:20] %>% abs() %>% plot() # first 20 coefficients
(all_pca$sdev^2/sum(all_pca$sdev^2)) %>% plot() # eigenvalues for the eigenfunctions
nonzero_all = which(pcamod_all$beta != 0)

# Cross validate combined model with MFCCs and formants
s <- split(1:400, nscv_log$id)
all_acc <- c()
predresp_all <- predlink_all <- pcpreds_all <- matrix(NA, nrow=100, ncol=4)
nonzero_cv_all <- betas_cv_all <- pcs_cv_all <- list()
print('Cross validating model with formants and MFCCs:')
for (i in 1:4){
  print(paste('Fold', i))
  testind <- s[[i]]
  
  train_fm = array(NA, dim=c(40, 300, 44))
  train_fm[,,1:4] <- nscv_cv_aligned[[i]]$train_aligned_formant
  train_fm[,,5:44] <- nscv_cv_aligned[[i]]$train_aligned_mfcc
  
  # Centre MFCC 1
  train_fm[,,5] <- apply(train_fm[,,5], 2, function(x) x-mean(x))
  
  # Stack training MFCC and formant curves into a single row
  trainstacked <- matrix(NA, nrow=300, ncol=40*44)
  for(j in 1:44){
    trainstacked[,(40*(j-1)+1):(40*j)] <- train_fm[,,j] %>% t()
  }
  
  test_fm = array(NA, dim=c(40, 100, 44))
  test_fm[,,1:4] <- nscv_cv_aligned[[i]]$test_aligned_formant
  test_fm[,,5:44] <- nscv_cv_aligned[[i]]$test_aligned_mfcc
  # Centre MFCC 1
  test_fm[,,5] <- apply(test_fm[,,5], 2, function(x) x-mean(x))
  
  # Stack testing MFCC curves into a single row
  teststacked <- matrix(NA, nrow=100, ncol=40*44)
  for(j in 1:44){
    teststacked[,(40*(j-1)+1):(40*j)] <- test_fm[,,j] %>% t()
  }
  
  train <- trainstacked
  test <- teststacked
  
  pcatrain <- prcomp(train)
  
  pccv <- cv.glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                    alpha=1, family=binomial())
  pcamod <- glmnet(pcatrain$x, nscv_log$accent[-testind]=='S', 
                   alpha=1, family=binomial(), lambda=pccv$lambda.min)
  
  testpcscores <- predict(pcatrain, test)
  predresp_all[,i] <- predict(pcamod, testpcscores, type='response')
  predlink_all[,i] <- predict(pcamod, testpcscores, type='link')
  pcpreds_all[,i] <- ifelse(predresp_all[,i] > 0.5, 1, 0)
  all_acc[i] <- sum(pcpreds_all[,i]==(nscv_log$accent[testind]=='S'))/
    length(testind)
  
  # Save chosen FPCs, and weighted loadings from each fold
  nonzero_cv_all[[i]] <- which(pcamod$beta != 0)
  pcs_cv_all[[i]] <- array(NA, dim=c(40,44,length(nonzero_cv_all[[i]])))
  for (j in 1:length(nonzero_cv_all[[i]])){
    pcs_cv_all[[i]][,,j] <- pcamod$beta[nonzero_cv_all[[i]][j]] * 
      pcatrain$rotation[,nonzero_cv_all[[i]][j]] %>% 
      matrix(nrow=40, ncol=44, byrow=F)
  }
  
}
print('Done')
all_acc
mean(all_acc) # 92.75
range(predlink_all) # -21, 12
hist(predlink_all)

# Confusion matrix for combined model predictions
predicted_com = as.vector(pcpreds_all)
table(predicted_com, nscv_log$accent)

# Look at weighted sum of FPCs in each fold
sum_pcs_cv_all <- lapply(pcs_cv_all, function(x) {apply(x, 1:2, sum)})
par(mfrow=c(2,2))
for (i in 1:length(pcs_cv_all)){
  image(sum_pcs_cv_all[[i]])
}


# Geographical variation model ####
# set up BNC data 
nspeaker <- bnc_log %>% 
  group_by(id) %>% 
  summarise(n=n(),
            soundnums=list(soundnum)) %>% 
  arrange(desc(n))

# Formants for BNC
bnc_f2 <- list(f2f=t(bnc_aligned_formant[,,2]),
               f2mid = bnc_aligned_formant[20,,2],
               # id=rep('p5', times=nrow(bnc_log)),
               id = bnc_log$id,
               time=matrix(seq(0, 1, length.out=40), ncol=40, nrow=nrow(bnc_log), byrow=T))

# Predictions for each sound using functional formant model
bncfpreds <- predict(formant_model, bnc_f2, type='response') # FLM predictions
bncflink <- predict(formant_model, bnc_f2, type='link')
bncfterms <- predict(formant_model, bnc_f2, type='terms')
range(as.numeric(bncflink)) # -238, 218


# MFCCs for BNC
bncmfccstacked <- matrix(NA, nrow=nrow(bnc_log), ncol=40*40)
# Centre MFCC 1
bncmfccstacked[,1:40] <- apply(bnc_aligned_mfcc[,,1], 2, function(x) x-mean(x)) %>% 
  t()
for(i in 2:40){
  bncmfccstacked[,(40*(i-1)+1):(40*i)] <- bnc_aligned_mfcc[,,i] %>% t()
}

# Project on FPCs to get the scores
bncpcs <- predict(mfccpm, bncmfccstacked)

# Predictions for each sound
bncppreds <- predict(pcamodfull, bncpcs, type='response') # PLR predicted probs
bncplink <- predict(pcamodfull, bncpcs, type='link')

# checking the range and distribution of linear predictors and probabilities
# between the NSCV cross validated predictions and BNC 
hist(bncplink)
hist(predlinkplr) #cross validated NSCV linear predictors
hist(bncppreds)
range(bncplink) #-13, 14 

table(bncppreds>0.5, bncfpreds>0.5) # agreement between FLM and PLR model predictions
# around 73% agreement


# BNC predictions with combined model
bnc_all_fm <- array(NA, dim=c(40, nrow(bnc_log), 44))
bnc_all_fm[,,1:4] <- bnc_aligned_formant
bnc_all_fm[,,5:44] <- bnc_aligned_mfcc

# centre first MFCC
bnc_all_fm[,,5] <- apply(bnc_all_fm[,,5], 2, function(x) x-mean(x))

bnc_allstacked <-  matrix(NA, nrow=nrow(bnc_log), ncol=40*44)
for(i in 1:44){
  bnc_allstacked[,(40*(i-1)+1):(40*i)] <- bnc_all_fm[,,i] %>% t()
}

bnc_pcscores = predict(all_pca, bnc_allstacked)
bncallpreds <- predict(pcamod_all, bnc_pcscores, type='response') # combined model predicted probs
bncalllink <- predict(pcamod_all, bnc_pcscores, type='link')


# Compute predictions for each speaker using average sound
nspeakerpreds <- bnc_log %>% 
  mutate(ppred = as.numeric(bncppreds),
         fpred = as.numeric(bncfpreds),
         plink = as.numeric(bncplink),
         flink = as.numeric(bncflink),
         combined_pred = as.numeric(bncallpreds),
         combined_link = as.numeric(bncalllink)) %>% 
  group_by(id) %>% 
  summarise(n=n(),
            soundnums=list(soundnum),
            avgppred = plogis(mean(plink)),
            avgplink = mean(plink),
            avgfpred = plogis(mean(flink)),
            avgflink = mean(flink),
            avg_combinedpred = plogis(mean(combined_link)),
            avg_link = mean(combined_link)) %>% 
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

# Combined model predictions
location_combined <- gam(avg_combinedpred ~ s(lon, lat, k=40, bs='so', xt=list(bnd=list(bn))),
                         knots=knotsfew, data=bnclocpred[lf.inside==1,], 
                         family = betar())
summary(location_combined)
