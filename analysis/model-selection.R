# Model selection for formants ####
source('analysis/setup.R')
source('analysis/functions.R')
source('analysis/models.R')

shorter = c(10, 20, 30) # time indices to use for 3 point formant models

# Midpoint formant models ####

# midpoint of f1
midf1re = formula((accent=='S') ~ f1mid+s(id, bs='re'))
af1re = check_model(midf1re, align=T)
bf1re = check_model(midf1re, align=F)
plot(af1re$model)

# midpoint of f2
midf2re = formula((accent=='S') ~ f2mid+s(id, bs='re'))
af2re = check_model(midf2re, align=T)
bf2re = check_model(midf2re, align=F)
plot(af2re$model)

# midpoint of f1 and f2
midf1f2re = formula((accent=='S') ~ f2mid + f1mid + s(id, bs='re'))
af1f2re = check_model(midf1f2re, align=T)
bf1f2re = check_model(midf1f2re, align=F)
plot(af1f2re$model)

# 3 point models ####

# F2 at 3 time points
shorterf2re = formula((accent=='S') ~ s(timeshorter, by=f2shorter, k=length(shorter), bs='cr')+s(id, bs='re'))
af2re2 = check_model(formula=shorterf2re, align=T)
bf2re2 = check_model(formula=shorterf2re, align=F)
plot(af2re2$model)

# F1 and F2 - 3 time points
shorterf1f2re = formula((accent=='S') ~ s(timeshorter, by=f2shorter, k=length(shorter), bs='cr') +
                          s(timeshorter, by=f1shorter, k=length(shorter), bs='cr') + s(id, bs='re'))
af1f2re2 = check_model(shorterf1f2re, align=T, optimizer=c('outer', 'bfgs'))
bf1f2re2 = check_model(shorterf1f2re, align=F, optimizer=c('outer', 'bfgs'))
plot(af1f2re2$model)

# mean(f1)+gap at 3 points
shortfmlre = formula((accent=='S') ~ f1+s(timeshorter, by=gapshorter, k=length(shorter), bs='cr')+s(id, bs='re'))
a5re = check_model(shortfmlre, align=T)
b5re = check_model(shortfmlre, align=F)

plot(a5$model, ylim=c(-0.03, 0))
abline(h=0)
plot(a5re$model)


# Functional models ####

# only f1 functional
fmlf1re = formula((accent=='S') ~ s(time, by=f1f, k=15, bs='cr')+s(id, bs='re'))
af1re = check_model(fmlf1re, align=T)
bf1re = check_model(fmlf1re, align=F)

gam.check(af1re$model)
plot(af1re$model)

# only f2 functional
fml8re = formula((accent=='S') ~ s(time, by=f2f, k=15, bs='cr')+s(id, bs='re'))
a8re = check_model(fml8re, align=T, optimizer=c('outer', 'bfgs'))
b8re = check_model(fml8re, align=F, optimizer=c('outer', 'bfgs'))

gam.check(a8re$model)
plot(a8re$model)
abline(h=0)

# only gap functional term
fml2re = formula((accent=='S') ~ s(time, by=gap, bs='cr', k=15)+s(id, bs='re'))
a2re = check_model(fml2re, align=T)
b2re = check_model(fml2re, align=F)

gam.check(a2re$model)
plot(a2re$model)
abline(h=0)


# mean(f1) and functional gap covariate
fml4re = formula((accent=='S') ~ f1+s(time, by=gap, bs='cr', k=15)+s(id, bs='re'))
a4re = check_model(fml4re, align=T, optimizer=c('outer', 'bfgs'))
b4re = check_model(fml4re, align=F, optimizer=c('outer', 'bfgs'))

summary(a4re$model)
gam.check(a4re$model)
plot(a4re$model)

# two functional terms for f1 and f2
fml7re = formula((accent=='S') ~ s(time, by=f1f, k=15, bs='cr')+
                   s(time, by=f2f, k=15, bs='cr')+s(id, bs='re'))
a7re = check_model(fml7re, align=T, optimizer=c('outer', 'bfgs'))
b7re = check_model(fml7re, align=F, optimizer=c('outer', 'bfgs'))

gam.check(a7re$model)

# Lasso random effects model ####
# glmmlasso

# Model with PC scores with full dataset
lassoredf = mfccpm$x %>% as.data.frame()
lassoredf$id =factor(nscv_log$id)
lassoredf$accent = (nscv_log$accent=='S')

# specify starting values for the very first fit
lambda <- seq(150, 0, by=-5)
Delta.start<-as.matrix(t(rep(0,4+401)))
Q.start<-0.1  
n = names(lassoredf)
f <- as.formula(paste("accent ~ 1 + ", paste(n[!n %in% c("accent", 'id')], collapse = " + ")))
AIC_vec <- c()
for(j in 1:length(lambda)){
  print(paste("Iteration ", j,sep=""))
  
  gtest <- glmmLasso(f, data = lassoredf, 
                     rnd = list(id=~1),  
                     family = binomial(), 
                     lambda=lambda[j], switch.NR=F,final.re=TRUE,
                     control = list(start=Delta.start[j,],q_start=Q.start[j]))  
  
  print(colnames(gtest$Deltamatrix)[2:7][gtest$Deltamatrix[gtest$conv.step,2:7]!=0])
  AIC_vec[j]<-gtest$aic
  Delta.start<-rbind(Delta.start,gtest$Deltamatrix[gtest$conv.step,])
  Q.start<-c(Q.start,gtest$Q_long[[gtest$conv.step+1]])
}

# Final mixed effect LASSO model trained on full dataset with optimal lambda
opt = which.min(AIC_vec)
nscvfinal <- glmmLasso(f, data = lassoredf, 
                       rnd = list(id=~1),  
                       family = binomial(), 
                       lambda=lambda[opt],
                       switch.NR=F,final.re=TRUE,
                       control = list(start=Delta.start[opt,],q_start=Q.start[opt]))  

summary(nscvfinal)

# cross validate
s <- split(1:400, nscv_log$id)
lre_acc <- c()
lambda <- seq(100, 0, by=-5)

predresp_lre <- predlink_lre <- pcpreds_lre <- matrix(NA, nrow=100, ncol=4)
nonzero_cv_lre <- pcs_cv_lre <-  betas_cv_lre <- list()
AIC_cv_lre = matrix(0, nrow=4, ncol=length(lambda))

print('Cross validating LASSO mixed effects MFCC model:')
for (i in 1:4){
  print(paste('Fold', i))
  testind <- s[[i]]
  
  train_fm = array(NA, dim=c(40, 300, 40))
  train_fm <- nscv_cv_aligned[[i]]$train_aligned_mfcc
  
  # Centre MFCC 1
  train_fm[,,1] <- apply(train_fm[,,1], 2, function(x) x-mean(x))
  
  # Stack training MFCC and formant curves into a single row
  trainstacked <- matrix(NA, nrow=300, ncol=40*40)
  for(j in 1:40){
    trainstacked[,(40*(j-1)+1):(40*j)] <- train_fm[,,j] %>% t()
  }
  
  test_fm <- nscv_cv_aligned[[i]]$test_aligned_mfcc
  # Centre MFCC 1
  test_fm[,,1] <- apply(test_fm[,,1], 2, function(x) x-mean(x))
  
  # Stack testing MFCC curves into a single row
  teststacked <- matrix(NA, nrow=100, ncol=40*40)
  for(j in 1:40){
    teststacked[,(40*(j-1)+1):(40*j)] <- test_fm[,,j] %>% t()
  }
  
  train <- trainstacked
  test <- teststacked
  pcatrain <- prcomp(train)
  
  traindf = pcatrain$x %>% as.data.frame()
  traindf$id =factor(nscv_log$id[-testind])
  traindf$accent = (nscv_log$accent=='S')[-testind] %>% as.integer()
  
  Delta.start <- as.matrix(t(rep(0,3+301)))
  Q.start<-0.1  
  n = names(traindf)
  f <- as.formula(paste("accent ~ 1 + ", paste(n[!n %in% c("accent", 'id')], collapse = " + ")))
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    
    mod_cv <- glmmLasso(f, data = traindf, 
                        rnd = list(id=~1),  
                        family = binomial(), 
                        lambda=lambda[j], switch.NR=F, final.re=TRUE,
                        control = list(start=Delta.start[j,], q_start=Q.start[j]))  
    
    AIC_cv_lre[i, j] <- mod_cv$aic
    Delta.start = rbind(Delta.start, mod_cv$Deltamatrix[mod_cv$conv.step,])
    Q.start <- c(Q.start, mod_cv$Q_long[[mod_cv$conv.step+1]])
  }
  
  optimal_lambda_cv = which.min(AIC_cv_lre[i,])
  cv_final_model <- glmmLasso(f, data = traindf, 
                              rnd = list(id=~1),  
                              family = binomial(), 
                              lambda=lambda[optimal_lambda_cv],
                              switch.NR=F,final.re=TRUE,
                              control = list(start=Delta.start[optimal_lambda_cv,], q_start=Q.start[optimal_lambda_cv]))  
  
  testpcscores <- predict(pcatrain, test) %>% as.data.frame()
  
  # need some dummy values for id and accent in test data 
  testpcscores$id = factor(unique(traindf$id)[3], levels=levels(factor(traindf$id)))
  testpcscores$accent = F
  
  predresp_lre[,i] <- predict(cv_final_model, testpcscores, type='response')
  predlink_lre[,i] <- predict(cv_final_model, testpcscores, type='link')
  pcpreds_lre[,i] <- ifelse(predresp_lre[,i] > 0.5, 1, 0)
  lre_acc[i] <- sum(pcpreds_lre[,i]==(nscv_log$accent[testind]=='S'))/length(testind)
  
  # Save chosen FPCs, and weighted loadings from each fold
  nonzero_cv_lre[[i]] <- which(cv_final_model$coefficients != 0)
  pcs_cv_lre[[i]] <- array(NA, dim=c(40, 40, length(nonzero_cv_lre[[i]])))
  for (j in 1:length(nonzero_cv_lre[[i]])){
    pcs_cv_lre[[i]][,,j] <- cv_final_model$coefficients[nonzero_cv_lre[[i]][j]] * 
      pcatrain$rotation[,nonzero_cv_lre[[i]][j]] %>% 
      matrix(nrow=40, ncol=40, byrow=F)
  }
}
print('Done')
lre_acc
mean(lre_acc) # 90.25
range(qlogis(predlink_lre), na.rm=T) # -36, 36
hist(qlogis(predlink_lre))
sum(sign(qlogis(predlink_lre))*sign(predlinkplr))/400

# compare mfccs in each fold
par(mfrow=c(2,2))
for (i in 1:4) {
  sumcontrib_i = apply(pcs_cv_lre[[i]], 1:2, sum)
  image(sumcontrib_i, main=paste('Fold', i, ':', (length(nonzero_cv_lre[[i]])-1), 'FPC chosen'))
}


# resynthesise sounds
nonzero_nscvfinal = which(nscvfinal$coefficients != 0)

betas_lre <- matrix(data=nscvfinal$coefficients[nonzero_nscvfinal][-1], 
                    nrow=length(nonzero_nscvfinal)-1, ncol=400,
                    byrow = F)

contrib_nscvfinal <- array(NA, dim=c(40,40,length(nonzero_nscvfinal)-1))
for (i in 2:length(nonzero_nscvfinal)){
  contrib_nscvfinal[,,i-1] <- nscvfinal$coefficients[nonzero_nscvfinal[i]] * 
    mfccpm$rotation[,nonzero_nscvfinal[i]] %>% 
    matrix(nrow=40, ncol=40, byrow=F)
}
sumcontrib_nscvfinal <- apply(contrib_nscvfinal, 1:2, sum)
image(sumcontrib_nscvfinal)

# resynthesise sounds - rename the matrix used in resynthesis
# sumcontrib_original=sumcontrib # MFCC matrix from lasso model in the paper
# sumcontrib = sumcontrib_nscvfinal
# sumcontrib = sumcontrib_original
perturbsound(infile = 'sounds/blast-S.wav',
             outfile = 'output/lasso-re-blast-StoN.wav',
             starttime = 0.150204, endtime = 0.407075,
             NtoS = F)
