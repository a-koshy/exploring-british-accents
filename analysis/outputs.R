# Figures in the paper #

source('analysis/models.R')

# Mel scale plot
pdf('output/figs/melhz.pdf', width=4.5, height=3.5)
par(mgp=c(1.5, 0.5, 0))
x <- 0:20000
plot(x, 2595*log10(1+(x/700)), type='l', bty='n',
     xlab='Hz', ylab='Mel', main='Mel versus Hz scale')
dev.off()

## Functional linear model with formants ####

#### NSCV EDA for F1-F4 ####
# 1: plot of each formant average against each other
# 2: plot each formant curves in N and S accents
formantavgs = apply(nscv_smooth_formant$formant_norm, c(2, 3), mean) %>% 
    as.data.frame()
colnames(formantavgs) = c('F1', 'F2', 'F3', 'F4')
formantavgs$accent = factor(nscv_log$accent)

pdf('output/figs/formant-avg-eda.pdf', width=5, height=5)
pairs(formantavgs[,1:4], 
      col=cols[formantavgs$accent],
      pch=20)
dev.off()

# For the table
formantavgs %>% 
    group_by(accent) %>% 
    summarise(avgf1 = mean(F1),
              avgf2 = mean(F2),
              avgf3 = mean(F3),
              avgf4 = mean(F4))

pdf('output/figs/formant-curve-eda.pdf', width=10, height=7)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 2.5, 1,0.5), mfrow=c(2,2))
for (i in 1:4){
    matplot(seq(0, 1, length.out=40), nscv_smooth_formant$formant_norm[,,i],
            type='l', lty=1,
            col=cols[as.factor(nscv_log$accent)], bty='n',
            xlab='Rescaled time', ylab=paste('Formant', i))
}
dev.off()

# Smoothed and aligned F2 curves
f2 <- nscv_aligned_formant[,,2]
pdf('output/figs/f2.pdf', width=5, height=3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 2.5, 1,0.5))
matplot(seq(0, 1, length.out=40), f2, type='l', lty=1,
        col=cols[as.factor(nscv_log$accent)], bty='n',
        xlab='Rescaled time', ylab=expression('F'[2]))
legend(0.1, 2500, legend=c('South', 'North'), col=opcols[2:1],
       lty=1, bty='n')
dev.off()

# GIF of NSCV F1 and F2 trajectories
print('Creating NSCV formant GIF:')
nscvgif <- makegif(100, 4, nscv_smooth_formant$formant_norm, nscv_log)
anim_save(filename='output/supplement/nscv.gif', nscvgif)


# Beta_1 curve
d <- plot(formant_model)
lims <- range(d[[1]]$fit+d[[1]]$se, d[[1]]$fit-d[[1]]$se)
pdf('output/figs/betahat.pdf', width = 5, height = 3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(seq(0, 1, length.out = 100), d[[1]]$fit*40, type='l',
     xlab='Rescaled time',ylab=expression(hat(beta)[1](t)), bty='n',
     ylim=lims*40, cex=0.5)
lines(seq(0, 1, length.out = 100), (d[[1]]$fit+d[[1]]$se)*40, lty=2)
lines(seq(0, 1, length.out = 100), (d[[1]]$fit-d[[1]]$se)*40, lty=2)
abline(h=0, lty=3, col='slategrey')
title(expression(hat(beta)[1](t)), line=-1)
dev.off()

# ROC curve
table(nscv_log$accent=='S', formant_model_cv$cv_preds)
roc <- simple_roc(nscv_log$accent=='S', scores=as.vector(formant_model_cv$cv_preds))
pdf('output/figs/flm-roc.pdf', width=4.5, height=3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(roc$FPR, roc$TPR, type='l', bty='n', xlab='False positive rate',
     ylab='True positive rate', main='ROC curve')
grid <- seq(0, 1, length.out = 40)
lines(grid, grid, lty=3)
points(8/200, 195/200, col='red', pch=20) #use confusion matrix to calculate FPR and TPR
text(.35, .97, 'Threshold = 0.5', pos=1, col='red')
dev.off()


# PLR model ####

# Scree plot from FPCA
pdf('output/figs/mfcc-screeplot.pdf', width=6, height=4)
plot(mfccpm$sdev[1:25]^2, type='b', main='Scree plot of first 25 eigenvalues',
     xlab='Component number', ylab='Eigenvalue', bty='n')
dev.off()


#### Plot FPCs chosen in each CV model ####
pdf('output/figs/cv-mfcc-images.pdf', height=8, width=12)
par(oma=c(3, 1, 1, 4), mfrow=c(2,2))
sum_pcs_cv <- lapply(pcs_cv, function(x) {apply(x, 1:2, sum)})
for (i in 1:length(pcs_cv)){
    image.plot(y=1:40, z=sum_pcs_cv[[i]],
               xlab='Rescaled time', ylab=ifelse(i==1 | i==3, 'MFCCs', ''), 
               main=paste0('CV Fold ', i, ': Accuracy ', pcacc[i]),
               col=hcl.colors(20, 'RdYlGn', rev=TRUE),
               zlim=range(sum_pcs_cv))
}
dev.off()

# MFCC accent perturbation matrix as image
pdf('output/figs/perturb.pdf', height = 4, width=4.5)
image.plot(y=1:40, z=sumcontrib, xlab='Rescaled time', ylab='MFCCs',
           main='Weighted sum of FPCs',
           col = hcl.colors(20, 'RdYlGn', rev=TRUE))
dev.off()

# First 9 MFCC components overlaid the corresponding components from 4 CV folds
pdf('output/figs/perturb-first9-cv.pdf', height=4, width=8)
par(mar=c(2,4,4,2))
sum_pcs_cv <- lapply(pcs_cv, function(x) {apply(x, 1:2, sum)})
h <- 0.3
plot(sumcontrib[,1], lty=1, type='n', col=1, xlim=c(1,360),
     bty='n', ylim=range(sum_pcs_cv)+c(0,0.03), xaxt='n',
     main='First 9 MFCC components', xlab='', ylab='')
text(20, h, 'MFCC\n1', col=1, cex=0.8)

for (j in 1:4){
    for(i in 1:9){
        x = (40*(i-1) +1)
        lines(x:(x+39), sum_pcs_cv[[j]][, i], col='grey')
    }
}

for(i in 1:9){
    x <- (40*(i-1) + 1)
    lines(x:(x+39), sumcontrib[, i],
          col=1, lwd=1.2)
    text(40*i-20, h, paste0('MFCC\n', i), col=, cex=0.8)
}

abline(h=0, lty=3)
abline(v=seq(41, 321, by=40), lty=2, col='lightgrey')
dev.off()

# ROC curve for PLR model
rocplr <- simple_roc(nscv_log$accent=='S', scores=as.vector(predrespplr))
pdf('output/figs/plr-roc.pdf', width=4.5, height=3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(rocplr$FPR, rocplr$TPR, type='l', bty='n', xlab='False positive rate',
     ylab='True positive rate', main='ROC curve for the MFCC model')
grid <- seq(0, 1, length.out = 40)
lines(grid, grid, lty=3)
points(11/200, 192/200, col='red', pch=20) #use confusion matrix to calculate FPR and TPR at threshold 0.5
text(.22, .95, 'Threshold = 0.5', pos=1, col='red')
dev.off()

# plot of beta values
pdf('output/figs/plrmodel-25coefs.pdf', width=4.5, height=3)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(pcamodfull$beta[1:25], type='h',
     bty='n', xlab='Coefficient',
     ylab='', col=(pcamodfull$beta[1:25]!=0)+1,
     main='Coefficients of the first 25 FPCs')
points(pcamodfull$beta[1:25], col=(pcamodfull$beta[1:25]!=0)+1, pch=19)
abline(h=0, lty=3)
dev.off()


# Perturb vowels ####
print('Perturbing vowels')
perturbsound(infile = 'sounds/class-N.wav', 
             outfile = 'output/supplement/class-NtoS.wav',
             starttime = 0.101085, endtime = 0.336657,
             NtoS = TRUE)

perturbsound(infile = 'sounds/blast-S.wav', 
             outfile = 'output/supplement/blast-StoN.wav',
             starttime = 0.150204, endtime = 0.407075,
             NtoS = FALSE)


# Geographical variation model ####
# Bubble plot of sounds at locations
avgd <- bnc_log %>% 
    group_by(placenamecleaned) %>%
    summarise(duration=mean(duration),
              nobs=n()) %>% 
    inner_join(locations, by=c('placenamecleaned'='location.names'))
datainside <- point.in.polygon(avgd$lon, avgd$lat,
                               bn$lon, bn$lat)

pdf('output/figs/obsnum.pdf', width=5, height=6)
par(mar=c(0,0,0,0), mfrow=c(1,3), mgp=c(1.5,0.5,0))
ggplot()+
    geom_polygon(aes(x=lon, y=lat), fill=NA, color='grey50', data=bn)+
    coord_fixed(1.3)+
    geom_point(aes(x=lon, y=lat, size=nobs), alpha=0.3, col=cols[1],
               data=avgd[datainside==1,])+
    theme_void()+
    xlab(label='')+ylab(label='')+
    labs(size='Sounds', title='Number of sounds by location')+
    scale_size(range = c(1, 10), name="Sounds", breaks=c(10,50,100,200))
dev.off()

# Bubble plot of speakers at locations
avgsp <- bnc_log %>% 
    group_by(placenamecleaned, id) %>%
    summarise(duration=mean(duration),
              n_each=n()) %>% 
    group_by(placenamecleaned) %>% 
    summarise(nid =n()) %>% 
    inner_join(locations, by=c('placenamecleaned'='location.names'))
avgsp.inside <- point.in.polygon(avgsp$lon, avgsp$lat,
                                 bn$lon, bn$lat)

pdf('output/figs/idnum.pdf', width=5, height=6)
par(mar=c(0,0,0,0), mfrow=c(1,3), mgp=c(1.5,0.5,0))
ggplot()+
    geom_polygon(aes(x=lon, y=lat), fill=NA, color='grey50', data=bn)+
    coord_fixed(1.3)+
    geom_point(aes(x=lon, y=lat, size=nid), alpha=0.3, col=cols[1],
               data=avgsp[datainside==1,])+
    theme_void()+
    xlab(label='')+ylab(label='')+
    labs(size='Speakers', title='Number of speakers by location')+
    scale_size(range = c(1, 10), name="Speakers", breaks=c(1, 5, 10))
dev.off()

# Boxplot of BNC predicted probs for 50 speakers
bncprobs <- bnc_log %>%
    mutate(fpred = as.numeric(bncfpreds),
           ppred = as.numeric(bncppreds),
           flink = as.numeric(bncflink),
           plink = as.numeric(bncplink))

hist(bncprobs$ppred)
hist(bncprobs$fpred)

bncboxplot <- bncprobs %>% 
inner_join(bncprobs %>% 
               group_by(id) %>% 
               summarise(nsound=n()) %>% 
               slice_max(order_by = nsound, n=50), by='id') %>% 
    select(id, placenamecleaned, fpred, ppred, nsound)
head(bncboxplot)
range(bncboxplot$nsound)

pdf('output/figs/predsboxplot.pdf', width=10, height=4)
boxplot(ppred~factor(id), data=bncboxplot, main='MFCC model predictions by id', 
        las=2, xlab='')
boxplot(fpred~factor(id), data=bncboxplot, main='Formant model predictions by id',
        las=2, xlab='')
dev.off()

# Soap film maps
plotsoapfilm(locationplr, name='', title='PLR')
plotsoapfilm(locationflm, name='', title='FLR')
plotsoapfilm(location_combined, name='', title='Combined Model')

# Table of social class for appendix
bnc_log %>%
    group_by(soc, id) %>% 
    summarise(num = n()) %>% 
    group_by(soc) %>% 
    summarise(numsoc=n())

# Histogram of sound durations
pdf('output/figs/durationhist.pdf', width=10, height=4)
par(mfrow=c(1,2))
hist(bnc_log$duration, xlab='Duration (s)', main='BNC sound lengths',
     xlim=c(0.2, 1), ylab='Observations')
hist(nscv_log$duration, xlim=c(0.2, 1), xlab='Duration (s)', 
     main='NSCV sound lengths', ylab='Observations')
dev.off()

# Histogram of sounds per speaker
nsound <- bnc_log %>% 
    group_by(id) %>% 
    summarise(n=n())
pdf('output/figs/sounds-per-speaker.pdf', width=5, height=4)
hist(nsound$n, xlab='Number of sounds per speaker', ylab='Speakers', 
     main='')
dev.off()

# Histogram of number of locations per speaker
bnclocs <- bnc_log %>% group_by(id, placenamecleaned) %>% summarise(n=n()) %>% 
    group_by(id) %>% summarise(nlocs=n())
hist(bnclocs$nlocs)
sum(bnclocs$nlocs==1)/nrow(bnclocs) # 88% have only one location

# How big a geographical range is spanned by multiple locations of each speaker?
# Bound by the dist between opposite vertices of the rectangle containing locations
bncdist <- bnc_log %>% 
    inner_join(locations, by=c('placenamecleaned'='location.names')) %>%
    group_by(id) %>% 
    summarise(maxlat = max(lat),
              minlat = min(lat),
              maxlon = max(lon),
              minlon = min(lon),
              maxdist = earth.dist(maxlon, maxlat, minlon, minlat)) %>% 
    arrange(desc(maxdist))

sum(bncdist$maxdist<=10)/nrow(bncdist) # 94% within 10km

# Comparing models and where they disagree ####
# combined model with formant model: 96.25%
table(predicted_com, predicted_flm)
# combined with mfcc model: 89%
table(predicted_com, predicted_plr)
# formant with mfcc model: 93.75%
table(predicted_flm, predicted_plr)

# comparing formant and mfcc models
disagree = which(predicted_plr != ifelse(predicted_flm=='N', 0, 1)) 
table(nscv_log$id[disagree]) # observations mostly from p2 and p3
table(nscv_log$accent[disagree]) # 15 N and 10 S
# compare the predictions with true accent
predicted_plr[disagree]
predicted_flm[disagree]
nscv_log$accent[disagree]

# get a sample of 50 observations where models agree
agree = sample(c(1:400)[-disagree], size=50)

pdf('output/figs/compare-disagreeing-formants-nscv.pdf', width=9, height=4)
par(mfrow=c(1,2))
tseq = seq(0, 1, length.out=40)
matplot(tseq, nscv_aligned_formant[,agree,2], type='l', lty=1,
        col=cols[factor(nscv_log$accent[agree])],
        ylim=c(700, 2000), bty='n', xlab='Aligned time', ylab='',
        main='NSCV F2 curves where models agree')
matplot(tseq, nscv_aligned_formant[,disagree,2], type='l', lty=1,
        col=cols[factor(nscv_log$accent[disagree])],
        ylim=c(700, 2000), bty='n', xlab='Aligned time', ylab='',
        main='NSCV F2 curves where models disagree')
dev.off()

pdf('output/figs/compare-disagreeing-formants-bnc.pdf', width=9, height=4)
par(mfrow=c(1,2))
bncdisagree = which((as.numeric(bncfpreds) < 0.5 & as.numeric(bncppreds) > 0.5) | (as.numeric(bncfpreds) >0.5 & as.numeric(bncppreds) < 0.5))

bncdis = sample(bncdisagree, size=50)
bncagree = sample(1:nrow(bnc_log)[-bncdisagree], size=50)
tseq = seq(0, 1, length.out=40)
matplot(tseq, bnc_aligned_formant[,bncagree,2], type='l', lty=1,
        col=cols[factor(bncppreds[bncagree]>0.5)],
        ylim=c(700, 2400), bty='n', xlab='Aligned time', ylab='',
        main='BNC F2 curves where models agree')
matplot(tseq, bnc_aligned_formant[,bncdis,2], type='l', lty=1, 
        col=cols[factor(bncppreds[bncdis]>0.5)],
        ylim=c(700, 2400), bty='n', xlab='Aligned time', ylab='',
        main='BNC F2 where models disagree')
dev.off()

pdf('output/figs/compare-link-predictions.pdf', width=8, height=4)
par(mfrow=c(1,2))
plot(bncplink, bncflink, pch=c(20,3)[factor(sign(as.numeric(bncplink))==sign(as.numeric(bncflink)))],
     col=cols[4:3][factor(sign(as.numeric(bncplink))==sign(as.numeric(bncflink)))],
     bty='n', main='Link predictions on BNC vowels', xlab='MFCC model',
     ylab='Formant model',
     ylim=c(-300, 400), xlim=c(-15, 15))

nscvlinkflm = as.numeric(formant_model_cv$cv_link)
nscvlinkplr = as.numeric(predlinkplr)
plot(nscvlinkplr, nscvlinkflm, pch=c(20,3)[factor(sign((nscvlinkflm))==sign((nscvlinkplr)))],
     col=cols[4:3][factor(sign(nscvlinkflm)==sign(nscvlinkplr))],
     bty='n', main='Link predictions on NSCV vowels', xlab='MFCC model', 
     ylab='Formant model',
     ylim=c(-300, 400), xlim=c(-15, 15))
dev.off()

# Perturbing vowels with MFCC model and predicting with formant model
set.seed(123)
testsounds = nscv_log %>% 
  group_by(accent) %>% 
  slice_sample(n=2) %>% 
  mutate(file= sapply(strsplit(filename, "sound/"), "[", 2))

compare_perturbed = list()
for (i in 1:4){
  compare_perturbed[[i]] = perturbsound_comparison(infile=testsounds$file[i],
                                                   starttime=0,
                                                   endtime=testsounds$duration[i],
                                                   f_model=formant_model,
                                                   sumcontrib = sumcontrib)
}


pdf('output/figs/comparing-perturbed-vowels-link.pdf')
r = range(lapply(compare_perturbed, function(x) range(x$formant_preds)))
plot(-12:12, compare_perturbed[[1]]$formant_preds, type='b', ylim=r,
     col=opcols[2], bty='n',
     xlab='Linear predictor from MFCC model', ylab='Linear predictor from formant model',
     main='Comparing predictions for perturbed vowels')
lines(-12:12, compare_perturbed[[2]]$formant_preds, type='b', col=opcols[2])
lines(-12:12, compare_perturbed[[3]]$formant_preds, type='b', col=opcols[1], pch=2)
lines(-12:12, compare_perturbed[[4]]$formant_preds, type='b', col=opcols[1], pch=2)
abline(h=0, v=0, lty=2)
dev.off()

pdf('output/figs/comparing-perturbed-vowels-prob.pdf')
plot(plogis(-12:12), plogis(compare_perturbed[[1]]$formant_preds), type='b', ylim=c(0,1),
     col=cols[2], bty='n',
     xlab='Predicted probability from MFCC model', ylab='Predicted probability from formant model',
     main='Comparing predictions for perturbed vowels')
lines(plogis(-12:12), plogis(compare_perturbed[[2]]$formant_preds), type='b', col=cols[2])
lines(plogis(-12:12), plogis(compare_perturbed[[3]]$formant_preds), type='b', col=cols[1], pch=2)
lines(plogis(-12:12), plogis(compare_perturbed[[4]]$formant_preds), type='b', col=cols[1], pch=2)
dev.off()
