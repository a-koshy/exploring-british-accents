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

# NSCV gap between F2 and F1
f1f2gap <- (nscv_smooth_formant$formant[,,2]-
                nscv_smooth_formant$formant[,,1])
pdf('output/figs/f1f2gap.pdf', width=5, height=3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 2.5, 1,0.5))
matplot(seq(0, 1, length.out=40), f1f2gap, type='l', lty=1,
        col=cols[as.factor(nscv_log$accent)], bty='n',
        xlab='Rescaled time', ylab=expression('F'[2]-'F'[1]))
legend(0.1, 2000, legend=c('South', 'North'), col=opcols[2:1],
       lty=1, bty='n')
dev.off()

# GIF of NSCV F1 and F2 trajectories
print('Creating NSCV formant GIF:')
nscvgif <- makegif(100, 4, nscv_smooth_formant$formant_norm, nscv_log)
anim_save(filename='output/supplement/nscv.gif', nscvgif)


# Beta_2 curve
d <- plot(gapmodel)
lims <- range(d[[1]]$fit+d[[1]]$se, d[[1]]$fit-d[[1]]$se)
pdf('output/figs/betahat.pdf', width = 5, height = 3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(seq(0, 1, length.out = 100), d[[1]]$fit*40, type='l',
     xlab='Rescaled time',ylab=expression(hat(beta)[2](t)), bty='n',
     ylim=lims*40, cex=0.5)
lines(seq(0, 1, length.out = 100), (d[[1]]$fit+d[[1]]$se)*40, lty=2)
lines(seq(0, 1, length.out = 100), (d[[1]]$fit-d[[1]]$se)*40, lty=2)
abline(h=0, lty=3, col='slategrey')
title(expression(hat(beta)[2](t)), line=-1)
dev.off()

# ROC curve
roc <- simple_roc(classgap$accent=='S', scores=as.vector(predrespflm))
pdf('output/figs/flm-roc.pdf', width=4.5, height=3.5)
par(mgp=c(1.5,0.5,0),mar=c(2.5, 3, 1,0.5))
plot(roc$FPR, roc$TPR, type='l', bty='n', xlab='False positive rate',
     ylab='True positive rate', main='ROC curve')
grid <- seq(0, 1, length.out = 40)
lines(grid, grid, lty=3)
points(.165, .97, col='red', pch=20) #use confusion matrix to calculate FPR and TPR
text(.35, .97, 'Threshold = 0.5', pos=1, col='red')
dev.off()


# PLR model ####

# Scree plot from FPCA
pdf('output/figs/mfcc-screeplot.pdf', width=6, height=4)
plot(mfccpm$sdev[1:25]^2, type='b', main='Scree plot of first 25 eigenvalues',
     xlab='Component number', ylab='Eigenvalue', bty='n')
dev.off()

# PC 1 and PC 2 scores
pdf('output/figs/mfccpc12.pdf', width=6, height=4)
par(mgp=c(1.5,0.5,0), mar=c(2.6, 2.4, 1.5, 1))
plot(mfccpm$x[,1:2], col=cols[factor(nscv_log$accent)],
     bty='n', pch=20, xlab='FPC1', ylab='FPC2')
title('FPC scores for Northern and Southern vowels')
legend(x=10, y=7, legend=c('North', 'South'), pch=19, col = opcols[1:2], 
       bty='n')
dev.off()

# First 9 MFCC components of PC 2
pdf('output/figs/pc2.pdf', width=8, height=4)
par(mar=c(2,4,4,2))
plot(mfccpm$rotation[1:40,2], lty=1, type='l', col=1, xlim=c(1,360),
     bty='n', ylim=c(-0.14, 0.17), xaxt='n',
     main='First 9 MFCC components of PC 2', xlab='', ylab='')
text(20, 0.14, 'MFCC\n1', col=1, cex=0.8)
for(i in 2:9){
    x <- (40*(i-1) + 1)
    # print(x)
    lines(x:(x+39), mfccpm$rotation[x:(39+x), 2],
          col=i)
    text(40*i-20, 0.14, paste0('MFCC\n', i), col=i, cex=0.8)
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
points(.01, .98, col='red', pch=20) #use confusion matrix to calculate FPR and TPR at threshold 0.5
text(.2, .98, 'Threshold = 0.5', pos=1, col='red')
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

# MFCC accent perturbation matrix as image
pdf('output/figs/perturb.pdf', height = 4, width=4.5)
image.plot(y=1:40, z=sumcontrib, xlab='Rescaled time', ylab='MFCCs',
           main='Weighted sum of FPCs',
           col = hcl.colors(20, 'YlOrRd', rev=TRUE))
dev.off()

# First 9 MFCC components
pdf('output/figs/perturb-first9.pdf', height=4, width=8)
par(mar=c(2,4,4,2))
h <- 0.24
plot(sumcontrib[,1], lty=1, type='l', col=1, xlim=c(1,360),
     bty='n', ylim=range(sumcontrib[,1:9])+c(0,0.03), xaxt='n',
     main='First 9 MFCC components', xlab='', ylab='')
text(20, h, 'MFCC\n1', col=1, cex=0.8)
for(i in 2:9){
    x <- (40*(i-1) + 1)
    lines(x:(x+39), sumcontrib[, i],
          col=i)
    text(40*i-20, h, paste0('MFCC\n', i), col=i, cex=0.8)
}
abline(h=0, lty=3)
abline(v=seq(41, 321, by=40), lty=2, col='lightgrey')
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

# Mean vs SD of BNC predicted probs
bncmeansd <- bncprobs %>% 
    inner_join(bncprobs %>% 
                   group_by(id) %>% 
                   summarise(nsound=n()) %>% 
                   dplyr::filter(nsound>=2), by='id') %>% 
    group_by(id) %>% 
    summarise(meanf=mean(fpred),
              sdf = sd(fpred)/n(),
              meanp = mean(ppred),
              sdp = sd(ppred)/n())
head(bncmeansd) 

pdf('output/figs/preds-meansd.pdf', width=5, height=4)
plot(sdf~meanf, data=bncmeansd, main='Mean vs SD of formant model predictions')
plot(sdp~meanp, data=bncmeansd, main='Mean vs SD of MFCC model predictions')
dev.off()

# Soap film maps
plotsoapfilm(locationplr, name='aligned-avglink-beta', title='PLR')
plotsoapfilm(locationflm, name='aligned-avglink-beta', title='FLM')


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

