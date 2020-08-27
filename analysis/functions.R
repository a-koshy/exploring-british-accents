# Preprocessing functions ####

# Create trimmed wav files using timestamps because ahocoder only processes the whole file
trim_wav <- function(infile, outfile, start, dur){
    system(paste("sox", infile, outfile, 'trim', start, dur))
}

sound_to_mfcc <- function(wavfile, ccord = 39, lframe = 80, n = 2000) {
    
    # temp files to store binary output for f0 and mfccs
    file_f0 <- tempfile('f0')
    file_fcc <- tempfile('fcc')
    file.create(file_f0)
    file.create(file_fcc)
    
    ret <- system(
        paste("ahocoder16_64",
              wavfile,
              file_f0,
              file_fcc,
              paste0("--CCORD=", ccord),
              paste0("--LFRAME=", lframe),
              sep = " "
        ),
        intern = TRUE
    )
    
    if (length(attr(ret, "status")) != 0) { # if there is an error
        f0 <- NULL
        ret1 <- system(paste("ahocoder16_64",
                             wavfile,
                             file_f0,
                             file_fcc,
                             paste0("--CCORD=", ccord),
                             paste0("--LFRAME=", lframe),
                             "--CCMETH=2",
                             sep = " "
        ),
        intern = TRUE
        )
        
        if (length(attr(ret1, "status")) != 0) { # if there is still an error
            return(invisible(list(mfcc = NULL, f0 = NULL)))
        }
    } else {
        ## reading f0
        to.read <- file(file_f0, open="rb")
        f0 <- readBin(to.read, numeric(), size = 4, n = n)
        close(to.read)
    }
    
    ## reading mfcc
    to.read <- file(file_fcc, "rb")
    tmp.read <- readBin(to.read, numeric(), size = 4, n = (ccord + 1) * n)
    close(to.read)
    mfcc <- t(matrix(tmp.read, (ccord + 1), length(tmp.read) / (ccord + 1)))
    rm(tmp.read)
    
    return(invisible(mget(c("mfcc", "f0"))))
}

# Input: list of MFCC matrices, 40 coefficients each but different durations
# Parameters: number of time points to measure smooth curves at.
# Output: array with dimensions nobs*timelength*40, vector of durations
smooth_mfccs <- function(raw_mfccs, resolution=40){
    nobs <- length(raw_mfccs) # number of observations
    nc <- ncol(raw_mfccs[[1]]) # number of MFCCs for each sound
    sound.lengths <- rep(NA, nobs)
    
    # If normalised time points not specified, take max of observed lengths
    if(is.na(resolution) ){
        resolution <- sapply(raw_mfccs, nrow) %>% 
            max()
    }
    
    # create grid on 0-1
    time_grid <- seq(0, 1, len=resolution)
    mfccs_grid <- array(NA, c(resolution, nobs, nc))
    
    # create progress bar
    print('Smoothing the functions')
    pb <- txtProgressBar(min = 0, max = nobs, style = 3)
    # For each observation, and for each column with a different MFCC
    # Fit a cubic spline in 0-1 using F0 vector and LOOCV
    # Store MFCC fitted values along the 0-1 grid, and original length of MFCCs
    for (k in 1:nobs){
        setTxtProgressBar(pb, k)
        for (i in 1:nc){ 
            fit <- smooth.spline(x=seq(0, 1, len=nrow(raw_mfccs[[k]])),
                                 y=raw_mfccs[[k]][, i], cv=TRUE)
            mfccs_grid[, k, i] <- predict(object=fit, x=time_grid)$y
        }
        sound.lengths[k] <- nrow(raw_mfccs[[k]])
    }
    close(pb)
    
    return(list(mfcc_norm=mfccs_grid, durations=sound.lengths))
}

# Smoothing formants with cubic splines, similar to MFCC smoothing
smooth_form <- function(raw_formants, resolution=40){
    nobs <- length(raw_formants) # number of observations
    nc <- ncol(raw_formants[[1]]) # number of formants for each sound
    sound.lengths <- rep(NA, nobs)
    
    # If normalised time points not specified, take max of observed lengths
    if( is.na(resolution) ){
        resolution <- sapply(raw_formants, nrow) %>% 
            max()
    }
    
    # create grid on 0-1
    time_grid <- seq(0, 1, len=resolution)
    form_grid <- array(NA, c(resolution, nobs, nc))
    
    # create progress bar
    print('Smoothing the functions')
    pb <- txtProgressBar(min = 0, max = nobs, style = 3)
    # For each observation, and for each column with a different formant
    # Fit a cubic spline in 0-1 using F0 vector and LOOCV
    # Store fitted values along the 0-1 grid, and original length of formants
    for (k in 1:nobs){
        setTxtProgressBar(pb, k)
        for (i in 1:nc){ 
            fit <- smooth.spline(x=seq(0, 1, len=nrow(raw_formants[[k]])),
                                 y=raw_formants[[k]][, i], cv=TRUE)
            form_grid[, k, i] <- predict(object=fit, x=time_grid)$y
        }
        sound.lengths[k] <- nrow(raw_formants[[k]])
    }
    close(pb)
    return(list(formant_norm=form_grid, durations=sound.lengths))
}

# Robust smoothing of formants with loess
smooth_robust <- function(raw_formants, resolution=40, r=0.4, degree=1){
    nobs <- length(raw_formants) # number of observations
    nc <- ncol(raw_formants[[1]]) # number of formants for each sound
    sound.lengths <- rep(NA, nobs)
    
    # If normalised time points not specified, take max of observed lengths
    if( is.na(resolution) ){
        resolution <- sapply(raw_formants, nrow) %>% 
            max()
    }
    
    # create grid on 0-1
    time_grid <- seq(0, 1, len=resolution)
    form_grid <- array(NA, c(resolution, nobs, nc))
    
    # create progress bar
    print('Smoothing the functions')
    pb <- txtProgressBar(min = 0, max = nobs, style = 3)
    # For each observation, and for each formant curve
    # Smooth with loess: span r and degree 2 or 1
    # Store fitted values along the 0-1 grid, save original length of formants
    for (k in 1:nobs){
        x <- seq(0,1,length.out = lengths(raw_formants)[k]/nc)
        setTxtProgressBar(pb, k)
        for (i in 1:nc){
            fit <- loess(raw_formants[[k]][,i]~x, span = r, degree = degree,
                         family = 'symmetric')
            form_grid[, k, i] <- predict(fit, time_grid)
        }
        sound.lengths[k] <- nrow(raw_formants[[k]])
    }
    close(pb)
    return(list(formant_norm=form_grid, durations=sound.lengths))
}


# Register MFCC or formant matrices according to given warpings
# Input: array of warpings, time-normalised curves to warp, index used to align
# Output: array of registered curves same dim as before,
align_matrix <- function(warpings, smooth.curves){
    aligned.sounds <- array(dim = dim(smooth.curves))
    nc <- dim(smooth.curves)[3] #number of MFCCs/formants curves per sound
    nobs <- dim(smooth.curves)[2]
    time.grid <- seq(0, 1, length.out = nrow(smooth.curves))
    # Interpolate each MFCC/formant curve using the corresponding warping function
    for (i in 1:nobs){ # each observation
        for (j in (1:nc)){ # every MFCC/formant curve
            aligned.sounds[, i, j] <- approx(time.grid, 
                                             smooth.curves[, i, j],
                                             xout=warpings[, i],
                                             rule=2)$y
        }
    }
    return(aligned.sounds)
}


get_inverse_warp2 <- function(warping.function, x.out){
    len <- length(warping.function)
    sapply(x.out, function (y) uniroot((function (x) approx(seq(0, 1, len=len),
                                                            warping.function, xout=x)$y  - y),
                                       lower = 0, upper = 1)[[1]], simplify='array')
}

mfcc_preprocess <- function(sound.data, sound.f0, cf=1, resolution=NA, ...)
{
    
    n <- length(sound.data)
    nc <- dim(sound.data[[1]])[2]
    
    sounds.lengths <- rep(NA, n)
    
    # time rescaling
    if( is.na(resolution) ){
        x <- max(lengths(sound.f0))
    } else {
        x  <-  resolution
    }
    eq.grid <- seq(0, 1, len=x)
    Data.grid <- array(NA, c(n, x, nc))
    f0.grid <- matrix(NA, n, x)
    
    # create progress bar
    
    print('Smoothing the functions')
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (k in 1:n){
        setTxtProgressBar(pb, k)
        for (i in 1:nc){ 
            #print(paste(k, i, sep='-'))
            fit <- smooth.spline(seq(0, 1, len=dim(sound.data[[k]])[1]),
                                 sound.data[[k]][, i], cv=TRUE)
            Data.grid[k, , i] <- predict(fit, eq.grid)$y
        }
        fit <- smooth.spline(seq(0, 1, len=dim(sound.data[[k]])[1]), sound.f0[[k]], cv=TRUE)
        f0.grid[k, ] <- predict(fit, eq.grid)$y
        sounds.lengths[k] <- length(sound.f0[[k]])
    }
    close(pb)
    
    abs <- seq(0, 1, len=x)
    
    FV <- align_fPCA(t(Data.grid[, , cf]), abs, ...)
    
    aligned.sounds <- array(NA, dim(Data.grid))
    aligned.f0 <- matrix(NA, n, x)
    for (i in 1:n){
        for (j in 1:nc){
            aligned.sounds[i, , j] <- approx(abs, Data.grid[i, , j], xout=FV$gam[, i], rule=2)$y
        }
    }
    
    for (i in 1:n){
        aligned.f0[i, ] <- approx(abs, f0.grid[i, ], xout=FV$gam[, i], rule=2)$y
    }
    
    warping.functions <- t(FV$gam)
    out <- list(mfcc=aligned.sounds, f0=aligned.f0, sounds.lengths=sounds.lengths, 
                warping.functions=warping.functions)
    return(invisible(out))
}

# Calculate distance in kilometers between two locations
earth.dist <- function (long1, lat1, long2, lat2){
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- long1 * rad
    b1 <- lat2 * rad
    b2 <- long2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 6378.145
    d <- R * c
    return(d)
}

# Animating F1-F2 curves ####
makegif <- function(nframes, duration, formantmatrix, log){
    tlen <- 40
    N <- dim(formantmatrix)[2]
    nobs <- sample(1:N, N)
    f1 <- formantmatrix[1:tlen, nobs, 1] %>% t() %>%  
        as.data.frame() %>% 
        mutate(id=log$id[nobs] %>% factor(),
               accent=log$accent[nobs] %>% factor(),
               word=log$word[nobs] %>% factor(),
               mic=log$mic[nobs] %>% factor(),
               obs=nobs %>% factor()) %>%
        gather(key=t, value=f1, -(41:45))
    
    f2 <- formantmatrix[1:tlen, nobs, 2] %>% t() %>%  
        as.data.frame() %>% 
        mutate(id=log$id[nobs] %>% factor(),
               accent=log$accent[nobs] %>% factor(),
               word=log$word[nobs] %>% factor(),
               mic=log$mic[nobs] %>% factor(),
               obs=nobs %>% factor()) %>%
        gather(key=t, value=f2, -(41:45))
    
    f1f2 <- f1 %>%
        inner_join(f2, by=c('id', 'obs', 't', 'accent', 'word', 'mic')) %>% 
        mutate(t=as.numeric(substr(t, 2, 5)),
               id=paste('Speaker', substr(id, 2, 2)) %>% factor())

    p <- ggplot(f1f2, aes(f1, f2, group=obs, col=accent))+
        geom_line(alpha=0.1)+
        geom_point(alpha=0.7)+
        transition_reveal(t)+
        theme_minimal()+
        scale_x_continuous(breaks=seq(0,1400, by=400),
                           limits=c(0,1400))+
        labs(title='Formant trajectories: time = {frame}',
             subtitle = "Northern and Southern vowels",
             x= 'F1 (Hz)', y='F2 (Hz)')+
        scale_colour_hue(labels=c('North', 'South'))
    
    anim <- animate(p, nframes=nframes, duration=duration)
    return(anim)
}

# Creating sounds ####
# MFCC to sound helper. length 40 - increase this to 100 or 200 or more?
# F0 is set to be constant pitch
mfcc_to_sound <- function(soundlist, output.file=NULL, lframe=80) {
    wavfile <- tempfile('sound', fileext='.wav')
    file_f0 <- tempfile('f0')
    file_fcc <- tempfile('fcc')
    
    zz <- file(file_f0, "wb")
    writeBin(soundlist$f0, zz, size = 4)
    close(zz)
    
    zz <- file(file_fcc, "wb")
    writeBin(as.vector(t(soundlist$mfcc)), zz, size = 4)
    close(zz)
    
    system(
        paste('ahodecoder16_64', file_f0, file_fcc, wavfile, 
              paste0("--LFRAME=", lframe),
              sep=" ")
    )
    
    sound <- readWave(wavfile)
    if(is.null(output.file))
    {
        return(invisible(sound))
    } else { ## output to file
        writeWave(sound, filename=output.file, extensible=FALSE)
        return(invisible(sound))
    }
}

createsound <- function(mfccmatrix, length=100, file=NULL, 
                        f0=NULL){
    if(length != nrow(mfccmatrix)){
        print('Smoothing MFCCs')
        a <- smooth_mfccs(list(mfccmatrix), resolution = length)$mfcc_norm[,,]
    } else {a <- mfccmatrix}

    if(is.null(f0)){
        f0 = rep(5, times=length)
    }
    
    if (is.null(file)){
        file=paste0('sounds/', deparse(substitute(mfccmatrix)), '.wav')
    }
    mfcc_to_sound(soundlist=list(mfcc=a, f0=f0),
                  output.file = file)
}

# Perturb sound towards S/N accent using PLR model
perturbsound <- function(infile, outfile, starttime, endtime, seqlength=5,
                         NtoS=T){
    # Extract MFCCs from infile
    insound <- sound_to_mfcc(infile)
    inmfcc <- insound$mfcc
    inf0 <- insound$f0
    # Identify MFCC section corresponding to vowel
    start <- floor(starttime*200)
    end <- ceiling(endtime*200)
    # Rescale perturbation matrix accordingly
    vowelcontrib <- smooth_mfccs(list(sumcontrib), resolution=end-start+1)
    # Pad with zeroes
    addmfcc <- matrix(0, nrow=nrow(inmfcc), ncol=40)
    addmfcc[start:end,] <- vowelcontrib$mfcc_norm
    
    # Create MFCC matrix with all sounds
    if (isTRUE(NtoS)) {add <- 1} else {add <- -1}
    alteredseq <- inmfcc
    alteredf0 <- inf0
    # Add a gap of silence between words
    if(!exists('sil')|!exists('f0sil')){
        silence <- sound_to_mfcc('sounds/silence.wav')
        sil <- silence$mfcc[1:50,]
        f0sil <- silence$f0[1:50]
    }
    for (i in 2:seqlength){
        alteredseq <- rbind(alteredseq, sil, inmfcc+(i*add*addmfcc))
        alteredf0 <- c(alteredf0, f0sil, inf0)
    }
    # Create sound outfile
    createsound(mfccmatrix=alteredseq, f0=alteredf0, file=outfile, 
                length=nrow(alteredseq))
    return(print('Success!'))
}

# Plotting and diagnostics ####

# ROC curve
simple_roc <- function(labels, scores){
    labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(TPR=cumsum(labels)/sum(labels), 
               FPR=cumsum(!labels)/sum(!labels), 
               labels)
}

# Soap film map point estimate and standard error plots
plotsoapfilm <- function(locationmodel, name='soapfilm', title='FLM'){
    # plot response scale
    p <- plot(locationmodel, n2=150, asp=1.3, too.far=0.5)
    xscale <- p[[1]]$xscale
    yscale <- p[[1]]$yscale
    fit <- p[[1]]$fit
    raw <- p[[1]]$raw
    
    # Plot standard errors
    pse <- plot(locationmodel, scheme=3, se=1.96, n2=150, asp=1.3, too.far=50)
    pxscale <- pse[[1]]$x
    pyscale <- pse[[1]]$y
    
    # linear predictors and SEs for M and F
    pfit1 <- matrix(pse[[1]]$fit, nrow=150, ncol=150)
    pfit2 <- matrix(pse[[1]]$se, nrow=150, ncol=150)
    
    zlimrange <- range(plogis(pfit1-pfit2), plogis(pfit1+pfit2), na.rm=T)
    
    pdf(paste0('output/figs/', title, name, '-map.pdf'))
    par(mar=c(0,0,2,0), mgp=c(1.5,0.5,0))
    image.plot(x=xscale, y=yscale, z=plogis(fit), asp=1.3,
               col=terrain.colors(100),
               xlab='', ylab='', 
               zlim=zlimrange, axes=F, bty='n',
               main=paste0('Probability of Southern vowel: ', title))
    lines(bn, col='grey50')
    contour(x=xscale, y=yscale, z=plogis(fit), add=T, labcex = 0.7, nlevel=7)
    points(raw$lon, raw$lat, pch=3)
    dev.off()
    
    pdf(paste0('output/figs/', title, name, '-se-maps.pdf'), width=4, height=4)
    par(mar=c(0,0,1.5,0), mfrow=c(1,2), mgp=c(1.5,0.5,0))
    image(x=pxscale, y=pyscale, z=plogis(pfit1-pfit2), asp=1.3, zlim=zlimrange,
          col=terrain.colors(100),
          xlab='', ylab='', axes=F)
    title('2.5th percentile', line=-2)
    lines(bn, col='grey50')
    contour(x=pxscale, y=pyscale, z=plogis(pfit1-pfit2), add=T, labcex = 0.7,
            nlevel=7)
    points(raw$lon, raw$lat, pch=3, cex=0.7)
    
    image(x=pxscale, y=pyscale, z=plogis(pfit1+pfit2), asp=1.3, zlim=zlimrange,
          col=terrain.colors(100),
          xlab='', ylab='', axes=F)
    title('97.5th percentile', line=-2)
    lines(bn, col='grey50')
    contour(x=pxscale, y=pyscale, z=plogis(pfit1+pfit2), add=T, labcex = 0.7,
            nlevel=7)
    points(raw$lon, raw$lat, pch=3, cex=0.7)
    mtext(paste('Confidence interval:', title), side = 3, line = -1.5, outer = TRUE)
    dev.off()
}
