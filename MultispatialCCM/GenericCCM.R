require(multispatialCCM)
require(ggplot2)
require(reshape2)

EmbeddingDimension <- function(df, maxE = 10){
    ## df is a data frame with each column containing one time series
    ## the column names of df will be used for embedding dimension plotting
    ## maxE is the maximum embedding dimension considered
    ## This should be based on the amount of data available

    Edf <- as.data.frame(matrix(NA, maxE-1, ncol(df)))
    names(Edf) <- names(df)

    ## Embedding dimension can't be lower than 2
    ## Nested loops here: consider vectorizing for speed?
    for(E in 2:maxE){
        for(i in 1:ncol(df)){
            Edf[E-1, i] <- SSR_pred_boot(A = df[,i], E = E, predstep = 1, tau = 1)$rho
        }
    }

    ## plot
    Edf <- data.frame(E = 2:maxE, Edf)
    Emelt <- melt(Edf, id.vars = 'E')
    names(Emelt) <- c('E', 'Species', 'rho')
    g <- ggplot(Emelt, aes(E, rho, colour = Species)) +
        geom_line() + theme_bw()
    dev.new()
    plot(g)
}

CheckNonlinearity <- function(df, Es, timesteps = 100){
    ## df is a data frame with each column containing a time series
    ## column names will be used for plotting
    ## Es is a vector of optimal embedding dimensions
    ## timesteps is the number of prediction steps to use for signal checking

    ## error checking
    if(ncol(df) != length(Es)){
        stop('There must be exactly as many time series as embedding dimensions!')
    }
    
    Signaldf <- as.data.frame(matrix(NA, timesteps, ncol(df)))
    names(Signaldf) <- names(df)

    for(i in 1:ncol(df)){
        Signaldf[,i] <- SSR_check_signal(A = df[,i],
                                         E = Es[i],
                                         tau = 1,
                                         predsteplist = 1:timesteps)$predatout$rho
    }

    ## plotting
    Signaldf <- data.frame(TimeStep = 1:timesteps, Signaldf)
    Signalmelt <- melt(Signaldf, id = 'TimeStep')
    names(Signalmelt) <- c('TimeStep', 'Species', 'rho')
    g <- ggplot(Signalmelt, aes(TimeStep, rho, colour = Species)) +
        geom_line() + theme_bw()
    dev.new()
    plot(g)
}

RunCCM <- function(df, Es, iterations = 100, filepath = NULL){
    ## iterations is the number of bootstrapping iterations for calculating
    ## standard deviations. 1000 is a good number, but 100 is much faster.

    ## for convenience over efficiency, we're growing this data frame.
    ## This is not an efficiency roadblock since the CCM takes so much longer.
    CCMdf <- data.frame(Species.A = character(0),
                        Species.B = character(0),
                        LibraryLen = numeric(0),
                        rho = numeric(0),
                        sdevhi = numeric(0),
                        sdevlo = numeric(0))
    nComparisons <- ncol(df)*(ncol(df) - 1)
    CCMsig <- as.data.frame(matrix(NA, nComparisons, 3))
    names(CCMsig) <- c('Species.A', 'Species.B', 'pval')
    k <- 1
    for(i in 1:(ncol(df)-1)){
        for(j in (i+1):ncol(df)){
            print(paste(names(df)[i], names(df)[j], sep = ', '))
            ccm.tmp1 <- CCM_boot(df[,i], df[,j], Es[i],
                                 tau = 1, iterations = iterations)
            ccm.tmp2 <- CCM_boot(df[,j], df[,i], Es[j],
                                 tau = 1, iterations = iterations)
            ccm.sig <- ccmtest(ccm.tmp1, ccm.tmp2)
            ccmlen1 <- length(ccm.tmp1$rho)
            ccmlen2 <- length(ccm.tmp2$rho)
            tmpdf <- data.frame(Species.A = c(rep(names(df)[i], ccmlen1),
                                    rep(names(df)[j], ccmlen2)),
                                Species.B = c(rep(names(df)[j], ccmlen1),
                                    rep(names(df)[i], ccmlen2)),
                                LibraryLen = c(ccm.tmp1$Lobs, ccm.tmp2$Lobs),
                                rho = c(ccm.tmp1$rho, ccm.tmp2$rho),
                                sdevhi = c(ccm.tmp1$rho + ccm.tmp1$sdevrho,
                                    ccm.tmp2$rho + ccm.tmp2$sdevrho),
                                sdevlo = c(ccm.tmp1$rho - ccm.tmp1$sdevrho,
                                    ccm.tmp2$rho - ccm.tmp2$sdevrho))
            CCMdf <- rbind(CCMdf, tmpdf)
            CCMsig$Species.A[k:(k+1)] <- c(names(df)[i], names(df)[j])
            CCMsig$Species.B[k:(k+1)] <- c(names(df)[j], names(df)[i])
            CCMsig$pval[k:(k+1)] <- c(ccm.sig['pval_a_cause_b'],
                                      ccm.sig['pval_b_cause_a'])
            k <- k + 2
        }
    }

    ## plotting: facet by species
    g <- ggplot(CCMdf, aes(LibraryLen, rho, fill = Species.B)) +
        theme_bw() +
        geom_ribbon(aes(ymin = sdevlo, ymax = sdevhi), alpha = 0.3) +
        geom_line(aes(colour = Species.B)) + facet_grid(. ~ Species.A) +
        geom_hline(aes(yintercept = 0))
        
    ## write to file, if path given
    if(!is.na(filepath)){
        pdf(paste0(filepath, '.pdf'))
        plot(g)
        dev.off()

        write.table(CCMsig, paste0(filepath, '-pval.txt'),
                    row.names = FALSE, quote = FALSE)
    } else {
        ## plot to R device
        dev.new()
        plot(g)
        return(CCMsig)
    }
}
