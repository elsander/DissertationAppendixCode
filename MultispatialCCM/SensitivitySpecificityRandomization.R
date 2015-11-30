## sensitivity = true positives/(true positives + false negatives)
## = ability to correctly identify existing interactions
## specificity = true negatives/(true negatives + false positives)
## = ability to correctly identify absent interactions

sensitivity <- function(testmat, tatoosh){
    truePos <- sum(tatoosh & testmat)
    falseNeg <- sum(tatoosh & !testmat)
    return(truePos/(truePos + falseNeg))
}

specificity <- function(testmat, tatoosh){
    trueNeg <- sum(!tatoosh & !testmat)
    falsePos <- sum(!tatoosh & testmat)
    return(trueNeg/(trueNeg + falsePos))
}

ccmRandomization <- function(tatoosh, ccm, reps = 10000000){
    ccmSensitivity <- sensitivity(ccm, tatoosh)
    ccmSpecificity <- specificity(ccm, tatoosh)
    print(paste('ccm sensitivity:', ccmSensitivity))
    print(paste('ccm specificity:', ccmSpecificity))
    if(reps %% 10000) print('# of reps will be rounded up to the nearest 10,000')

    sampCCM <- rep(list(ccm), 10000)
    S <- length(ccm)
    nMoreSensitive <- nMoreSpecific <- 0
    for(i in 1:ceiling(reps/10000)){
        ##since it's reshuffling the original vector over
        ##and over, the only bias this will introduce will be a result
        ##of the randomizer used
        sampCCM <- lapply(sampCCM, sample, size = S,
                            replace = FALSE)
        sensitivities <- unlist(lapply(sampCCM, sensitivity, tatoosh = tatoosh))
        specificities <- unlist(lapply(sampCCM, specificity, tatoosh = tatoosh))
        nMoreSensitive <- nMoreSensitive + sum(sensitivities >=
                                               ccmSensitivity)
        nMoreSpecific <- nMoreSpecific + sum(specificities >=
                                             ccmSpecificity)
    }
    sensitivityPval <- nMoreSensitive/(ceiling(reps/10000)*10000)
    specificityPval <- nMoreSpecific/(ceiling(reps/10000)*10000)
    return(list(sensitivityPval = sensitivityPval, specificityPval = specificityPval))
}

tatoosh <- rbind(c(1,1,1,1,0),
                 c(1,1,1,1,0),
                 c(1,1,1,0,0),
                 c(1,1,1,1,1),
                 c(1,1,1,0,1))
ccm <- rbind(c(0,0,0,1,0),
             c(0,0,0,0,0),
             c(0,0,0,1,1),
             c(1,0,1,0,0),
             c(1,0,0,1,0))
colnames(ccm) <- rownames(ccm) <- c('Mc', 'Mt', 'Pp', 'Sc', 'Cv')
colnames(tatoosh) <- rownames(tatoosh) <- c('Mc', 'Mt', 'Pp', 'Sc', 'Cv')

## those are the networks in matrix form, but for the randomization
## it's convenient to flatten them into vectors
tatoosh <- as.vector(tatoosh)
ccm <- as.vector(ccm)
##set.seed(7415623)
pvals <- ccmRandomization(tatoosh, ccm)
