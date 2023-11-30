##Here's a much simpler implementation of genetic drift


##Set some parameters
N <- 50 #set a population size
p <- .7 #set a starting freqeuncy for allele A
gens = 100 # set how many generations of drift we want to simulate
nReps = 100 #how many replications to run

#define drift function
drift <- function(N, p, gens){
    ps <- numeric()
    P <- p # make a copy of p, so that we don't reset the global version
    for (i in 1:gens){
        ps <- c(ps, P)
        p.prime <- sum(rbinom(N, 1, P))/N #here's all the action
        P <- p.prime
    }
    return(ps)
}

#call, the function we just made a bunch of times
finalFreqs <- numeric()
for (i in 1:nReps){
    finalFreq <- drift(N, p, gens)[gens]
    finalFreqs <- c(finalFreqs, finalFreq)
}

fixed <- length(finalFreqs[finalFreqs == 1])
print(paste("Number of runs that ended with fixation = ", fixed))


