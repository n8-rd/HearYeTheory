#A simple quantitative genetic model. Assumes infitite sites, that is, that 
#the genome is big enough that two mutations won't happen at the same site.
#Also assumes no recombination.

##Set some model parameters
N <- 500 #set a population size
mu <- 0.05 #mutation rate per individuals, per generation
sigma <- 0.1 #variance of QTL effects
#O <- 1 #optimal phenotype
#s <- 3 #strength of selection (smaller is stronger)
P <- rep(0, N) #initialize phenotypes
W <- numeric(length=N) #intialize vector of fitnesses
gens <- 1000 #how many generations to simulate

#here's polynomial fitness function that generates a multi-optima topology
fitness <- function(x){
    w <- -x^4 + 8.2*x^3 - 22.2*x^2 + 23*x
    w[w<0.2] <- 0.2
    return(w)
}

evolve <- function(){
    Pmeans <- numeric(length=gens)
    for (i in 1:gens){
        #do a round of QTL mutation
        for (n in 1:N){ #for each individual
            mut <- rbinom(1,1,mu) #roll the ol' mutation dice
            if (mut==1){ #if a mutation happened
                qtl.effect <- rnorm(1,0,sigma) #get an allele effect
                P[n] <- P[n] + qtl.effect #add it to phenotype value
            }
        }
        #do a round of selection and population regulation
        #Use our gnarly, multi-peak fitness function    
        x <- seq(0,5,0.01)
        Z <- fitness(x)
        scale <- max(Z)
        W <- fitness(P) / scale
        P <- sample(P,N,replace=T,prob=W) 
        Pmeans[i] <- mean(P)
    }
    return(Pmeans)
}

par(mfrow=c(3,3))
for (i in 1:9){
    Pmeans <- evolve()
    plot(Pmeans ~ seq(1,gens), type='l', col='magenta', lwd=3, xlab='generation', ylab='Mean Phenoype value', main='Rugged fitness topology')
}


