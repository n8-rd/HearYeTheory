#A simple quantitative genetic model. Assumes infitite sites, that is, that 
#the genome is big enough that two mutations won't happen at the same site.
#Also assumes no recombination.

##Set some model parameters
N <- 500 #set a population size
mu <- 0.05 #mutation rate per individuals, per generation
sigma <- 0.1 #variance of QTL effects
O <- 1 #optimal phenotype
s <- 3 #strength of selection (smaller is stronger)
P <- rep(0, N) #initialize phenotypes
W <- numeric(length=N) #intialize vector of fitnesses
gens <- 1000 #how many generations to simulate

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
    #this is vectorized, which makes computation much faster
    scale <- dnorm(O,O,s)
    W <- dnorm(P,O,s)/scale
    P <- sample(P,N,replace=T,prob=W)

    Pmeans[i] <- mean(P)
}

plot(Pmeans ~ seq(1,gens), type='l', col='green', lwd=3, xlab='generation', ylab='Mean Phenoype value', main='Simple Quantitative Genetic Model | O=1.0')
