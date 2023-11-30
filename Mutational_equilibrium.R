#A simple model of mutational equilibrium at a di-allelic locus

##Initialize a gene pool
N <- 500 #set a population size
p <- .2 #set a starting freqeuncy for allele A
q <- 1-p #set a starting freqeuncy for allele a
As <- p * N #use those frequencies and the pop size to determine the actuall number of alleles
as <- q * N #same for the number of little a alleles
mu.aA <- 1e-4 #set rate of mutation from a to A
mu.Aa <- 1e-4 # set rate of mutation from A to a
gens = 100 # set how many generations of mutation we want to simulate
gen <- seq(1:gens) #this will help with plotting

#We can simulate a mutational allele equilibrium with a simple algorithm
#Each generation, we have some As mutate to as and vice versa, depending on our mutation rates
freqs <- numeric(length=gens)
for (i in 1:gens){
    #these next two lines are where the magic happens
    #we use R to sample from a Bernoulli random variable. This is like a coin toss where we can specify an arbitrary 
    #probability for each of the two outcomes. Here our outcomes are (0) not to mutate, (1) to mutate.
    # In each generation, for each copy of the a allele, we sample from a Bernoulli distribution with a probability
    # equal to the mutation rate for a to A. Then we sum up the number of mutations that actually happened. We do the
    #same thing for the A alleles. Then its just a matter of adding and subtracting and keeping track.
    new.A <- sum(rbinom(as,as,mu.aA)) 
    new.a <- sum(rbinom(As,As,mu.Aa))
    As <- As + new.A - new.a
    as <- as - new.A + new.a
    new.p <- As/N
    freqs[i] <- new.p
    }
    
par(mfrow=c(1,2))
plot(gen, freqs, type='l', col='green', lwd=3, xlab='generation', ylab='frequency of allele A', main='Mutational equlibrium | mu.Aa / mu.aA = 1')

##Initialize a gene pool
N <- 500 #set a population size
p <- .2 #set a starting freqeuncy for allele A
q <- 1-p #set a starting freqeuncy for allele a
As <- p * N #use those frequencies and the pop size to determine the actuall number of alleles
as <- q * N #same for the number of little a alleles
mu.aA <- 4e-4 #set rate of mutation from a to A
mu.Aa <- 1e-4 # set rate of mutation from A to a
gens = 100 # set how many generations of mutation we want to simulate
gen <- seq(1:gens) #this will help with plotting

#We can simulate a mutational allele equilibrium with a simple algorithm
#Each generation, we have some As mutate to as and vice versa, depending on our mutation rates
freqs <- numeric(length=gens)
for (i in 1:gens){
    #these next two lines are where the magic happens
    #we use R to sample from a Bernoulli random variable. This is like a coin toss where we can specify an arbitrary 
    #probability for each of the two outcomes. Here our outcomes are (0) not to mutate, (1) to mutate.
    # In each generation, for each copy of the a allele, we sample from a Bernoulli distribution with a probability
    # equal to the mutation rate for a to A. Then we sum up the number of mutations that actually happened. We do the
    #same thing for the A alleles. Then its just a matter of adding and subtracting and keeping track.
    new.A <- sum(rbinom(as,as,mu.aA)) 
    new.a <- sum(rbinom(As,As,mu.Aa))
    As <- As + new.A - new.a
    as <- as - new.A + new.a
    new.p <- As/N
    freqs[i] <- new.p
    }
    
plot(gen, freqs, type='l', col='green', lwd=3, xlab='generation', ylab='frequency of allele A', main='Mutational equlibrium | mu.Aa = 4* mu.aA')
