##Initialize a gene pool
N <- 40 #set a population size
p <- .2 #set a starting freqeuncy for allele A
q <- 1-p #set a starting freqeuncy for allele a
As <- rep('A',N*p) #use those frequencies and the pop size to determine the actuall number of alleles
as <- rep('a',N*q)
alleles <- c(As, as)
allele.pool <- sample(alleles,40,replace=FALSE)
gens = 100 # set how many generations of drift we want to simulate

#We can simulate a generation of genetic drift with a simple algorithm
#Follow the following three steps N times (where N is the number of haploid individuals):
# 1. Chose an allele at random from the N alleles in the parental generation
# 2. Make a copy of that allele
# 3. Place that copy in the new generation
#That's it. But note that we are assuming non-overlapping generations and fixed population size

##Now a function that will do 'gens' numbers of rounds of reproduction and drift
drift <- function(N, allele.pool, gens){
  freqs <- vector(length=gens) #make an empty container to store p at each generation
  for (i in 1:gens){ #loop over generations
    p <- length(which(allele.pool=='A')) / N #calculalte p
    freqs[i] <- p #stick p in the container
    next.generation <- character() #make an empty container to hold next generations's alelles
    for (j in 1:N){ #loop over individuals in the next generation
      pick <- sample(allele.pool,1,replace=TRUE) #each picks a parental allele at random from the pool
      next.generation <- c(next.generation, pick) #these get added to the next gen. allele container
    }
    allele.pool <- next.generation #the complete pool for this generation, becomes the source for the next one
  }
  return(freqs) #at the end, we want a record of the frequency of A at each generation
}

#call, the function we just made and plot the results
#we'll do this a bunch of times
par(mfrow=c(2,4))
for (i in 1:8){
    genX <- drift(N,allele.pool,gens)
    gen <- seq(1,gens,1)
    plot(gen, genX, type='l', col='magenta', lwd=3, xlab='generation', ylab='frequency of allele A', main='Drift in action. N=100')
    }
