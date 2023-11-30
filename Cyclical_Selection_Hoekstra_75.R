##This is the Hoekstra 1975 model of course-grained cyclical selection
##As per normal, we've got a di-allelic locus with alleles a and A.
##A diploid species with this locus lives in a fluctating two-seasonal evironment, such that alternating generations are exposed to season 1 and season 2

#initialize allele frequencies
#p1 <- 0.2 #freq of allele A. anything will do
#q1 <- 1-p1 #freq of allele a

#frequencies of the genotypes are according to Hardy Weinberg
#fAA1 <- p1^2
#fAa1 <- 2*p1*q1
#faa1 <- q1^2

#The fitnesses of genotypes in season 1:
wAA1 <- 0.8 ##fitness of AA genotype in season 1
wAa1 <- 1.0 ##fitness of Aa genotype in season 1
waa1 <- 1.05 ##fitness of aa genotype is season 1

#The frequencies of genotypes after selection on generation 1 in season 1:
#w.bar1 <- fAA1*wAA1 + fAa1*wAa1 + faa1*waa1
#fAA1.prime <- (fAA1*wAA1)/w.bar1
#fAa1.prime <- (fAa1*wAa1)/w.bar1
#faa1.prime <- (faa1*waa1)/w.bar1
#These forms should start looking farmiliar

#OK. Now for generation 2 in season 2
#These are the starting frequencies of alleles in the new gamete pool
#p2 <- fAA1.prime + 0.5*fAa1.prime
#q2 <- 1 - p2
#And the new genotype freqs are:
#fAA2 <- p2^2
#fAa2 <- 2*p2*q2
#faa2 <- q2^2

#The fitnesses in season 2:
wAA2 <- 1.2
wAa2 <- 1.0
waa2 <- 0.9

#The frequencies of genotype after selection on generation 2 in season 2:
#w.bar2 <- fAA2*wAA2 + fAa2*wAa2 + faa2*waa2
#fAA2.prime <- (fAA2*wAA2)/w.bar2
#fAa2.prime <- (fAa2*wAa2)/w.bar2
#faa2.prime <- (faa2*waa2)/w.bar2

#Those are the pieces, let's put them together
###########################################################################
#Let's look to see how p changes as a function of p. 

#We need a function that will calculate the expected change in p after one
#two-genertation cycle of environmental fluctuation

cyclical.selection <- function(p1){
    #starting freqs
    q1 <- 1-p1
    fAA1 <- p1^2
    fAa1 <- 2*p1*q1
    faa1 <- q1^2
    
    #freqs after selection in generation 1 in season 1
    w.bar1 <- fAA1*wAA1 + fAa1*wAa1 + faa1*waa1
    fAA1.prime <- (fAA1*wAA1)/w.bar1
    fAa1.prime <- (fAa1*wAa1)/w.bar1
    faa1.prime <- (faa1*waa1)/w.bar1

    #Now for generation 2 in season 2
    #These are the starting frequencies of alleles in the new gamete pool
    p2 <- fAA1.prime + 0.5*fAa1.prime
    q2 <- 1 - p2
    #And the new genotype freqs are:
    fAA2 <- p2^2
    fAa2 <- 2*p2*q2
    faa2 <- q2^2
    #The frequencies of genotype after selection on generation 2 in season 2:
    w.bar2 <- fAA2*wAA2 + fAa2*wAa2 + faa2*waa2
    fAA2.prime <- (fAA2*wAA2)/w.bar2
    fAa2.prime <- (fAa2*wAa2)/w.bar2
    faa2.prime <- (faa2*waa2)/w.bar2
    #and the freq of A
    p3 <- fAA2.prime + 0.5*fAa2.prime
    delta <- p3 - p1
    return(delta)
}

fAs <- seq(0,1,0.01) #a range of freqs to look at
deltas <- numeric(length=length(fAs))
for (i in 1:length(fAs)){
    p <- fAs[i]
    d <- cyclical.selection(p)
    deltas[i] <- d
}

par(mfrow=c(1,2))
plot(deltas~fAs, type='l', lwd=3, col='purple', ylab='Change in p', xlab='p', main='Cyclical Selection')
abline(h=0, lty=2, col='gray', lwd=2)

#We saw this kind of plot first when we talked about competition selection. 
#Basically, every time the line crosses zero, that is an equlibrium point. 
#In other words,those are situations
#in the which p does not change. But there are two kinds of equilibria: stable and unstable.
#What are the equilibria? Which are stable? Which are unstable?
##################################################################################
#Next let's actually watch p change over several generations
#cook time
gens = 100 # set how many generations we want to simulate
p1 = 0.9 #set and initial frequency for p


##Now a function that will do 'gens' numbers of rounds of competition selection
evolve.it <- function(p1, gens){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        delta <- cyclical.selection(p1)
        p.prime <- p1 + delta
        freqs[i] <- p.prime
        p1 <- p.prime
        }
    return(freqs)
    }
    
#plot the results
genX <- evolve.it(p1,gens)
gen <- seq(1,gens,1)
plot(gen, genX, type='l', col='magenta', lwd=3, xlab='two-generation cycle', ylab='frequency of allele A', main='Starting freqs = 0.9, 0.7, 0.2', ylim=c(0,1))
#Add a couple traces for different starting alleles freqs
p1 = 0.7
genY <- evolve.it(p1,gens)
lines(gen,genY,col='coral', lwd=3, lty=2)
p1 = 0.2
genY <- evolve.it(p1,gens)
lines(gen,genY,col='pink', lwd=3, lty=3)
