##Initialize populations and set viabilities or fitness values
fA <- 0.2 #set a starting freqeuncy for allele A
fa <- 1-fA #set a starting freqeuncy for allele a

##fitness values for the haploid selection model
##But wait! There is pleiotropy
#so for phenotype 1 the fitnesses are:
wA1 <- 1.3
wa1 <- 1.0
#and for phenotype 2 the fitnesses are:
wA2 <- 0.9
wa2 <- 1.0
gens = 100 # set how many generations we want to simulate

##A haploid selection function
haploid.selection <- function(fA){
    fa <- 1 -fA
    #the simplest thing is to just take the product of each allele frequency and fitness value
    #w depends on survival probability and reproductive rate
    fAprime <- (wA1*wA2*fA)/((wA1*wA2*fA)+(wa1*wa2*fa))
    return(fAprime)
    }
    

##A function that will do 'gens' numbers of rounds haploid selection
evol.hap <- function(fA, gens){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        fA.prime <- haploid.selection(fA)
        freqs[i] <- fA.prime
        fA <- fA.prime
        }
    return(freqs)
    }

#make some plots
par(mfrow=c(1,2))
genX <- evol.hap(fA,gens)
gen <- seq(1,gens,1)
plot(gen, genX, type='l', col='magenta', lwd=3, xlab='generation', ylab='frequency of allele A', main='Selection w/ Pleiotropy. wA1=1.3, wa1=1, wA2=0.9, wa2=1')
#Now let's change the allele effects on the second phenotype
#and for phenotype 2 the fitnesses are:
wA2 <- 0.7
wa2 <- 1.0
genX <- evol.hap(fA,gens)
plot(gen, genX, type='l', col='red', lwd=3, xlab='generation', ylab='frequency of allele A', main='Selection w/ Pleiotropy. wA1=1.3, wa1=1, wA2=0.7, wa2=1')
