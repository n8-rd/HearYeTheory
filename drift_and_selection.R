##Initialize populations and set viabilities or fitness values
N <- 50 #set a population size
fA <- 0.1#1/N #set a starting freqeuncy for allele A
fa <- 1-fA #set a starting freqeuncy for allele a
##viabilities for diploid selection model
vAA <- 1
vAa <- 0.95
vaa <- 0.9
##fitness values for the haploid selection model
wA <- 1
wa <- 0.99
gens = 500 # set how many generations we want to simulate

##A haploid selection function
haploid.selection <- function(fA,wA,wa){
    fa <- 1 -fA
    #the simplest thing is to just take the product of each allele frequency and fitness value
    #w depends on survival probability and reproductive rate
    fAprime <- (wA*fA)/((wA*fA)+(wa*fa))
    #an equivalent expression:
    #fAnext <- fA / (fA + (wa/wA)*fa)
    #another equivalent expression:
    #wa/wA = 1-s
    #fAnext <- fA / (fA+((1-s)*fa))
    return(fAprime)
    }
    
##A diploid viability selection function 
viability.selection <- function(vAA, vAa, vaa, fA){
    fa = 1 - fA
    #all the action happens in this one line:
    fAprime <- ((vAA*fA^2) + (vAa*fA*fa)) / ((vAA*fA^2) + (2*vAa*fA*fa) + (vaa*fa^2))
    #the viability coefficients determine the fraction of individuals with each genotype that 
    #suvive to adulthood and have the chance of reproducing.
    #Let's break it down. With fA and fa in the adult populaiton, 
    #the frequency of genotypes in the newborns of the next generation will be fAA=fA^2, fAa=2*fA*fa, faa=fa^2
    #That's just Hardy-Weinberg
    #The Number o those babies that survice to adulthood is given by the viabilities:
    #So AA adults = N*vAA*fA^2, Aa adults = N*vAa*2*fA*fa, aa adults = N*vaa*fa^2. Yeah?
    #We get the genotype frequencies by dividing by the number of surviving adults.
    #That is just N * v-bar, the average viability, which is (vAA*fA^2 + 2*vAa*fA*fa + vaa*fa^2) 
    #Sooo, Adult genotype freqs after selections are f'AA = vAA*fA^2/v-bar, f'Aa = 2*fA*fa/v-bar, f'aa = vaa*faa^2
    #Random mating doesn't change that, so the f'A of the next gen is f'AA + 1/2*f'Aa = (vAA*fA^2 + vAa*fA*fa) / v-bar
    return(fAprime)
    }

##A function that will do 'gens' numbers of rounds of drift and haploid selection
evol.hap <- function(N, fA, gens, wA, wa){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        fA.prime <- haploid.selection(fA, wA, wa)
        NAs <- rbinom(N, 1, fA.prime)
        fAnext <- sum(NAs) / N
        fA <- fAnext
        freqs[i] <- fA
        }
    return(freqs)
    }

##Now a function that will do 'gens' numbers of rounds of drift and viability selection
evol.dip.viability <- function(N, fA, gens, vAA, vAa, vaa){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        fA.prime <- viability.selection(vAA, vAa, vaa, fA)
        number.gametes <- 2*N
        gametes <- rbinom(number.gametes, 1, fA.prime)
        fAnext <- sum(gametes) / (2*N)
        fA <- fAnext
        freqs[i] <- fA
        }
    return(freqs)
    }
    
#call the functions we just made and plot the results
#we'll do this a bunch of times
par(mfrow=c(2,4))
for (i in 1:8){
    genX <- evol.dip.viability(N,fA,gens,vAA,vAa,vaa)
    #genX <- evol.hap(N, fA, gens, wA, wa)
    gen <- seq(1,gens,1)
    plot(gen, genX, type='l', col='magenta', lwd=2, xlab='generation', ylab='frequency of allele A', main='Drift and Selection. N=100')
    }
