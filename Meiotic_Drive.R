##This is a modification of the script we used to model
#positive selection in diploids. We make allele A have negative fitness
#effects. We also make it do meiotic drive.

##Initialize populations, set viabilities, set k, the drive parameter
fA <- 0.2 #set a starting freqeuncy for allele A
fa <- 1-fA #set a starting freqeuncy for allele a
##viabilities for diploid selection model
vAA <- 0.9 #make A deleterious
vAa <- 0.9 #A is also completely dominant
vaa <- 1.0
gens = 100 # set how many generations we want to simulate
##meiotic drive
k <- 0.5 #set this to 1/2 for the special case in which there is no drive

##Here's out drive function. In normal meiosis, in heterozygotes each
##allele has an equal chance of making it into a zygote, and then going
#on to be involved in a fertilization event. So if we focus on
#the A allele, it has a 1/2 chance of making it in to a zygote. 
#Here we let A affect that chance, changing it from 1/2 to k.
meiotic.drive <- function(fA, k){
    fa <- 1 - fA
    #get genotype freqs assuming Hardy-Weinberg
    fAA <- fA^2
    fAa <- 2*fA*fa
    faa <- fa^2
    #meiotic drive only kicks in for heterozygotes
    fAprime <- fA + fAa*(k-1/2)
    return(fAprime)
    }

    
##A diploid viability selection function 
viability.selection <- function(vAA, vAa, vaa, fA){
    fa = 1 - fA
    #all the action happens in this one line:
    fAprime <- ((vAA*fA^2) + (vAa*fA*fa)) / ((vAA*fA^2) + (2*vAa*fA*fa) + (vaa*fa^2))
    return(fAprime)
    }


##Now a function that will do 'gens' numbers of rounds of drive and viability selection
evol.dip.viability <- function(fA, gens, vAA, vAa, vaa){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        fA.gametes <- meiotic.drive(fA,k)
        fA.prime <- viability.selection(vAA, vAa, vaa, fA.gametes)
        freqs[i] <- fA.prime
        fA <- fA.prime
        }
    return(freqs)
    }
    
#call the functions we just made and plot the results
par(mfrow=c(1,2))
genX <- evol.dip.viability(fA,gens,vAA,vAa,vaa)
gen <- seq(1,gens,1)
plot(gen, genX, type='l', col='magenta', lwd=3, xlab='generation', ylab='frequency of allele A', main='Drive vs. Selection. vAA, vAa = 0.9 | k = 0.5', ylim=c(0,1))
k <- 0.6
genY <- evol.dip.viability(fA,gens,vAA,vAa,vaa)
plot(gen, genY, type='l', col='green', lwd=3, xlab='generation', ylab='frequency of allele A', main='Drive vs. Selection. vAA, vAa =0.9 | 0.6')
