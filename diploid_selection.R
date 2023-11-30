##Initialize populations and set viabilities or fitness values
fA <- 0.9 #set a starting freqeuncy for allele A
fa <- 1-fA #set a starting freqeuncy for allele a
##viabilities for diploid selection model
vAA <- 0.9
vAa <- 0.96
vaa <- 0.82
gens = 100 # set how many generations we want to simulate

    
##A diploid viability selection function 
#We make this it's own function, and this call this function inside of the 
#next function to make it easier to extend the model later to account for genetic drift
viability.selection <- function(vAA, vAa, vaa, fA){
    fa = 1 - fA
    #all the action happens in this one line:
    fAprime <- ((vAA*fA^2) + (vAa*fA*fa)) / ((vAA*fA^2) + (2*vAa*fA*fa) + (vaa*fa^2))
    return(fAprime)
    }


##Now a function that will do 'gens' numbers of rounds of drift and viability selection
evol.dip.viability <- function(fA, gens, vAA, vAa, vaa){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        fA.prime <- viability.selection(vAA, vAa, vaa, fA)
        freqs[i] <- fA.prime
        fA <- fA.prime
        }
    return(freqs)
    }
    
#call one of the functions we just made and plot the results
#par(mfrow=c(1,2))
genX <- evol.dip.viability(fA,gens,vAA,vAa,vaa)
gen <- seq(1,gens,1)
plot(gen, genX, type='l', col='magenta', lwd=2, xlab='generation', ylab='frequency of allele A', main='Diploid Selection. vAA=0.9 | vAa=0.96 | vaa=0.8', ylim=c(0,1))
fA <- 0.1
fa <- 1-fA
genY <- evol.dip.viability(fA,gens,vAA,vAa,vaa)
lines(gen,genY,col='green', lwd=2)
#plot(gen, genY, type='l', col='green', lwd=2, xlab='generation', ylab='frequency of allele A', main='Diploid Selection. vAA=0.9 | vAa=0.96 | vaa=0.8')
