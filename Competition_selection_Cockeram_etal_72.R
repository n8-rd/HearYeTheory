##This is the Cockeram et al. 1972 model of selection on competition traits
##di-allelic locus with alleles a and A. With diploids there are nine different
##competitive interactions between genotypes. Each has its own fitness outcome.

##The average competitive fitness of each genotype is frequency dependent
#So instead of fixed individual-level viabilities like this:
vAA <- 0.9
vAa <- 0.96
vaa <- 0.82
#We have fixed competive pair viabilities like this:
w22 <- 0.9 ##this is fitness of an AA individual competing with another AA
w21 <- 1 #AA against Aa
w20 <- 1 #AA against aa
w12 <- 1 #Aa against AA
w11 <- 1 #Aa against Aa
w10 <- 1 #Aa against aa
w02 <- 1 #aa against AA
w01 <- 1 #aa against Aa
w00 <- 0.9 # aa against aa
#Then, when we map those interaction fitneses to individual genotype,
#assuming Hardy-Weinberg we get average fitnesses that look like this:
fA <- 0.4 #just something to start with
fa <- 1-fA
wAA <- (fA^2)*w22 + 2*fA*fa*w21 + (fa^2)*w20
wAa <- (fA^2)*w12 + 2*fA*fa*w11 + (fa^2)*w10
waa <- (fA^2)*w02 + 2*fA*fa*w01 + (fa^2)*w00
#then the overall average fitness of the population is
#That's the model
###########################################################################
#Let's look to see how fA changes as a function of A. 
#If there was no frequency dependence, then this would be a flat line
fAs <- seq(0,1,0.01) #a range of freqs to look at

#We need a function that will calculate the expected change in fA after one
#generation of competition

#Our old diploid viability selection function will help:
viability.selection <- function(vAA, vAa, vaa, fA){
    fa = 1 - fA
    #all the action happens in this one line:
    fAprime <- ((vAA*fA^2) + (vAa*fA*fa)) / ((vAA*fA^2) + (2*vAa*fA*fa) + (vaa*fa^2))
    return(fAprime)
    }

#Our competition fitness function just needs to roll in the frequency-dependent
#fitnesses of each genotype. And let's have it return the difference between fA and fA
#prime instead of just the new freq of A

competition.selection <- function(fA){
    fa <- 1-fA
    wAA <- (fA^2)*w22 + 2*fA*fa*w21 + (fa^2)*w20
    wAa <- (fA^2)*w12 + 2*fA*fa*w11 + (fa^2)*w10
    waa <- (fA^2)*w02 + 2*fA*fa*w01 + (fa^2)*w00
    fA.prime <- (wAA*fA^2 + wAa*fA*fa) / (wAA*fA^2 + 2*wAa*fA*fa + waa*fa^2)
    delta <- fA.prime - fA
    return(delta)
}

deltas <- numeric(length=length(fAs))
for (i in 1:length(fAs)){
    p <- fAs[i]
    d <- competition.selection(p)
    deltas[i] <- d
}

par(mfrow=c(1,2))
plot(deltas~fAs, type='l', lwd=3, col='green', ylab='Change in fA', xlab='fA', main='Selection on Competition')
abline(h=0, lty=2, col='gray', lwd=2)

#Cool. These kinds of plots might take some getting used to. Basically, every time
#the line crosses zero, that is an equlibrium point. In other words,those are situations
#in the which fA does not change. But there are two kinds of equilibria: stable and unstable.
#What are the equilibria? Which are stable? Which are unstable?
##################################################################################
#Next let's actually, watch fA change over several generations
#cook time
gens = 100 # set how many generations we want to simulate
fA = 0.1 #set and initial frequency for fA


##Now a function that will do 'gens' numbers of rounds of competition selection
evolve.it <- function(fA, gens){
    freqs <- vector(length=gens)
    for (i in 1:gens){
        delta <- competition.selection(fA)
        fA.prime <- fA + delta
        freqs[i] <- fA.prime
        fA <- fA.prime
        }
    return(freqs)
    }
    
#plot the results
genX <- evolve.it(fA,gens)
gen <- seq(1,gens,1)
plot(gen, genX, type='l', col='magenta', lwd=3, xlab='generation', ylab='frequency of allele A', main='fA0 = 0.1', ylim=c(0,1))

#Cool! Now check this out. With different wij value, we can can some very different dynamics. 
#For example, we can have systems with more stable equilibria. Play around with the script and 
#see if you can find parameters that do this!
