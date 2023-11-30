#Here's a simple matching-alleles model
#it's a foil against the Leonard77 G4G model

k = 0.4 #the cost for the parasite of non-matching
c = 0.3 #the cost for the host of parasite matching

#pathogen genotype fitness values on each host genotype
#here we make parasite fitness directly proportional to the current fitness 
#of hosts
V11 = 1
V21 = 1 - k
V12 = 1 - k
V22 = 1

#host genotype fitess values with each pathogen genotype
W11 = 1 - c 
W21 = 1
W12 = 1
W22 = 1 - c

#starting allele frequencies in hosts and pathogens
fA1 <- 0.5
fA2 <- 1 - fA1
fC1 <- 0.3
fC2 <- 1 - fC1

#how many generations to evolve
gens=300

##A haploid selection function of matching alleles coevolution
#Extends our basic haploid selection function
haploid.coevolution <- function(fA2,fC2){
    fA1 <- 1 - fA2
    fC1 <- 1 - fC2
    #new pathogen genotype frequencies
    fC2prime <- (fC2*fA1*V21 + fC2*fA2*V22) / (fC1*fA1*V11 + fC1*fA2*V12 + fC2*fA1*V21 + fC2*fA2*V22)
    #new host genotype frequencies
    fA2prime <- (fA2*fC1*W21 + fA2*fC2*W22) / (fA1*fC1*W11 + fA1*fC2*W12 + fA2*fC1*W21 + fA2*fC2*W22)
    return(c(fA2prime,fC2prime))
    } 


#Iterate over the generations and log the frequencies of host resistace and
#pathogen virulence alleles
x <- seq(1,gens)
fA2s <- rep(0,gens)
fC2s <- rep(0,gens)

for (i in 1:gens){
    fA2s[i] <- fA2
    fC2s[i] <- fC2
    gent1 <- haploid.coevolution(fA2, fC2)
    fA2 <- gent1[1]
    fC2 <- gent1[2]
}

#Make a plot
plot(fC2s ~ x, type='l', lwd=3, lty=2, col='purple', xlab='generation', ylab='frequency of alleles', main='Matching Alleles Model of Coevolution', ylim=c(-0.2,1.1))
lines(fA2s~x, lwd=3, col='green')
legend(1, -.05, legend=c("Resistance", "Virulence"),
       col=c("green", "purple"), lty=1:2, cex=1.0)


