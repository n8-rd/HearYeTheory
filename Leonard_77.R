#Leonard 1977

#host has resistance locus A
#A1 hosts are suceptible and A2 hosts are resistant
#two pathogen genotypes C1, C2
#Vij is the fitness of pathogen i on host j
#Wij is the fitness of host i attacked by pathogen j
#On resistant host A2, fitness of narrowly virulent pathogen C1 cut by factor t
#factor t is thus the effectiveness of resistance
#C2 is the broadly resistant pathogen, counter-adapted to resistant host 
#C2 genotype comes with fitness penalty k, i.e., vost of broad virulene
#k cuts C2 fitness on all host genotypes
#optionally factor a gives C2 a species boost on the resistant host
#this is instead of just having two ks
#hosts with the resistance allele pay a fitness penalty c, i.e., the cost of resistanc
#all hosts infected by C1 or C2
#host fitness reduced proportional to pathogen strain relative fitness, weighted by s
#which is the disease severity factor

#Set model parameters
k=0.2 #broad virulence fitness penalty
t=0.3 #effectiveness of resistance 
c=0.1 #host resistance cost
a=0.05 #broad resistance advantage on resistance host
s=0.4 #disease severity

#pathogen genotype fitness values on each host genotype
V11 = 1
V21 = 1 - k
V12 = 1 - t
V22 = 1 - k + a

#host genotype fitess values with each pathogen genotype
#Note that host fitness is inversly proportional to within-host pathogen fitness
#this is the classical trade-off model of virulence evolution
W11 = 1 - s*V11
W21 = 1 - c - s*V12
W12 = 1 - s*V21
W22 = 1 - c - s*V22

#starting allele frequencies in hosts and pathogens
fA1 <- 0.5
fA2 <- 1 - fA1
fC1 <- 0.3
fC2 <- 1 - fC1

#how many generations to evolve
gens=500

##A haploid selection function of G4G coevolution
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
plot(fC2s ~ x, type='l', lwd=3, lty=2, col='purple', xlab='generation', ylab='frequency of alleles', main='Leonard 1977 Model of Coevolution', ylim=c(0,1))
lines(fA2s~x, lwd=3, col='green')
legend(1, 0.97, legend=c("Resistance", "Virulence"),
       col=c("green", "purple"), lty=1:2, cex=1.0)


