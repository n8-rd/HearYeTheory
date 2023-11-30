#Rausher 84 model of habitat (or host) preference evolution

#mated females put offspring in one or two habitats
#offspring mature in natal habitat, then mate randomly with respect
#to their origin
#repeat

#For easy recall, I've put the starting value for each parameter
#in square brackets.

#proportion of eggs a mother puts in each habitat is governed by
#a di-allelic locus with alleles a and A
P1 = 0.9 #proportion of eggs put in habitat 1 by AA mothers [0.9]
P2 = 0.7 #proprtion of eggs put into habitat 1 by Aa mothers [0.7]
P3 = 0.5 #proportion of eggs put into habitat 1 by aa mothers [0.5]

#within a habitat all genotypes have equal viability
#but they may differ in fecundity. One logic for this is that mothers 
#that prefer a rarer habitat might have a lowe egg laying rate, 
#and total clutch size. Or maybe one habitat is sparser, or induces
#non fatal behavioral changes...

F1 = 1.3 #fecundity of AA moms [1.3]
F2 = 1.0 #fecundity of Aa moms [1.0]
F3 = 0.9 #fecundity of aa moms [0.9]

#soft selection; pop. regulation independent in each habitat
c = 0.5 # fraction of individuals in mating pool from habitat 1 [0.5]

#G1 is frequency of AA
#G2 is frequency of Aa
#G3 is frequency of aa
#p is frequency of A
#q is frequency of a, which is 1-p

selection <- function(p){
    q = 1-p
    G1 = p^2
    G2 = 2*p*q
    G3 = q^2
    A = G1*P1*F1 + 0.5*G2*P2*F2 #fraction of A alleles in h1 from AA mothers + fraction of A alleles in h1 from Aa mothers
    B = G1*(1-P1)*F1 + 0.5*G2*(1-P2)*F2 #fraction of A alleles in h2 from AA mothers + fraction of A alleles in h2 from Aa mothers
    M = 0.5*G2*P2*F2 + G3*P3*F3 #fraction of a alleles h1 from Aa mothers + fraction of a alleles in h1 from aa mothers
    N = 0.5*G2*(1-P2)*F2 + G3*(1-P3)*F3 #fraction of a alleles in h2 from Aa mothers + fraction of a alleles in h2 from aa mothers
    Ti = A+M #total A and a fractions in h1
    Tii = B+N #total A and a fractions in h2
    g1.prime = p*((c*A/Ti) + ((1-c)*B/Tii)) # = p*p
    g2.prime = q*((c*A/Ti) + ((1-c)*B/Tii)) + p*((c*M/Ti) + ((1-c)*N/Tii)) # = q*p + p*q
    g3.prime = q*((c*M)/Ti + ((1-c)*N/Tii)) # = q*q
    #p.prime = 0.5*p + 0.5*((c*A/Ti) + ((1-c)*B)/Tii)
    p.prime = g1.prime + 0.5*g2.prime
    delta <- p.prime - p
    return(delta)
}

ps <- seq(0,1,0.01)
deltas <- numeric()

for (p in ps){
    delta <- selection(p)
    deltas <- c(deltas, delta)
}

plot(deltas ~ ps, type='l', lwd=3, col='green', xlab='p', ylab='change in p', main='Rausher 84 Model of Habitat Preference Evolution')
abline(h=0, col='gray', lwd=3, lty=2)
