#Decay of D over generations

#Specify some parameters
c = 0.00142 # recombination rate
D.t0 = 1# initial value of linkage coefficient
gens <- seq(1,1000) #run it for 1000 generations

#Now do the calculations. Note that we could 
#do this by looping over every element in gens, but
#that would be inefficient and slow.
#Instead we can take advantage of the fact that 
#almost everything in R is a vector. This let's us 
#*Vectorize* calculations. Which mean that we can tell
#R to use some highly-optimized linear algebra libraries
#Vectorized code also tend to be simpler:

Ds <- (1-c)^gens * D.t0

#That's it!

#Now Plot it.
plot(Ds ~ gens, type='l', col='magenta', main='Decay of D over generations', lwd=2, xlab='generations', ylab='D')


