library(ggplot2)
library(dplyr)

scaleFUN <- function(x) sprintf("%.5f", x)

d <- read.csv('flatFlux.csv')
d <- d[-1,] #an ugly hack, drops the first line in the log file which tends to be junk
gd <- d %>% group_by(V) %>% summarise(winner=mean(winner))
p <- ggplot(d, aes(x=V, y=winner)) + 
geom_point() + 
geom_bar(data=gd, stat="identity", alpha=.3) +

scale_x_continuous(trans='log10') +
theme_minimal() + 
labs(x="Rate of Env. Change (bigger is slower)", y="Prop. of simulation in which L1 wins", title="Survival of the Flattest!")

print(p)


