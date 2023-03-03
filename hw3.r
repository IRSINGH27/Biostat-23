# Exercise 3: Probability and probability distributions
# a) Discrete probability distributions
# b) continuous probability distributions
# c) check if my data follows normal distribution
# d) fitting distributions
# e) TODO

# --------------------------------------
# a) Discrete probability distributions
# --------------------------------------

# binomial distribution
# Example 1: On campus, 80% of the individuals got vaccinated. If we randomly 
# sample 5 individuals, what is the probability that the number of 
# following individuals got vaccinated? 
# Exactly 2
# Fewer than 2
# 2 and above

# first find out the values for the parameters
# number of trails (n)=5
# probability of success (p)=0.8
# probability of failure (1-p) = 0.2

# for exactly 2, that is, for two successes (k=2)
# we can use the probability mass function formula of binomial distribution
factorial(5)/(factorial(2)*factorial(3))*(0.8^2)*(0.2^3)
# alternatively we can use the density function of binomial distribution available in R 
?dbinom
dbinom(2, 5, 0.8)

# for k fewer than 2, P(X<=1) = P(0) + P(1)
dbinom(0, 5, 0.8)+dbinom(1, 5, 0.8)
# alternatively we can use cumulative distribution function to get P(X<=1)
?pbinom
pbinom(1, 5, 0.8)

# for k = 2 or more 
sum(dbinom(2:5, 5, 0.8))
# alternatively use distribution function to compute left tail P(X<=1) and take one minus of it (to get the other side of the distribution)
1-pbinom(1, 5, 0.8)

#########################
# Negative binomial distribution

# Example: An exploratory study to find new species in an island (20% chance of success of finding new species), 
# what is the probability of finding 2 new species out of 6 total species spotted?

# find out the parameters required for dnbinom and input values:
?dnbinom

# number of failure (x): 4
# number of success (size): 2
# prob (p) = 0.2

# compute the probability for number of failures until the xth successes
dbinom(x=4, size=2, prob=0.2)

# extra comparison
# compute the probability for number of failures until the 1st success
dbinom(x=5, size=1, prob=0.2)
# this is equivalent to the geometric distribution (which deals for at least one success)
?dgeom
dgeom(x=5, prob=0.2)

############################
# Poisson distribution

# Example: On average a genome acquired 3 mutations per megabase (MB), 
# what is the probability that a randomly selected genomic region of 
# size 1 MB has at least one mutation?

?ppois
1-ppois(q=0, lambda=3)

# plot to see occurrence of mutations ranging from 0 to 20
plot(0:20, dpois(x=0:20, 3), type="o", xlab="number of occurences", ylab = "P(X=x)")
###########################

# --------------------------------------
# b) Continuous probability distribution
# --------------------------------------

# Normal distribution

df<-read.csv("normal_skewed_distr.csv")
?hist # to get help page
hist(df$geneA, freq=FALSE, xlab="expression of geneA", ylab="Probability of x", main="expression of geneA")
lines(density(df$geneA))

# what is the probability of finding expression value exactly 50
?dnorm
dnorm(50, mean=mean(df$geneA), sd=sd(df$geneA))

# what it the probability of finding expression value >=70
1-pnorm(69, mean=mean(df$geneA), sd=sd(df$geneA))

# what is the probability of finding expression values between >20 and <80
pnorm(80, mean=mean(df$geneA), sd=sd(df$geneA))-pnorm(20, mean=mean(df$geneA), sd=sd(df$geneA))

###########################

# ------------------------------------------------
# (c) check if my data follows normal distribution
# ------------------------------------------------
# Example 1
# randomly select some values that follows normal distribution, with mean 10 and sd 1. 
d1 <- rnorm(500, 10, 1)

# check the assumptions for normal distributions
# plot this histogram or boxplot to see if the data is symmetrically distributed
hist(d1, freq=FALSE)
lines(density(d1))

boxplot(d1) # check if the data is spread equally 

# check mean and median are located closely
summary(d1)

# Q-Q plot to see the quantile distribution
?qqnorm
qqnorm(d1) # produces a normal QQ plot of the variable
qqline(d1, col = "red", lwd = 2) # adds a reference line

# goodness of fit test
ks.test(d1, "pnorm", mean=mean(d1), sd=sd(d1)) # Kolmogorov–Smirnov test
shapiro.test(d1) # Shapiro-Wilk test (for sample size 3 to 5000)
##########

ks.test(rnorm(3000, mean = 5, sd = 3), "pnorm", mean=5, sd=3)
shapiro.test(rnorm(3000, mean = 5, sd = 3))



# Example 2
# now randomly select values using lognormal and check whether it matches the normal distribution
?rlnorm
d2 <- rlnorm(500)

# check the assumptions for normal distributions
hist(d2, freq=FALSE)
lines(density(d2))

boxplot(d2) 

summary(d2)

qqnorm(d2)
qqline(d2, col = "red", lwd = 2)

ks.test(d2, "pnorm", mean=mean(d2), sd=sd(d2)) # Kolmogorov–Smirnov test
shapiro.test(d2) # Shapiro-Wilk test (for smaller sample size)
###############

# ---------------------------------------------------
# (d) fitting distribution using fitdistrplus package
# ---------------------------------------------------
library(fitdistrplus)

# if the above package is not installed on your sytem, then install using the following command
install.packages("fitdistrplus")

# first explore our data using histogram, density plot and cummulative density function
plotdist(d1, histo = TRUE, demp = TRUE)

# function to fit distribution
?fitdist

hist(d1)
mean(d1)
sd(d1)


f1<-fitdist(d1, "norm")
plot(f1)
summary(f1)
gofstat(f1)

f3<-fitdist(d1, "gamma")
summary(f3)
###################

# now do it for d2
plotdist(d2, histo = TRUE, demp = TRUE)

# function to fit distribution
f2<-fitdist(d2, "lnorm")
plot(f2)
summary(f2)
gofstat(f2)
# --------
# e) TODO 
# --------

# An RNA sequence of length 22 is generated, with the probability of purines equals 0.7. 
# What is the probability that our sequence contains:
# Exactly 14 purines
# Fewer than 14 purines
# 10 or more purines

# Find out which distribution better explains the geneA, geneB and geneC expression values?