# Exercise 2: Standard deviation, Estimates and CI

# Part 1
# ------
# read the input data file
df<-read.csv("normal_skewed_distr.csv")

head(df)

# plot histogram
hist(df$geneA, xlab="expression of geneA", main="Histogram of geneA")

# measure of central tendency
# compute mean
m<-mean(df$geneA)
# add mean value to the plot
abline(v=m, col="red", lw=3)

# measure of spread using standard deviation
s<-sd(df$geneA)

# sd and coverage/intervals
# using 68–95–99.7 rule find the interval for 68%, 95% and 99.7% 
# for 68% - mean plus/minus one standard deviation
c(m-s, m+s)
# higlight those ranges in the plot
abline(v=c(m-s, m+s), lty=2, lw=2, col="blue")
# check proporation our data points covered within this range
lt=m-s # lower limit of the range
ut=m+s # upper limit of the range
# number of data points within that range divided by total number of data points 
length(df$geneA[df$geneA>=lt & df$geneA<=ut])/length(df$geneA)*100

# for 95% - mean plus/minus two standard deviation
c(m-2*s, m+2*s) 
# higlight those ranges in the plot
abline(v=c(m-2*s, m+2*s), lty=2, lw=2, col="blue")
# check proporation our data points covered within this range
lt=m-2*s # lower limit of the range
ut=m+2*s # upper limit of the range
# number of data points within that range divided by total number of data points
length(df$geneA[df$geneA>=lt & df$geneA<=ut])/length(df$geneA)*100

# for 99% - mean plus/minus three standard deviation
c(m-3*s, m+3*s)
# higlight those ranges in the plot
abline(v=c(m-3*s, m+3*s), lty=2, lw=2, col="blue")
# check proporation our data points covered within this range
lt=m-3*s # lower limit of the range
ut=m+3*s # upper limit of the range
# number of data points within that range divided by total number of data points
length(df$geneA[df$geneA>=lt & df$geneA<=ut])/length(df$geneA)*100

# z-score - to find how many standard deviation away from the mean for a given observation
zscore<-(df$geneA-m)/s

# view actual and z-score value as two columns
cbind(df$geneA, zscore)

# plot the distribution of zscore
hist(zscore, xlab="zscore of expression of geneA")
#####################################
# TODO:
# 1) Chebyshev's theorem - non-normal distribution
# check what proportion of values covered under 2 and 3 standard deviation for geneB and geneC
# distribution to see if it follows the Chebyshev's theorem 
# (i.e, whether 2 sd covers minimum 75% of data and 3 sd covers 88% of data)
#####################################

# Part 2
# ------
# draw random samples from geneA to check if the distribution of sample means follow normal distribution

# define a variable to collect mean values
geneA_sm <- c()
# using the below function to randomly select 10 data points from geneA dataset for 100 times
for(i in c(1:100)){ # iterations from 1 to 100
  sam<-sample(df$geneA, 10) # use sample function to randomly select 10 data points
  m_sam <- mean(sam) # compute the mean of sampled values
  geneA_sm <- c(geneA_sm, m_sam) # append mean value to the list
}

# plot the distribution of sample means
hist(geneA_sm, main="distribution of samples means of gene A")

# mean of sample means
mean(geneA_sm)

# standard deviation of sample means
sd(geneA_sm)

# standard error of the mean
sem<-sd(geneA_sm)/sqrt(length(geneA_sm))
sem

# compute 95% confidence interval
# first compute marginal error using z-value for 95% interval
me <- 1.96*sem
me

lower_limit = mean(geneA_sm) - me
upper_limit = mean(geneA_sm) + me
print(c(lower_limit, upper_limit))

# highlight the confidence interval in the plot
hist(geneA_sm, main="distribution of samples means of gene A")
abline(v=mean(geneA_sm), col="red", lw=2)
abline(v=c(lower_limit, upper_limit), lw=2)
#####################################
# TODO
# 2) Central limit theorem
# draw 100 samples of size 10 from the geneB, geneC distributions separately and plot the distribution of means
# for each gene to see if they're normally distributed

# 3) uncertainty of population mean estimation and sample size relation
# draw samples from geneA distribution with sizes ranging from 10 to 200 and check the relation
# between sample size and standard error of the mean. 