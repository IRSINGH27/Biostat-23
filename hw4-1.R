# E4 Hypothesis testing
# parametric and non-parametric tests

# parametric test for continuous data (following normal distribution)
# a) One sample Z-test
# ====================

mu<-70.8 # population mean
sd<-8.63 # population standard deviation
xbar<-78.6 # sample mean
n=14 # sample size

# state the hypothesis
# Ho = mu == xbar
# H1 = mu != xbar

# compute test statistic
z = (xbar-mu)/(sd/sqrt(n))
z

# find critical values given alpha = 0.05
alpha=0.05

# plot a standard normal distribution to understand the position of critical values
curve(dnorm, -4, 4, xlab="z", ylab="Probability (z)")

z.HA.1 <- qnorm(alpha/2, lower.tail = FALSE)
z.HA.1
abline(v=z.HA.1, col="blue")
abline(v=-z.HA.1, col="blue")

z.HA.2 <- qnorm(alpha)
abline(v=z.HA.2, col="green")
z.HA.2

z.HA.3 <- qnorm(alpha, lower.tail = FALSE)
abline(v=z.HA.3, col="orange")

print(c(z.HA.1, z.HA.2, z.HA.3))

# check if my z statistic is in the rejection area
abline(v=z, col="red", lty=2)
z
print(z>z.HA.1 || z < -z.HA.1) # two tail
print(z < z.HA.2) # lower end of the tail
print(z > z.HA.3) # upper end of the tail

# p-value approach
upper <- pnorm(abs(z), lower.tail = FALSE)
lower <- pnorm(-abs(z), lower.tail = TRUE)
p.z.1 <- upper + lower
print(p.z.1)
# since p-value of this test is lower than significance level (alpha) of 5%, we can reject the null hypothesis 
# and accept the alternate that the two means are different.

p.z.2 <- pnorm(z, lower.tail = TRUE)
print(p.z.2)
# since p-value of this test is higher than significance level (alpha) of 5%, we can't reject the null hypothesis. 

p.z.3 <- pnorm(z, lower.tail = FALSE)
print(p.z.3)
# since p-value of this test is lower than significance level (alpha) of 5%, we can reject the null hypothesis. 
# =================== #

# b) One sample t-test
# ====================
mu0 <- 70.8 # population mean
sd <- 8.63 # population standard deviation
alpha=0.05 # significance level

# here we're using the population sd mainly to randomly draw some sample data
sample_x<-rnorm(10, mean=70.8, sd=sd) # random sampling of values with the given mean and sd
# ======

# calculate t test statistic given x and population mean
tstat <- (mean(sample_x)-mu0)/(sd(sample_x)/sqrt(length(sample_x)))
tstat

# find the critical value for two tailed, given alpha = 0.05
?qt
HA1 <- qt(p=alpha/2, df=length(sample_x)-1, lower.tail=FALSE)
HA1
tstat

# check if tstat is in the critical range/rejection area
print(tstat > HA1 || tstat < -HA1) 

# tstat is not in the rejection area, so we can't reject the null hypothesis. 
# that is, the mean of the random sample is almost equal to the population mean

# p-value approach for two-sided test (in other words, for the given tstat what is the probability seeing this extreme or more)
pval = pt(abs(tstat), df=length(sample_x)-1, lower.tail = FALSE) + pt(-abs(tstat), df=length(sample_x)-1, lower.tail=TRUE)
pval
# the pvalue is larger than the significance level (alpha) of 5%, so we can't reject the null hypotheis. 

# now do that above steps thorugh the R built-in t.test function
?t.test
# for the t-test we give the sample data (x) to x paramters, and population mean to mu parameter, and
# alpha 0.05 to alpha
t.test(x = sample_x, mu = mu0, alpha = 0.05, alternative = 'two.sided')
# check if the statistic and the p-value obtained through this test matches with the above. 

# test for one sided
t.test(x = sample_x, mu = mu0, alpha = 0.05, alternative = 'greater')
# =================== #

# c) two sample t-test
# =====================
# example 1 (with equal variance)
gene1<-rnorm(100, mean=50, sd=1)
gene2<-rnorm(100, mean=50, sd=1)

# check if the data follows normal distribution
qqnorm(gene1)
qqline(gene1)
qqnorm(gene2)
qqline(gene2)

# check normality using ks test
# here we're comparing observed with the theoretical normal distribution calculated from observed mean and sd
ks.test(gene1, "pnorm", mean(gene1), sd(gene1)) 
ks.test(gene2, "pnorm", mean(gene2), sd(gene2))
# check for normality using Shapiro test
shapiro.test(gene1)
shapiro.test(gene2)

boxplot(gene1, gene2)
var.test(gene1, gene2) # F-test to check if the variance from two distribution are equal

# given that the variance are same for the two distribution, let's run t-test with equal variance parameter
t.test(gene1, gene2, var.equal = TRUE)


# example 2 (with equal variance)
gene1<-rnorm(100, mean=50, sd=1)
gene2<-rnorm(100, mean=60, sd=1)

boxplot(gene1, gene2)
var.test(gene1, gene2)

t.test(gene1, gene2, var.equal = TRUE)

# example 3 (with unequal variance)
gene1<-rnorm(100, mean=50, sd=1)
gene2<-rnorm(100, mean=60, sd=2)

boxplot(gene1, gene2)
var.test(gene1, gene2)

# given that the variance of the two distributions are difference, we will not provide variance equal to true parameter
t.test(gene1, gene2)

# example 4 (paired test)
t.test(gene1, gene2, paired = TRUE)

# ================ #

# non-parametric tests
# --------------------
# one sample test
# ===============
x<-c(0.25, 0.5, 0.81, 0.82, 0.84, 0.95)
mu = 0.4

# compute W = number of values above mu
b=sum(x>mu)
b
# change of observing a value of W under the null that is at least as extreme (W ≤ 1 or W ≥ 5) using binomial distribution
# two tailed
2*pbinom(b-1, 6, p=0.5, lower.tail = F)
# one tailed
pbinom(b-1, 6, p=0.5, lower.tail = F)

# use Wilcox.test for sign test
?wilcox.test
wilcox.test(x, mu=mu)
# note the p-value obtained from the above function will be different from above binomial
# because, the R wilcox function uses different function to compute p-value

# two sample test
# ==============
# x - same as above
y<-c(0.01, 0.41, 0.43, 0.45)
z<-c(0.04, 0.05, 0.21, 0.8)

# visualize two sample distribution x and y
boxplot(x, y)
 
# compare x and y using Wilcox rank sum test
wilcox.test(x, y)
wilcox.test(y, x)
# Note:
# wilcox function by default uses first sample to compute the W statistics (or sum of ranks),
# so you will see the different W value if you give x, y or y, x in the above function
# but p-value will be the same. 

# compare x and y using Wilcox rank sum test
boxplot(x, z)
wilcox.test(x, z)
wilcox.test(z, x)

# report the confidence interval around the difference in location
wilcox.test(z, x, conf.int = T)

# x and y with some tie values
x<-c(0.25, 0.5, 0.81, 0.82, 0.84, 0.95)
y<-c(0.01, 0.41, 0.43, 0.5)
wilcox.test(y, x)

# paired test (or two sample independent test)
x<-c(0.25, 0.5, 0.81, 0.82, 0.84, 0.95)
y<-c(0.01, 0.41, 0.43, 0.5, 0.6, 0.8)
wilcox.test(x, y, paired = TRUE)
# here the p-value is lower than the significance level of 5%, so we can reject the null hypothesis
# and accept the alternate hypothesis that the location of two distrubitions are different

# you may noticed in the above wilcox.test result, instead of W, test statics V is listed output. 
# the V stastics is just sum of positive rank values, see below on how to get this value
# V = sum of positive ranks
diff<-c(x-y) # first compute the difference from mu for each point
diff<-diff[diff!=0] # remove differences equal to zero 
diff.rank <- rank(abs(diff)) # find out the ranks of the absolute difference
diff.rank.sign<-diff.rank * sign(diff) # assign the signs to the ranks
sum(diff.rank.sign[diff.rank.sign>0]) # sum of positive ranks