# hypothesis testing for categorical data

# goodness-of-fit test
# ====================

#Example 1
# check whether the bird count during a week is normally distributed
# Mon, Tues, Wed, Thur, Fri
c(23, 18, 17, 19, 23)

sum(c(23, 18, 17, 19, 23))
# if we assume that the counts are uniformly distributed, then the probability for each category should be total_number divided by number of categories
# in our case, it will be 100/5, probability = 0.2. 
# If we multiply probability with total counts, we get an expected count for each category
0.2*100

# compute chi-square statistic sum((0bs-exp)**2/exp)
sum = (((23-20)**2)/20) + (((18-20)**2)/20) + (((17-20)**2)/20) + (((19-20)**2)/20) + (((23-20)**2)/20)
sum

# compute p-value using cdf function for upper tail 
# degress of freedom = number of categories - number of parameters estimated - 1
pchisq(sum, df=4,  lower.tail=FALSE)

# use R built-in function
chisq.test(c(23, 18, 17, 19, 23), p=c(0.2, 0.2, 0.2, 0.2, 0.2))

# view all values computed by chisq.test
tmp<-chisq.test(c(23, 18, 17, 19, 23), p=c(0.2, 0.2, 0.2, 0.2, 0.2))
names(tmp)
tmp$observed
tmp$expected # expected value calculated by chisq.test for the given probabilities
# -----------

#Example 2
#Two alleles related to seed shape: B and b  
#Crossing of B and b = BB, Bb, bb with probabilities 0.25, 0.5, 0.25, respectively
#Observed frequency: BB (BB + Bb): 5474, bb: 1850
#Expected probabilities: BB: 0.75, bb: 0.25

OBB=5474
Obb=1850
EBB = 0.75*(OBB+Obb)
Ebb = 0.25*(OBB+Obb)
print(matrix(c(OBB, EBB, Obb, Ebb), ncol=2,  nrow=2))

# chi-square statistic
sum = (((OBB-EBB)**2)/EBB) + (((Obb-Ebb)**2)/Ebb)
 
# degree of freedom: 2-1
df=1

# compute p-value for the right-tail
pchisq(sum, df=df,  lower.tail=FALSE)

# run above with the R built-in function
?chisq.test
chisq.test(c(OBB,Obb), p=c(0.75, 0.25))
# ------------

# Example 3
# Assume two alleles A and a with frequency p and 1 − p, respectively. 
# Under random mating and survival, we expect offspring genotypes AA, Aa and aa with probabilities p2, 2p(1 − p) and (1 − p)2, respectively. (Hardy–Weinberg equilibrium) 
# we simulate this for 100 randomly selected individuals, together with expected counts based on p = 0.25 

p=0.25
pAA = p*p
pAa = 2*p*(1-p)
paa = (1-p)**2
print(c(pAA, pAa, paa))

chisq.test(c(10,50,40), p=c(pAA, pAa, paa))
# --------

# Example 4
# same as above but, we estimate probability of A from the observed data
p = (2 * 10 + 1 * 50)/200
p
pAA = p*p
pAa = 2*p*(1-p)
paa = (1-p)**2
print(c(pAA, pAa, paa))

# don't run the below test alone, because the below test doesn't know that we have computed one parameter
# chisq.test(c(10,50,40), p=c(pAA, pAa, paa)) 

# we have to compute the p-value with df = 1 (number of categories - number of parameters test -1; in our case, it will be 3-1-1)
out<-chisq.test(c(10,50,40), p=c(pAA, pAa, paa)) # get the chisquare statistics from this test
out$parameter <-c(df=1) # edit the degree of freedom value
out$p.value = pchisq(out$statistic, df=out$parameter, lower.tail = FALSE) # recompute p-value with the new df value
out
# --------

# Example 5
# see slide 10 for more details
# estimate probability using binomial distribution for the number of trails from 0 to 4
for (i in seq(0, 4)){print(dbinom(i, 4, 0.5))}

# use the probabilities for the goodness-of-fit test
chisq.test(c(10, 41, 70, 57, 22), p=c(0.0625, 0.2500, 0.3750, 0.2500, 0.0625)) 


# Tests of independence
# ---------------------
# using chisquare and fisher's exact test
# Example 1:
# see slide 14 for experiment details
matrix(c(3,45,38,7,5,2),2, byrow=TRUE)

chisq.test(matrix(c(3,45,38,7,5,2),2, byrow=TRUE))

fisher.test(matrix(c(3,45,38,7,5,2),2, byrow=TRUE))

# Example 2: with larger expected counts
# see slide 15 for experiment details
chisq.test(matrix(c(3,27,8,7,23,30),2, byrow=TRUE))

fisher.test(matrix(c(3,27,8,7,23,30),2, byrow=TRUE))

# Example 3: fisher's exact test
# see slide 16 for experiment details
# check if the number of oncogenes is enriched in chromsome 1 as compared to the rest of the genome
fisher.test(matrix(c(100,2000,300,6000),2, byrow=TRUE))

fisher.test(matrix(c(300,500,3000,6000),2, byrow=TRUE))
