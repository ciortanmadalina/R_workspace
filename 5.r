d<- data.frame( count = c(8,5,6), cat = c('a', 'b', 'c'))
order(d$count)
d[c(1),]
d[order(d$count, decreasing = T),]
d[order(-d$count, decreasing = T),]
sort(d$count)

# uniform distribution
dunif(0.5) #“d” returns the height of the probability density function
runif(3)
# normal distribution
dnorm(0.5)
rnorm(7, mean = 3, sd = 1) #“r” returns randomly generated numbers 

#“p” returns the cumulative density function
#“q” returns the inverse cumulative density function (quantiles)

dnorm(2)
dnorm(3, mean =2) # donne la probabilité d'obtenir x pour mean (0 default ), sf=1

rnorm(3, mean =2) # distribution poisson
rpois(10, 3.4)
x <-rnorm(10, mean =2)
hist(x)


s <- seq(-10, 20, 1)

hist(x, breaks = seq(-10, 10, 1))

hist(rnorm(100, mean =0), breaks = seq(-10, 10, 0.1))
hist(rnorm(10000, mean =0), breaks = seq(-10, 10, 0.1))


a <- rnorm(20, mean = 3, sd = 1)
b <- rnorm(20, mean = 3, sd = 1)
a-b
plot(a - b)
hist(a, breaks = s)

quantile(a)
quantile(a, c(0, 0.5))
quantile(b)

plot(quantile(a), quantile(b))

# Quantile Quantile plot
qqplot(a, b) 

qqplot(rnorm(100), rnorm(100), xlim =c(-5,5), ylim = c(-5,5)) 

duplicated(col1, col2) #donne les elements communs


#remove after. 
rownames(counts) <- gsub('[.][0-9]+', rownames(counts))

fpkm <- function(col) {
  1000 * 1e6 *col/size/sum(col)
  }

#ENCFF154GIE.bed.gz ENCFF682FCZ.bed.gz

setwd("D:\\workspace\\bioinformatics")


f1 <- read.table("ENCFF154GIE.bed.gz" , sep = '\t')
f2 <- read.table("ENCFF682FCZ.bed.gz" , sep = '\t')

head(f1)

"
chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
The 9 additional optional BED fields are:

name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
shade	 	 	 	 	 	 	 	 	 
score in range  	≤ 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	≥ 945
strand - Defines the strand - either '+' or '-'.


5th: integer score for display
7th: fold-change
8th: -log10pvalue
9th: -log10qvalue
10th: relative summit position to peak start

"
head(f1)
#Load the peaks in R and with usable column names 
colnames(f1) <- c('chrom', 'start', 'end', 'peak', 'score', 'strand', 'foldChange', 'logPValue', 'logQValue', 'pos')
colnames(f2) <- c('chrom', 'start', 'end', 'peak', 'score', 'strand', 'foldChange', 'logPValue', 'logQValue', 'pos')

f1$pValue <- 10 ** - f1$logPValue
f2$pValue <- 10 ** - f2$logPValue
#What is the lowest and highest pvalue for replicate 1? 

min(f1$pValue)
max(f1$pValue)

#Sort the peaks based on their Fold Change. What are the Ids of top 10 peaks?
sorted <- f1[order(-f1$foldChange),]
top10Peaks <- sorted[1:10,]$peak
top10Peaks

#Draw the Pvalue distribution for each replicate (histogram) 
hist(f1$pValue,  breaks = seq(0, 0.01, 0.0001), col='red')

pnorm(0) #probabilite d'etre superior à 0
dnorm(0) # probabilite d'observer 0

pnorm(0, lower.tail = T)

pnorm(62.5, mean = 60, lower.tail = F)
pnorm(62.5, mean = 60, lower.tail = T)

hist(rnorm(1000, mean = 60), nclass = 100)
pnorm(62, mean = 60, lower.tail = F)
pnorm(rnorm(100, mean=60), mean = 60, lower.tail = F) #pvaleur

hist(pnorm(rnorm(100, mean=60), mean = 60, lower.tail = F) , breaks = seq(0,1,0.1))

#Add the expected Pvalue on the plots
pValue <- function(x){
  pnorm(x, mean=mean(x), sd = sd(x), lower.tail = F)
}

f1$expectedPValue <- pValue(f1$pValue)
hist(f1$expectedPValue)
min(f1$expectedPValue)
mean(f1$expectedPValue)
qqplot(f1$pValue,f1$expectedPValue )
