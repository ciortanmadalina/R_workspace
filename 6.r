
setwd("C:\\workspace\\bioinformatics")
d1 <- data.frame(id = seq(1,10), score1 = rnorm(10))
d2 <- data.frame(id = seq(1,10), score2 = rnorm(10))

head(d1)
d3 <- cbind(d1, d2)
d2b <- d2[order(d2$score2),]
d3 <- cbind(d1, d2b)

head(d3)
d4 <- merge(d1,d2) # merge sur rownames
rownames(d4)
d4 <- merge(d1,d2b)
head(d4)
a<- rnorm(100)
summary(a)
plot(density(a))
hist(a, nclass = 20)
plot(density(a, bw = 0.1))
a<- rnorm(10000)
b<- rnorm(10000, mean= 0.2)
plot(density(a, bw = 0.5), col = 'red', lwd = 1 )
lines(density(b, bw = 0.5), col = 'blue', lwd = 3)

#correlation : cor, score divise par 100
#https://www.encodeproject.org/experiments/ENCSR000DFM/
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM999397

f1 <- read.table("tp6\\dnaME\\GSM999397_hg19_wgEncodeHaibMethylArrayH1hescHaibSitesRep1.bed.gz" )
f2 <- read.table("tp6\\dnaME\\GSM999397_hg19_wgEncodeHaibMethylArrayH1hescHaibSitesRep2.bed.gz" )
head(f1)
#score = beta value
score<- f1$V5
summary(score)
summary(f2$V5)

beta.value <- f1$V5 /100
hist(beta.value, nclass = 100)
plot(density(beta.value), lwd = 3)

#Compare the replicate (scatter plot, correlation) 
par(mfrow= c(1,2))

par(mfrow= c(1,1))
table(f1$V4 == f2$V4) #check order

plot(x = f1$V5 /100 , y = f2$V5 /100)
cor(x = f1$V5 /100 , y = f2$V5 /100)

b3 <- sample(beta.value, length(beta.value))
cor(beta.value, b3 )
plot(beta.value, b3 )
#Extract the beta values for RRBS data 
r1 <- read.table("tp6\\bedMethyl\\ENCFF001TMC.bed.gz",skip = 1 )
r2 <- read.table("tp6\\bedMethyl\\ENCFF001TMD.bed.gz",skip = 1 )

nrow(r1)
nrow(r2)
summary(r1$V5/1000) #on divise par 1000
#merge by chromosome + positions V1 + V2

plot(density(r1$V5/1000))
table(r1$V5/1000 > 0.7) # nombre des zone methylées

head(r1)

colnames(r1)
rBeta1 <- r1$V5 /100
head(rBeta1)
#Describe the beta values distribution 
summary(rBeta1)

head(r1)


k1 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep1.bed.gz", stringsAsFactors = FALSE )
k2 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep2.bed.gz" , stringsAsFactors = FALSE)
k4<- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562HaibSitesRep4.bed.gz" , stringsAsFactors = FALSE)

cor(k1$V5 /100, k2$V5 /100)
plot(density(k1$V5 /100), lwd = 3)

table(k1$V4 == k4$V4) 

commun <- intersect(k1$V4, k4$V4)

k1 <- k1[k1$V4 %in% commun,]
k1< k1[order(k1$V4),]
k4 <- k4[k4$V4 %in% commun,]
k4< k4[order(k4$V4),]
k1$V4 <- factor(k1$V4)
k4$V4 <- factor(k4$V4)
table(k1$V4 == k4$V4)

# or read table as 
k4<- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562HaibSitesRep4.bed.gz", stringsAsFactors = FALSE )

cor(k1$V5, k4$V5)


k1 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep1.bed.gz", row.names = 4 )
k2 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep2.bed.gz" , row.names = 4)
k4<- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562HaibSitesRep4.bed.gz" , row.names = 4)


k1$betak1 <- k1$V5/100
k4$betak4 <- k4$V5/100

x<- merge(k1, k4, by = 0)#merge by row names
nrow(k1)
nrow(k4)

head(x)
cor(x$betak1, x$betak4)

plot(x$betak1, x$betak4)



k1 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep1.bed.gz")
k2 <- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562SitesRep2.bed.gz")
k4<- read.table("tp6\\k562\\GSM999412_hg19_wgEncodeHaibMethylArrayK562HaibSitesRep4.bed.gz")

x<- merge(k1, k4, by = 'V4')#merge by v4
x<- merge(k1, k4, by = 'V4')#merge by row names

