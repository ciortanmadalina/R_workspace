"
a <-1
class(a)
class('x')
a <- c(1,2)
length(a)
a <- c(1,2)
a <- 1
length(a)
class(c(1,2))
#vector is also numeric type
#en R tout e valeur, numeric = vector
b<- 1
b[1]
# b<-1 == b<- c(1)
b[1][1][1]
#list + data frame heterogen. Vectors autoconvert data
b<- list(1, '2', TRUE)
length(b)
class(b)
b[[2]]
#acceder aux elements [[]]
b[[1]] +1

b[[3]] && FALSE
!!!FALSE
a<-2
(a<4) &&(a>1)



"
install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)

str_to_upper("test")
v <- c('a', 'b')
v2<-sapply(v, str_to_upper)
v2

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library('Biostrings')

setwd("D:\\workspace\\bioinformatics")
#http://hgdownload.soe.ucsc.edu/downloads.html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
#cat chr21.fa chr22.fa > genome.fa
fa <- readDNAStringSet('genome.fa')

chr21_1 <-fa$chr21
chr22_1 <- fa$chr22

names(fa)
class(fa)
class(fa[[1]])

chr21 <- fa[[1]]
chr22 <- fa[[2]]

head(chr21)
?str_count

chr21_counts <- str_count(chr21, 'A')
chr21_counts

a_counts <- str_count(chr21, 'A')
c_counts <- str_count(chr21, 'C')
g_counts <- str_count(chr21, 'G')
t_counts <- str_count(chr21, 'T')
a_counts + c_counts + g_counts + t_counts


total <- str_count(chr21, c('A', 'C', 'G', 'T'))
c_g_total <- total[2] + total[3]
c_g_total
c_g_percentage <- (c_g_total * 100) / sum(total)
c_g_percentage

total <- str_count(chr22_1, c('A', 'C', 'G', 'T'))
total

t <- read.table('ref', header = T, comment.char = '')[,-1]

c_g_percentage <- function(v){
  total <- str_count(v, c('A', 'C', 'G', 'T'))
  c_g_total <- total[2] + total[3]
  c_g_percentage <- (c_g_total * 100) / sum(total)
  c_g_percentage
}

c_g_percentage(chr21_1)
c_g_percentage(chr22_1)

head(t)

annotation_21 <- t[t$chrom == 'chr21', ]
annotation_21 <- annotation_21[annotation_21$strand == '+', ]
annotation_21 <- annotation_21[grep('NM', annotation_21$name),]

annotation_21_ <- t[t$chrom == 'chr21', ]
annotation_21_ <- annotation_21_[annotation_21_$strand == '-', ]
annotation_21_ <- annotation_21_[grep('NM', annotation_21_$name),]



annotation_22 <- t[t$chrom == 'chr22', ]
annotation_22 <- annotation_22[annotation_21$strand == '+', ]
annotation_22 <- annotation_22[grep('NM', annotation_22$name),]

which(t$chrom == 'chr21')


nrow(annotation_21)
head(annotation_21)
DNAStringSet(chr21_1, start= c(1, 300), end =c(10, 400))

start_21_plus <- DNAStringSet(chr21_1, start= annotation_21$cdsStart + 1, end = annotation_21$cdsStart + 3)
start_21_moins <- DNAStringSet(chr21_1, start= annotation_21_$cdsStart -2, end = annotation_21_$cdsStart  )

fist.codon <- reverseComplement(start_21_moins)
fist.codon
table(fist.codon)
table(start_21_moins)
start_21_moins <- start_21_moins[1:5]

start_21_moins
reverse_start_21 <- sapply(start_21_moins, revert)

table(start_21_moins)

reverseComplement(DNAString('ATT'))

#la frequence du 2eme codon

head(reverse_start_21)
head(start_21_moins)
revert <- function(v){
  reverseComplement(v)
}

?reverseComplement
reverse_start_21
table(start_21)

d <- annotation_21$cdsEnd - annotation_21$cdsStart
table(d>0)
#shows frequency
table(start_21)

start_21 <- str_locate_all( chr21_1, 'ATG')
start_21
#for - strand use reverseComplement
reverseComplement(start_21)

str_locate_all( 'testtest', 'te')

?str_locate_all

pos <- gregexpr('ATG', chr21_1)
pos

s