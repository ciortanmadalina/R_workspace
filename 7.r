setwd("C:\\workspace\\bioinformatics")

?pbinom
x <- 1:10
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

source("https://bioconductor.org/biocLite.R")
biocLite()

#3.2 R version
library('Biostrings')

seq <- "CTACGATCATTCG"
seq
table(seq)
s <- DNAString("CTACGATCATTCG")
length(s)

s <- DNAString('ATC')

d <- DNAString("TTGAAAA-CTC-N")
d

?DNAStringSet
x0 <- c("CTACCAGTAT", "TTGA", "TACCTAGAG")
x0[1]
s <- DNAStringSet(x0)
s[[1]]
s1 <- DNAString('AGTTTC')
reverse(s1)
reverseComplement(s1)
s
reverseComplement(s)
reverseComplement(s[[1]])
reverseComplement(s)[[1]]

alphabet(s1)
letterFrequency(s1, c('A', 'C', 'T'))
letterFrequency(s1, c('A', 'C', 'T'), as.prob = TRUE)
letterFrequency(s1, c('A', 'C', 'T'))/length(s1)

s <- readDNAStringSet('chr21.fa')
s
chr21 <- s[[1]]
chr21
length(chr21)
alphabetFrequency(chr21 )
alphabetFrequency(chr21, as.prob = TRUE )

t <-letterFrequency(chr21, c('A', 'C', 'G', 'T'))
t/sum(t) #expected frequency

p <- t/sum(t)
p
#"ATCCG"

prob <- function (s, p) {
  result <- 1
  for( i in 1 : length(s))
    result <- result * p[as.character(s[i])]
  as.numeric(result)
}
prob(DNAString("ATCCG"), p)

prob(DNAString("AA"), p)
0.25 ^2
pBernoulli <- c('A' = 0.25, 'C' = 0.25, 'G' = 0.25, 'T' = 0.25)

prob(DNAString("AA"), pBernoulli)

#5.Compute the oligomers (7nt) frequency in the chr21
oligonucleotideFrequency(chr21, width = 7)
#k nombre de occurances de oligo k

#6. on veut > x, pas inferieur
k<-10
pValue <- pbinom(k - 1 , 10000, 0.0625, lower.tail = F ) # on calcule tout ce qui est superior à k

pbinom(1040 - 1 , 10000, 0.0625, lower.tail = F )
# lower.tail ; R calcule strictement superieur, nous on veut superior ou egal, on calcule k - 1
pValue
?pbinom
# ou 
1 - pbinom(1040 + 1 , 10000, 0.0625)

#frequence attendu = bernoulli, frequence observe = comptage


##B - "h19_refSeqGenes_tp7.gz"

genes <- read.table( "h19_refSeqGenes_tp7.gz", header = TRUE, comment.char = '')
head(genes)

#data frame avec: 

d <- data.frame(id = genes$name, chr = genes$chrom, start = genes$txStart, end = genes$txEnd, strand = genes$strand)
head(d)
a <- d[d$chr == "chr21",]
head(a)
#handle strand -, end is start

plus <- a[a$strand == '+', ]
minus <- a[a$strand == '-', ]

promoteur_plus <- data.frame(chr = 'chr21', i = plus$start - 2000, j = plus$start + 500, strand = '+')
promoteur_minus <- data.frame(chr = 'chr21', i = plus$end - 500, j = plus$end + 2000, strand = '-')


extractSeq<-function(promoteur_input){
  #subseq(chr21, as.numeric(promoteur$i), as.numeric(promoteur$j))
  promoteur_input$i
}
apply(promoteur_plus[1:5,], 1 ,extractSeq)

for (i in 1:nrow(promoteur_plus)){
  start <- promoteur_plus[i, "i"]
  start <- as.numeric(start)
  start
  
  end <- promoteur_plus[i, "j"]
  end <- as.numeric(end)
  
  seq <- subseq(chr21, start, end)
  t <- oligonucleotideFrequency(seq, 7)
}


seq_promoteur_plus <- 

count_s1 <- oligonucleotideFrequency(s1, 3)
promoteur_plus

prob(DNAString('AAA'), p)

sum(oligonucleotideFrequency(s1, 1)) #l - k + 1

#faire le reverse complement pour les sequences sur le brain -

pbinom(66 - 1, 2500, prob(DNAString('AAA'), p), lower.tail = FALSE)

data.frame(list(t1,t2))


library('Biostrings')
seq <- DNAString("CTACGGATCACT")
letterFrequency(seq, c('A', 'C', 'G', 'T'))
letterFrequency(seq, c('A', 'C', 'G', 'T'), as.prob = TRUE)

setwd("C:\\workspace\\bioinformatics")
chr21 <-readDNAStringSet('chr21.fa')[[1]]
#Nucleotide frequencies as %
alphabetFrequency(chr21,as.prob = TRUE )
#Nucleotide frequencies as counts
alphabetFrequency(chr21)
#Number of missing values
alphabetFrequency(chr21)['N']
#Frequency of A/C/G/T nucleotides only
t <-letterFrequency(chr21, c('A', 'C', 'G', 'T'))
t/sum(t) #expected frequency



word <- DNAString("ATCCG")

prob <- function (s, p) {
  result <- 1
  for( i in 1 : length(s))
    result <- result * p[as.character(s[i])]
  as.numeric(result)
}
#Bernoulli model probability
pBernoulli <- c('A' = 0.25, 'C' = 0.25, 'G' = 0.25, 'T' = 0.25)
prob(word, pBernoulli)
#Nucleotide frequency probability
pNucleotide <- t/sum(t)
prob(word, pNucleotide)
#As counts (print just head)
head(oligonucleotideFrequency(chr21, 7))
#As % (print just head)
head(oligonucleotideFrequency(chr21, 7,as.prob = TRUE ))

alphabetFrequency(chr21)

x <- head(oligonucleotideFrequency(chr21, 7))

df<- as.data.frame(oligonucleotideFrequency(chr21, 7))
colnames(df) <- c('counts')
df$name <-rownames(df)
df$obsFreq <- oligonucleotideFrequency(chr21, 7,as.prob = TRUE)
df
#In bernoulli model we have the same probability for all nucleotides, so this column will be constant
#it is enough to calculate the bernoulli probability for 1 oligomer and also its pValue
sample <- DNAString("AAAAAAA")
oligoBernoulliProb <- prob(sample, pBernoulli)
df$bernoulliProb <-rep(oligoBernoulliProb, nrow(df))

observedFrequency <- alphabetFrequency(chr21,as.prob = TRUE)

calcPValue <- function(df, observedFrequency, size){
  seq <- DNAString(df['name'])
  p<-prob(seq, observedFrequency)
  k <-length(seq)
  pValue <- pbinom(k - 1 , size, p, lower.tail = F ) 
}

head(df)


df$observedPValue <-apply(df,1, calcPValue, observedFrequency, nrow(df))


k<-length(sample)
pValue <- pbinom(k - 1 , nrow(df), oligoBernoulliProb, lower.tail = F ) 
df$pValueBernoulli <-rep(oligoBernoulliProb, nrow(df))



genes <- read.table( "h19_refSeqGenes_tp7.gz", header = TRUE, comment.char = '')
genes


d <- data.frame(id = genes$name, chr = genes$chrom, start = genes$txStart, end = genes$txEnd, strand = genes$strand)
head(d)
a21 <- d[d$chr == "chr21",]

head(a21)

plus <- a21[a21$strand == '+', ]
minus <- a21[a21$strand == '-', ]

promoteur_plus <- data.frame(chr = 'chr21', i = plus$start - 2000, j = plus$start + 500, strand = '+')
promoteur_minus <- data.frame(chr = 'chr21', i = plus$end - 500, j = plus$end + 2000, strand = '-')

pp <- head(promoteur_plus)
pp

extractSeq <- function(promoteur){
  if ( as.character(promoteur['strand']) == '+') {
    seq <- subseq(chr21, as.numeric(promoteur['i']), as.numeric(promoteur['j']))
    paste(seq, collapse="") #keep string seq because it is atomic
  } else {
    seq <-reverseComplement(subseq(chr21, as.numeric(promoteur['i']), as.numeric(promoteur['j'])))
    paste(seq, collapse="") #keep string seq because it is atomic
  }
}


promoteur_plus$seq <- apply(promoteur_plus,1, extractSeq)
promoteur_minus$seq <- apply(promoteur_minus,1, extractSeq)

head(promoteur_plus)
head(promoteur_minus)


head(promoteur_plus)
#for (i in 1:nrow(promoteur_plus)){
for (i in 1:3){
  seqString <- promoteur_plus[i, "seq"]
  
  seq <- DNAString(seqString)
  t <- oligonucleotideFrequency(seq, 7)
  t
}

t <- oligonucleotideFrequency(DNAString("ATTATTCAGTT"), 3)
t + t

#As aggregating pValues for all 7nt oligos for all promoters requires 
#a well defined function, the most significant oligo is also the most encountered one
#this method thus returns the counts structure which can easily be aggregated as sum for
#all promoters by strand
calculateStats <-function (sequences, size){
  result <- calculateStatsForSequence(sequences[1])
  for (i in 2:size){
    counts <- calculateStatsForSequence(sequences[1])
    result <- result + counts
  }
  result
}

calculateStatsForSequence <-function (sequence){
  seq <-DNAString(sequence)
  counts <- oligonucleotideFrequency(seq, 7)
  observedFrequency <- oligonucleotideFrequency(seq, 7, as.prob = TRUE)
  random7NString <- DNAString("AAAAAAA")
  oligoBernoulliProb <- prob(random7NString, pBernoulli)
  df<- as.data.frame(oligonucleotideFrequency(seq, 7))
  df$name <-rownames(df)
  #observedPValue <-apply(df,1, calcPValue, observedFrequency, nrow(df))
  counts
}

r <- calculateStatsForSequence("ATTATTCAGTTATTATTCAGTT")
mat <- as.matrix(r)
head(mat)
head(apply(mat,2,sort))
sort.list(mat[,1])


mostSignificantOligos <-function (allCounts){
  df<- as.data.frame(allCounts)
  colnames(df) <- c('counts')
  df$name <-rownames(df)
  df<-df[order(df$counts, decreasing = TRUE), ]
  head(df) 
}

mostSignificantOligos(allCountsPlus)

allCountsPlus <- calculateStats(promoteur_plus$seq, nrow(promoteur_plus))

setA<-c("a", "b", "c", "d", "e")
setB<-c("d", "e", "f", "g")

union(setA,setB)

setA<-c("a", "b")
setB<-c(1,2)

m <- data.frame(value = c(1,2))
row.names(m) <- c("a", "b")
m["a",][1] = m["a",][1] + 8
m
colnames(m)
