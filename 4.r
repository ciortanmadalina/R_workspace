setwd("D:\\workspace\\bioinformatics")
f <- c('red', 'blue', 'red', 'red')
a <- factor(f)
a
levels(a)
b <- as.vector(a)
levels(b)
levels(a) <- c('y', 'n')

as.vector(a)
a[1]
a[1]=='n'
as.factor(f)

x <- c(1, 3, 2, 4, 2)
sapply(x, sum)
sum(x)
df <- data.frame(a = x, b = x +1)
df
sapply(df, sum)
sum(sapply(df, sum))
#equivalent to
apply(df, 2, sum)
apply(df, 1, sum)
?sapply

z <- list(c(1,2,2), c(3,6,8,9))
z
lapply(z, sum)
sapply(z, sum)
simplify2array(lapply(z, sum))


x <- c(1, 3, 2, 4, 2)
f <- c('red', 'blue', 'red', 'red', 'green')
d<- data.frame(value = x, color =f)
d$color

d$color =='red'
d$color[d$color =='red']
length(d$color[d$color =='red'])
d$value[d$color =='red']
d[d$color =='red',]
d[d$color =='red', 'value']
sum(d[d$color =='red', 'value'])

tapply(d$value, d$color, sum) 
aggregate(d$value, list(col =d$color), sum) 

x <- c('a1', 'b1', 'c1')
grep('1', x)
sub('1', '_2', x) #substitute 1 par _2 
#substiturer tous les chiffres
sub('[1-4]', '_x', x) 
sub('[0-9]', '_x', x) 

x <- c('a16', 'b18', 'c19')
sub('[0-9]+', '_x', x) 
"
http://www.ensembl.org/biomart/martview/9b92821532164ae36e99b02fe5e8603a
https://gdc-portal.nci.nih.gov/

mart_export.txt
B
maplot 1 vs 2 samples FPKM

hugo name HGNC

"
t <- read.table('mart_export.txt', header = T, comment.char = '')[,-1]
setwd("D:\\workspace\\bioinformatics\\BIOLF449\\INFOF434\\input")
t <- read.table('gdc_manifest_20161026_135408.txt', header = T)
head(t)
t$id
t$filename


readAllCountFiles <- function(filenames) {
  allCounts <- read.table(as.character(filenames[1]), header = F)
  for(i in 2:length(filenames)) {
    sample <- read.table(as.character(filenames[i]), header = F)
    allCounts<- rbind(allCounts, sample)
    #print(paste(i , ' total file length :' , nrow(allCounts), ' file len ' , nrow(sample)))
  }
  allCounts
}

allCounts <- readAllCountFiles(t$filename)


