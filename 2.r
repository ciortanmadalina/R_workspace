#execute script
source("D:\\workspace\\bioinformatics\\BIOLF449\\1.r")

setwd("D:\\workspace\\bioinformatics\\BIOLF449")
"
function.name <- function(arg1, arg2)
{}
function.name('val', 'val1')

add <- function(a,b = 0)
{
  print('hi')
  a <- a+ 2
  #returns last line
  a + b
  
}
add(1)
add(add(1,2),5)
add(a = 3)

t <- read.table('GSE17487_series_matrix.txt', header = TRUE, comment.char = '!', row.names = 1)
sum(t[,1])
colSums(t)
rowMeans(t)
hist(t[,3])
hist(t[,2])
head(t)

#make sum on the second dimension = columns. First dimension = lines
apply(t, 2 , sum)

apply(t, 1 , sum)

d1 <- t[1:20,]
apply(d1, 1 , sum)
v <- c(1,2,5,7)
mean(v)
sum(v)/length(v)

mean2 <- function(v){
sum(v)/length(v)
}

apply(d1, 2 , mean2)
x <- t[,2]
y <- t[,3]
plot(x,y)
#la difference entre x et y
plot(x - y)
plot(t[,1] -t[,4], ylim = c(-10,10), ylab ='diff', col = 'blue')
col1 <- t[,1]
hist(col1, xlab = 'expression level', nclass =60, col = 'red')
hist(col1 , xlab = 'expression level', nclass =60, col = 'red', xlim = c(0,16))
apply(t, 2, hist)
x
pdf('test.pdf')
hist(x)
dev.off()

"
#delete first column
#t <- read.table('ref', header = T, comment.char = '', stringsAsFactors = F)[,-1]
t <- read.table('ref', header = T, comment.char = '')[,-1]
head(t)

d1 <- t[1:20,]
d1[1, ]
newdata <- subset(t, name = 'NM_032291')
newdata
t['name']
XM_

?grep
v <- c('A', 'B', 'A')
grep('A',v)
v[grep('A',v)]

v <- c('A12', 'BA', 'A122')
v[grep('A',v)]

names <- t[,'name']
t[c(1,2)] #select 2 lines
length(names)
names[grep('NM_',names)]

length(names[grep('NR_',names)])
length(names[grep('NM_',names)])

nrow(t)
t[grep('NM', t$name),]
nrow(t[grep('NM', t$name),])
d <-t[1:20,1:3] 
d
grep('NR',d$name)
d[grep('NR',d$name),]

#https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.what_is_the_difference_between
#https://genome.ucsc.edu/cgi-bin/hgTables
