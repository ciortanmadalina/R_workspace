setwd("D:\\workspace\\bioinformatics\\BIOLF449")
#source("1.r")
print(getwd())
#?read
"
x <- c(1, 2, 3); x
is.numeric(x)
as.character(x)
b <- 7
a <- b
?c
#arrays start at 1 : x[1], negative indexes return vector without given element
?vector
vector(length = 10)
a <- c(1,2)
a
a +a
a*a -> multiplie cell by cell
a <-c(2, 3)
a * a
length(a)
TRUE FALSE
> 1:4
[1] 1 2 3 4
d <- c(FALSE, TRUE)
seq(1, 10, 2) # 2 = sequence step
a[c(1,2)]
a <- seq(1,10, 2)
a[a[c(1,2)]]
b <- c(1,5)
a[b]
m <- matrix(1:10)
m <- matrix(1:4, 2, 2)
m[1,]
m[,1]

m <-matrix(1:10, 2, 5)
m <-matrix(1:10, 2, 5, byrow = TRUE)

m <-matrix(1:10, nrow = 2, ncol = 5, byrow = TRUE)
m
dim(m)

dim(m) <- c(5,2)
a <- 1:10
dim(a) <- c(2,5)
a + a
a * a
(a +a )^2
(a+ a )**2
a
a[1:2, 2:3]
a[1:2, c(1,5)][,2][2]
a<- c(1, 'hi');a
list(1, TRUE)[[1]]
df <- data.frame(col1 = 1:10, col2=10:1)
df[2, 'col2']
as.matrix(df)
d <- data.frame(ID_REF = c('x1', 'x2', 'x3'), GSM45 = c(1.4, 4.6, 6.3), GSM46 = c(1.4, 4.6, 6.3))

d[2,]
d[2, 'GSM45']
colnames(d)
rownames(d)

colnames(d) <- c('id', 'n1', 'n1')
rownames(d) <- d[,1]
rownames(d) <- d[,'id']
d['x1',]
m <-matrix(1:10, 2, 5)
m
?sum
sum(m[,2])

colSums(m)
"
?read.table

t <- read.table("GSE17487_series_matrix.txt", header = TRUE, comment.char = "!", row.names = 1)
t[1,2]
head(t, 10)
m <- as.matrix(t)
t['10339650',]
colnames(t) <- c('WT1', 'WT2','WT3', 'S1', 'S2', 'S3')
head(t)

wt <- t[, 1:3]
head(wt)
z <- wt +1
head(z)
dim(t)
ncol(t)
nrow(t)

rowSums(t)

rowMeans(t)
colSums(t)

hist(t[,3])
