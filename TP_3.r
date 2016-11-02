setwd("D:\\workspace\\bioinformatics")
pdf("histograms.pdf")
# load table data in variable t
t <- read.table('ref', header = T, comment.char = '')[,-1]

coding <- t[grep('NM', t$name),]

#Generate a table of the transcripts length (from TSS = Transcription start site 
#to TES = transcription end site)

# First approach: simply make the difference between columns
deltaApproach1 <-  coding$txEnd -coding$txStart

# Second approach : use apply function
# Numeric vector values are converted to string so we need to reconvert them to 
# integer in order to perform the operation
deltaLength <- function(v){
  as.integer(v['txEnd']) - as.integer(v['txStart'])
}
deltaApproach2 <- apply(coding, 1 , deltaLength)

# Check both appoaches produce equal results
all (deltaApproach1 == deltaApproach2)

# Plot an histogram of gene length distribution
hist(deltaApproach1, nclass =100, col = 'red', main ="Histogram of gene length distribution",
     xlab = 'Gene Length', ylab = 'Frequency')

# Plot an histogram of gene length distribution
hist(deltaApproach1, nclass =300, col = 'red', main ="Histogram of gene length distribution xlim = 200000",
     xlab = 'Gene Length', ylab = 'Frequency', xlim = c(0,200000))

#Select all the long non coding transcript (length ≥ 2000bp)
nonCoding <- t[grep('NR_',t$name),]

nonCodingLengths <- nonCoding$txEnd -nonCoding$txStart

deltaNonCoding <- nonCodingLengths[nonCodingLengths[] > 20000]
deltaNonCoding2 <- subset(deltaNonCoding, deltaNonCoding[] > 20000)

# Check both appoaches produce equal results
all (deltaNonCoding == deltaNonCoding2)

#Plot an histogram of long non coding gene length distribution
hist(deltaNonCoding, nclass =100, col = 'blue', main ="Histogram of long non coding gene length distribution",
     xlab = 'Gene Length', ylab = 'Frequency')

#Compare the coding / long non coding transcript distribution


hist(deltaApproach1, nclass =100, col = 'red', main ="Coding / long non coding transcript distribution",
     xlab = 'Gene Length', ylab = 'Frequency')
hist(deltaNonCoding, add=T, nclass =100, col = 'blue')

legend("topright", c("Long non coding lengths", "Coding lengths"), fill=c("blue", "red"))

dev.off()

source("D:\\workspace\\bioinformatics\\TP_3.r")

