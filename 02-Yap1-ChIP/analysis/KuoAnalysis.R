library(limma)
library(Ringo)
setwd("C:/Users/Zoheb/Desktop/HE LAB/kuo_data/")

#Read in files
file.1.gpr <- c("GSM397447.gpr")
file.2.gpr <- c("GSM397448.gpr")
file.3.gpr <- c("GSM397449.gpr")

#Get reads
reads1 <- read.maimages(files = file.1.gpr, source = 'genepix')
reads2 <- read.maimages(files = file.2.gpr, source = 'genepix')
reads3 <- read.maimages(files = file.3.gpr, source = 'genepix')

#Background Correction
RG1 <- backgroundCorrect(reads1)
RG2 <- backgroundCorrect(reads2)
RG3 <- backgroundCorrect(reads3)

par(mar=c(0.01,0.01,0.01,0.01), bg="black")
image(RG1, 1, channel="red", dim1="Column", dim2="Row", mycols=c("sienna","darkred","orangered"))
X <- preprocess(RG1[RG1$genes$ControlType==0,], method="nimblegen", idColumn="ProbeName")


print(typeof(RG1))

#Get MA values, reverse sign since dyes are given in reversed order
MA1 <- normalizeWithinArrays(RG1)
MA1$M <- (-1)*MA1$M
MA2 <- normalizeWithinArrays(RG2)
MA2$M <- (-1)*MA2$M
MA3 <- normalizeWithinArrays(RG3)
MA3$M <- (-1)*MA3$M

#Get rows where read count is > 2x for immunoprecipitated due to log2 transformation
rows.of.interest1 <- which(MA1$M > 1)
rows.of.interest2 <- which(MA2$M > 1)
rows.of.interest3 <- which(MA3$M > 1)

#Find overlapping rows/genes
overlap <- intersect(intersect(rows.of.interest1,rows.of.interest2),rows.of.interest3)
gene.hits <- MA1$genes$Name[overlap]


plot(as.factor(gene.hits),MA1$M[overlap])
print(gene.hits)
print(MA1$M[overlap])
print(MA2$M[overlap])
print(MA1$A[overlap])
