setwd("C:/Users/Zoheb/Desktop/HE LAB/YAP1 Transcriptome Files")

# As per the E-MTAB-4457.sdrf.txt file:
# US92003687_252794010003_S01_GE2_107_Sep09_1_2.gpr is Rep1.gpr
# US92003687_252794010003_S01_GE2_107_Sep09_2_2.gpr is Rep2.gpr
# Note: Filtered SDRF is found under /data/

# A-MEXP-2402.adf.txt is GeneConversion.txt


library(marray)
library(data.table)

# Read in data
rawData <- read.GenePix(fnames = c("Rep1.gpr","Rep2.gpr"), path = ("."), name.Gf = "F532 Median",
             name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",
             name.W ="Flags", layout = NULL, gnames = NULL, targets = NULL,
             notes = NULL, skip=NULL, sep = "	", quote = "\"", DEBUG=FALSE)

# Perform LOESS normalization
data.norm <- maNorm(rawData, norm = "loess")


# Note since Cy5/Cy3 are flipped in rep1, the M value after normalization needs 
# to be multiplied by -1 for the correct log2 fold-change because M = log2(Cy5/Cy3) in maNorm()
rep1.fold.change <- (-1)*data.norm@maM[,1]
rep2.fold.change <- data.norm@maM[,2]

# Get reporter names associated with fold-change
reporter.names <- data.norm@maGnames@maInfo$Name

# Create fold-change dataframe
expression.values <- data.frame(reporter.names, rep1.fold.change, rep2.fold.change)


# Make vector containing gene names corresponding to reporter names
conversionTable <- fread("GeneConversion.txt", skip=15)
gene.names <- vector(mode = "character")
for(name in expression.values$reporter.names){
  index <- (which(grepl(name, conversionTable$`Reporter Name`)))[1]
  gene.names <- c(gene.names, conversionTable$`Reporter Database Entry[genolevures]`[index])
}


# Create gene name based fold-change dataframe
expression.values <- data.frame(gene.names, rep1.fold.change, rep2.fold.change)

# Remove rows that have both replicates as NA
expression.values <- expression.values[!(is.na(expression.values$rep1.fold.change)) 
                                       | !(is.na(expression.values$rep2.fold.change)),]

# Remove Controls
gene.expression.total <- expression.values[!(expression.values$gene.names==""),]


# For each individual gene, flatten the rows to get mean and variance
mean.log2.expression.rep1 <- vector(mode= "numeric")
variance.rep1 <- vector(mode="numeric")
mean.log2.expression.rep2 <- vector(mode= "numeric")
variance.rep2 <- vector(mode="numeric")
gene.name <- vector(mode = "character")

for(gene in unique(gene.expression.total$gene.names)){
  temp.df <- gene.expression.total[which(gene.expression.total$gene.names == gene), ]
  mean1 <- mean(rm.na(temp.df$rep1.fold.change))
  var1 <- var(rm.na(temp.df$rep1.fold.change))
  mean2 <- mean(rm.na(temp.df$rep2.fold.change))
  var2 <- var(rm.na(temp.df$rep2.fold.change))
  
  mean.log2.expression.rep1 <- c(mean.log2.expression.rep1, mean1)
  variance.rep1 <- c(variance.rep1, var1)
  mean.log2.expression.rep2 <- c(mean.log2.expression.rep2, mean2)
  variance.rep2 <- c(variance.rep2, var2)
  gene.name <- c(gene.name, gene)
}

# Final data frame with gene names, expression, and variance
transcriptome.data <- data.frame(gene.name, mean.log2.expression.rep1, variance.rep1, 
                        mean.log2.expression.rep2, variance.rep2)

# Clean up to remove rows with no gene name, fix column names
transcriptome.data <- transcriptome.data[!(is.na(transcriptome.data$gene.name)),]
names(transcriptome.data) <- c("Gene", "Mean log2 fold change Rep1", "Variance Rep1", 
                               "Mean log2 fold change Rep2", "Variance Rep2")

# Export
write.table(transcriptome.data, "./TranscriptomeData.txt", sep="\t", row.names = FALSE)