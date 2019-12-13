library(marray)
library(data.table)

# Cy5 is red, Cy3 is green

# Rep 1
# mena_wt.gpr is cy5 (wild type)
# mena_sy22.gpr is cy3 (yap1skn7 double deletion)

# Rep2
# mena_sy.gpr is cy5 (yap1skn7 double deletion)
# mena_wt22.gpr is cy3 (wild type)

# Read in data
rawData <- read.GenePix(fnames = c("mena_wt22.gpr", "mena_sy.gpr"), path = ("."))

maPlot(rawData)
