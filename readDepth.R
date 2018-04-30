readDepth <- list.files("~/Desktop/hello")

setwd("~/Desktop/hello")

meanCoverage <- vector()
sample_ID <- vector()

for(i in 1:length(readDepth)){
  temp <- read.table(file = readDepth[i], header = FALSE)
  meanCoverage[i] <- mean(temp$V3)
  sample_ID[i] <- str_split_fixed(readDepth[i], pattern = "_", n = 2)[,1]
  
}

data.frame(sample_ID, meanCoverage)
min(meanCoverage)
mean(meanCoverage)
max(meanCoverage)
