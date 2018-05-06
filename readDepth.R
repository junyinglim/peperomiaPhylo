library(stringr)

readDepth <- list.files("~/Desktop/depthCoverage_1st run")

setwd("~/Desktop/depthCoverage_1st run")

meanCoverage <- vector()
sample_ID <- vector()
medianCoverage <- vector()

for(i in 1:length(readDepth)){
  temp <- read.table(file = readDepth[i], header = FALSE)
  meanCoverage[i] <- mean(temp$V3)
  medianCoverage[i] <- median(temp$V3)
  sample_ID[i] <- str_split_fixed(readDepth[i], pattern = "_", n = 2)[,1]
}

mean(medianCoverage)

quantile(temp$V)
data.frame(sample_ID, meanCoverage)
min(meanCoverage)
mean(meanCoverage)
max(meanCoverage)
hist(meanCoverage)
median(meanCoverage)
