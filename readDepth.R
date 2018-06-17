library(stringr)

main.dir <- "~/Dropbox/Projects/2015/Peperomia/data/chloroplast_bwa_readdepth/"
output.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
readDepth <- list.files(main.dir)

meanCoverage <- vector()
sample_ID <- vector()
medianCoverage <- vector()

for(i in 1:length(readDepth)){
  temp <- read.table(file = file.path(main.dir, readDepth[i]), header = FALSE)
  meanCoverage[i] <- mean(temp$V3)
  medianCoverage[i] <- median(temp$V3)
  sample_ID[i] <- str_split_fixed(readDepth[i], pattern = "_", n = 2)[,1]
}

summary <- data.frame(sample_ID, meanCoverage = signif(meanCoverage, 3), medianCoverage)

write.csv(summary, file.path(output.dir, "readDepthSummary.csv"), row.names = FALSE)


# What percent of reads was chloroplast