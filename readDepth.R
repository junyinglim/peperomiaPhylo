library(stringr)

main.dir <- "~/Dropbox/Projects/2015/Peperomia/data/chloroplast_bwa_readdepth/"
output.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
readDepth <- list.files(main.dir)

meanCoverage <- vector()
sample_ID <- vector()
medianCoverage <- vector()
completeness <- vector()

for(i in 1:length(readDepth)){
  temp <- read.table(file = file.path(main.dir, readDepth[i]), header = FALSE)
  meanCoverage[i] <- mean(temp$V3)
  medianCoverage[i] <- median(temp$V3)
  completeness[i] <- sum(temp$V3 == 0) / nrow(temp)
  sample_ID[i] <- str_split_fixed(readDepth[i], pattern = "_", n = 2)[,1]
}

depthsummary <- data.frame(sample_ID, meanCoverage = signif(meanCoverage, 3), medianCoverage)

write.csv(summary, file.path(output.dir, "readDepthSummary.csv"), row.names = FALSE)

# Calculate percentage ambiguous
library(seqinr)

seq.dir <- "~/Dropbox/Projects/2015/Peperomia/data/chloroplast_bwa_assembled/"

seqfiles <- list.files(seq.dir)
assembledFastas<- seqfiles[grep(seqfiles, pattern = "_assembled.fasta")]


sample_ID <- vector()
propAmbiguous <- vector()
for(i in 1:length(assembledFastas)){

  temp <- read.dna(file.path(seq.dir, assembledFastas[i]), format = "fasta", as.character = TRUE)
  propAmbiguous[i] <- sum(tolower(temp[1,]) == "n") / length(temp)
  sample_ID[i] <- str_split_fixed(assembledFastas[i], pattern = "_", n = 2)[,1]
}
summary2 <- data.frame(sample_ID, signif(propAmbiguous, 3))

summaryFinal <- merge(depthsummary, summary2)

write.csv(summaryFinal, file.path(output.dir, "readDepthSummary.csv"), row.names = FALSE)

# What percent of reads was chloroplast