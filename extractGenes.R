## Extract genes

# Packages
library(ape)
data_dir <- "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"
assembled_dir <- file.path(data_dir, "assembled")
output_dir <- file.path(data_dir, "alignments")

# Extract out from alignments
alignPositions <- read.table(file.path(assembled_dir, "blatOutput.fsl"),
                             col.names = c("matches", "misMatches", "repMatches", "nCount",
                                           "qNumInsert", "qBaseInsert","tNumInsert","tBaseInsert","strand", "qName", "qSize", "qStart", "qEnd","tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"), stringsAsFactors = FALSE)

alignPositions <- subset(alignPositions, strand == "+")

# Read alignment
rawalignment <- read.dna(file.path(assembled_dir, "pep_seq_aligned.fasta"), format = "fasta", as.matrix = TRUE)

# Extract subalignments that map back to gene sequences from the reference 
geneNames <- unique(alignPositions$qName)
extractedSeq <- list()
for(i in geneNames){
  temp = subset(alignPositions, qName == i)
  starts_ends <- vector()
  for(j in 1:nrow(temp)){
    if(temp$blockCount[j] > 1){
      starts <- as.numeric(unlist(strsplit(temp$tStarts, split = ",")))
      ends <- starts + as.numeric(unlist(strsplit(temp$blockSizes, split = ",")))
      for(k in 1:length(starts)){
        starts_ends <- c(starts_ends, starts[k]:ends[k])
      }
    } else {
      starts_ends <- c(starts_ends, temp$tStart[j]:temp$tEnd[j])
    }
  }
  extractedSeq[[i]] <- rawalignment[,starts_ends]
}

# Rename the reference sequence name
for(i in 1:length(extractedSeq)){
  rownames(extractedSeq[[i]])[grep(rownames(extractedSeq[[i]]), pattern = "KT223569")] <- "Piper_kadsura"  
}

# Export as phylip files
for(i in 1:length(extractedSeq)){
  write.dna(extractedSeq[[i]],
            file = file.path(output_dir, paste0(names(extractedSeq)[i], "_align.fasta")),
            format = "fasta", colsep = "")
}