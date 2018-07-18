## Extract protein-coding and non-coding sequences from BLAT output

# Packages
library(ape)
data_dir <- "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"
assembled_dir <- file.path(data_dir, "chloroplast_bwa_assembled")
output_dir <- file.path(data_dir, "bwa_alignments")

# Extract out from alignments
blatcol <- c("matches", "misMatches", "repMatches", "nCount",
            "qNumInsert", "qBaseInsert","tNumInsert","tBaseInsert","strand", "qName", "qSize", "qStart", "qEnd","tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")
codingRegionPositions <- read.table(file.path(assembled_dir, "blatOutput.fsl"),
                             col.names = blatcol, stringsAsFactors = FALSE)
noncodingRegionPositions <- read.table(file.path(assembled_dir, "blatNonCodingOutput.fsl"),
                                       col.names = blatcol, stringsAsFactors = FALSE)

codingRegionPositions <- subset(codingRegionPositions, strand == "+")
noncodingRegionPositions <- subset(noncodingRegionPositions, strand == "+")

# Read alignment
rawalignment <- read.dna(file.path(assembled_dir, "pep_seq_aligned.fasta"), format = "fasta", as.matrix = TRUE)

# Extract subalignments that map back to gene sequences from the reference 
extractLoci <- function(blat, alignment){
  # Extract loci from an alignment based on blat output
  # Arguments:
  #    blat = blat output as data.frame object
  #    alignment = DNAbin object
  # Returns:
  #    list of DNAbin objects for each locus
  geneNames <- unique(blat$qName)
  extractedSeq <- list()
  for(i in geneNames){
    temp = subset(blat, qName == i)
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
    extractedSeq[[i]] <- alignment[,starts_ends]
  }
  return(extractedSeq)
}

codingRegions <- extractLoci(blat = codingRegionPositions, alignment = rawalignment)
noncodingRegions <- extractLoci(blat = noncodingRegionPositions, alignment = rawalignment)

# for(i in geneNames){
#   temp = subset(alignPositions, qName == i)
#   starts_ends <- vector()
#   for(j in 1:nrow(temp)){
#     if(temp$blockCount[j] > 1){
#       starts <- as.numeric(unlist(strsplit(temp$tStarts, split = ",")))
#       ends <- starts + as.numeric(unlist(strsplit(temp$blockSizes, split = ",")))
#       for(k in 1:length(starts)){
#         starts_ends <- c(starts_ends, starts[k]:ends[k])
#       }
#     } else {
#       starts_ends <- c(starts_ends, temp$tStart[j]:temp$tEnd[j])
#     }
#   }
#   extractedSeq[[i]] <- rawalignment[,starts_ends]
# }

# Rename the reference sequence name
for(i in 1:length(codingRegions)){
  rownames(codingRegions[[i]])[grep(rownames(codingRegions[[i]]), pattern = "KT223569")] <- "Piper_kadsura"  
}

for(i in 1:length(noncodingRegions)){
  rownames(noncodingRegions[[i]])[grep(rownames(noncodingRegions[[i]]), pattern = "KT223569")] <- "Piper_kadsura"  
}

# Export as phylip files
for(i in 1:length(codingRegions)){
  write.dna(codingRegions[[i]],
            file = file.path(output_dir, coding, paste0(names(codingRegions)[i], "_align.fasta")),
            format = "fasta", colsep = "")
}

for(i in 1:length(noncodingRegions)){
  write.dna(noncodingRegions[[i]],
            file = file.path(output_dir, "noncoding", paste0(names(noncodingRegions)[i], "_align.fasta")),
            format = "fasta", colsep = "")
}