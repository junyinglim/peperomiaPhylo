# Peperomia of the Pacific

Molecular phylogenetic study of Peperomia (Piperaceae) in the Pacific

<img src="https://github.com/junyinglim/peperomiaPhylo/blob/master/peperomiaImage.jpg" width="300">

**Citation:** Lim, J.Y., Marshall, C.R., Zimmer, E.A. & Wagner, W.L. (in prep) Parallel radiations of peppers (Peperomia, Piperaceae) on the Hawaiian Islands suggest role of ecological release in island biotic assembly.


## Summary of pipeline
**Step 1: Reference mapping to Piper chloroplast genome (`referenceAssembly.sh`)**
* Adapter trimming using `Trimmomatic`
* Index the reference genome using `bowtie2` (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)
* Mapping to reference genome using `bowtie2`
* Evaluating coverage of chloroplast genomes
<!---
### Step 1B: Guided de novo assembly of chloroplast genomes (`denovoAssembly.sh`)
* Use GetOrganelle's (https://github.com/Kinggerm/GetOrganelle) reference mapping (`bowtie2`) coupled with read extension algorithms to filter out organelle reads
1. To use SPAdes as part of GetOrganelle pipeline, you will need to copy `spades.py` into `/usr/local/bin`
2. However, I'm not sure I want to use GetOrganelle to do the de novo assembly; tried several times but it does not seem to produce the assembly graphs
-->

**Step 2: Extract out chloroplast genes using reference annotations**
* Multiple sequence alignment using MAFFT, and then BLAT reference sequences of protein-coding genes back to the alignment (`extractAnnotation.py`)
* Using BLAT annotations, individual locus alignments are parsed out into separate alignment files (`extractLoci.R`)
* Concatenate and generate a partition file (`concatenateLoci.py`)

**Step 3: Phylogeny estimation using `IQ-TREE`**
* Step was done in CIPRES (https://www.phylo.org/portal2/)
  
**Step 4: Divergence time estimation of Piperales (`fossilBD.rev`)**
* Bayesian relaxed clock (uncorrelated log-normal) and fossilized birth-death models
* Includes 11 fossil tip constraints


**Step 5: Divergence time estimation using non-parametric rate smoothing**
* Summarize the maximum a posteriori and maximum clade credibility tree from posterior distribution of trees (`summarizeFossilBD.Rev`)
* Non-parametric rate smoothing implemented in `TreePL` (`runTreePL.py`)
* Randomly samples Peperomia-Piper divergence times from the posterior distribution of fossilized birth death trees (see Step 4)
* Used as a fixed root constraint in `TreePL`

## Plotting scripts
* `plotPhylo.R`: Plots phylogeny and maps biogeographic regions onto phylogeny
* `plotPacificMap.R`: Recenters worldmap over the Pacific, and reprojects in the Mollweide (equal-area) projection (https://en.wikipedia.org/wiki/Mollweide_projection)


<!---
## Ideas:
* Find ORF vs. using annotation codes to extract taxa
* De novo (reference-free) assembly using `SOAPdenovo` / `quickassembler` / `spades` / `velvet`
    * Align with those obtained through reference mapping to see if they are fairly congruent
    * Indels often masked when mapping to reference
    * Use `BLAST` to identify and extract chloroplast contigs
* Consensus call assuming diploid genomes to deal with heterozygosity
* Mapping to putative single-copy nuclear genes
* Use `bwa` or `novoalign` instead for mapping assemblies

## Reference-based assembly notes:
* Increasing the minimum and maximum fragment size in `bowtie`, from default (0 - 500) to (200 - 700) did not increase coverage.
* Most gaps are concentrated in non-coding regions.

* Local presets are less stringent (they don't assume the whole read aligns)
    * 'Soft-clipping' results in better alignments at ends of the genome (assembler does not consider a circular genome)
    * Reduces the number of ambiguities by half (from around 20,000 to 10,000 bases)
    * Potentially more error? Observed many "local" specific disagreements with the end-to-end assemblies (and the reference), but perhaps that is to be expected as the end-to-end requires the entire read to match perfectly

* Also impt to remember that this is based on the reference and is thus biased in an unknown way
-->