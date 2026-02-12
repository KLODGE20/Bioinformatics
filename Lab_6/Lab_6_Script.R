library(Biostrings)
#BiocManager::install("seqnir")
library(seqinr)
#upload packages
Elephant1 <- readDNAStringSet("Elephant1.fasta")
Elephant2 <- readDNAStringSet("Elephant2.fasta")
Elephant3 <- readDNAStringSet("Elephant3.fasta")
Elephant4 <- readDNAStringSet("Elephant4.fasta")
Elephant5 <- readDNAStringSet("Elephant5.fasta")
#upload files into R
sequences <- c(Elephant1, Elephant2, Elephant3, Elephant4, Elephant5)
sequences
#combine the files
ElephantAlignment <- msa(sequences)
library(msa)
print(ElephantAlignment, show="complete")
#Put sequences in a multiple sequence alignment use code below
Ealign <- msaMuscle(sequences)
print(Ealign, show="complete")
Estring <- readDNAStringSet(Ealign)
GC <- as(Ealign, "DNAStringSet")
letterFrequency(GC, "-")
letterFrequency(GC, "GC")
letterFrequency(GC, letters="CG", as.prob = TRUE)
#Convert into stringset and tabulate variables.
Eseqnir <- msaConvert(Ealign, type = "seqinr::alignment")
as.matrix(dist.alignment(Eseqnir, "identity"))
names(sequences) <- c("Elephant1", "Elephant2", "Elephant3", "Elephant4", "Elephant5")
#This is how closely related they are

read.fasta("Elephant1.fasta")
E1DNA <- readDNAStringSet("Elephant1.fasta")
E1AA <- Biostrings::translate(E1DNA)
E1AA
#Change DNA to amino acids



BiocManager::install("phagorn")
if (!require("BiocManager", quietly = TRUE))
install.packages("phangorn")  
library("phangorn")
inst
