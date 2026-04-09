# install.packages

library(Biostrings)
library(GenomicAlignments)
library(UniprotR)
library(protti)
library(r3dmol)

#Upload alignment
alignment = readDNAStringSet("metazoa_alignment.gene.fasta")

#Prepare alignment to translate the homo sapiens sequence
HS = alignment$"Homo_sapiens"
refineHS <- gsub("-","",HS)
Homo <- as(refineHS, "DNAStringSet") 

#Translate into protein sequence 
HSAA <- Biostrings::translate(Homo)

#Make into a fasta file
Biostrings::writeXStringSet(HSAA, "HS_AA.fasta")

#save your uniprotr accession numbers to a character ser
Accession <- read.csv("MidAccessions.txt")

#Gives you info on the proteins based on accession numbers
Goresult <- GetProteinGOInfo(Accession)
plot <- PlotGoInfo(Goresult)

# provides 3D structures of protein
PD <- fetch_alphafold_prediction(Accession)
