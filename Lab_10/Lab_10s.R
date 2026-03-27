# install.packages("protti")
# install.packages("UniprotR")
# install.packages("r3dmol")

library(Biostrings)
library(GenomicAlignments)
library(UniprotR)
library(protti)
library(r3dmol)

read.fasta("Elephant1.fasta")
E1DNA <- readDNAStringSet("Elephant1.fasta")
E1AA <- Biostrings::translate(E1DNA)
# translate DNA into amino acid sequence

Biostrings::writeXStringSet(E1AA, "Elephant_AA.fasta")
# make the sequence into a fasta file

Accessions <- read.csv("Lab_10_#.txt", header = FALSE)
# save your uniprotr accession numbers to a character set

GOResult <- GetProteinGOInfo(Accessions)
PlotGoInfo(GOResult)
# Gives you info on the proteins based on accession numbers

PlotGOAll(GOObj = GOResult, Top = 10, directorypath = getwd(), width = 8, height = 5)
# Pushes GO terms to github

PathR <- GetPathology_Biotech(Accessions)
Get.diseases(PathR)
# Gives you in on the diseases/pathologies associated with the gene

Fetch1 <- fetch_uniprot(Accessions$V1)
Fetch2 <- fetch_pdb(Accessions$V1)
# Provides information on things such as protein structure

Proteinpic <- fetch_alphafold_prediction(Accessions$V1)
# Provides 3D structure of protein