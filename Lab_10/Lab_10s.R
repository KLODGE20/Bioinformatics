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

Accessions <- read.csv("Lab_10_#.txt")
# save your uniprotr accession numbers to a character set

GOResult <- GetProteinGOInfo(Accessions)
PlotGoInfo(GOResult)
# Gives you info on the proteins based on accession numbers

PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5)
# Pushes GO terms to github

GetPathology_Biotech(Accessions)
Get.diseases(Accessions)
# Gives you in on the diseases/pathologies associated with the gene

fetch_uniprot(Accessions)
fetch_pdb(Accessions)
# Provides information on things such as protein structure

fetch_alphafold_prediction(Accessions)
# Provides 3D structure of protein