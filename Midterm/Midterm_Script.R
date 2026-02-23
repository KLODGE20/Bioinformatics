#install.packages("Biostrings")
library(Biostrings)
#BiocManager::install("seqnir")
library(seqinr)
#install.packages("msa")
library(msa)
#BiocManager::install("phagorn")
library("phangorn")
#Upload packages to your script


HSsequences <- readDNAStringSet("sequences.fasta")
#Upload sequences into your script

HSAlignment <- msa(HSsequences)
print(HSAlignment, show="complete")
#Create an alignment between the sequences


HSalign <- msaMuscle(HSsequences)
print(HSalign, show="complete")
HSaligned <- as(HSalign, "DNAStringSet")
#transform MSA back into a DNA String Set

consensus <- consensusString(HSalign)
print(consensus)
#find consensus sequence "AACTCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATAAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAACTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAATTCATGTCATAGGAAGGGG"


GC <- sum(letterFrequency(HSaligned, letters = "CG", as.prob = TRUE))
print(GC)
#Calculate GC Content 10.31308


HSseqnir <- msaConvert(HSalign, type = "seqinr::alignment")
as.matrix(dist.alignment(HSseqnir, "identity"))
#Find how distantly aligned the sequences are.
#In most to least distant it is Homo_sapiens_6, Homo_sapiens_10, and Homo_sapiens_4
#The close identity scores of most of the individuals reveals that this is a good alignment.

#Question 5 There are three individuals that vary from each other 4, 6, and 10 due to point mutations.


#Question 6 for Homo_Sapiens 10 and 4 it is the hbb gene for beta globin, Accession # LC121775.1
#Question 6 for Homo_Sapiens 6 it is the HBB mutant hemoglobin beta chain, Accession # AY356351.1

read.fasta("seqdump-2.fasta")
hs6DNA <- readDNAStringSet("seqdump-2.fasta")
hs6AA <- Biostrings::translate(hs6DNA)
#Uploaded Homo_sapiens_6 into a new file so it would be easy to translate into amino acids.


writeXStringSet(hs6AA, "Homo_sapien_6_protein.fasta")
#upload protein into a fasta

#Question 8 the protein with the best match to Homo_sapiens_6 is hemoglobin subunit beta, accession # KAI2558340.1.

# There are two diseases associated with the HBB gene, sickle cell anemia and beta thalassemia. This patient has beta Thalassemia. 