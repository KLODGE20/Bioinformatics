#!/bin/bash
raxml-ng --all --msa /Users/kristenlodge/Documents/GitHub/Bioinformatics/Midterm_2/metazoa_alignment.5k_modified.fasta --model GTR --prefix /Users/kristenlodge/Documents/GitHub/Bioinformatics/Midterm_2/metazoa_alignment.5k_1 --seed 41557 --bs-metric tbe --tree rand{10} --bs-trees 100
