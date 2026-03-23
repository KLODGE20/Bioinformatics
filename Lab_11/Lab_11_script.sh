#!/bin/bash
raxml-ng --all --msa /Users/kristenlodge/Downloads/primate.phy --model GTR --prefix /Users/kristenlodge/Downloads/primate_1 --seed 331404 --bs-metric tbe --tree rand{1} --bs-trees 100
