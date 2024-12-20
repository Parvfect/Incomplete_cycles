Modelling the effects of capping v no capping for different coupling rates in base level dna storage synthesis.

naive_model/ has the code I use to simulate different parts of the pipeline, from a simple synthesis model, writing it to files, recovering files and clustering. The important files are as follows - 

1. creating_synthesis_file.py
    > This simulates a single strand repeat for different coupling rates and capping flags. In a later stage, I plan to do it for different strands as well. 

I then pass those files to simulate using PBSIM - (https://github.com/pfaucon/PBSIM-PacBio-Simulator). I use the default configuration for the first model end to end. Badread can be used instead - but very slow in comparision (https://github.com/rrwick/Badread/tree/main)

pbsim --data-type CLR --depth {repeats-per-strand} --model_qc data/model_qc_clr sample/{your-sample}.fasta


2. testing_pbsim_alignment.ipynb
    > This reads the files outputted by PBSIM and attempts to recover the information using clustering - implemented in aligned_clustering.py. Most of the code is a straight copy of - https://github.com/MLI-lab/noisy_dna_data_storage/blob/master/LSH_clustering.ipynb.
    > I've changed the code a bit to run 1000 strands together, each with varying coupling rates, and capping flag on and off. I then isolate using the clusters and merge clusters based on top-k candidates (more than 0.3 matching the strand) using consensus. Seems to work really well.
    They basically cluster into groups of (5-15) and then run multiple sequence alignment in each cluster using MUSCLE (https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein), the defacto way to do it for natural bio purposes. Taking consensus of the top k candidates (right now checking against the original strand). Leads to a pretty high accuracy (on one test - need to repeat).


My goal is to make an educated guess of effect of synthesis parameters (capping/no-capping and coupling rate) on sequencing recovery - so I am not attempting to optimise recovery (at least at first).

Currently putting my code on the cluster to run for a lot of repeats to make the recovery results more consistent and a good guess. Looks like this right now, 

![output](https://github.com/user-attachments/assets/f78aa32f-70ea-4bd1-bd3b-5b2b0a33e7b4)

X axis is coupling rate and y axis is recovery percentage - github seems to black it out for some reason. Redoing the results - merging top k cluster candidates by consesus gives a 90 % + recovery percentage.
