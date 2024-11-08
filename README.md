Modelling the effects of capping v no capping for different coupling rates in base level dna storage synthesis.

naive_model/ has the code I use to simulate different parts of the pipeline, from a simple synthesis model, writing it to files, recovering files and clustering. The important files are as follows - 

1. creating_synthesis_file.py
    > This simulates a single strand repeat for different coupling rates and capping flags. In a later stage, I plan to do it for different strands as well. 

I then pass those files to simulate using PBSIM - (https://github.com/pfaucon/PBSIM-PacBio-Simulator). I use the default configuration for the first model end to end. Badread can be used instead - but very slow in comparision (https://github.com/rrwick/Badread/tree/main)

pbsim --data-type CLR --depth {repeats-per-strand} --model_qc data/model_qc_clr sample/{your-sample}.fasta


2. testing_pbsim_alignment.ipynb
    > This reads the files outputted by PBSIM and attempts to recover the information using clustering - implemented in aligned_clustering.py. Most of the code is a straight copy of - https://github.com/MLI-lab/noisy_dna_data_storage/blob/master/LSH_clustering.ipynb. 
    They basically cluster into groups of (5-15) and then run multiple sequence alignment in each cluster using MUSCLE (https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein), the defacto way to do it for natural bio purposes.
    Ideally, you take consensus post, but for now I'm taking the best candidate from all the clusters. This will change.


My goal is to make an educated guess of effect of synthesis parameters (capping/no-capping and coupling rate) on sequencing recovery - so I am not attempting to optimise recovery (at least at first).

Currently putting my code on the cluster to run for a lot of repeats to make the recovery results more consistent and a good guess. Looks like this right now, 

![image.png](https://prod-files-secure.s3.us-west-2.amazonaws.com/35c4b3df-13a8-46ea-beb6-726a908fd7bb/d6bf7cb7-fb54-45a1-b6b8-b8235145f144/image.png)
