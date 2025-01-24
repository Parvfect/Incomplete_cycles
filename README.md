Modelling the effects of capping v no capping for different coupling rates in base level dna storage synthesis.

Model is in v2. End to end pipeline automated does not exist. The synthesis files are created by synthesis_pipeline.py, where you can change the coupling_rates, sim_repeats, strand repeats and length. Currently only works with random strand, but will need functionality to run with user drawn strands.
The synthesis simulation yields a few run files and the synthesised pool in fasta. This pool is run through Badread, using the required coverage, which outputs a fastq file. This fastq file is used by cluster_reads.py (with command line support) that generates possbile strand candidates by aligned clustering. These candidates should be tested against the original strands to see the recovery percentages.

The requirements are in v2 as requirements.txt. Additionally, you will need to install MUSCLE - https://www.drive5.com/muscle/, and make the code in aligned_clustering.py (line 22 - 26) point to where your MUSCLE package is stored. Will be changed to get from the command line in the next version.

More description of the modelling specifics in - https://illustrious-club-b62.notion.site/Incomplete-Cycles-Experiment-Description-154cd662ddb880088bf8d85cd6309cba
