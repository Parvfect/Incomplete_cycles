#!/bin/bash
#PBS -lwalltime=20:00:00
#PBS -lselect=1:ncpus=1:mem=60gb

# 2. Clone Incomplete cycles and run it
git clone https://github.com/Parvfect/Incomplete_cycles.git
cd Incomplete_cycles
cd v2
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
python cluster_reads.py --reads_filepath /rds/general/user/pa1123/home/clustering_data/strand_pool.fasta --output_filepath /rds/general/user/pa1123/home/clustering_data/ --hpc --min_cluster_length 70 --max_cluster_length 180
