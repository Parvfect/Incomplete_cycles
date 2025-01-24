#!/bin/bash
#PBS -lwalltime=20:00:00
#PBS -lselect=1:ncpus=1:mem=20gb

module load miniforge/3
miniforge-setup
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/pa1123/home/anaconda3/envs/pytorch_gpu
git clone https://github.com/rrwick/Badread.git
pip3 install ./Badread
badread simulate --reference /rds/general/user/pa1123/home/strands_for_badread/synthesized_0.8_million.fasta --quantity 2x | gzip > /rds/general/user/pa1123/home/strands_for_badread/reads_0.8_million_more_depth_3.fastq.gz
