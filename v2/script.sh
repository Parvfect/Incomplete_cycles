#!/bin/bash
#PBS -lwalltime=10:00:00
#PBS -lselect=1:ncpus=1:mem=100gb

# 2. Clone Incomplete cycles and run it
git clone https://github.com/Parvfect/Incomplete_cycles.git
cd Incomplete_cycles
cd v2
pip install -r requirements.txt
python synthesis_pipeline.py