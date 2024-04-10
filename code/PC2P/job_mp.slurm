#!/bin/bash
#SBATCH --job-name=scl1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=32
#SBATCH --output=/home/yshebron/scl-bioinfo/test.out

# Activate your environment, source command is required
module load miniconda3
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate scl

# Insert your code here.
# TODO: Force code to run multiprocess while ray is being debugged (or while pympi is being implemented)
cd ~/scl-bioinfo
srun -u python code/PC2P/PC2P.py code/PC2P/Yeast/Integrated_CYC_ppirel.txt code/PC2P/Results/Integrated -p

# [Optional] You can run other commands here

# Sample training code format
# srun bash -c 'CUDA_VISIBLE_DEVICES=0 python train.py --workers 1 --train_data data_lmdb_release/training --valid_data data_lmdb_release/validation --select_data MJ-ST --batch_ratio 0.5-0.5 --Transformation None --FeatureExtraction RCNN --SequenceModeling None --Prediction Attn'
