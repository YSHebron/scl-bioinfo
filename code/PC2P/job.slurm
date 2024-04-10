#!/bin/bash
#SBATCH --job-name=scl1
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=/home/yshebron/scl-bioinfo/test.out

# Activate your environment, source command is required
module load miniconda3
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate scl

# Insert your code here.
cd /home/yshebron/scl-bioinfo
echo TEST
# srun python /home/yenzy_urson_hebron/scl-bioinfo/code/PC2P/PC2P.py /home/yenzy_urson_hebron/scl-bioinfo/code/PC2P/Yeast/Integrated_CYC_ppirel.txt /home/yenzy_urson_hebron/scl-bioinfo/code/PC2P/Results/Integrated -p

# [Optional] You can run other commands here

# Sample training code format
# srun bash -c 'CUDA_VISIBLE_DEVICES=0 python train.py --workers 1 --train_data data_lmdb_release/training --valid_data data_lmdb_release/validation --select_data MJ-ST --batch_ratio 0.5-0.5 --Transformation None --FeatureExtraction RCNN --SequenceModeling None --Prediction Attn'
