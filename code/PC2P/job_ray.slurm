#!/bin/bash
#SBATCH --job-name=scl1
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output=/home/yshebron/scl-bioinfo/test.out

# Activate your environment, source command is required
module load miniconda3
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate scl-ray

# Getting the node names
nodes=$(scontrol show hostnames "$SLURM_JOB_NODELIST")
nodes_array=($nodes)

head_node=${nodes_array[0]}
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" hostname --ip-address)

# if we detect a space character in the head node IP, we'll
# convert it to an ipv4 address. This step is optional.
if [[ "$head_node_ip" == *" "* ]]; then
IFS=' ' read -ra ADDR <<<"$head_node_ip"
if [[ ${#ADDR[0]} -gt 16 ]]; then
  head_node_ip=${ADDR[1]}
else
  head_node_ip=${ADDR[0]}
fi
echo "IPV6 address detected. We split the IPV4 address as $head_node_ip"
fi

port=6379
ip_head=$head_node_ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "Starting HEAD at $head_node"
srun --nodes=1 --ntasks=1 -w "$head_node" \
    ray start --head --node-ip-address="$head_node_ip" --port=$port \
    --num-cpus "${SLURM_CPUS_PER_TASK}" --block &

# optional, though may be useful in certain versions of Ray < 1.0.
sleep 10

# number of nodes other than the head node
worker_num=$((SLURM_JOB_NUM_NODES - 1))

for ((i = 1; i <= worker_num; i++)); do
    node_i=${nodes_array[$i]}
    echo "Starting WORKER $i at $node_i"
    srun --nodes=1 --ntasks=1 -w "$node_i" \
        ray start --address "$ip_head" \
        --num-cpus "${SLURM_CPUS_PER_TASK}" --num-gpus "${SLURM_GPUS_PER_TASK}" --block &
    sleep 5
done

# Ensure that when our program starts, it will connect to the active ray cluster
sleep 10

# Insert your code here.
cd ~/scl-bioinfo
# python -u code/PC2P/PC2P.py code/PC2P/Yeast/FilteredPPINs/Collins_SGD_perprotein_weighted.txt code/PC2P/Results/Dummy -p
python -u code/PC2P/PC2P.py code/PC2P/Yeast/FilteredPPINs/Collins_CYC_perprotein_weighted.txt code/PC2P/Results/Dummy -p
# python -u code/PC2P/PC2P.py code/PC2P/Yeast/Integrated_CYC_ppirel.txt code/PC2P/Results/Integrated -p

# [Optional] You can run other commands here

# Sample training code format
# srun bash -c 'CUDA_VISIBLE_DEVICES=0 python train.py --workers 1 --train_data data_lmdb_release/training --valid_data data_lmdb_release/validation --select_data MJ-ST --batch_ratio 0.5-0.5 --Transformation None --FeatureExtraction RCNN --SequenceModeling None --Prediction Attn'
