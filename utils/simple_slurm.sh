#!/bin/bash 
#testing slurm environment on cluster
# set the number of nodes
#SBATCH --nodes=1
# set max wallclock time
#SBATCH --time=00:00:10 
# allocate memory in MB
#SBATCH --mem=100 
# exlude a node
#SBATCH --exclude=leak59,leak58

sleep 1
echo "running on node: $HOSTNAME"
