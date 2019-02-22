#!/bin/bash
#SBATCH --job-name=deapTest
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:09:00
#SBATCH --partition=debug
  
cd /localdisk/thierry/output/deap_test
echo "nodes=2, ntasks=32, cpus-per-task=1"
export OMP_NUM_THREADS=1
/usr/bin/python2 -m scoop /localdisk/thierry/tc_install/py/tests/try_parallel.py -p 1000 -g 10
# end of run.sh
