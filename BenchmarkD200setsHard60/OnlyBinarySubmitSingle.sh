#!/bin/bash
#BATCH --job-name="exampleFromHermesMarch26egcasFile.%j.%N.out"
#SBATCH --output="exampleFromHermesMarch26egcasFile.%j.%N.out"
#SBATCH --partition gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --gres=gpu:p100:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eamaccar@aggies.ncat.edu
#SBATCH --export=ALL
#SBATCH -t 10:00:00
#SBATCH -A nct104

module purge
# module load gnutools
export MODULEPATH="/share/apps/compute/modulefiles/applications:${MODULEPATH}"
module load parallel/20120122

module load pgi
module load cuda/8.0
module load gnu
#export PGI_ACC_TIME=1
export PATH="/share/apps/compute/parallel/bin:${PATH}"
export MODULEPATH="/share/apps/compute/modulefiles:$MODULEPATH"
#module load gnu/6.2.0
#module load gcc
#module load gnu
# /opt/packages/cuda/10.0/bin/nvprof ./casForPGImodulesWorkingMatEliSep11.o
#This job runs with 3 nodes, 24 cores per node for a total of 82 cores.
#ibrun in verbose mode will give binding detail
# time ./tasser
perl /home/maccarth/I-TASSER5.1/I-TASSERmod/runI-TASSER.pl -libdir /home/maccarth/I-TASSER5.1/ITLIB  -seqname seq.fasta -datadir   /home/maccarth/I-TASSER5.1/BenchmarkD200setsHard60  -runstyle gnuparallel -hours 10 -outdir /home/maccarth/I-TASSER5.1/BenchmarkD200setsHard60
