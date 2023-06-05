# OpenMP offload of ITASSER:  
I-TASSER (Iterative Threading ASSEmbly Refinement) is a hierarchical approach to protein structure prediction and structure-based function annotation. The I-TASSER suite predicts protein structures through four main steps. These include threading template identification, iterative structure assembly simulation, model selection, and refinement, and the final step being structure-based function annotation. The structure folding and reassembling stage is conducted by replica-exchange Monte Carlo simulations. For details on I-TASSER for protein structure prediction, please refer to https://zhanggroup.org/I-TASSER/
The OpenMP offload version of I-TASSER uses OpenMP to offload replica-exchange Monte Carlo regions of the I-TASSER pipeline to the device for optimal performance. 
Please follow these steps to compile and run OpenMP offload version of ITASSER: 
1. Please ensure the NVHPC compiler  (or XL or GCC, refer to single-src branch)  is loaded: eg. module load nvhpc/(version).
3. Please locate the Makefile (Makefile_xl or Makefile_gcc from single-src branch) and make or make -j 8
4. This process should generate cas, the binarey file for the OpenMP I-TASSER source code. Please put this binary file in the /I-TASSERmod directory and run the standalone version of I-TASSER using the perl script, runI-TASSER.pl. 
5. If input files for cas already exist from a previous I-TASSER run (as in the directory BenchmarkD200setsHard60/ here), then I-TASSER could be run with the srun/jsrun line by simply putting the binary cas in the directory BenchmarkD200setsHard60/ : 
   a. srun -n1 --gpus=1 ./cas (for OLCF Wombat cluster) or 
   b. jsrun -n1 -c1 -a1 -g1 --smpiargs=off ./cas (for OLCF Summit cluster)
