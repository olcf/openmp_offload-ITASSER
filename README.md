OpenMP offload of ITASSER:
I-TASSER (Iterative Threading ASSEmbly Refinement) is a hierarchical approach to protein structure prediction and structure-based function annotation. The I-TASSER suite predicts protein structures through four main steps. These include threading template identification, iterative structure assembly simulation, model selection, and refinement, and the final step being structure-based function annotation. The structure folding and reassembling stage is conducted by replica-exchange Monte Carlo simulations. For details on I-TASSER for protein structure prediction, please refer to https://zhanggroup.org/I-TASSER/ The OpenMP offload version of I-TASSER uses OpenMP to offload replica-exchange Monte Carlo regions of the I-TASSER pipeline to the device for optimal performance. Please follow these steps to compile and run OpenMP offload version of ITASSER:

Please ensure any of these 3 compilers are loaded: NVHPC, GCC and XL (eg. module load nvhpc)

Please locate the Makefile for the loaded compilers respective amd make

This process should generate cas, the binarey file for the OpenMP I-TASSER source code.

If input files for cas already exist from a previous I-TASSER run (as in the directory BenchmarkD200setsHard60/ here), then I-TASSER could be run as below after copying cas to the BenchmarkD200setsHard60/ directory

a. eg: srun -n1 --gpus=1 ./cas (for OLCF Wombat cluster) or

b. eg: jsrun -n1 -c1 -a1 -g1 --smpiargs=off ./cas (for OLCF Summit cluster)
