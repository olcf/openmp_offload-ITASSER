#!/bin/bash
cd /scratch/maccarth/34469527/maccarth/ITseq.fasta
cp -f /home/maccarth/I-TASSER5.1/bin/pair99 ./pair
./pair /
sync
sleep 1
cp -f pair.3  /home/maccarth/I-TASSER5.1/BenchmarkD200setsHard60/pair3.dat
cp -f pair.1  /home/maccarth/I-TASSER5.1/BenchmarkD200setsHard60/pair1.dat
