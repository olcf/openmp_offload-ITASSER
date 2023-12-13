#!/bin/bash
read -r -d '' CMD_LST << EOF
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_1A
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_2A
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_3A
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_4A
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_1M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_2M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_3M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_4M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_5M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_6M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_7M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_8M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_9M
/home/maccarth/I-TASSER5.1/I-TASSER5.1/BenchmarkD200setsHard60/seq.fastasim_10M
EOF
echo "$CMD_LST"|parallel -k
