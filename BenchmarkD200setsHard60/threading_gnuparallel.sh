#!/bin/bash
read -r -d '' CMD_LST << EOF
/home/maccarth/I-TASSER5.1/BenchmarkD200setsHard60/runpair.sh
EOF
echo "$CMD_LST"|parallel -k
