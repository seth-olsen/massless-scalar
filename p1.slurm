#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8
#SBATCH --time=1:00:00

./p1 -outfile d210 -ic_Dsq 230.0

./p1 -outfile d215 -ic_Dsq 235.0

./p1 -outfile d210 -ic_Dsq 240.0

./p1 -outfile d225 -ic_Dsq 245.0

