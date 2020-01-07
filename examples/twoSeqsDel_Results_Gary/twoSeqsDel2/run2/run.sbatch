#!/bin/bash
#SBATCH -o output.o
java -jar ../../../../../bondalign/StatAlign.jar -seed=3948933 -verbose=true -plugin:structal[printRmsd] -mcmc=200k,200m,100k,300 -log=E,T,P,A,B,S,R,C,L,N 2_globins.fasta 2_globins.coor> output.o

