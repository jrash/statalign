#!/bin/bash
#SBATCH -o output.o
java -jar ../../../../bondalign/StatAlign.jar -seed=3948933 -verbose=true -plugin:structal[localSigma,epsilon=0,printRmsd] -mcmc=50k,500m,100k -log=E,T,P,A,B,S,R,C,L,N 5_globins.fasta 5_globins.coor > output.o

