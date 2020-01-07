#!/bin/bash
#SBATCH -o output.o
java -jar ../../../../bondalign/StatAlign.jar -verbose=true -plugin:structal[epsilon=0,printRmsd] -mcmc=50k,50m,10k -log=E,T,P,A,B,S,R,C,L,N 5_globins.fasta 5_globins.coor > output.o

