#!/bin/bash
#SBATCH -o output.o
java -jar ../../../../../bondalign/StatAlign.jar -verbose=true -plugin:structal[printRmsd] -mcmc=20m,200m,100k,300 -log=E,T,P,A,B,S,R,C,L,N 2_globins.fasta 2_globins.coor > output.o

