#!/bin/bash
#SBATCH -o output.o
java -jar ../../../../structalign/StatAlign.jar -seed=-434563 -verbose=true -plugin:structal[localSigma,epsilon=0,printRmsd] -mcmc=50k,500m,100k -log=E,T,P,A,B,S,R,C,L,N 5_globins.fasta 5_globins.coor > output.o

