#!/bin/bash

for instance in X-n237-k14.vrp X-n242-k48.vrp X-n256-k16.vrp X-n261-k13.vrp X-n266-k58.vrp X-n275-k28.vrp X-n284-k15.vrp X-n294-k50.vrp X-n303-k21.vrp ; do

#for instance in X-n200-k36.vrp X-n214-k11.vrp X-n223-k34.vrp X-n237-k14.vrp X-n242-k48.vrp X-n256-k16.vrp X-n261-k13.vrp X-n266-k58.vrp X-n275-k28.vrp X-n284-k15.vrp X-n294-k50.vrp X-n303-k21.vrp ; do

#for instance in X-n101-k25.vrp X-n106-k14.vrp X-n110-k13.vrp X-n115-k10.vrp X-n120-k6.vrp X-n125-k30.vrp X-n129-k18.vrp X-n134-k13.vrp X-n139-k10.vrp X-n143-k7.vrp X-n148-k46.vrp X-n153-k22.vrp X-n157-k13.vrp X-n162-k11.vrp X-n167-k10.vrp X-n172-k51.vrp X-n176-k26.vrp X-n181-k23.vrp X-n186-k15.vrp X-n190-k8.vrp X-n195-k51.vrp; do

#for instance in X-n129-k18.vrp ; do

for crossover in 0 1 2 3 4 5 6 7 ; do

#for seed in 0 1 2 3 4 5 6 7 8 9 ; do

for seed in 5 6 7 8 9 ; do

./hgs ../Instances/CVRP/$instance $crossover mySolution.sol -seed $seed  -t 60

done
done
done

