#!/bin/bash

# Paralelo  para a flag -O0 - esse é serial

for i in 1 2 3 4 5 6 7 8;
do
        echo "Estou na Thread $i no programa sem flags (-O0)"
        
        gcc -fopenmp calormp.c -O0 -lm
       	export OMP_NUM_THREADS=$i 
        (time ./a.out) 2> tempo_O0_T$i.out
        rm *.dat
	rm a.out
done

# agora para a melhor flag
for i in 1 2 3 4 5 6 7 8;
do
        echo "Estou na flag O1 na thread $i"
        
        gcc -fopenmp calormp.c  -O1 -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mavx -mtune=native -march=native  -lm

        export OMP_NUM_THREADS=$i
        (time ./a.out) 2> tempo_O1flags_T$i.out
        mkdir flagsO1_T$i
	mv *.dat flagsO1_T$i
	rm a.out
done

