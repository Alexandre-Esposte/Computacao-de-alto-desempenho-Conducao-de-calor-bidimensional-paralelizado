#!/bin/bash

# Paralelo  para a flag -O0 - esse é serial

for i in 3 4 5 6 8;
do
        echo "Estou com $i processos no programa sem flags (-O0)"
        
        mpicc calormpi.c -O0 -lm
	ulimit -v unlimited
        (time mpirun -np $i ./a.out) 2> tempo_O0_Proc$i.out
        rm *.dat
	rm a.out
done

# agora para a melhor flag
for i in 1 2 3 4 5 6 8;
do
        echo "Estou na flag O1 com $i processos"
        
        mpicc calormpi.c -O1 -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mavx -mtune=native -march=native  -lm

	ulimit -v unlimited
        (time mpirun -np $i ./a.out) 2> tempo_O1_proc$i.out
	rm *.dat
	rm a.out
done

