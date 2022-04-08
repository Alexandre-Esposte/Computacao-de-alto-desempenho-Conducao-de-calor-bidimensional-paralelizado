#!/bin/bash


# agora para a melhor flag
for i in 1 2 3 4 5 6 7 8;
do
        echo "Estou na flag O1 com 2 processos na thread $i"
        
        mpicc -fopenmp calorhibrido.c  -O1 -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mavx -mtune=native -march=native  -lm

        export OMP_NUM_THREADS=$i
        (time mpirun -np 2 ./a.out) 2> tempo_hibrido_T$i.out
        mkdir dados_hibrido_T$i
	mv *.dat dados_hibrido_T$i
	rm a.out
done

