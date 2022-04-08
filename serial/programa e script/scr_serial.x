#!/bin/bash

# Para o serial primeiro fazendo o profiling aqui

echo "Estou na flag O0"
mkdir O0
gcc -pg -g calorserial.c -O0 -lm
mv a.out O0
(time ./O0/a.out) 2> tempoO0.out
mv *.dat O0

# agora para as outras flags Ox
for i in 3 2 1;
do
	echo "Estou na flag O$i"
	mkdir O$i
	gcc calorserial.c -O$i -lm
	mv a.out O$i
	(time ./O$i/a.out) 2> tempoO$i.out
	mv *.dat O$i
done

# agora para as outras flags
for i in 3 2 1 0;
do
        echo "Estou na flag O$i com mais flags"
        mkdir flagsO$i
        gcc calorserial.c -O$i -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mavx -mtune=native -march=native  -lm
        mv a.out flagsO$i
        (time ./flagsO$i/a.out) 2> tempo_flagsO$i.out
	mv *.dat flagsO$i
done

