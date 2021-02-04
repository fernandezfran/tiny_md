#!/bin/bash

### Configuración del trabajo

### Nombre de la tarea
#SBATCH --job-name=omp256

### Cantidad de nodos a usar
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
#SBATCH --ntasks-per-node=28
### Cores por proceso (para MPI+OpenMP)
#SBATCH --cpus-per-task=1

### Para que corra solo
#SBATCH --exclusive
#---------------------------------------------------------------------------------
make clean
make N=256 > out0256.txt
for i in {1..28}
do
    echo "threads: ${i}" >> out0256.txt
    #OMP_PROC_BIND=true OMP_NUM_THREADS=$i perf stat -r 30 ./tiny_md >> out0256.txt
    OMP_NUM_THREADS=$i perf stat -r 30 ./tiny_md >> out0256.txt
done
