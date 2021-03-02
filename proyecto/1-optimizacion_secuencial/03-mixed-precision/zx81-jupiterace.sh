#!/bin/bash

### Configuración del trabajo

### Nombre de la tarea
#SBATCH --job-name=md-mixed

### Cantidad de nodos a usar
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
#SBATCH --ntasks-per-node=1
### Cores por proceso (para MPI+OpenMP)
#SBATCH --cpus-per-task=1

### Para que corra solo
#SBATCH --exclusive
#--------------------------------------------------------------------------------------

### Script que se ejecuta al arrancar el trabajo
cat flaglist.txt | python3 exploracion.py
