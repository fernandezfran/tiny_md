/*
 * Tiny Molecular Dynamics
 *
 * Unidades: Lennard-Jones
 *
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>  // printf(), fprintf()
#include <math.h>   // sqrt(), cbrt(), pow(), M_PI
#include <stdlib.h> // rand()
#include <time.h>   // time(NULL)
#include <omp.h>    // omp_get_wtime()

#ifndef N // número de particulas (debe ser un 4m^3 para el cristal inicial)
#define N 864
#endif

#ifndef SEED // rand SEED para las velocidades
#define SEED (time(NULL))
#endif

#ifndef T0 // isoterma
#define T0 2.0
#endif

#ifndef Rhoi // densidad inicial
#define Rhoi 1.2
#endif

#ifndef rcut // radio de corte
#define rcut 2.5
#endif

#ifndef dt // paso temporal ~ 1.6 fs para el Ar
#define dt 0.005
#endif

#ifndef teq // pasos de equilibracion
#define teq 500
#endif

#ifndef trun // trun - teq: pasos de medicion
#define trun 2000
#endif

#ifndef tmes // cada cuantos pasos se mide
#define tmes 10
#endif

#define ecut (4.0*(pow(rcut,-12) - pow(rcut,-6)))

#endif
