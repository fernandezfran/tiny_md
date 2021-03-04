#include "parameters.h"
#include "init.h"
#include "lj_force.h"
#include "md_step.h"

#include <stdio.h>  // printf(), fprintf()
#include <math.h>   // sqrt(), cbrt(), pow(), M_PI
#include <stdlib.h> // rand()
#include <time.h>   // time(NULL)
#include <omp.h>    // omp_get_wtime()

int main(){ // tiny_md : ecuación de estado
    FILE *file_xyz, *file_thermo;
    file_xyz = fopen("trajectory.xyz", "w");
    file_thermo = fopen("thermo.log", "w");
    double start = 0.0, elapsed = 0.0;
    double Ekin, Epot, Temp, Pres;       // variables macroscopicas
    double Rho, V, box_size, tail, Etail, Ptail;
    double *rxyz, *vxyz, *fxyz;          // variables microscopicas

    rxyz = (double *)malloc(3*N*sizeof(double));
    vxyz = (double *)malloc(3*N*sizeof(double));
    fxyz = (double *)malloc(3*N*sizeof(double));

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", teq);
    printf("# Pasos de medición:         %d\n", trun - teq);
    printf("# (mediciones cada %d pasos)\n", tmes);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    fprintf(file_thermo, "# t Temp Pres Epot Etot\n");

    srand(SEED);
    double t = 0.0, sf;
    double Rhob;
    Rho = Rhoi;
    init_pos(rxyz, Rho);
    start = omp_get_wtime();
    for (int m = 0; m < 9; m++){
        Rhob   = Rho;
        Rho    = Rhoi - 0.1*(double)m;
        V = (double)N / Rho;
        box_size = cbrt(V);
        tail   = 16.0*M_PI*Rho*((2.0/3.0)*pow(rcut,-9) - pow(rcut,-3))/3.0;
        Etail  = tail*(double)N;
        Ptail  = tail*Rho;

        int i = 0;
        sf = cbrt(Rhob/Rho);
        for (int k = 0; k < 3*N; k++){ // reescaleo posiciones a nueva densidad
            rxyz[k] *=sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, V, box_size);
    
        for (i = 1; i < teq; i++){ // loop de equilibracion
            md_step(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, V, box_size);
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0; 
        for (i = teq; i < trun; i++){ // loop de medicion

            md_step(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, V, box_size);
            
            if (i % tmes == 0){
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;
                
                fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot+Ekin);
                fprintf(file_xyz, "%d\n\n", N);
                for (int k = 0; k < N; k++){
                    fprintf(file_xyz, "Ar %e %e %e\n", rxyz[k], rxyz[N + k], rxyz[2*N + k]);
                }
            }

            t += dt;

        }
        printf("%f\t%f\t%f\t%f\n", Rho, V, epotm/(double)mes, presm/(double)mes);
    }

    elapsed = omp_get_wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t*1.6);
    printf("# ns/day = %f\n", (1.6e-6*t)/elapsed*86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}
