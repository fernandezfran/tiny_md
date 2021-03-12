#include "parameters.h"
#include "lj_force.h"

#include <omp.h>    // #pragma omp parallel


static double minimum_image(double cordi, const double box_size){ 
    // imagen más cercana
    if (cordi <= -0.5*box_size)
        cordi += box_size;
    else if (cordi > 0.5*box_size)
        cordi -= box_size;
    return cordi;
}


void forces(const double *rxyz, double *fxyz, double *epot, double *pres,
            const double *temp, const double rho, const double V, const double L){
    // calcula las fuerzas LJ (12-6)
   
    for (int k = 0; k < 3*N; k++){
        fxyz[k] = 0.0;
    }
    double pres_vir = 0.0;
    double rcut2 = rcut*rcut;
    double pot = 0.0;
    double rij2;
    double xi[3], xj[3], rij[3];

    #pragma omp parallel                          \
        default(shared) private(rij2, xi, xj, rij) \
        reduction(+:pot,pres_vir)
    {   // manual schedule
        int threads = omp_get_num_threads();    // cantidad de hilos
        int tid     = omp_get_thread_num();     // id del hilo en el que estoy

        int work = (N + threads - 1) / threads; // divido el trabajo segun n_hilos
        int from = tid * work;                  // desde donde empieza el hilo id
        int to   = (tid + 1) * work;            // hasta donde va el hilo id
        if (to > N) to = N;                     // para no pasarme de N

        for (int i = from; i < to; i++){
            
            for (int k = 0; k < 3; k++) xi[k] = rxyz[k*N + i];
        
            for (int j = 0; j < N; j++){

                if (i == j) continue;

                for (int k = 0; k < 3; k++) xj[k] = rxyz[k*N + j];

                // distancia mínima entre r_i y r_j
                for (int k = 0; k < 3; k++) rij[k] = xi[k] - xj[k];
                for (int k = 0; k < 3; k++) rij[k] = minimum_image(rij[k], L);
                
                rij2 = 0.0;
                for (int k = 0; k < 3; k++) rij2 += rij[k]*rij[k];
                

                if (rij2 <= rcut2){
                    double r2inv = 1.0/rij2;
                    double r6inv = r2inv*r2inv*r2inv;

                    double fr = 24.0*r2inv*r6inv*(2.0*r6inv - 1.0);

                    for (int k = 0; k < 3; k++){
                        fxyz[k*N + i] += fr*rij[k];
                    }

                    pot += 4.0*r6inv*(r6inv - 1.0) - ecut;
                    pres_vir += fr*rij2;
                }
            }
        }
    }
    *epot = 0.5*pot;
    pres_vir /= (3.0 * V);
    *pres = *temp * rho + 0.5*pres_vir;
}
