#include "parameters.h"
#include "verlet.h"


void v_firststep(double *rxyz, double *vxyz, const double *fxyz){
    // actualizo posiciones y mitad de las velocidades 
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 3; j++){
            vxyz[j*N + i] += 0.5 * fxyz[j*N + i] * dt;
            rxyz[j*N + i] += vxyz[j*N + i] * dt;
        }
    }
}


void v_secondstep(double *vxyz, const double *fxyz, double *ekin, double *temp){
    // actualizo la otra mitad de las velocidades
    double sumv2 = 0.0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 3; j++){
            vxyz[j*N + i] += 0.5 * fxyz[j*N + i] * dt;
            sumv2 += vxyz[j*N + i] * vxyz[j*N + i];
        }
    }
    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
