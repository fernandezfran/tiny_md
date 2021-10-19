#include "init.h"
#include "parameters.h"

#include <math.h>   // sqrt()
#include <stdlib.h> // rand()

void init_pos(double *rxyz, const double rho) {
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                rxyz[idx + 0] = i * a;         // x
                rxyz[N + idx + 0] = j * a;     // y
                rxyz[2 * N + idx + 0] = k * a; // z
                                               // del mismo átomo
                rxyz[idx + 1] = (i + 0.5) * a;
                rxyz[N + idx + 1] = (j + 0.5) * a;
                rxyz[2 * N + idx + 1] = k * a;

                rxyz[idx + 2] = (i + 0.5) * a;
                rxyz[N + idx + 2] = j * a;
                rxyz[2 * N + idx + 2] = (k + 0.5) * a;

                rxyz[idx + 3] = i * a;
                rxyz[N + idx + 3] = (j + 0.5) * a;
                rxyz[2 * N + idx + 3] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}

void init_vel(double *vxyz, double *temp, double *ekin) {
    // inicialización de velocidades aleatorias

    double sf, sumv2 = 0.0;
    double sumv[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            vxyz[j * N + i] = rand() / (double)RAND_MAX - 0.5;
            sumv[j] += vxyz[j * N + i];
            sumv2 += vxyz[j * N + i] * vxyz[j * N + i];
        }
    }

    for (int j = 0; j < 3; j++)
        sumv[j] /= (double)N;

    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) {     // elimina la velocidad del centro de masa
        for (int j = 0; j < 3; j++) { // y ajusta la temperatura
            vxyz[j * N + i] = sf * (vxyz[j * N + i] - sumv[j]);
        }
    }
}
