#include "rescale.h"
#include "parameters.h"

#include <math.h> // sqrt()

void rescale(double *vxyz, const double *temp) {
    // reescaleo de velocidades
    double sf = sqrt(T0 / *temp);
    for (int i = 0; i < 3 * N; i++) {
        vxyz[i] *= sf;
    }
}
