#include "parameters.h"
#include "pbc.h"


void pbc(double *rxyz, const double box_size){ 
    // condiciones periodicas de contorno coordenadas entre [0,L)

    for (int i = 0; i < N; i++){
        for (int j = 0; j < 3; j++){
            if (rxyz[j*N + i] <= 0.0)
                rxyz[j*N + i] += box_size;
            if (rxyz[j*N + i] > box_size)
                rxyz[j*N + i] -= box_size;
        }
    }
}
