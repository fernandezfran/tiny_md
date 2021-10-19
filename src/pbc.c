#include "pbc.h"
#include "parameters.h"

void pbc(double *rxyz, const double box_size) {
    // condiciones periodicas de contorno coordenadas entre [0,L)

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            while (rxyz[j * N + i] <= 0.0)
                rxyz[j * N + i] += box_size;
            while (rxyz[j * N + i] > box_size)
                rxyz[j * N + i] -= box_size;
        }
    }
}
