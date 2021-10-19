#include "tiny_md.h"

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

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        vxyz[i] = rand() / (double)RAND_MAX - 0.5;
        sumvx += vxyz[i];
        sumv2 += vxyz[i] * vxyz[i];

        vxyz[N + i] = rand() / (double)RAND_MAX - 0.5;
        sumvy += vxyz[N + i];
        sumv2 += vxyz[N + i] * vxyz[N + i];

        vxyz[2 * N + i] = rand() / (double)RAND_MAX - 0.5;
        sumvz += vxyz[2 * N + i];
        sumv2 += vxyz[2 * N + i] * vxyz[2 * N + i];
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
                                  // y ajusta la temperatura
        vxyz[i] = (vxyz[i] - sumvx) * sf;
        vxyz[N + i] = (vxyz[N + i] - sumvy) * sf;
        vxyz[2 * N + i] = (vxyz[2 * N + i] - sumvz) * sf;
    }
}

static double pbc(double cordi, const double cell_length) {
    // condiciones periodicas de contorno coordenadas entre [-L/2,L/2)
    // e imagen más cercana

    if (cordi <= -0.5 * cell_length)
        cordi += cell_length;
    else if (cordi > 0.5 * cell_length)
        cordi -= cell_length;
    return cordi;
}

void forces(const double *rxyz, double *fxyz, double *epot, double *pres,
            const double *temp, const double rho, const double V,
            const double L) {
    // calcula las fuerzas LJ (12-6)

    for (int k = 0; k < 3 * N; k++) {
        fxyz[k] = 0.0;
    }
    double pres_vir = 0.0;
    double rcut2 = rcut * rcut;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        double xi = rxyz[i];
        double yi = rxyz[N + i];
        double zi = rxyz[2 * N + i];

        for (int j = i + 1; j < N; j++) {

            double xj = rxyz[j];
            double yj = rxyz[N + j];
            double zj = rxyz[2 * N + j];

            // distancia mínima entre r_i y r_j
            double rx = xi - xj;
            rx = pbc(rx, L);
            double ry = yi - yj;
            ry = pbc(ry, L);
            double rz = zi - zj;
            rz = pbc(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz;

            if (rij2 <= rcut2) {
                double r2inv = 1.0 / rij2;
                double r6inv = r2inv * r2inv * r2inv;

                double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                fxyz[i] += fr * rx;
                fxyz[N + i] += fr * ry;
                fxyz[2 * N + i] += fr * rz;

                fxyz[j] -= fr * rx;
                fxyz[N + j] -= fr * ry;
                fxyz[2 * N + j] -= fr * rz;

                *epot += 4.0 * r6inv * (r6inv - 1.0) - ecut;
                pres_vir += fr * rij2;
            }
        }
    }
    pres_vir /= (3.0 * V);
    *pres = *temp * rho + pres_vir;
}

void velocity_verlet(double *rxyz, double *vxyz, double *fxyz, double *epot,
                     double *ekin, double *pres, double *temp, const double rho,
                     const double V, const double L) {

    for (int i = 0; i < N; i++) { // actualizo posiciones
        vxyz[i] += 0.5 * fxyz[i] * dt;
        vxyz[N + i] += 0.5 * fxyz[N + i] * dt;
        vxyz[2 * N + i] += 0.5 * fxyz[2 * N + i] * dt;

        rxyz[i] += vxyz[i] * dt;
        rxyz[N + i] += vxyz[N + i] * dt;
        rxyz[2 * N + i] += vxyz[2 * N + i] * dt;

        rxyz[i] = pbc(rxyz[i], L);
        rxyz[N + i] = pbc(rxyz[N + i], L);
        rxyz[2 * N + i] = pbc(rxyz[2 * N + i], L);
    }

    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        vxyz[i] += 0.5 * fxyz[i] * dt;
        vxyz[N + i] += 0.5 * fxyz[N + i] * dt;
        vxyz[2 * N + i] += 0.5 * fxyz[2 * N + i] * dt;

        sumv2 += vxyz[i] * vxyz[i] + vxyz[N + i] * vxyz[N + i] +
                 vxyz[2 * N + i] * vxyz[2 * N + i];
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}

int main() {
    FILE *file_xyz, *file_thermo;
    file_xyz = fopen("trajectory.xyz", "w");
    file_thermo = fopen("thermo.log", "w");
    double start = 0.0, elapsed = 0.0;
    double Ekin, Epot, Temp, Pres; // variables macroscopicas
    double Rho, cell_V, cell_L, tail, Etail, Ptail;
    double *rxyz, *vxyz, *fxyz; // variables microscopicas

    rxyz = (double *)malloc(3 * N * sizeof(double));
    vxyz = (double *)malloc(3 * N * sizeof(double));
    fxyz = (double *)malloc(3 * N * sizeof(double));

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
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = Rhoi - 0.1 * (double)m;
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0 * M_PI * Rho *
               ((2.0 / 3.0) * pow(rcut, -9) - pow(rcut, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < 3 * N;
             k++) { // reescaleo posiciones a nueva densidad
            rxyz[k] *= sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < teq; i++) { // loop de equilibracion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho,
                            cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0;
        for (i = teq; i < trun; i++) { // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho,
                            cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }

            if (i % tmes == 0) {
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;

                fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot,
                        Epot + Ekin);
                fprintf(file_xyz, "%d\n\n", N);
                for (int k = 0; k < N; k++) {
                    fprintf(file_xyz, "Ar %e %e %e\n", rxyz[k], rxyz[N + k],
                            rxyz[2 * N + k]);
                }
            }

            t += dt;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (double)mes,
               presm / (double)mes);
    }

    elapsed = omp_get_wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}
