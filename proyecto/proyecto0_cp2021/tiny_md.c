#include "parameters.h"
#include "core.h"


int main(){
    FILE *file_xyz, *file_thermo;
    file_xyz = fopen("trajectory.xyz", "w");
    file_thermo = fopen("thermo.log", "w");
    double start = 0.0, elapsed = 0.0;
    double Ekin, Epot, Temp, Pres;       // variables macroscopicas
    double Rho, cell_V, cell_L, tail, Etail, Ptail;
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
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail   = 16.0*M_PI*Rho*((2.0/3.0)*pow(rcut,-9) - pow(rcut,-3))/3.0;
        Etail  = tail*(double)N;
        Ptail  = tail*Rho;

        int i = 0;
        sf = cbrt(Rhob/Rho);
        for (int k = 0; k < 3*N; k++){ // reescaleo posiciones a nueva densidad
            rxyz[k] *=sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);
    
        for (i = 1; i < teq; i++){ // loop de equilibracion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);
             
            sf = sqrt(T0/Temp);
            for (int k = 0; k < 3*N; k++){ // reescaleo de velocidades
                vxyz[k] *=sf;
            }
            
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0; 
        for (i = teq; i < trun; i++){ // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);
            
            sf = sqrt(T0/Temp);
            for (int k = 0; k < 3*N; k++){ // reescaleo de velocidades
                vxyz[k] *=sf;
            }

            if (i % tmes == 0){
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;
                
                fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot+Ekin);
                fprintf(file_xyz, "%d\n\n", N);
                for (int k = 0; k < 3*N; k+=3){
                    fprintf(file_xyz, "Ar %e %e %e\n", rxyz[k + 0], rxyz[k + 1], rxyz[k + 2]);
                }
            }

            t += dt;

        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm/(double)mes, presm/(double)mes);
    }

    elapsed = omp_get_wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t*1.6);
    printf("# ns/day = %f\n", (1.6e-6*t)/elapsed*86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}
