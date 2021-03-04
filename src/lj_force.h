/*  Calculo de las fuerzas entre las particulas para un potencial interatómico
 *  de Lennard Jones 12-6
 */
#ifndef LJ_FORCE_H
#define LJ_FORCE_H

void forces(const double *rxyz, double *fxyz, double *epot, double *pres,
            const double *temp, const double rho, const double V, const double L);

#endif
