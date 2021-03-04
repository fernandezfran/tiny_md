/*  Algoritmo de evolución temporal Velocity Verlet, dividido en dos funciones
 *  correspondientes al primer paso y al segundo
 */
#ifndef VERLET_H
#define VERLET_H

void v_firststep(double *rxyz, double *vxyz, const double *fxyz);
void v_secondstep(double *vxyz, const double *fxyz, double *ekin, double *temp);

#endif
