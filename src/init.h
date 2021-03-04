/*  Inicialización del sistema: posiciones en un cristal FCC y velocidades
 *  aleatorias distribuidas uniformemente entre -T0/2 y T0/2
 */
#ifndef INIT_H
#define INIT_H

void init_pos(double *rxyz, const double rho);
void init_vel(double *vxyz, double *temp, double *ekin);

#endif
