/*  Paso de dinámica molecular teniendo en cuenta:
 *      + primer paso del algoritmo de Velocity Verlet
 *      + condiciones periódicas de contorno
 *      + cálculo de fuerzas LJ 12-6
 *      + segundo paso del algoritmo de Velocity Verlet
 *      + reescaleo de velocidades para controlar la temperatura del sistema
 */
#ifndef MD_STEP_H
#define MD_STEP_H

void md_step(double *rxyz, double *vxyz, double *fxyz, double *epot,
             double *ekin, double *pres, double *temp, const double rho,
             const double V, const double L);

#endif
