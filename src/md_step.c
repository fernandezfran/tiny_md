#include "parameters.h"
#include "lj_force.h"
#include "pbc.h"
#include "verlet.h"
#include "rescale.h"
#include "md_step.h"


void md_step(double *rxyz, double *vxyz, double *fxyz, double *epot,
             double *ekin, double *pres, double *temp, const double rho,
             const double V, const double L){
    // un paso de dinámica molecular
    v_firststep(rxyz, vxyz, fxyz);
    pbc(rxyz, L);
    forces(rxyz, fxyz, epot, pres, temp, rho, V, L);
    v_secondstep(vxyz, fxyz, ekin, temp);
    rescale(vxyz, temp);
}
