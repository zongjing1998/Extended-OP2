#ifndef _FLUX_2D_H_
#define _FLUX_2D_H_

#define FLUXNND 1
#define FLUXROE 2
#define FLUXLLF 3
#define FLUXKIN 4
#define FLUXOSH 5

void FluxSel(
    unsigned int sel,
    const double *cvL, const double *cvL0,
    const double *cvR, const double *cvR0, const double *norm, double *flux);
void FluxLLF(const double *cvL, const double *cvL0,
	     const double *cvR, const double *cvR0, const double *norm, double *flux);
void FluxRoe(const double *qL, const double *qR,
	     double hL, double hR, const double *norm, double *diss);
#endif