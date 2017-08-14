#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_simulation_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	C2,
	centr,
	V2,
	depot,
	KA,
	CL;

	V2 = par_ptr[0];
	KA = par_ptr[1];
	CL = par_ptr[2];

	depot = __zzStateVar__[0];
	centr = __zzStateVar__[1];

	C2 = centr / V2;
	__DDtStateVar__[0] = InfusionRate[0] + - KA * depot;
	__DDtStateVar__[1] = InfusionRate[1] + KA * depot - CL * C2;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_simulation_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	C2,
	centr,
	V2,
	depot,
	KA,
	CL;

	V2 = par_ptr[0];
	KA = par_ptr[1];
	CL = par_ptr[2];

	depot = __zzStateVar__[0];
	centr = __zzStateVar__[1];

	C2 = centr / V2;

	lhs[0]=C2;
}
