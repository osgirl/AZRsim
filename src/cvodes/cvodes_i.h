/* Common header file for model C code function and cvodesIQRinterface.c */

/* ParamData (contains pointer to parameter values passed to integrator) */
typedef struct {
  double *parametervector;
} ParamData;

/* CVODES-IQR Tools interface related flags */
#define DOFLAG_DDT         0          /* Returns derivatives for integration */
#define DOFLAG_VARREAC     1          /* Returns variable and reaction values for given time and states */
#define DOFLAG_EVENTS      2          /* Evaluates event triggers / root finding */
#define DOFLAG_EVENTASSIGN 3          /* Determines event assignments */
