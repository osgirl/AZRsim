#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#include "cvodes/cvodes_i.h"
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_dense.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"

#define Ith(v,i) NV_Ith_S(v,i-1)
typedef int model_func_type(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
static int  f(double time, N_Vector u, N_Vector udot, void *f_data);
static int  g(double time, N_Vector y, double *gout, void *g_data);
static void addVec2Mat(double *matrix, double *rowvector, int row, int nrows, int ncols);
static void freeMem(int freeCVODE);
static void errorMsg(char *text);
model_func_type *model_func;
double *timesimvector;
double *initialconditions;
double *parametervector;
double *resultXdot;
double *statevalues;
double *variablevalues;
double *variablevec;
double *reactionvalues;
double *reactionvec;
double *eventdataold;
int    *eventvec;
N_Vector u;
void    *cvode_mem;

SEXP cvodes_simdosingtable(SEXP r_model_func_ptr,
                      SEXP r_model_elements_nr,
                      SEXP r_timevector,
                      SEXP r_ICs,
                      SEXP r_paramvalues,
                      SEXP r_opt_method_stiff,
                      SEXP r_opt_abstol,
                      SEXP r_opt_reltol,
                      SEXP r_opt_minstep,
                      SEXP r_opt_maxstep,
                      SEXP r_opt_initstep,
                      SEXP r_opt_maxnumsteps,
                      SEXP r_opt_maxerrtestfails,
                      SEXP r_opt_maxorder_stiff,
                      SEXP r_opt_maxorder_nonstiff,
                      SEXP r_opt_maxconvfails,
                      SEXP r_opt_maxnonlineariter,
                      SEXP r_opt_xdotcalc,
                      SEXP r_verbose,
                      SEXP r_dosing_effect_start_times,
                      SEXP r_simtime_lengths
) {
  int     NRSTATES              = INTEGER(r_model_elements_nr)[0];
  int     NRPARAMETERS          = INTEGER(r_model_elements_nr)[1];
  int     NRVARIABLES           = INTEGER(r_model_elements_nr)[2];
  int     NRREACTIONS           = INTEGER(r_model_elements_nr)[3];
  int     NREVENTS              = INTEGER(r_model_elements_nr)[4];
  int    opt_method_stiff      = INTEGER(r_opt_method_stiff)[0];
  double opt_abstol            = REAL(r_opt_abstol)[0];
  double opt_reltol            = REAL(r_opt_reltol)[0];
  double opt_minstep           = REAL(r_opt_minstep)[0];
  double opt_maxstep           = REAL(r_opt_maxstep)[0];
  double opt_initstep          = REAL(r_opt_initstep)[0];
  int    opt_maxnumsteps       = INTEGER(r_opt_maxnumsteps)[0];
  int    opt_maxerrtestfails   = INTEGER(r_opt_maxerrtestfails)[0];
  int    opt_maxorder_stiff    = INTEGER(r_opt_maxorder_stiff)[0];
  int    opt_maxorder_nonstiff = INTEGER(r_opt_maxorder_nonstiff)[0];
  int    opt_maxconvfails      = INTEGER(r_opt_maxconvfails)[0];
  int    opt_maxnonlineariter  = INTEGER(r_opt_maxnonlineariter)[0];
  int    opt_xdotcalc          = INTEGER(r_opt_xdotcalc)[0];
  int    verbose               = INTEGER(r_verbose)[0];
  int     NRTIMESTEPS;
  int initialSimtime = 0;
  int initialOutputRow = 0;
  int parameterListIterator = 0;
  int start_times=INTEGER(r_dosing_effect_start_times)[0];
  int numberOfRowsOfOutput = LENGTH(r_timevector);
  int output_nrow=REAL(r_simtime_lengths)[0];
  SEXP output;
  PROTECT (output = allocMatrix(REALSXP, numberOfRowsOfOutput, 1+INTEGER(r_model_elements_nr)[0]+INTEGER(r_model_elements_nr)[2]+INTEGER(r_model_elements_nr)[3]));
  if (verbose) {
    Rprintf("--------------------------------------\n");
    Rprintf("AZR Tools CVODES Integrator Interface!\n");
    Rprintf("--------------------------------------\n");
    Rprintf("\n");
    Rprintf("Model Information:\n");
    Rprintf("------------------\n");
    Rprintf("NRSTATES:                %d\n",NRSTATES);
    Rprintf("NRPARAMETERS:            %d\n",NRPARAMETERS);
    Rprintf("NRVARIABLES:             %d\n",NRVARIABLES);
    Rprintf("NRREACTIONS:             %d\n",NRREACTIONS);
    Rprintf("NREVENTS:                %d\n",NREVENTS);
    Rprintf("\n");
    Rprintf("General settings:\n");
    Rprintf("-----------------------\n");
    Rprintf("opt_xdotcalc:            %d\n",opt_xdotcalc);
    Rprintf("verbose:                 %d\n",opt_xdotcalc);
    Rprintf("\n");
    Rprintf("CVODES options:\n");
    Rprintf("---------------\n");
    Rprintf("opt_method_stiff:        %d\n",opt_method_stiff);
    Rprintf("opt_abstol:              %f\n",opt_abstol);
    Rprintf("opt_reltol:              %f\n",opt_reltol);
    Rprintf("opt_minstep:             %f\n",opt_minstep);
    Rprintf("opt_maxstep:             %f\n",opt_maxstep);
    Rprintf("opt_initstep:            %f\n",opt_initstep);
    Rprintf("opt_maxnumsteps:         %d\n",opt_maxnumsteps);
    Rprintf("opt_maxerrtestfails:     %d\n",opt_maxerrtestfails);
    Rprintf("opt_maxorder_stiff:      %d\n",opt_maxorder_stiff);
    Rprintf("opt_maxorder_nonstiff:   %d\n",opt_maxorder_nonstiff);
    Rprintf("opt_maxconvfails:        %d\n",opt_maxconvfails);
    Rprintf("opt_maxnonlineariter:    %d\n",opt_maxnonlineariter);
    Rprintf("\n");
  }
  for(int p=0;p<start_times;p++)
  {
    timesimvector                 = NULL;
    initialconditions             = NULL;
    parametervector               = NULL;
    resultXdot                    = NULL;
    statevalues                   = NULL;
    variablevalues                = NULL;
    variablevec                   = NULL;
    reactionvalues                = NULL;
    reactionvec                   = NULL;
    eventdataold                  = NULL;
    eventvec                      = NULL;
    u                             = NULL;
    cvode_mem                     = NULL;
    model_func                    = (model_func_type *) R_ExternalPtrAddr(r_model_func_ptr);
    NRTIMESTEPS           = REAL(r_simtime_lengths)[p];
    timesimvector                 = (double *) malloc(NRTIMESTEPS*sizeof(double));
    initialconditions             = (double *) malloc(NRSTATES*sizeof(double));
    parametervector               = (double *) malloc(NRPARAMETERS*sizeof(double));
    for(int k=0,tt=initialSimtime; k<REAL(r_simtime_lengths)[p]; k++,tt++) timesimvector[k] = REAL(r_timevector)[tt];
    initialSimtime += NRTIMESTEPS;
    for(int k=0; k<NRSTATES; k++) initialconditions[k] = REAL(r_ICs)[k];
    for(int k=0,tt=parameterListIterator; k<NRPARAMETERS; k++,tt++) parametervector[k] = REAL(r_paramvalues)[tt];
    parameterListIterator += NRPARAMETERS;
    ParamData modeldata;
    modeldata.parametervector = parametervector;
    if (opt_xdotcalc) {
      SEXP outputIntegrator;
      resultXdot = (double *) malloc(NRSTATES*sizeof(double));
      model_func(timesimvector[0], initialconditions, resultXdot, &modeldata, DOFLAG_DDT, NULL, NULL, NULL, NULL);
      PROTECT (outputIntegrator = allocVector(REALSXP, NRSTATES));
      for (int k=0; k<NRSTATES; k++) REAL(outputIntegrator)[k] = resultXdot[k];
      UNPROTECT(1);
      freeMem(0);
      return(outputIntegrator);
    }
    statevalues       = (double *) malloc(NRTIMESTEPS*NRSTATES*sizeof(double));
    if (NRVARIABLES > 0) {
      variablevalues  = (double *) malloc(NRTIMESTEPS*NRVARIABLES*sizeof(double));
      variablevec     = (double *) malloc(NRVARIABLES*sizeof(double));
    }
    if (NRREACTIONS > 0) {
      reactionvalues  = (double *) malloc(NRTIMESTEPS*NRREACTIONS*sizeof(double));
      reactionvec     = (double *) malloc(NRREACTIONS*sizeof(double));
    }
    if (NREVENTS > 0) {
      eventdataold      = (double *) malloc(NREVENTS*sizeof(double));
      eventvec          = (int *)    malloc(NREVENTS*sizeof(int));
    }
    u = N_VMake_Serial(NRSTATES,initialconditions);
    if (opt_method_stiff == 1) {
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      CVodeSetMaxOrd(cvode_mem, opt_maxorder_stiff);
    }
    else {
      cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
      CVodeSetMaxOrd(cvode_mem, opt_maxorder_nonstiff);
    }
    CVodeMalloc(cvode_mem, f, timesimvector[0], u, CV_SS, opt_reltol, &opt_abstol);
    if (NREVENTS > 0) {
      CVodeRootInit(cvode_mem, NREVENTS, g, &modeldata);
    }
    CVodeSetFdata(cvode_mem, &modeldata);
    if (opt_minstep > 0) CVodeSetMinStep(cvode_mem, opt_minstep);
    if (opt_maxstep > 0) CVodeSetMaxStep(cvode_mem, opt_maxstep);
    if (opt_maxnumsteps > 0) CVodeSetMaxNumSteps(cvode_mem, opt_maxnumsteps);
    CVodeSetMaxErrTestFails(cvode_mem, opt_maxerrtestfails);
    CVodeSetMaxConvFails(cvode_mem, opt_maxconvfails);
    CVodeSetInitStep(cvode_mem, opt_initstep);
    CVodeSetMaxNonlinIters(cvode_mem, opt_maxnonlineariter);
    CVDense(cvode_mem,NRSTATES);
    addVec2Mat(statevalues,initialconditions,0,NRTIMESTEPS,NRSTATES);
    model_func(timesimvector[0], initialconditions, NULL, &modeldata, DOFLAG_VARREAC, variablevec, reactionvec,NULL,NULL);
    if (NRVARIABLES > 0) addVec2Mat(variablevalues,variablevec,0,NRTIMESTEPS,NRVARIABLES);
    if (NRREACTIONS > 0) addVec2Mat(reactionvalues,reactionvec,0,NRTIMESTEPS,NRREACTIONS);
    if (NREVENTS > 0) {
      model_func(timesimvector[0], initialconditions, NULL, &modeldata, DOFLAG_EVENTS, NULL, NULL, eventdataold, NULL);
    }
    int     ktime           = 1;
    double  tendstep        = timesimvector[ktime];
    double  treturn         = 0;
    int     flag_CVODE      = 0;
    char    stringbuffer[256];
    double *statevec        = NULL;
    double  eventcorrecttest;
    while(1) {
      flag_CVODE = CVode(cvode_mem, tendstep, u, &treturn, CV_NORMAL);
      if (flag_CVODE < 0) {
        if (flag_CVODE == CV_TOO_MUCH_WORK) errorMsg("CVODE Error: CV_TOO_MUCH_WORK");
        else if (flag_CVODE == CV_TOO_MUCH_ACC) errorMsg("CVODE Error: CV_TOO_MUCH_ACC");
        else if (flag_CVODE == CV_ERR_FAILURE || flag_CVODE == CV_CONV_FAILURE) errorMsg("CVODE Error: CV_ERR_FAILURE");
        else {
          sprintf(stringbuffer, "CVODE Error Flag: %d",flag_CVODE);
          errorMsg(stringbuffer);
        }
      }
      statevec = NV_DATA_S(u);
      if (flag_CVODE == CV_ROOT_RETURN) {
        CVodeGetRootInfo(cvode_mem, eventvec);
        eventcorrecttest = 0;
        model_func(treturn, statevec, &eventcorrecttest, &modeldata, DOFLAG_EVENTASSIGN, NULL, NULL, eventdataold, eventvec);
        flag_CVODE = CVodeReInit(cvode_mem, f, treturn, u, CV_SS, opt_reltol, &opt_abstol);
      }
      if (tendstep == treturn) {
        addVec2Mat(statevalues,statevec,ktime,NRTIMESTEPS,NRSTATES);
        model_func(tendstep, statevec, NULL, &modeldata, DOFLAG_VARREAC, variablevec, reactionvec, NULL, NULL);
        if (NRVARIABLES > 0) addVec2Mat(variablevalues,variablevec,ktime,NRTIMESTEPS,NRVARIABLES);
        if (NRREACTIONS > 0) addVec2Mat(reactionvalues,reactionvec,ktime,NRTIMESTEPS,NRREACTIONS);
        ktime = ktime+1;
        if (ktime < NRTIMESTEPS)
          tendstep = timesimvector[ktime];
        else
          break;
      }
      if (NREVENTS > 0) {
        model_func(treturn, statevec, NULL, &modeldata, DOFLAG_EVENTS, NULL, NULL, eventdataold, NULL);
      }
    }
  for (int krow=0,output_krow=initialOutputRow; krow<NRTIMESTEPS; krow++,output_krow++) {
    REAL(output)[output_krow] = timesimvector[krow];
    for (int kcol=0; kcol<NRSTATES; kcol++) {
      REAL(output)[output_krow + numberOfRowsOfOutput*(kcol+1)] = statevalues[krow + NRTIMESTEPS*kcol];
    }
    for (int kcol=0; kcol<NRVARIABLES; kcol++) {
      REAL(output)[output_krow + numberOfRowsOfOutput*(kcol+NRSTATES+1)] = variablevalues[krow + NRTIMESTEPS*kcol];
    }
    for (int kcol=0; kcol<NRREACTIONS; kcol++) {
      REAL(output)[output_krow + numberOfRowsOfOutput*(kcol+NRSTATES+NRVARIABLES+1)] = reactionvalues[krow + NRTIMESTEPS*kcol];
    }
  }
    output_nrow = initialOutputRow+NRTIMESTEPS-1;
    initialOutputRow += NRTIMESTEPS;
    int temp=0;
      for(int tt=1;tt<NRSTATES+1;tt++)
      {
        REAL(r_ICs)[temp]= REAL(output)[output_nrow + (numberOfRowsOfOutput*(tt))];
        temp++;
      }
    freeMem(1);
  }
  UNPROTECT(1);
  return(output);
}

static int f(double time, N_Vector u, N_Vector udot, void *f_data)
{
  double *statevec, *DDTvector;
  ParamData *paramdataPtr;
  paramdataPtr = (ParamData*) f_data;
  statevec = NV_DATA_S(u);
  DDTvector = NV_DATA_S(udot);
  model_func(time, statevec, DDTvector, paramdataPtr, DOFLAG_DDT, NULL, NULL, NULL, NULL);
  return(0);
}


static int g(double time, N_Vector y, double *gout, void *g_data)
{
  double *statevec;
  ParamData *paramdataPtr;
  paramdataPtr = (ParamData*) g_data;
  statevec = NV_DATA_S(y);
  model_func(time, statevec, NULL, paramdataPtr, DOFLAG_EVENTS, NULL, NULL, gout, NULL);
  return(0);
}


static void addVec2Mat(double *matrix, double *rowvector, int row, int nrrows, int nrcols)
{
  int k;
  for (k=0;k<nrcols;k++) {
    matrix[row+k*nrrows] = rowvector[k];
  }
}


static void freeMem(int freeCVODE)
{
  if(timesimvector != NULL)     free(timesimvector);
  if(initialconditions != NULL) free(initialconditions);
  if(parametervector != NULL)   free(parametervector);
  if(resultXdot != NULL)        free(resultXdot);
  if (freeCVODE) {
    if(statevalues != NULL)       free(statevalues);
    if(variablevalues != NULL)    free(variablevalues);
    if(variablevec != NULL)       free(variablevec);
    if(reactionvalues != NULL)    free(reactionvalues);
    if(reactionvec != NULL)       free(reactionvec);
    if(eventdataold != NULL)      free(eventdataold);
    if(eventvec != NULL)          free(eventvec);
    N_VDestroy_Serial(u);
    CVodeFree(&cvode_mem);
  }
}


static void errorMsg(char *text)
{
  freeMem(1);
  error(text);
}


