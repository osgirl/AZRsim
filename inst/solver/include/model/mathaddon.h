#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#define pi 3.141592653589793

double minAZR(int nargs, ...)
{
    int k;
    double minimum;
    double checkmin;
    va_list vararg;
    va_start(vararg, nargs);
    minimum = va_arg(vararg, double);
    for (k=1; k<nargs; k++ ) {
        checkmin = va_arg(vararg, double);
        if (checkmin < minimum) {
            minimum = checkmin;
        }
    }
    va_end (vararg);
    return minimum;
}

double maxAZR(int nargs, ...)
{
    int k;
    double maximum;
    double checkmax;
    va_list vararg;
    va_start(vararg, nargs);
    maximum = va_arg(vararg, double);
    for (k=1; k<nargs; k++ ) {
        checkmax = va_arg(vararg, double);
        if (checkmax > maximum) {
            maximum = checkmax;
        }
    }
    va_end (vararg);
    return maximum;
}

double indexmaxAZR(int nargs, ...)
{
    int k;
    double maximum;
    double indexmaximum;
    double checkmax;
    va_list vararg;
    va_start(vararg, nargs);
    maximum = va_arg(vararg, double);
    indexmaximum = 1;
    for (k=1; k<nargs; k++ ) {
        checkmax = va_arg(vararg, double);
        if (checkmax > maximum) {
            maximum = checkmax;
            indexmaximum = (double)k+1;
        }
    }
    va_end (vararg);
    return indexmaximum;
}

double sign(double a)
{
    if (a<0) return -1;
    else if (a==0) return 0;
    else if (a>0) return 1;
    else return a;
}

double absAZR(double a)
{
    if (a < 0) return -a;
    else return a;
}

double mod(double a, double b)
{
    if (b == 0) return a;
    else if (a == b) return 0;
    return sign(b)*absAZR(a-floor(a/b)*b);
}

double nthrootAZR(double a, double n)
{
    return pow(a,1.0/n);
}

double andAZR(int nargs, ...)
{
    int k;
    va_list vararg;
    va_start(vararg, nargs);
    for (k=0; k<nargs; k++ ) {
        if (va_arg(vararg, double) == 0) {
            va_end (vararg);
            return 0;
        }
    }
    va_end (vararg);
    return 1;
}

double orAZR(int nargs, ...)
{
    int k;
    va_list vararg;
    va_start(vararg, nargs);
    for (k=0; k<nargs; k++ ) {
        if (va_arg(vararg, double) != 0) {
            va_end (vararg);
            return 1;
        }
    }
    va_end (vararg);
    return 0;
}

double piecewiseAZR(int nargs, ...)
{
    int k, oddnumber;
    double *data;
    double returnvalue;
    va_list vararg;
    va_start(vararg, nargs);
    data = (double *)calloc(nargs,sizeof(double));
    for (k=0; k<nargs; k++ ) {
        data[k] = va_arg(vararg, double);
    }
    oddnumber = nargs % 2;
    for (k=0;k<nargs-oddnumber;k+=2) {
        if (data[k+1] != 0) {
            returnvalue = data[k];
            free(data); 
            return returnvalue;
        }
    }
    if (oddnumber != 0) {
            returnvalue = data[nargs-1];
            free(data); 
            return returnvalue;
    }
    else {
        error("A piecewise statement is wrongly defined - missing (but needed) default value.");
        free(data); 
        return 0; 
    }
}

double gt(double a, double b)
{
    if (a > b) return 1;
    else return 0;
}

double ge(double a, double b)
{
    if (a >= b) return 1;
    else return 0;
}

double lt(double a, double b)
{
    if (a < b) return 1;
    else return 0;
}

double le(double a, double b)
{
    if (a <= b) return 1;
    else return 0;
}

double eq(double a, double b)
{
    if (a*b < 0) return 0;
    else if ( absAZR(a-b) < 10.0*FLT_MIN ) return 1;
    else return 0;
}

