#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

void spline_pchip_set( int n, double x[], double f[], double d[] );
void spline_pchip_val( int n, double x[], double f[], double d[],  int ne, double xe, double *feptr );
int i4_max( int i1, int i2 );
int chfev( double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe, double *feptr, int next[] );
double r8_min( double x, double y );
double r8_max( double x, double y );
double pchst( double arg1, double arg2 );


/****************************************************************************/

double interpcsAZR(int nargs, int n, double xe, ...)

/****************************************************************************

  Purpose:

    interpcsAZR is the calling function to evaluate a piecewise cubic spline
    function and calculte the Hermite interpolant.


  Modified:

    1 October 2007

  Author:

    Basti Bergdahl,
    Department of Applied Microbiology,
    Lund Insitute of Technology.

  Calling syntax:

    interpcsAZR(NARGS, N, XE, X1,...,XN,F1,...,FN)

  Parameters:

    Input, int NARGS, the total number of input arguments (not counting NARGS itself).
    The number should be even, otherwise an error will occure.

    Input, int N, the number of data points (rows in the table).  N must be at least 2.

    Input, double XE, point at which the function is to be evaluated.

    Input, double X1,...,XN, the strictly increasing independent variable values.

    Input, double F1,...,FN, dependent variable values to be interpolated.  This
    routine is designed for monotonic data, but it will work for any F-array.
    It will force extrema at points where monotonicity switches direction.

    Output, double FE, the value of the cubic Hermite function at XE.
*/
{
    va_list vararg;
    int k, oddnumber;
    int ne = 1;
    double *data, *x, *f, *d, *feptr;
    double fe;

    /* Check if table data has been given */
    if (nargs == 2)
    {
        error("interpcsAZR error: Input incorrect. No data for the table has been supplied.");
        return 0;
    }

    /* Check if odd or even number of input arguments */
    oddnumber = nargs % 2;
    if (oddnumber != 0)
    {
        error("interpcsAZR error: Input incorrect. The number of inputs should be even.");
        return 0;
    }

    /* If number of input arguments in odd, check if
     * the X and Y vectors are of the same size */
    if ((nargs-2)/2 != n)
    {
        error("interpcsAZR error: Input incorrect. The X and Y vectors are either of different lengths oran incorrect number for the data points has been given.");
        return 0;
    }

    /* Read in all variable input arguments in double array */
    va_start(vararg, xe);
    data = (double *)calloc((n*2),sizeof(double));
    for (k=0; k<nargs-2; k++ )
    {
        data[k] = va_arg(vararg, double);
    }

    /* Create the X and Y arrays making up the table */
    x = (double *)calloc(n,sizeof(double));
    f = (double *)calloc(n,sizeof(double));
    for (k=0; k<n; k++)
    {
        x[k] = data[k];
        f[k] = data[k+n];
    }
    free(data); /* free temporary array */

    /* handle off limit values */
    if (xe < x[0]) return(f[0]);
    if (xe > x[n-1]) return(f[n-1]);

    /* Calculate the derivatives for the piecewise cubic Hermite interpolant */
    d = (double *)calloc(n,sizeof(double));
    spline_pchip_set(n, x, f, d);

    /* Evaluate the cubic Hermite function at x */
    feptr = &fe;
    spline_pchip_val(n, x, f, d, ne, xe, feptr);

    return(fe);
}



/****************************************************************************/

void spline_pchip_set( int n, double x[], double f[], double d[] )

/****************************************************************************

  Purpose:

    SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.

  Discussion:

    This routine computes what would normally be called a Hermite
    interpolant.  However, the user is only required to supply function
    values, not derivative values as well.  This routine computes
    "suitable" derivative values, so that the resulting Hermite interpolant
    has desirable shape and monotonicity properties.

    The interpolant will have an extremum at each point where
    monotonicity switches direction.

    The resulting piecewise cubic Hermite function may be evaluated
    by SPLINE_PCHIP_VAL..

    This routine was originally called "PCHIM".

  Modified:

    14 August 2005

  Author:

    Fred Fritsch,
    Mathematics and Statistics Division,
    Lawrence Livermore National Laboratory.

    C++ translation by John Burkardt.

  Reference:

    Fred Fritsch, Ralph Carlson,
    Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 2, April 1980, pages 238-246.

    Fred Fritsch, Judy Butland,
    A Method for Constructing Local Monotone Piecewise
    Cubic Interpolants,
    SIAM Journal on Scientific and Statistical Computing,
    Volume 5, Number 2, 1984, pages 300-304.

  Parameters:

    Input, int N, the number of data points.  N must be at least 2.

    Input, double X[N], the strictly increasing independent
    variable values.

    Input, double F[N], dependent variable values to be interpolated.  This
    routine is designed for monotonic data, but it will work for any F-array.
    It will force extrema at points where monotonicity switches direction.

    Output, double D[N], the derivative values at the
    data points.  If the data are monotonic, these values will determine
    a monotone cubic Hermite function.
*/
{
  double del1;
  double del2;
  double dmax;
  double dmin;
  double drat1;
  double drat2;
  double dsave;
  double h1;
  double h2;
  double hsum;
  double hsumt3;
  int i;
  int ierr;
  int nless1;
  double temp;
  double w1;
  double w2;

  /*  Check the arguments. */

  if ( n < 2 )
  {
    ierr = -1;
    error("SPLINE_PCHIP_SET - Fatal error! Number of data points less than 2.");
  }

  for ( i = 1; i < n; i++ )
  {
    if ( x[i] <= x[i-1] )
    {
      ierr = -3;
      error("SPLINE_PCHIP_SET - Fatal error! X array not strictly increasing.");
    }
  }

  ierr = 0;
  nless1 = n - 1;
  h1 = x[1] - x[0];
  del1 = ( f[1] - f[0] ) / h1;
  dsave = del1;

/*  Special case N=2, use linear interpolation. */

  if ( n == 2 )
  {
    d[0] = del1;
    d[n-1] = del1;
    return;
  }

/*  Normal case, 3 <= N. */

  h2 = x[2] - x[1];
  del2 = ( f[2] - f[1] ) / h2;

/*  Set D(1) via non-centered three point formula, adjusted to be
    shape preserving. */

  hsum = h1 + h2;
  w1 = ( h1 + hsum ) / hsum;
  w2 = -h1 / hsum;
  d[0] = w1 * del1 + w2 * del2;

  if ( pchst ( d[0], del1 ) <= 0.0 )
  {
    d[0] = 0.0;
  }

/*  Need do this check only if monotonicity switches. */

  else if ( pchst ( del1, del2 ) < 0.0 )
  {
     dmax = 3.0 * del1;

     if ( fabs ( dmax ) < fabs ( d[0] ) )
     {
       d[0] = dmax;
     }

  }

/*  Loop through interior points. */

  for ( i = 2; i <= nless1; i++ )
  {
    if ( 2 < i )
    {
      h1 = h2;
      h2 = x[i] - x[i-1];
      hsum = h1 + h2;
      del1 = del2;
      del2 = ( f[i] - f[i-1] ) / h2;
    }

/*  Set D(I)=0 unless data are strictly monotonic. */

    d[i-1] = 0.0;

    temp = pchst ( del1, del2 );

    if ( temp < 0.0 )
    {
      ierr = ierr + 1;
      dsave = del2;
    }

/*  Count number of changes in direction of monotonicity. */

    else if ( temp == 0.0 )
    {
      if ( del2 != 0.0 )
      {
        if ( pchst ( dsave, del2 ) < 0.0 )
        {
          ierr = ierr + 1;
        }
        dsave = del2;
      }
    }

/*  Use Brodlie modification of Butland formula. */

    else
    {
      hsumt3 = 3.0 * hsum;
      w1 = ( hsum + h1 ) / hsumt3;
      w2 = ( hsum + h2 ) / hsumt3;
      dmax = r8_max ( fabs ( del1 ), fabs ( del2 ) );
      dmin = r8_min ( fabs ( del1 ), fabs ( del2 ) );
      drat1 = del1 / dmax;
      drat2 = del2 / dmax;
      d[i-1] = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }

/*  Set D(N) via non-centered three point formula, adjusted to be
    shape preserving. */

  w1 = -h2 / hsum;
  w2 = ( h2 + hsum ) / hsum;
  d[n-1] = w1 * del1 + w2 * del2;

  if ( pchst ( d[n-1], del2 ) <= 0.0 )
  {
    d[n-1] = 0.0;
  }
  else if ( pchst ( del1, del2 ) < 0.0 )
  {

/*  Need do this check only if monotonicity switches. */

    dmax = 3.0 * del2;

    if ( fabs ( dmax ) < abs ( d[n-1] ) )
    {
      d[n-1] = dmax;
    }

  }
  return;
}

/****************************************************************************/

void spline_pchip_val( int n, double x[], double f[], double d[],
  int ne, double xe, double *feptr )

/****************************************************************************

  Purpose:

    SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.

  Description:

    This routine may be used by itself for Hermite interpolation, or as an
    evaluator for SPLINE_PCHIP_SET.

    This routine evaluates the cubic Hermite function at the points XE.

    Most of the coding between the call to CHFEV and the end of
    the IR loop could be eliminated if it were permissible to
    assume that XE is ordered relative to X.

    CHFEV does not assume that X1 is less than X2.  Thus, it would
    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
    decreasing X array by simply running the IR loop backwards
    and reversing the order of appropriate tests.

    The present code has a minor bug, which I have decided is not
    worth the effort that would be required to fix it.
    If XE contains points in [X(N-1),X(N)], followed by points less than
    X(N-1), followed by points greater than X(N), the extrapolation points
    will be counted (at least) twice in the total returned in IERR.

    The evaluation will be most efficient if the elements of XE are
    increasing relative to X; that is, for all J <= K,
      X(I) <= XE(J)
    implies
      X(I) <= XE(K).

    If any of the XE are outside the interval [X(1),X(N)],
    values are extrapolated from the nearest extreme cubic,
    and a warning error is returned.

    This routine was originally named "PCHFE".

  Modified:

    14 August 2005

  Author:

    Fred Fritsch,
    Mathematics and Statistics Division,
    Lawrence Livermore National Laboratory.

    C++ translation by John Burkardt.

  Reference:

    Fred Fritsch, Ralph Carlson,
    Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 2, April 1980, pages 238-246.

  Parameters:

    Input, int N, the number of data points.  N must be at least 2.

    Input, double X[N], the strictly increasing independent
    variable values.

    Input, double F[N], the function values.

    Input, double D[N], the derivative values.

    Input, int NE, the number of evaluation points.

    Input, double XE, point at which the function is to
    be evaluated.

    Output, double FE, the value of the cubic Hermite
    function at XE.
*/
{
  int i;
  int ierc;
  int ierr;
  int ir;
  int j;
  int j_first;
  int j_new;
  int j_save;
  int next[2];
  int nj;

/*  Check arguments. */

  if ( n < 2 )
  {
    ierr = -1;
    error("SPLINE_PCHIP_VAL - Fatal error! Number of data points less than 2.");
  }

  for ( i = 1; i < n; i++ )
  {
    if ( x[i] <= x[i-1] )
    {
      ierr = -3;
      error("SPLINE_PCHIP_VAL - Fatal error! X array not strictly increasing.");
    }
  }

  if ( ne < 1 )
  {
    ierr = -4;
    error("SPLINE_PCHIP_VAL - Fatal error! Number of evaluation points less than 1.");
  }

  ierr = 0;

/*  Loop over intervals.
  The interval index is IL = IR-1.
  The interval is X(IL) <= X < X(IR).
*/
  j_first = 1;
  ir = 2;

  for ( ; ; )
  {
/*
  Skip out of the loop if have processed all evaluation points.
*/
    if ( ne < j_first )
    {
      break;
    }
/*
  Locate all points in the interval.
*/
    j_save = ne + 1;

    for ( j = j_first; j <= ne; j++ )
    {
      if ( x[ir-1] <= xe )
      {
        j_save = j;
        if ( ir == n )
        {
          j_save = ne + 1;
        }
        break;
      }
    }
/*
  Have located first point beyond interval.
*/
    j = j_save;

    nj = j - j_first;
/*
  Skip evaluation if no points in interval.
*/
    if ( nj != 0 )
    {
/*
  Evaluate cubic at XE(J_FIRST:J-1).
*/
      ierc = chfev ( x[ir-2], x[ir-1], f[ir-2], f[ir-1], d[ir-2], d[ir-1],
        nj, xe, feptr, next );

      if ( ierc < 0 )
      {
        ierr = -5;
        error("SPLINE_PCHIP_VAL - Fatal error! Error return from CHFEV.");
      }
/*
  In the current set of XE points, there are NEXT(2) to the right of X(IR).
*/
      if ( next[1] != 0 )
      {
        if ( ir < n )
        {
          ierr = -5;
          error("SPLINE_PCHIP_VAL - Fatal error! IR < N.");
        }
/*
  These are actually extrapolation points.
*/
        ierr = ierr + next[1];

      }
/*
  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
*/
      if ( next[0] != 0 )
      {
/*
  These are actually extrapolation points.
*/
        if ( ir <= 2 )
        {
          ierr = ierr + next[0];
        }
        else
        {
          j_new = -1;

          for ( i = j_first; i <= j-1; i++ )
          {
            if ( xe < x[ir-2] )
            {
              j_new = i;
              break;
            }
          }

          if ( j_new == -1 )
          {
            ierr = -5;
            error("SPLINE_PCHIP_VAL - Fatal error!Could not bracket the data point.");
          }
/*
  Reset J.  This will be the new J_FIRST.
*/
          j = j_new;
/*
  Now find out how far to back up in the X array.
*/
          for ( i = 1; i <= ir-1; i++ )
          {
            if ( xe < x[i-1] )
            {
              break;
            }
          }
/*
  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .

  Reset IR, recognizing that it will be incremented before cycling.
*/
          ir = i4_max ( 1, i-1 );
        }
      }

      j_first = j;
    }

    ir = ir + 1;

    if ( n < ir )
    {
      break;
    }

  }

  return;
}

/****************************************************************************/

int chfev( double x1, double x2, double f1, double f2, double d1, double d2,
  int ne, double xe, double *feptr, int next[] )

/****************************************************************************

  Purpose:

    CHFEV evaluates a cubic polynomial given in Hermite form.

  Discussion:

    This routine evaluates a cubic polynomial given in Hermite form at an
    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
    be useful directly as an evaluator for a piecewise cubic
    Hermite function in applications, such as graphing, where
    the interval is known in advance.

    The cubic polynomial is determined by function values
    F1, F2 and derivatives D1, D2 on the interval [X1,X2].

  Modified:

    12 August 2005

  Author:

    Fred Fritsch,
    Mathematics and Statistics Division,
    Lawrence Livermore National Laboratory.

    C++ translation by John Burkardt.

  Reference:

    Fred Fritsch, Ralph Carlson,
    Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 2, April 1980, pages 238-246.

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input, double X1, X2, the endpoints of the interval of
    definition of the cubic.  X1 and X2 must be distinct.

    Input, double F1, F2, the values of the function at X1 and
    X2, respectively.

    Input, double D1, D2, the derivative values at X1 and
    X2, respectively.

    Input, int NE, the number of evaluation points.

    Input, double XE, the point at which the function is to
    be evaluated.  If the value of XE is outside the interval
    [X1,X2], a warning error is returned in NEXT.

    Output, double FE, the value of the cubic function
    at the point XE.

    Output, int NEXT[2], indicates the number of extrapolation points:
    NEXT[0] = number of evaluation points to the left of interval.
    NEXT[1] = number of evaluation points to the right of interval.

    Output, int CHFEV, error flag.
    0, no errors.
    -1, NE < 1.
    -2, X1 == X2.
*/
{
  double c2;
  double c3;
  double del1;
  double del2;
  double delta;
  double h;
  int ierr;
  double x;
  double fe;
  double xma;
  double xmi;

  if ( ne < 1 )
  {
    ierr = -1;
    error("CHFEV - Fatal error!\n");
    error("  Number of evaluation points is less than 1.\n");
    printf("  NE = %d\n", ne);
    return ierr;
  }

  h = x2 - x1;

  if ( h == 0.0 )
  {
    ierr = -2;
    error("\n");
    error("CHFEV - Fatal error!\n");
    error("  The interval [X1,X2] is of zero length.\n");
    return ierr;
  }
/*
  Initialize.
*/
  ierr = 0;
  next[0] = 0;
  next[1] = 0;
  xmi = r8_min ( 0.0, h );
  xma = r8_max ( 0.0, h );
/*
  Compute cubic coefficients expanded about X1.
*/
  delta = ( f2 - f1 ) / h;
  del1 = ( d1 - delta ) / h;
  del2 = ( d2 - delta ) / h;
  c2 = -( del1 + del1 + del2 );
  c3 = ( del1 + del2 ) / h;
    x = xe - x1;
    fe = f1 + x * ( d1 + x * ( c2 + x * c3 ) );
    *feptr = fe;
/*
  Count the extrapolation points.
*/
    if ( x < xmi )
    {
      next[0] = next[0] + 1;
    }

    if ( xma < x )
    {
      next[1] = next[1] + 1;
    }


  return 0;
}

/****************************************************************************/

double pchst( double arg1, double arg2 )

/****************************************************************************

  Purpose:

    PCHST: PCHIP sign-testing routine.

  Discussion:

    This routine essentially computes the sign of ARG1 * ARG2.

    The object is to do this without multiplying ARG1 * ARG2, to avoid
    possible over/underflow problems.

  Modified:

    12 August 2005

  Author:

    Fred Fritsch,
    Mathematics and Statistics Division,
    Lawrence Livermore National Laboratory.

    C++ translation by John Burkardt.

  Reference:

    Fred Fritsch, Ralph Carlson,
    Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 2, April 1980, pages 238-246.

  Parameters:

    Input, double ARG1, ARG2, two values to check.

    Output, double PCHST,
    -1.0, if ARG1 and ARG2 are of opposite sign.
     0.0, if either argument is zero.
    +1.0, if ARG1 and ARG2 are of the same sign.
*/
{
  double value = 0.0;

  if ( arg1 == 0.0 )
  {
    value = 0.0;
  }
  else if ( arg1 < 0.0 )
  {
    if ( arg2 < 0.0 )
    {
      value = 1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = -1.0;
    }
  }
  else if ( 0.0 < arg1 )
  {
    if ( arg2 < 0.0 )
    {
      value = -1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = 1.0;
    }
  }

  return(value);
}

/****************************************************************************/

double r8_max( double x, double y )

/****************************************************************************

  Purpose:

    R8_MAX returns the maximum of two R8's.

  Modified:

    10 January 2002

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  if ( y < x )
  {
    return(x);
  }
  else
  {
    return(y);
  }
}
/****************************************************************************/

double r8_min( double x, double y )

/****************************************************************************

  Purpose:

    R8_MIN returns the minimum of two R8's.

  Modified:

    09 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  if ( y < x )
  {
    return(y);
  }
  else
  {
    return(x);
  }
}


/****************************************************************************/

int i4_max ( int i1, int i2 )

/****************************************************************************

  Purpose:

    I4_MAX returns the maximum of two I4's.

  Modified:

    13 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.

*/
{
  if ( i2 < i1 )
  {
    return(i1);
  }
  else
  {
    return(i2);
  }

}
