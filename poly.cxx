// this is a direct translation of the NUMEICAL RECEIPES function into C++
// Arrays start with index 0 instead of 1
// the functions are wrapped in a c++ class.
// jorg kewisch 2002

// the following copy right notice applies only to the changes I made:

// Notice: This computer software was prepared by Brookhaven Science
// Associates, (the Contractor), under Federal Contract No. DE-AC02-98CH10886
// (the Contract) with the U.S. Department of Energy (DOE).  All rights in
// this computer software are reserved by DOE on behalf of the United States
// Government and the Contractor as provided in the Contract.  You are
// authorized to use this computer software solely for U.S. Governmental
// purposes.  This software is not to be released or distributed in any form
// or manner except as may be provided in your Software User Acknowledgment.
// NEITHER THE UNITED STATES GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
// THIS SOFTWARE.  This notice including this sentence must appear on any
// copies of this computer software.



#include <stdio.h>
#include <math.h>
#include "poly.hxx"

//	fit of a polynom of order n to m data points (x[],y[])
//	coeff is a m long array of double pointers pointing to 
//	m+1 long arrays of doubles.
//	a is a n long array of doubles for the result.

int polyInit(double * x, int m, double * coef[])
{
   // m = vector length

   // Construct the orthogonal polynomials.
   for(int i = 0; i<m; i++)
   for(int j = 0; j<=m; j++)
   coef[i][j] = 0.0;

   coef[0][0] = 1.0;
   for (int i = 1; i < m; i++)
   {
//       printf("polyInit %d\n", i);
      double w1   = 0.;
      double w2   = 0.;
      double den1 = 0.;
      double den2 = 0.;
      for (int j = 0; j < m; j++)
      {
   	 double pm1= 0.;
   	 double pm2= 0.;
	 for (int k = m - 1; k >= 0; k--)
	 {
	    pm1 = coef[i - 1][k] + x[j] * pm1;
	    if (i > 1) pm2 = coef[i - 2][k] + x[j] * pm2;
	 }
	 w1 += x[j] * pm1 * pm1;
	 w2 += x[j] * pm1 * pm2;
	 den1 += pm1 * pm1;
	 den2 += pm2 * pm2;
      }
      if(den1 == 0. )
      {
	fprintf(stderr,"polyInit: denominator 1 = %e\n", den1);
	return 1;
      }
      else  w1 /= den1;

      if(i > 1)
      {
         if(den2 == 0.)
         {
	   fprintf(stderr,"polyInit: denominator 2 = %e\n", den2);
	   return 1;
         }
	 else w2 /= den2;
	 coef[i][0] = -w1 * coef[i - 1][0] - w2 * coef[i - 2][0];
	 for (int j = 1; j < m; j++)
	 {
	    coef[i][j] = coef[i - 1][j - 1] - w1 * coef[i - 1][j] -
	       w2 * coef[i - 2][j];
	 }
      }
      else
      {
	 coef[i][0] = -w1 * coef[i - 1][0];
	 for (int j = 1; j < m; j++)
	 {
	    coef[i][j] = coef[i - 1][j - 1] - w1 * coef[i - 1][j];
	 }
      }
   }
   return 0;
}

// -------------------------------------------------------------------------

double polyFit(double * x, double * y, int m, double * coef[], double * a, int n)
{
   double * dPoly = new double[n];

   // Compute the coefficient on each polynomial.
   for (int i = 0; i < n; i++)
   {
//       printf("polyFit1 %d\n", i);
      double den1 = 0.;
      double w1 = 0.;
      for (int j = 0; j < m; j++)
      {
	 double pm1 = 0.0;
	 for (int k = m - 1; k >= 0; k--) {
	    pm1 = coef[i][k] + x[j] * pm1;
	 }
	 w1 += y[j] * pm1;
	 den1 += pm1 * pm1;
      }
      if (den1 != 0.0) dPoly[i] = w1 / den1;
      else             dPoly[i] =0.;
   }

   // Compute coefficients for the approximating polynomial.
   for (int i = 0; i < n; i++)
   {
//       printf("polyFit22 %d\n", i);
      double w1 = 0.0;
      for(int j = 0; j < n; j++) w1 += coef[j][i] * dPoly[j];
      a[i] = w1;
   }
   delete [] dPoly;

   // Compute the chi-square error.
   double chisq = 0.0;
   for (int i = 0; i < m; i++)
   {
//       printf("polyFit3 %d\n", i);
      double pm1 = 0.0;
      for (int j = n - 1; j >= 0; j--) pm1 = a[j] + x[i] * pm1;
      double diff = y[i] - pm1;
      chisq += diff * diff;
   }

   return chisq/m;
}




double polyValue(double x, double * a, int n)
{
	n--;
	double r = a[n]; 
	n--;
	while(n >= 0)
	{
		r *= x; r += a[n]; n--;
	}
	return r;
}
	
double polyPrime(double x, double * a, int n)
{
	n--;
	double r = n*a[n]; 
	n--;
	while(n >= 1)
	{
		r *= x; r += n*a[n]; n--;
	}
	return r;
}
	
double polyDoublePrime(double x, double * a, int n)
{
	n--;
	double r = n*(n-1)*a[n]; 
	n--;
	while(n >= 2)
	{
		r *= x; r += n*(n-1)*a[n]; n--;
	}
	return r;
}



	
int polyTrackFill(double ymax, int maxRows, double * xr, double * yr, double * yp,
	double * a, int n, double x0, double xf, int nrows, int first)
{
  double tolerance = 2.e-4*ymax;
  int nrows1 = nrows;
  for(int j=0; j<10; j++)
    {
      nrows = polyTrack(tolerance, maxRows, xr, yr, yp, a, n, x0, xf, nrows1, first);
      if(n < 3)  break;// all linear
      printf(" tol = %f, # of rows = %d\n", tolerance, nrows);
      if(nrows <= maxRows && nrows > maxRows*0.9) break;
      double f= double(maxRows*0.95)/double(nrows);
      tolerance /= f*f;
//       if(tolerance < 1.e-5*ymax) break;
    }
    return nrows;
}




int polyTrack(double tolerance, int maxRows, double * xr, double * yr, double * yp,
	double * a, int n, double x0, double xf, int nrows, int first)
{
	if(first && nrows < maxRows)
	{
		xr[nrows] = x0;
		yr[nrows] = polyValue(x0, a, n);
		if(yp) yp[nrows] = polyPrime(x0, a, n);
		nrows++;
	}

	if(n < 3) // all linear
	{
		xr[nrows] = xf;
		yr[nrows] = polyValue(xf, a, n);
		if(yp) yp[nrows] = polyPrime(xf, a, n);
		nrows++;
		return nrows;
	}

	double x =x0;
	while(x < xf)
	{
		double dp = polyDoublePrime(x, a, n);
		double dx;
		if(dp == 0.) dx = (xf - x0)/maxRows;
		else dx = sqrt(fabs(8.*tolerance/dp));
		printf("dp = %e dx = %e tol = %e\n", dp,dx,tolerance);
		x +=dx;
		if(x > xf) x = xf;


		if(nrows < maxRows)
		{
			xr[nrows]  = x;
			yr[nrows]  = polyValue(x, a, n);
			if(yp) yp[nrows] = polyPrime(x, a, n);
		}
		nrows++;
	}

	return nrows;
}

