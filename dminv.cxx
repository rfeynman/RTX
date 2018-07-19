//   14/11/86 611282107  member name  dminv    (petlib)      fortran77
//
//     .................................................................
//
//        subroutine dminv
//
//        purpose
//           invert a matrix
//
//        usage
//           double d = dminv(a,n)
//
//        description of parameters
//           a - input matrix, destroyed in computation and replaced by
//               resultant inverse.
//           n - order of matrix a
//           d - resultant determinant
//
//        remarks
//           matrix a must be a general matrix
//
//        subroutines and function subprograms required
//           none
//
//        method
//           the standard gauss-jordan method is used. the determinant
//           is also calculated. a determinant of zero indicates that
//           the matrix is singular.
//
//        translated into fortran77 11.86 jkl
//        translated into c++       05.98 jkl
//
//
//	copyright: The original version is from cern (I think).
//
//
//
//     .................................................................
//

#include <math.h>


double dminv(double *a, int n)
{
      int * l = new int[n];
      int * m = new int[n];
      double d=1.0;
// 	printf("matrix of %d\n",n);
// 	for(int i=0; i<n*n; i++) printf("%3.1lf  ", a[i]);
// 	printf("\n");

      int nk=-n;
      for(int k=0; k<n; k++)
      {
         nk=nk+n;
         l[k]=k;
         m[k]=k;
         int kk=nk+k;
         double biga=a[kk];
	 for(int j=k; j<n; j++)
	 {
            int iz=n*j;
	    for(int i=k; i<n; i++)
	    {
               int ij=iz+i;
               if(fabs(biga) < fabs(a[ij]))
	       {
                  biga=a[ij];
                  l[k]=i;
                  m[k]=j;
               }
	    }
	 }
         int j=l[k];
         if(j > k) {
            int ki=k-n;
	    for(int i=0; i<n; i++)
	    {
               ki=ki+n;
               double hold = -a[ki];
               int ji=ki-k+j;
               a[ki]=a[ji];
               a[ji] =hold;
	    }
         }
         int i=m[k];
         if(i > k) {
            int jp=n*i;
	    for(int j=0; j<n; j++)
	    {
               int jk=nk+j;
               int ji=jp+j;
               double hold= - a[jk];
               a[jk]=a[ji];
               a[ji] =hold;
	    }
         }
         if(biga == 0.0) {
	    delete [] l;
	    delete [] m;
            return 0.0;
         }
	 for(int i=0; i<n; i++)
	 {
            if(i != k) {
               int ik=nk+i;
               a[ik]=a[ik]/(-biga);
            }
	 }

	 for(i=0; i<n; i++)
	 {
            int ik=nk+i;
            double hold=a[ik];
            int ij=i-n;
	    for(int j=0; j<n; j++)
	    {
               ij=ij+n;
               if(i != k && j != k) {
                  int kj=ij-i+k;
                  a[ij]=hold*a[kj]+a[ij];
               }
	    }
	 }
         int kj=k-n;
	 for(j=0; j<n; j++)
	 {
            kj=kj+n;
            if(j != k) a[kj]=a[kj]/biga;
	 }
         d=d*biga;
         a[kk]=1.0/biga;
      }

      for(int k=n-2; k>=0; k--)
      {
         int i=l[k];
         if(i > k) {
            int jq=n*k;
            int jr=n*i;
	    for(int j=0; j<n; j++)
	    {
               int jk=jq+j;
               double hold=a[jk];
               int ji=jr+j;
               a[jk]=-a[ji];
               a[ji] =hold;
	    }
         }
         int j=m[k];
         if(j > k) {
            int ki=k-n;
	    for(int i=0; i<n; i++)
	    {
               ki=ki+n;
               double hold=a[ki];
               int ji=ki-k+j;
               a[ki]=-a[ji];
               a[ji] =hold;
	    }
         }
      }
      delete [] l;
      delete [] m;
      return d;
}
