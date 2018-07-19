// read pandira output files
//  jorg kewisch 2005


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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


void realft(double *data,int n,int isign);

main(int argc, char ** argv)
{
	FILE * fp = fopen(argv[1], "r");
	if(!fp)
	{
		fprintf(stderr, "cant open %s for read\n", argv[1]);
		exit(-1);
	}

	FILE * fg = fopen("_plot.agr", "w");
	if(!fg)
	{
		fprintf(stderr, "cant open %s for write\n", argv[1]);
		exit(-1);
	}

	double rmin, rmax, zmin, zmax;
	int rstep, zstep;

	fscanf(fp, "%lf%lf%d", &rmin, &rmax, &rstep);
	fscanf(fp, "%lf%lf%d", &zmin, &zmax, &zstep);

	rstep++; zstep++;

	double br[zstep][rstep];
	double bz[zstep][rstep];

	int k=0;
	for(int z=0; z< zstep; z++)
	for(int r=0; r< rstep; r++)
	{
		if(fscanf(fp, "%lf%lf", &(br[z][r]), &(bz[z][r])) != 2) fprintf(stderr, "read err z=%d, r=%d\n", z,r);
// 		fprintf(stderr, "%e  %e\n", bz[z][r], br[z][r]);
// 		printf( "%d  %e\n", k++, br[z][r]);
	}
	

#if 1

	int zbest=-1;
	double dbobbest=1e33;
	for(int z=0; z< zstep; z++)
	{
		double zz=zmin+z*(zmax-zmin)/(zstep-1);
		double mean =0.;
		double sigm =0.;
		int nm=0;
		for(int r=0; r< rstep; r++)
		{
			double rr=rmin+r*(rmax-rmin)/(rstep-1);
			if(rr < 0.7)
			{
				double bb=bz[z][r];
				mean += bb;
				sigm += bb*bb;
				nm++;
			}
		}
		mean /= nm;
		sigm /= nm;
		double dbob=sqrt((sigm-mean*mean)/(mean*mean));
// 		fprintf(fg, "%e %e\n", zz, dbob);
// 		fprintf(stderr, "z=%d zz=%f dbob %e\n", z, zz, dbob);
		if(dbob < dbobbest)
		{
			dbobbest=dbob;
			zbest=z;
		}
	}
	double zz=zmin+zbest*(zmax-zmin)/(zstep-1);
	fprintf(stderr, "z=%d zz=%f dbob %e\n", zbest, zz, dbobbest);
	for(int r=0; r< rstep; r++)
	{
		double rr=rmin+r*(rmax-rmin)/(rstep-1);
		if(rr < 0.7) fprintf(fg, "%e %e\n",rr, bz[zbest][r]);
	}
	fprintf(fg, "&\n");
#endif


#if 0
	for(int r=0; r< 30; r+=3)
	{
		double rr=rmin+r*(rmax-rmin)/(rstep-1);
		fprintf(stderr, "r=%f\n",rr);
		for(int z=0; z< zstep; z++)
		{
			double zz=zmin+z*(zmax-zmin)/(zstep-1);
			fprintf(fg, "%e %e\n",zz, bz[z][r]);
// 			fprintf(fg, "%e %e\n",zz, br[z][r]);
		}
		fprintf(fg, "&\n");
	}
#endif



// 	for(int i=0; i< ystep; i++)
// 		printf("%d %e\n",i, bz[0][i]);



#if 0
	const int nfft=2048;
	double fftdata[nfft];
	for(int i=0; i< nfft; i++) fftdata[i]=0.;


	int r=0;
	{
		double rr=rmin+r*(rmax-rmin)/(rstep-1);
		fprintf(stderr, "r=%f\n",rr);
		for(int z=0; z< zstep; z++)
		{
			double zz=zmin+z*(zmax-zmin)/(zstep-1);
// 			fprintf(fg, "%e %e\n", zz, bz[z][r]);
			if(nfft/2+z == nfft) break;
			fftdata[nfft/2+z] =bz[z][r];
			fftdata[nfft/2-z] =bz[z][r];

		}
// 		fprintf(fg, "&\n");
	}


	realft(fftdata-1, nfft, 1);

	for(int i=0; i< nfft; i++) fprintf(fg, "%f\n",  fftdata[i]);
// 	for(int i=0; i< nfft; i++) fprintf(fg, "%d %f\n", i, fftdata[i]);
#endif

}



