


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


// this program creates a particle distribution for parmela and writes it into a file "PART_RFQ.DST"
// The particle coordinated can then be read by parmela with the statement:
// INPUT 40 <number of particles>
// PART_RFQ.DST is the default name of the particle file, it can be changed in the LANL.INI file by setting
// Part_In_Dst=<file name>
// Remember that the file names follow DOS conventions

// preParmela  reads the parmela input file and looks for !@ lines. 
// The first line in the parmela input file  must be the "run" line.

// the following commands are recognized:
// @elypsoid <nmax> <lcathode> <radius> <a> <b> <c> <d> <e> <f> <xpmax>
//	make an elypsoid with <nmax> particles, length <lcathode> in degrees, <a> <b> <c> <d> <e> <f> dimensions, <xpmax> transerse divergence.: 
//
// @elypsoid2    <nmax> <lcathode> <r1> <r2> <a> <b> <xpmax>
// @beercan      <nmax> <lcathode> <r1> <xpmax>
// @beercan5     <nmax> <lcathode> <r1> <xpmax>, like @beercan and scale xpmax /= 0.32810;
// @cylinder     <nmax> <lcathode> <r1> <a> <xpmax>, like @beercan and do something funny to the edges
// @beercancone  <nmax> <lcathode> <r1> <xpmax> <head density> <tail density>  @beercan5 can where the density varies linearly
// @gaussian     <nmax> <lcathode> r1> <xpmax>
// @ellipsoid4   <nmax> <lcathode> <ra> <rb> <pa> <pb> <xpmax>
// @ellipsoid5   <nmax> <lcathode> <rm> <rrat> <pa> <pb> <xpmax>


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Solver.hxx"
#include "tools.hxx"
#include "ObjectiveFunction.hxx"

#ifdef WIN32
#include <crtdbg.h>
#endif

#include <unistd.h>

//  calculate  particle distribution for parmela
//  Jorg Kewisch, 3/2006
//
// input: !@dipole  <rfdata#> <lhalf> <lmag> <lfringe> <energy> <angle> <f1> <f2>
// where e1 = angle*f1...
//

void sobseq(int *n, double x[]);
inline double square(double x) { return x*x;}
inline double qube(double x) { return x*x*x;}

const double elmass=   .5109999;	// [MeV]    must agree with line 9 of ImpactT.in
const double elchar=   1.601917e-19;	// [C]
const double clight=   2.997924562e8;	// [m/s]

int hist[501];  
int histn[501];  
double histx[501];
double histxp[501];


//
//	two elipses: center at pa and pb on the x-axis, half axis ra, la  and rb, lb
//	conditions: -l < pa < +l,  -l < pb < +l
//	if l, pa, and pb is given then la=l+pa and lb=l-pb; pa is normally negative
//	calculate: xa and xb, where the two ellipses are connected by a itangetial straight line;
//

class Tangent : public UnconstrainedObjectiveFunction
{
  public:
	double xa, xb;
	double fa, fb;
	double fap, fbp;
	Tangent(double l, double _pa,  double _ra, double _pb,  double _rb);
    	~Tangent(){};
    	double eval(Vector v, int *nerror=NULL);
	double pa;
	double la;
	double ra;

	double pb;
	double lb;
	double rb;
	double diffRadius, diffSlope;
};




double Tangent::eval(Vector X, int *nerror)
{

    	double *x=X; xa=x[0]; xb = x[1];
	if(nerror) if( fabs(xa-pa) >= la || fabs(xb-pb) >= lb) { *nerror=1; return 1e33;} else *nerror =0;
	
	// the y value of the elipses
	fa = ra * sqrt( 1. - square(xa - pa)/square(la));
	fb = rb * sqrt( 1. - square(xb - pb)/square(lb));

	// the slope value of the elipses
	fap = -square(ra)/square(la) * (xa-pa)/ fa;
	fbp = -square(rb)/square(lb) * (xb-pb)/ fb;

	// the slope should be the same
	diffSlope = fap - fbp;
	// the line should go through both points
	diffRadius = fa + fap*(xb-xa)-fb;


	double r=diffRadius*diffRadius+diffSlope*diffSlope;
// 	printf("eval x(%15.7f %15.7f) f(%15.7f %15.7f) f'(%15.7e %15.7e) slope %15.7e connect %15.7e val %15.7e\n", xa, xb, fa, fb, fap, fbp, diffSlope, diffRadius, r);
    	return r;
}

Tangent::Tangent(double l, double _pa,  double _ra, double _pb,  double _rb)
{

	pa=_pa;
	ra=_ra;
	pb=_pb;
	rb=_rb;
	la=l+pa;
	lb=l-pb;

	diffRadius=1.e33;
	diffSlope=1.e33;
    	t=1;
    	xStart.setSize(2);
    	xStart[0]= pa ;
    	xStart[1]= pb ;
}


//    print 4 particles into the file, so that the beam is symmetric in radius and slope
void printOne(FILE * ff, double x, double xp, double y, double yp, double z, double wr)
{
	double scalexy=1.;
	fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  scalexy*x,  xp,   scalexy*y,  yp, z, wr);
}
//    print 4 particles into the file, so that the beam is symmetric in radius and slope
void printFour(FILE * ff, double x, double xp, double y, double yp, double z, double wr)
{
	double scalexy=1.;	// this used to be 100., don't know why
	fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  scalexy*x,  xp,   scalexy*y,  yp, z, wr);
	fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  scalexy*x,  xp,  -scalexy*y, -yp, z, wr);
	fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n", -scalexy*x, -xp,   scalexy*y,  yp, z, wr);
	fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n", -scalexy*x, -xp,  -scalexy*y, -yp, z, wr);
}


main(int argc, char ** argv)
{
	int pcount=0;
	char * home_directory = getenv("HOME");
	if( !home_directory)
	{
		fprintf(stderr, "$HOME is not defined\n");
		exit(-1);
	}

	int shapeagr=0;
	int atCathode=0;
	while( ! strcmp(argv[1], "-s")   || ! strcmp(argv[1], "-t" ))
	{
		if ( ! strcmp(argv[1], "-s")  ) shapeagr=1;
		if ( ! strcmp(argv[1], "-t")  ) atCathode=1;
		argc--; argv++;
	}


	if(argc != 2)
	{
		fprintf(stderr, "usage: %s [ -s ] [ -t ] <file>\n", argv[0]);
		fprintf(stderr, "option: -s   print shape.agr file\n");
		fprintf(stderr, "option: -t   print atCathode.txt file\n");
		exit(-1);
	}

	FILE * fi=fopen(argv[1], "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n",argv[1]);
		perror(argv[1]);
		exit(-1);
	}


	FILE * ff = fopen("PART_RFQ.DST", "w");
	if(!ff)
	{
		fprintf(stderr, "cant open PART_RFQ.DST\n");
		exit(-1);
	}
	
        FILE * fh = fopen("hist.agr", "w");
        if(!fh)
        {
                fprintf(stderr, "cant open hist.txt\n");
                exit(-1);
        }


        memset(hist,  0, 501*sizeof(int));
        memset(histn, 0, 501*sizeof(int));
        memset(histx, 0, 501*sizeof(double));
        memset(histxp,0, 501*sizeof(double));





	char line[256];
	fgets(line, 255, fi);
	for(int i=0; i<3; i++) line[i]=tolower(line[i]);
	if (strncmp(line, "run", 3)  )
	{
		fprintf(stderr, "first line does not start with: run\n%s\n", line);
		exit(-1);
	}
	int i1,i2,i3;
	double freq, z0, w0;
	if( sscanf(line+3,"%d%d%lf%lf%lf%d", &i1, &i2, &freq, &z0, &w0, &i3)  != 6)
	{
		fprintf(stderr, "could not read run line:\n%s\n", line);
		exit(-1);
	}
	double v0=clight*sqrt(2.*w0/elmass);	// [m/s]
	double vdfreq = v0/(freq*360.e6);
	z0 *= 0.01;				// [m]
	printf(" found freq=%e [MHz]  z0 = %e[m]  w0=%e[MeV]  v0=%e[m/s]  vdfreq=%e [m/degree]\n", freq, z0, w0, v0, vdfreq);





	while( fgets(line, 255, fi))
	{
		char comnd[100];
		if( sscanf(line, "%s", comnd) != 1)  continue;	// line empty?
// 		printf("%s :: %s \n", comnd,line);
		
		if( !strcmp(comnd,"!@elypsoid") )
		{
			int nmax=0;
			double lcathode, radius, a,b,c,d,e,f,xpmax;


			sscanf( line+10, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &nmax, &lcathode, &radius, &a, &b, &c, &d, &e, &f, &xpmax);
			printf("elypsoid nmax %d lcathode %f radius %f \n",  nmax, lcathode, radius);
			printf("elypsoid space a %e b %e c %e \n",   a, b, c);
			printf("elypsoid density d %e e %e f %e \n",   d, e, f);
			printf("elypsoid xpmax  %f  \n", xpmax);



			double A = b;
			double B = (c - a)/(2.* lcathode);
			double C = (c + a - 2*b)/(2.*lcathode*lcathode);
			printf("A B C %e %e %e\n", A, B, C);

			double D = e;
			double E = (f - d)/(2.* lcathode);
			double F = (f + d - 2*e)/(2.*lcathode*lcathode);
			printf("D E F %e %e %e\n", D, E, F);
			double maxc = d; if(e > maxc) maxc=e; if(f > maxc) maxc=f; 



			double rmax=0.;
			for(double z=-lcathode; z < lcathode; z += lcathode/100)
			{
				double rz= radius * sqrt( 1. - (z*z/(lcathode*lcathode)))* (A + B*z + C*z*z);
				if(rmax < rz) rmax = rz;
			}



			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=6;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = rmax * ran[1];
				double y = rmax * ran[2];
				double z = lcathode * (2.*ran[3]-1.);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
			
				double charge = (D + E*z + F*z*z);
				if(ran[6]*maxc  > charge) continue;
			
				double rz= radius * sqrt( 1. - (z*z/(lcathode*lcathode)))* (A + B*z + C*z*z);
				double r2max = rz*rz;
					
				double r2=x*x+y*y;
				if(r2max >= r2) 
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printFour(ff, x, xp, y, yp, z, wr);
					i++;

				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);




		}

		
		if( !strcmp(comnd,"!@elypsoid2") )
		{
			int nmax=0;
			double lcathode, r1, r2,  a,b,xpmax;


			FILE * fe = fopen("ellipsoid2.txt", "w");
			if(!fe)
			{
				fprintf(stderr, "cant open ellipsoid2.txt\n");
				exit(-1);
			}

			int h =sscanf( line+11, "%d%lf%lf%lf%lf%lf%lf", &nmax, &lcathode, &r1, &r2, &a, &b, &xpmax);
			if(h != 7) 
			{
				fprintf(stderr," not 8 parameters\n");
				exit(-1);
			}
			printf("elypsoid nmax %d lcathode %f radius1 %f radius2 %f \n",  nmax, lcathode, r1, r2);
			printf("elypsoid xpmax  %f  \n", xpmax);
			printf("elypsoid space a %e b %e \n",   a, b);


			a = fabs(a);
			b = fabs(b);
			r1 = fabs(r1);
			r2 = fabs(r2);
			lcathode = fabs(lcathode);
			double linner=2.*lcathode-a-b;
			if(linner < 0.) 
			{
				double shrink = 2.*lcathode/(a+b);
				a *= shrink;
				b *= shrink;
				printf("shrink dimensions by %e\n", shrink);;
				linner=0.;
			}

			printf("elypsoid space a %e b %e li %e l %e\n",   a, b, linner, lcathode);

			double rmax = r1 < r2 ? r2 : r1;



			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);

                        n=7;
                        for(int i=0; i< 5000; i++)
                                sobseq(&n, ran);
			n=5;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = rmax * ran[1];
				double y = rmax * ran[2];
				double z = lcathode * (2.*ran[3]-1.);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
			
			
				double r2max;  // square radius of this 
				if( z > lcathode - b)
				{
					double z1 = z - lcathode + b;
					r2max = r2*r2 * ( 1. - (z1*z1/(b*b)));
				}
				else if( z > -lcathode + a)
				{
					double rz = r1 + (r2 - r1)*(z + lcathode - a)/linner;
					r2max = rz*rz;
				}
				else
				{
					double z1 = - z - lcathode + a;
					r2max = r1*r1 * ( 1. - (z1*z1/(a*a)));
				}
					
					
				double r2=x*x+y*y;
				double rp2=square(x*xp+y*yp)/r2;
				if(r2 <= r2max) 
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);
                                        int iz = int( (z + lcathode)/(2.*lcathode)*500);
                                        histn[iz] ++;
                                        histx[iz] += r2;
                                        histxp[iz] += rp2; 

// 					printFour(ff, x, xp, y, yp, z, wr);
					double dz= z0 + z*vdfreq;	// distance to cathode
					x += xp*dz;			// track backwards
					y += yp*dz;
					printFour(ff, x, xp, y, yp, z, wr);
					i++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);




                        for(int i=0; i<=500; i++)
                        {
                                if(histn[i] == 0) histn[i]=8888888;
                                fprintf(fh, "%d %e\n", i, (histxp[i]/histn[i]) );
                        }






                        ff = fopen("shape.agr", "w");
                        if(!ff)
                        {
                                fprintf(stderr, "cant open shape.agr\n");
                                exit(-1);
                        }


                        for(double z= - lcathode; z <= lcathode; z += 0.1)
                        {


                                double r2max;  // isquare radius of this
                                if( z > lcathode - b)
                                {
                                        double z1 = z - lcathode + b;
                                        r2max = r2*r2 * ( 1. - (z1*z1/(b*b)));
                                }
                                else if( z > -lcathode + a)
                                {
                                        double rz = r1 + (r2 - r1)*(z + lcathode - a)/linner;
                                        r2max = rz*rz;
                                }
                                else
                                {
                                        double z1 = - z - lcathode + a;
                                        r2max = r1*r1 * ( 1. - (z1*z1/(a*a)));
                                }

                                fprintf(ff,"%f %e\n", z, sqrt(r2max));
                        }
			printf("ellipsoid2 done\n");
		}
















		
		if( !strcmp(comnd,"!@beercan") )
		{
			int nmax=0;
			double lcathode, r1,  xpmax;


			int h =sscanf( line+10, "%d%lf%lf%lf", &nmax, &lcathode, &r1, &xpmax);
			if(h != 4) 
			{
				fprintf(stderr," not 4 parameters\n");
				exit(-1);
			}
			printf("beercan nmax %d lcathode %f radius1 %f  xpmax  %f \n",  nmax, lcathode, r1, xpmax);


			r1 = fabs(r1);
			lcathode = fabs(lcathode);


			double r2max=r1*r1;




			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=5;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = r1 * ran[1];
				double y = r1 * ran[2];
				double z = lcathode * (2.*ran[3]-1.);
				int iz = int( (z + lcathode)/(2.*lcathode)*500);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
				double c = ran[6];
				double cmax;
				double r2=x*x+y*y;
				if(r2 < r2max)
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printFour(ff, x, xp, y, yp, z, wr);
					i++;
					hist[iz]++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);





			for(int i=0; i<=500; i++) fprintf(fh, "%d %d\n", i, hist[i]);
		}




//			with xpmax in eV
		
		if( !strcmp(comnd,"!@beercan5") )
		{
			int nmax=0;
			double lcathode, r1,  xpmax;


			int h =sscanf( line+10, "%d%lf%lf%lf", &nmax, &lcathode, &r1, &xpmax);
			if(h != 4) 
			{
				fprintf(stderr," not 4 parameters\n");
				exit(-1);
			}
			printf("beercan nmax %d lcathode %f radius1 %f  xpmax  %f \n",  nmax, lcathode, r1, xpmax);

			xpmax /= 0.32810;   //emperical to make number in eV, this seems to be the only difference to !@beercan

			r1 = fabs(r1);
			lcathode = fabs(lcathode);


			double r2max=r1*r1;




			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=5;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = r1 * ran[1];
				double y = r1 * ran[2];
				double z = lcathode * (2.*ran[3]-1.);
				int iz = int( (z + lcathode)/(2.*lcathode)*500);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
				double c = ran[6];
				double cmax;
				double r2=x*x+y*y;
				if(r2 < r2max)
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printFour(ff, x, xp, y, yp, z, wr);
					i++;
					hist[iz]++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);





			for(int i=0; i<=500; i++) fprintf(fh, "%d %d\n", i, hist[i]);
		}



//			with xpmax in eV
		
		if( !strcmp(comnd,"!@beercancone") )
		{
			int nmax=0;
			double lcathode, r1,  xpmax, headDensity, tailDensity;


			int h =sscanf( line+14, "%d%lf%lf%lf%lf%lf", &nmax, &lcathode, &r1, &xpmax, &headDensity, &tailDensity);
			if(h != 6) 
			{
				fprintf(stderr," not 6 parameters, but %d\n", h);
				exit(-1);
			}
			printf("beercan nmax %d lcathode %f radius1 %f  xpmax  %f  head density %f   tail density %f\n",  nmax, lcathode, r1, xpmax, headDensity, tailDensity);

			xpmax /= 0.32810;   //emperical to make number in eV, this seems to be the only difference to !@beercan

			r1 = fabs(r1);
			lcathode = fabs(lcathode);
			double maxDens = tailDensity > headDensity ? tailDensity : headDensity;
			headDensity /= maxDens;
			tailDensity /= maxDens;


			double r2max=r1*r1;




			// write them all
			pcount = nmax /1;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=5;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = r1       * (2.*ran[1]-1.);
				double xp =xpmax    * (2.*ran[5]-1.);
				double y = r1       * (2.*ran[2]-1.);
				double yp =xpmax    * (2.*ran[4]-1.);
				double z = lcathode * (2.*ran[3]-1.);

// 				double x = r1       * (2.*double(random())/RAND_MAX-1.);
// 				double xp =xpmax    * (2.*double(random())/RAND_MAX-1.);
// 				double y = r1       * (2.*double(random())/RAND_MAX-1.);
// 				double yp =xpmax    * (2.*double(random())/RAND_MAX-1.);
// 				double z = lcathode * (2.*double(random())/RAND_MAX-1.);


				double zrel = (z + lcathode)/(2.*lcathode);	// longitudinal position from 0 to 1
				int iz = int( zrel*500);
				double densz = headDensity+ (tailDensity-headDensity)*zrel; // wanted density at z
				double d = double(random())/RAND_MAX;
				double r2=x*x+y*y;
// 				printf("z %f denz = %f d = %f\n", zrel, densz, d);
				if(r2 < r2max && d < densz)
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printOne(ff, x, xp, y, yp, z, wr);
					hist[iz]++;
					i++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);





			for(int i=0; i<=500; i++) fprintf(fh, "%d %d\n", i, hist[i]);
		}



		
		if( !strcmp(comnd,"!@gaussian") )
		{
			int nmax=0;
			double lcathode, r1,  xpmax;


			int h =sscanf( line+10, "%d%lf%lf%lf", &nmax, &lcathode, &r1, &xpmax);
			if(h != 4) 
			{
				fprintf(stderr," not 4 parameters\n");
				exit(-1);
			}
			printf("gaussian nmax %d lcathode %f radius1 %f  xpmax  %f \n",  nmax, lcathode, r1, xpmax);


			r1 = fabs(r1);
			lcathode = fabs(lcathode);


			double r2max=r1*r1;




			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=6;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = r1 * ran[1];
				double y = r1 * ran[2];

				double z = lcathode * (2.*ran[3]-1.);
				int iz = int( (z + lcathode)/(2.*lcathode)*500);

				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
				double c = ran[6];
				double cmax= exp(-z*z/(lcathode*lcathode/2.))/sqrt(M_PI*lcathode);
				double r2=x*x+y*y;
				if(c  <  cmax)
				if(r2 < r2max)
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printFour(ff, x, xp, y, yp, z, wr);
					i++;
					hist[iz]++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);




			for(int i=0; i<=500; i++) fprintf(fh, "%d %d\n", i, hist[i]);

		}



		
		if( !strcmp(comnd,"!@cylinder") )
		{
			int nmax=0;
			double lcathode, r1,  a,  xpmax;


			int h =sscanf( line+10, "%d%lf%lf%lf%lf", &nmax, &lcathode, &r1, &a, &xpmax);
			if(h != 5) 
			{
				fprintf(stderr," not 5 parameters\n");
				exit(-1);
			}
			printf("cylinder nmax %d lcathode %f radius1 %f  xpmax  %f \n",  nmax, lcathode, r1, xpmax);
			printf("cylinder space a %e\n",   a);


			a = fabs(a);
			r1 = fabs(r1);
			lcathode = fabs(lcathode);

			if(a > 1.85*lcathode) a =  1.85*lcathode;
			double b = 2*lcathode-a;

			printf("cylinder space a %e b %e  l %e\n",   a, b, lcathode);
			double r2max=r1*r1;




			// write them all
			pcount = nmax /4;
// 			double w0=0.4e-6;
			double wt2=0.;

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=6;

			for(int i=0; i< pcount; )
			{
				sobseq(&n, ran);
				double x = r1 * ran[1];
				double y = r1 * ran[2];
				double z = lcathode * (2.*ran[3]-1.);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
				double c = ran[6];
				double cmax;
				double r2=x*x+y*y;
				int iz = int( (z + lcathode)/(2.*lcathode)*500);
				if(r2 > r2max) continue;
			
			
				if( z > lcathode-b)
				{
					double z1 = z - lcathode + b;
					cmax =  1. - (z1*z1/(b*b));
				}
				else
				{
					double z1 = - z - lcathode + a;
					cmax =  1. - (z1*z1/(a*a));
				}
					
					
				if(c <= cmax) 
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printFour(ff, x, xp, y, yp, z, wr);
					i++;
					hist[iz]++;
				}
			}
			wt2 /= pcount;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);




			for(int i=0; i<=500; i++) fprintf(fh, "%d %d\n", i, hist[i]);


		}

		
		if( !strcmp(comnd,"!@ellipsoid4") )
		{
			int nmax=0;
			double lcathode, ra, rb,  pa, pb,xpmax;


			int h =sscanf( line+12, "%d%lf%lf%lf%lf%lf%lf", &nmax, &lcathode, &ra, &rb, &pa, &pb, &xpmax);
			if(h != 7) 
			{
				fprintf(stderr," not 8 parameters\n");
				exit(-1);
			}
			printf("ellipsoid3 nmax %d lcathode %f radius1 %f radius2 %f pa %e pb %e xpmax  %f\n",  nmax, lcathode, ra, rb,   pa, pb, xpmax);
//  The bunch has the length of +- lcathode and is made out of two ellipses with the center of (pa,0) and (pb,0). |pa| < lcathode, |pb < lcathode.
// the half axis of the elipses are therefor (lcathoder+pa, ra) and (lcathode-pb, rb).
//  we find the line that is tangent to both ellipses . it intersects at (xa,fa) and (xb,fb)
// 	if(nerror) if( fabs(xa-pa) >= la || fabs(xb-pb) >= lb) { *nerror=1; return 1e33;} else *nerror =0;


			xpmax /= 0.32810;   //emperical to make number in eV
			ra = fabs(ra);
			rb = fabs(rb);
			double la=lcathode+pa;
			double lb=lcathode-pb;
			lcathode = fabs(lcathode);
			if(pa > pb) { double t=pa; pa=pb; pb=t; }
			if(pa < -lcathode) pa= -lcathode*0.99;
			if(pa >  lcathode) pa=  lcathode*0.99;
			if(pb < -lcathode) pb= -lcathode*0.99;
			if(pb >  lcathode) pb=  lcathode*0.99;
			printf("ellipsoid4 nmax %d lcathode %f radius1 %f radius2 %f pa %e pb %e xpmax  %f\n",  nmax, lcathode, ra, rb,   pa, pb, xpmax);


			printf("ellipsoid4 space a %e b %e l %e\n",   pa, pb,  lcathode);

    			double rhoStart=1e-0, rhoEnd=1e-7;
    			int niter=1000;

    			Tangent *of = new Tangent(lcathode,pa,ra,pb,rb);

    			CONDOR(rhoStart, rhoEnd, niter, of);

			double fa=of->fa;
			double fb=of->fb;
			double xa=of->xa;
			double xb=of->xb;
    			delete of;


			printf("ellipsoid4 xa %e xb %e \n",   xa, xb);
			printf("ellipsoid4 fa %e fb %e \n",   fa, fb);
			double rmax = ra < rb ? rb : ra;

                        FILE * fs = fopen("shape.agr", "w");
                        if(!fs)
                        {
                                fprintf(stderr, "cant open shape.agr\n");
                                exit(-1);
                        }


			for(double z= -lcathode; z <= xa; z += 0.1)
			{
				double z1 =  z - pa;
				double r2max = ra*ra * ( 1. - (z1*z1/(la*la)));
				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
			}
			fprintf(fs, "&\n");

			for(double z= xa; z <= xb; z += 0.1)
			{
				double rz = fa + (fb - fa)/(xb - xa)*(z-xa);
				double r2max = rz*rz;
				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
			}
			fprintf(fs, "&\n");

			for(double z= xb; z <= lcathode; z += 0.1)
			{
				double z1 = z - pb;
				double r2max = rb*rb * ( 1. - (z1*z1/(lb*lb)));
				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
			}
			fprintf(fs, "&\n");


// 			for(double z= -lcathode; z <= lcathode; z += 0.1)
// 			{
// 				double r2max;  // square radius of this particle
// 				if( z > xb)
// 				{
// 					double z1 = z - pb;
// 					r2max = rb*rb * ( 1. - (z1*z1/(lb*lb)));
// 				}
// 				else if( z > xa)
// 				{
// 					double rz = fa + (fb - fa)/(xb - xa)*(z-xa);
// 					r2max = rz*rz;
// 				}
// 				else
// 				{
// 					double z1 =  z - pa;
// 					r2max = ra*ra * ( 1. - (z1*z1/(la*la)));
// 				}
// 				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
// 			}
			fclose(fs);

			// write them all
			double wt2=0.;		// transverse energy

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=5;

			// increase counter only if a new particle is found inside 
			for(int i=0; i< nmax; )
			{
				sobseq(&n, ran);
				double x = rmax * (2.*ran[1]-1.);
				double y = rmax * (2.*ran[2]-1.);
				double z = lcathode * (2.*ran[3]-1.);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
			
			
				double r2max;  // square radius of this particle
				if( z > xb)
				{
					double z1 = z - pb;
					r2max = rb*rb * ( 1. - (z1*z1/(lb*lb)));
				}
				else if( z > xa)
				{
					double rz = fa + (fb - fa)/(xb - xa)*(z-xa);
					r2max = rz*rz;
				}
				else
				{
					double z1 =  z - pa;
					r2max = ra*ra * ( 1. - (z1*z1/(la*la)));
				}
					
					
				double r2=x*x+y*y;
				if(r2 <= r2max) 
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					double dz= z0 + z*vdfreq;	// distance to cathode
					x -= xp*dz;			// track backwards
					y -= yp*dz;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

					printOne(ff, x, xp, y, yp, z, wr);
// 					fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  100.*x,  xp,   100.*y,  yp, z, wr);
					i++;
				}
			}
			wt2 /= nmax;
			printf("nmax = %d  wt = %f eV\n", nmax, sqrt(wt2)*1e6);




		}



//  ellipsoid5
//  The bunch has the length of +- lcathode and is made out of two ellipses with the center of (pa,0) and (pb,0). |pa| < lcathode, |pb < lcathode.
// the half axis of the elipses are therefor (lcathoder+pa, ra) and (lcathode-pb, rb).
// In order to force a magnetization for magnetized beams we specify not ra and rb, but rm=3./2.*sqrt( int r ds)/(2*lcathode)  ) and rrat=rb/ra.
// ( if the beam is a pure ellipse rm is the maximum radius)
//  we find the line that is tangent to both ellipses . it intersects at (xa,fa) and (xb,fb)

		
		if( !strcmp(comnd,"!@ellipsoid5") )
		{
			int nmax=0;
			double lcathode, rm, rrat,  pa, pb,xpmax, ra, rb;


			FILE * fe = fopen("ellipsoid5.txt", "w");
			if(!fe)
			{
				fprintf(stderr, "cant open ellipsoid5.txt\n");
				exit(-1);
			}

			int h =sscanf( line+12, "%d%lf%lf%lf%lf%lf%lf", &nmax, &lcathode, &rm, &rrat, &pa, &pb, &xpmax);
			if(h != 7) 
			{
				fprintf(stderr," not 8 parameters\n");
				exit(-1);
			}
			printf("ellipsoid5 nmax %d lcathode %f radius %f rrat %f pa %e pb %e xpmax  %f\n",  nmax, lcathode, rm, rrat,   pa, pb, xpmax);
// 

			xpmax /= 0.32810;   //emperical to make number in eV
			rm = fabs(rm);
			rrat = fabs(rrat);
			lcathode = fabs(lcathode);

			double zmax = - z0/vdfreq;
			printf("zmax %f lcathode %f\n", zmax, lcathode);
			if(lcathode > zmax)
			{
				fprintf(stderr, " *********************  bunch is too long, change z0 ****************\n");
				fprintf(stderr, " *********************  bunch is too long, change z0 ****************\n");
				fprintf(stderr, " *********************  bunch is too long, change z0 ****************\n");
				fprintf(stderr, " *********************  bunch is too long, change z0 ****************\n");
				exit(-1);
			}





//			we assume first that ra = 1 and rb=rrat, later we scale.
			ra=1.;
			rb=rrat;


// 	printf(" found freq=%e [MHz]  z0 = %e[m]  w0=%e[MeV]  v0=%e[m/s]  vdfreq=%e [m/degree]\n", freq, z0, w0, v0, vdfreq);


			if(pa < -lcathode*0.99) pa= -lcathode*0.99;
			if(pa >  lcathode*0.99) pa=  lcathode*0.99;
			if(pb < -lcathode*0.99) pb= -lcathode*0.99;
			if(pb >  lcathode*0.99) pb=  lcathode*0.99;
			if(pa > pb) { double t=pa; pa=pb; pb=t; }

			double la=lcathode+pa;
			double lb=lcathode-pb;
			printf("ellipsoid5 nmax %d lcathode %f radius1 %f radius2 %f pa %e pb %e xpmax  %f\n",  nmax, lcathode, ra, rb,   pa, pb, xpmax);



    			double rhoStart=1e-0, rhoEnd=1e-7;
    			int niter=1000;

    			Tangent *of = new Tangent(lcathode,pa,ra,pb,rb);

    			CONDOR(rhoStart, rhoEnd, niter, of);

			double xa=of->xa;
			double xb=of->xb;
			double fa=of->fa;
			double fb=of->fb;
			double fap=of->fap;
			double fbp=of->fbp;
			double diffRadius=of->diffRadius;
			double diffSlope=of->diffSlope;
    			delete of;

			printf("ellipsoid5 xa %e xb %e \n",   xa, xb);
			printf("ellipsoid5 fa %e fb %e \n",   fa, fb);
			if( fabs(diffRadius)  > 1e-6 )
			{
				// the line should go through both points
				double correction = (fa + fap*(xb-xa))/fb;
				fb *= correction;
				rb *= correction;
				rrat *= correction;
				printf("ellipsoid5 fa %e fb %e rb %e correction\n",   fa, fb, rb);
			}

			double a= -lcathode-pa;
			double b= xa -pa;
			double rmm1 = b - a - (qube(b) - qube(a))/(3.*square(la));

			double rmm2 = (fa*fb + fa*fa + fb*fb) * (xb-xa) / 3.;


			a= xb-pb;
			b=lcathode -pb;
			double rmm3 = square(rrat) * ( b - a - (qube(b) - qube(a))/(3.*square(lb)) );

			double rmm = (rmm1 + rmm2 + rmm3) /(2.*lcathode);

			printf("ellipsoid5  rmm1 %e   rmm2 %e rmm3 %e  rmm %e\n", rmm1, rmm2, rmm3, rmm);

			double rm2 = 2./3.*rm*rm;
			double fact = sqrt(rm2/rmm);

			printf("ellipsoid5  rm %e   rmm %e fact %e\n", rm, rmm, fact);

			// sclale all r values
			ra *= fact;
			rb *= fact;
			fa *= fact;
			fb *= fact;
			double rmax = ra < rb ? rb : ra;


			printf("ellipsoid5  ra %e   rb %e fa %e  fb %e  rmax %e \n", ra, rb, fa, fb,  rmax);



// 	plot the shape
                        FILE * fs = fopen("shape.agr", "w");
                        if(!fs)
                        {
                                fprintf(stderr, "cant open shape.agr\n");
                                exit(-1);
                        }


			for(double z= -lcathode; z <= xa; z += 0.1)
			{
				double z1 =  z - pa;
				double r2max = ra*ra * ( 1. - (z1*z1/(la*la)));
				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
			}
			fprintf(fs, "&\n");

			fprintf(fs, "%e  %e\n", xa, fa);
			fprintf(fs, "%e  %e\n", xb, fb);
			fprintf(fs, "&\n");

			for(double z= xb; z <= lcathode; z += 0.1)
			{
				double z1 = z - pb;
				double r2max = rb*rb * ( 1. - (z1*z1/(lb*lb)));
				fprintf(fs, "%e  %e\n", z, sqrt(r2max));
			}
			fprintf(fs, "&\n");

			fclose(fs);

// write them all
			double wt2=0.;		// transverse energy
			double M=0.;		// magnetization

			double ran[7];
			int n=  -1;
			sobseq(&n, ran);
			n=5;

			// increase counter only if a new particle is found inside 
			for(int i=0; i< nmax; )
			{
				sobseq(&n, ran);
				double x = rmax * (2.*ran[1]-1.);
				double y = rmax * (2.*ran[2]-1.);
				double z = lcathode * (2.*ran[3]-1.);
				double xp =  (2.*ran[4]-1.)*xpmax;
				double yp =  (2.*ran[5]-1.)*xpmax;
			
			
				double r2max;  // square radius of this particle
				if( z > xb)
				{
					double z1 = z - pb;
					r2max = rb*rb * ( 1. - (z1*z1/(lb*lb)));
				}
				else if( z > xa)
				{
					double rz = fa + (fb - fa)/(xb - xa)*(z-xa);
					r2max = rz*rz;
				}
				else
				{
					double z1 =  z - pa;
					r2max = ra*ra * ( 1. - (z1*z1/(la*la)));
				}
					
					
				double r2=x*x+y*y;
				if(r2 <= r2max) 
				{
					double ww = sqrt(1.+xp*xp+yp*yp);
					xp /= ww;
					yp /= ww;
					double wr=w0*ww-w0;
					wt2 += w0*w0*ww*ww*(xp*xp+yp*yp);

// 					fprintf(fe, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  100.*x,  xp,   100.*y,  yp, z, wr);
					double dz= z0 + z*vdfreq;	// distance to cathode
					x -= xp*dz;			// track backwards
					y -= yp*dz;

					M += r2;
					printOne(ff, x, xp, y, yp, z, wr);
// 					fprintf(ff, "%17.9e  %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n",  100.*x,  xp,   100.*y,  yp, z, wr);
					i++;
				}
			}
			wt2 /= nmax;
			M /= nmax;
			printf("ellipsoid5 nmax = %d  wt = %f eV  M=%e\n", nmax, sqrt(wt2)*1e6, M);
			printf("ellipsoid5 done\n");




			fclose(fe);
		}






	}
	fclose(ff);
	fclose(fh);
}





