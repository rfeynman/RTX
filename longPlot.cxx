// make plots of the longitudinal phase space from parmela tape2 file for the grace (xmgr) plotting packege.
// usage  lonplot [l,x,y] <start element (typically 1)> <element> [<element>......]

//  Jorg Kewisch, BNL, 2004-2014
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
#include "pmla.hxx"
#include "xmgr_class.hxx"

const double energ = 55*1e6;    // energy in eV
const double clight = 2.997924562e8;		// in m/s
const double elmass = 0.51104;

inline double square(double x) { return x*x; }

void printmat(double xx[4][4])
{
	printf("\n");
	for(int i =0; i < 4; i++)
	{
		for(int j =0; j < 4; j++)
		{
			printf("%f  ", xx[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void mkdots(FILE * fp, int i)
{
	fprintf(fp, "@with g0\n");
	fprintf(fp, "@    s%d symbol 9\n",i);
	fprintf(fp, "@    s%d symbol size 0.060000\n",i);
	fprintf(fp, "@    s%d line type 0\n",i);
}


main(int argc, char ** argv)
{
	char cmd[200];
	char path[200];
	char ext[20];
	ext[0]=0;

	char cwd[200];
	getcwd(cwd, 200);

	int cutWhat=0;
	double cut=-1.; 
// 	printf("cut  (0.99)? "); scanf("%lf", &cut);
	if(argc < 4)
	{
		fprintf(stderr, "usage: %s [-c] what n0 n1 [n2 ...]  \n", argv[0]);
		exit(-1);
	}
		
	int doColor=1;
	if( ! strcmp(argv[1], "-c")  )
	{
		doColor = 0;
		argv++; argc --;
	}




	char elemfilename[100];
	strcpy(elemfilename,"tape2.t2");
	strcat(elemfilename,ext);
	char timefilename[100];
	strcpy(timefilename,"tape3.t3");
	strcat(timefilename,ext);
	Pmla *p = new Pmla(elemfilename, timefilename);

	int n1 = atoi(argv[2]);
	if(p->getElement(n1))
	{
		fprintf(stderr, "no element %d\n", n1);
		exit(-1);
	}
	//  sort the reference
	p->sortZ('7');
	double **scord0;
	int nbuf = p->nbuf;
	scord0 = new double*[nbuf];
	for(int i=0; i<nbuf; i++)
	{
		scord0[i] = new double[7];
		memcpy(scord0[i], p->scord[i], 7*sizeof(double));
	}


	for(int n =3; n < argc; n++)
	{
		int ne = atoi(argv[n]);
		if(ne >= p->nplots)
		{
			fprintf(stderr, "number %d larger than nplots = %d\n", ne, p->nplots);
			break;
		}
		if(p->getElement(ne))
		{
			fprintf(stderr, "no element %d\n", n1);
			exit(-1);
		}
		p->sortZ('7');

		if( argv[1][0] == 'L' )		// try to fit a thin lense 1st and 3rd harmonic cavity, does not work so far
		{
			p->Pmla::longFit(p->scord, nbuf);
		}


		double pavg = 0;
		double lavg = 0;
		double xmin= 1e33;
		double xmax=-1e33;
		double ymin= 1e33;
		double ymax=-1e33;
		double zmin= 1e33;
		double zmax=-1e33;
		double pmin= 1e33;
		double pmax=-1e33;
		double rmax=-1e33;
		for(int i=0; i<nbuf; i++)
		{
			if(p->scord[i][6] != scord0[i][6])
			{
				fprintf(stderr, "mismatch particle %d %e %e\n", i, p->scord[i][6], scord0[i][6]);
				exit(-1);
			}
			p->scord[i][6] = scord0[i][4];
			lavg += p->scord[i][4];
			pavg += p->scord[i][5];
			double r = sqrt(  square(p->scord[i][0]) + square(p->scord[i][2]));
			if(rmax < r) rmax = r;
			if(xmin > p->scord[i][0]) xmin = p->scord[i][0];
			if(xmax < p->scord[i][0]) xmax = p->scord[i][0];
			if(ymin > p->scord[i][2]) ymin = p->scord[i][2];
			if(ymax < p->scord[i][2]) ymax = p->scord[i][2];
			if(zmin > p->scord[i][4]) zmin = p->scord[i][4];
			if(zmax < p->scord[i][4]) zmax = p->scord[i][4];
			if(pmin > p->scord[i][5]) pmin = p->scord[i][5];
			if(pmax < p->scord[i][5]) pmax = p->scord[i][5];
		}
		lavg /= nbuf;
		pavg /= nbuf;

		p->sortZ('7');
		char filename[100];
		sprintf(filename, "longPhase%04d.agr", ne);
		FILE * fp = fopen(filename, "w");
		if(!fp)
		{
			perror(filename);
			exit(-1);
		}
        	fprintf(fp, "@with g0\n");
		sprintf(cmd, "%s, s=%f", cwd, p->zloc[ne]);
        	fprintf(fp, "@subtitle \"%s\"\n", cmd); 
		if( argv[1][0] == 'l' || argv[1][0] == 'L' ) 
		{
        		fprintf(fp, "@xaxis label \"RF Phase\"\n");
        		fprintf(fp, "@yaxis label \"Delta P/P0\"\n");
		}
		if( argv[1][0] == 'x' ) 
		{
        		fprintf(fp, "@xaxis label \"x'\"\n");
        		fprintf(fp, "@yaxis label \"x [m]\"\n");
		}
		if( argv[1][0] == 'y' )
		{
        		fprintf(fp, "@xaxis label \"y'\"\n");
        		fprintf(fp, "@yaxis label \"y [m]\"\n");
		}
		if( argv[1][0] == 'r' )
		{
        		fprintf(fp, "@xaxis label \"r'\"\n");
        		fprintf(fp, "@yaxis label \"r [m]\"\n");
		}
        	fprintf(fp, "@hardcopy device \"PNG\"\n");
        	fprintf(fp, "@device \"PNG\" DPI 1200\n");






		const int nhist=nbuf/30;
		double hist[nhist+1];
		double hmin[nhist+1];
		double hmax[nhist+1];
		for(int i=0; i<= nhist; i++)
		{
			hist[i]=0;
			hmin[i]= 1e33;
			hmax[i]=-1e33;
		}
		double dl;
		if( argv[1][0] == 'l' || argv[1][0] == 'L' ) dl = (zmax-zmin)/nhist;
		if( argv[1][0] == 'p' ) dl = (pmax-pmin)/nhist;
		if( argv[1][0] == 'x' ) dl = (xmax-xmin)/nhist;
		if( argv[1][0] == 'y' ) dl = (ymax-ymin)/nhist;
		if( argv[1][0] == 'r' ) dl = (rmax)/nhist;  // rmin=0
// 		if( argv[1][0] == 'l' || argv[1][0] == 'L' )
			printf("zmin %e %e %e \n", zmin, zmax, dl);
// 		if( argv[1][0] == 'x' )
			printf("xmin %e %e %e \n", xmin, xmax, dl);
// 		if( argv[1][0] == 'y' )
			printf("ymin %e %e %e \n", ymin, ymax, dl);
// 		if( argv[1][0] == 'r' )
			printf("rmin %e %e %e \n",   0., rmax, dl);
// 		if( argv[1][0] == 'p' )
			printf("pmin %e %e %e \n", pmin, pmax, dl);




		int nsets=0;
		mkdots(fp, nsets++);
		for(int i=0; i<nbuf; i++)
		{
			double x, y;
// 			if( argv[1][0] == 'l' || argv[1][0] == 'L' || argv[1][0] == 'p' ) {x = p->scord[i][4]-lavg; y = (p->scord[i][5]-pavg)/pavg; }
			if( argv[1][0] == 'l' || argv[1][0] == 'L' || argv[1][0] == 'p' ) {x = p->scord[i][4]-lavg; y = p->scord[i][5]; }
			if( argv[1][0] == 'x' ) { x = p->scord[i][0]; y =  p->scord[i][1]; }
			if( argv[1][0] == 'y' ) { x = p->scord[i][2]; y =  p->scord[i][3]; }
			if( argv[1][0] == 'r' )
			{
				double r =  sqrt( p->scord[i][0]*p->scord[i][0] + p->scord[i][2]*p->scord[i][2]  );
				double rp = sqrt( p->scord[i][0]*p->scord[i][1] + p->scord[i][2]*p->scord[i][3]  );
				x = r; y = rp;
			}
			fprintf(fp, "%e %e\n", x, y);
			if(doColor  &&  (((i+1)%(nbuf/20)) == 0)  )
			{
				fprintf(fp, "\n");
				mkdots(fp, nsets++);
			}
			int ihist;
			if( argv[1][0] == 'l' || argv[1][0] == 'L' ) ihist = (p->scord[i][4]- zmin)/dl;
			if( argv[1][0] == 'x' ) ihist = (p->scord[i][0]- xmin)/dl;
			if( argv[1][0] == 'y' ) ihist = (p->scord[i][2]- ymin)/dl;
			if( argv[1][0] == 'p' ) ihist = (p->scord[i][5]- pmin)/dl;
			if( argv[1][0] == 'r' )
			{
				double r =  sqrt( p->scord[i][0]*p->scord[i][0] + p->scord[i][2]*p->scord[i][2]  );
				ihist = r/dl;
			}

			if(ihist >= nhist)  ihist = nhist-1;
			if(ihist < 0) ihist=0;
			if(hmin[ihist] > y ) hmin[ihist] = y;
			if(hmax[ihist] < y ) hmax[ihist] = y;
			hist[ihist]++;
		}
		fclose(fp);
		strcpy(cmd, "gracebat ");
		strcat(cmd, filename);
		system(cmd);


		sprintf(filename, "longHist%c%04d.agr", argv[1][0], ne);
		FILE * fh = fopen(filename, "w");
		if(!fh)
		{
			perror(filename);
			exit(-1);
		}
        	fprintf(fp, "@with g0\n");
		sprintf(cmd, "%s, s=%f", cwd, p->zloc[ne]);
        	fprintf(fp, "@subtitle \"%s\"\n", cmd); 
        	fprintf(fp, "@yaxis label \"Density\"\n");
		if( argv[1][0] == 'l' || argv[1][0] == 'L' ) 
		{
        		fprintf(fp, "@xaxis label \"RF Phase\"\n");
		}
		if( argv[1][0] == 'x' ) 
		{
        		fprintf(fp, "@xaxis label \"x [m]\"\n");
		}
		if( argv[1][0] == 'y' )
		{
        		fprintf(fp, "@xaxis label \"y [m]\"\n");
		}
		if( argv[1][0] == 'r' )
		{
        		fprintf(fp, "@xaxis label \"r [m]\"\n");
		}
		if( argv[1][0] == 'p' )
		{
        		fprintf(fp, "@xaxis label \"P [MeV]\"\n");
		}
        	fprintf(fh, "@with g0\n");
        	fprintf(fh, "@subtitle \"%s\"\n", cwd); 
        	fprintf(fh, "@hardcopy device \"PNG\"\n");
        	fprintf(fh, "@device \"PNG\" DPI 1200\n");


		for(int h=0; h<=nhist; h++)
		{
			if(hist[h] == 0.) hist[h] =1e-3;
			if( argv[1][0] == 'l' || argv[1][0] == 'L' ) fprintf(fh, "%e %e\n", double(h)*dl+dl/2.-lavg,  hist[h]);
			if( argv[1][0] == 'x' ) fprintf(fh, "%e %e\n", double(h)*dl+dl/2.+ xmin, hist[h]);
			if( argv[1][0] == 'y' ) fprintf(fh, "%e %e\n", double(h)*dl+dl/2.+ ymin, hist[h]);
			if( argv[1][0] == 'r' ) fprintf(fh, "%e %e\n", double(h)*dl+dl/2.,       hist[h]/(double(h)*dl+dl/2.));
			if( argv[1][0] == 'p' ) fprintf(fh, "%e %e\n", double(h)*dl+dl/2.+ pmin, hist[h]);
		}
		fclose(fh);
		strcpy(cmd, "gracebat ");
		strcat(cmd, filename);
		system(cmd);




		sprintf(filename, "longSpread%04d.agr", ne);
		FILE * fs = fopen(filename, "w");
		if(!fs)
		{
			perror(filename);
			exit(-1);
		}
        	fprintf(fh, "@with g0\n");
        	fprintf(fh, "@subtitle \"%s\"\n", cwd); 

		for(int h=0; h<=nhist; h++)
		{
			if( argv[1][0] == 'l' || argv[1][0] == 'L' ) fprintf(fs, "%e %e\n", double(h)*dl+dl/2.-lavg,  hmax[h]-hmin[h]);
			if( argv[1][0] == 'x' ) fprintf(fs, "%e %e\n", double(h)*dl+dl/2.+ xmin, hmax[h]-hmin[h]);
			if( argv[1][0] == 'y' ) fprintf(fs, "%e %e\n", double(h)*dl+dl/2.+ ymin, hmax[h]-hmin[h]);
			if( argv[1][0] == 'r' ) fprintf(fs, "%e %e\n", double(h)*dl+dl/2.,       (hmax[h]-hmin[h])/(double(h)*dl+dl/2.));
		}
		fclose(fs);
		strcpy(cmd, "gracebat ");
		strcat(cmd, filename);
		system(cmd);



		printf("done %d\n", ne);


	}

}

