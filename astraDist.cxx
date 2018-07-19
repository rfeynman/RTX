// make plots from parmela tape2 file for the grace (xmgr) plotting packege.
//  Jorg Kewisch, BNL, 2004-2006
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

const double clight = 2.997924562e8;		// in m/s
const double elmass = 0.51104;





int main(int argc, char ** argv)
{
	char path[200];
	char ext[20];
	ext[0]=0;

	char cwd[200];
	getcwd(cwd, 200);

	int cutWhat=0;
	double cut=-1.; 
// 	printf("cut  (0.99)? "); scanf("%lf", &cut);
// 	if(argc < 4)
// 	{
// 		fprintf(stderr, "usage: %s what n0 n1 [n2 ...]  \n", argv[0]);
// 		exit(-1);
// 	}
		



	char elemfilename[100];
	strcpy(elemfilename,"tape2.t2");
	strcat(elemfilename,ext);
	char timefilename[100];
	strcpy(timefilename,"tape3.t3");
	strcat(timefilename,ext);
	Pmla *p = new Pmla(elemfilename, timefilename);

	int n1=-1;
	int n2=-1;

	for(int i = 0; i< 5; i++)
	{
		printf("zloc %d %e\n", i, p->zloc[i]);
		if(p->zloc[i] > 0.)
		{
			n2 = i; n1 = i-1;
			break;
		}
	}
	if(n1 < 0)
	{
		fprintf(stderr, "n2 not found\n");
		exit(-1);
	}
	//  sort the reference
	printf(" n1 =%d n2=%d\n", n1, n2);


	if(p->getElement(n1))
	{
		fprintf(stderr, "no element %d\n", n1);
		exit(-1);
	}
	//  sort the reference
	p->sortZ('7');
	double **scord0;
	int nbuf0 = p->nbuf;
	scord0 = new double*[nbuf0];
	for(int i=0; i<nbuf0; i++)
	{
		scord0[i] = new double[7];
		memcpy(scord0[i], p->scord[i], 7*sizeof(double));
	}
	double zloc0 = p->zloc[n1];


	if(p->getElement(2))
	{
		fprintf(stderr, "no element %d\n", 2);
		exit(-1);
	}
	p->sortZ('7');
	double **scord = p->scord;
	int nbuf = p->nbuf;
	if(nbuf0 != nbuf)
	{
		fprintf(stderr, "nbuf mismatch\n");
		exit(-1);
	}

	double savg = 0.;
	double smax = -1e33;
	double smin =  1e33;
	for(int i=0; i<nbuf; i++)
	{
		if(p->scord[i][6] != scord0[i][6])
		{
			fprintf(stderr, "mismatch particle %d %e %e\n", i, p->scord[i][6], scord0[i][6]);
			exit(-1);
		}
		savg += scord0[i][4];
		if(scord0[i][4] > smax) smax = scord0[i][4];
		if(scord0[i][4] < smin) smin = scord0[i][4];
	}
	savg /= nbuf;
	printf("savg %e %e %e\n", savg, smax, smin);
	double zloc = p->zloc[n2];
	double deriv = 1./(zloc-zloc0);
	printf("zloc %e %e %e\n", zloc0, zloc, deriv);

	char filename[100];
	sprintf(filename, "astra-init");
	FILE * fp = fopen(filename, "w");
	if(!fp)
	{
		perror(filename);
		exit(-1);
	}
	for(int i=0; i<nbuf; i++)
	{
		double pz = sqrt(scord[i][5]/(2.*elmass))*1e6;
		double px = scord[i][1]*pz;
		double py = scord[i][3]*pz;
		double t = (scord[i][4]-savg)*1e9/(p->freq*1e6*360.);
		
		fprintf(fp, "%15e %15e %15e %15e %15e %15e %15e -3e-5 1 -1\n", 
				scord0[i][0], scord0[i][2], 0., px, py, pz, t);
	}


	fclose(fp);


	FILE * fi = fopen("init_dist", "r");
	if(!fi)
	{
		perror("init_dist");
		exit(-1);
	}

	FILE * fo = fopen("init_dist_form", "w");
	if(!fo)
	{
		perror("init_dist_form");
		exit(-1);
	}
	double x,y,z,px,py,pz,t,c;
	int f1,f2;
	while( fscanf(fi, "%le%le%le%le%le%le%le%le%d%d", &x,&y,&z,&px,&py,&pz,&t,&c,&f1,&f2) == 10)
		fprintf(fo, "%15e %15e %15e %15e %15e %15e %15e %15e %4d %4d\n", 
				x,y,z,px,py,pz,t,c,f1,f2);
	fclose(fi);
	fclose(fo);

}

