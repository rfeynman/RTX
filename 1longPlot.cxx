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
		fprintf(stderr, "usage: %s what n0 n1 [n2 ...]  \n", argv[0]);
		exit(-1);
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
		if(ne >= p->nplots) break;
		printf("%d \n", ne);
		if(p->getElement(ne))
		{
			fprintf(stderr, "no element %d\n", n1);
			exit(-1);
		}
		p->sortZ('7');

		double pavg = 0;
		for(int i=0; i<nbuf; i++)
		{
			if(p->scord[i][6] != scord0[i][6])
			{
				fprintf(stderr, "mismatch particle %d %e %e\n", i, p->scord[i][6], scord0[i][6]);
				exit(-1);
			}
			p->scord[i][6] = scord0[i][4];
			pavg += p->scord[i][5];
		}
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


		sprintf(filename, "longHist%04d.agr", ne);
		FILE * fh = fopen(filename, "w");
		if(!fh)
		{
			perror(filename);
			exit(-1);
		}



		const int nhist=nbuf/30;
		int hist[nhist];
		for(int i=0; i< nhist; i++) hist[i]=0;
		double dl = (p->scord[nbuf-1][4] - p->scord[0][4])/nhist;
		int nsets=0;
		mkdots(fp, nsets++);
		for(int i=0; i<nbuf; i++)
		{
			if( argv[1][0] == 'l' ) fprintf(fp, "%e %e\n", p->scord[i][4], (p->scord[i][5]-pavg)/pavg);
			if( argv[1][0] == 'x' ) fprintf(fp, "%e %e\n", p->scord[i][0], p->scord[i][1]);
			if( argv[1][0] == 'y' ) fprintf(fp, "%e %e\n", p->scord[i][2], p->scord[i][3]);
			if( argv[1][0] == 'r' )
			{
				double r =  sqrt( p->scord[i][0]*p->scord[i][0] + p->scord[i][2]*p->scord[i][2]  );
				fprintf(fp, "%e %e\n", r, p->scord[i][0]*p->scord[i][1] + p->scord[i][2]*p->scord[i][3] );
			}
			if(((i+1)%(nbuf/20)) == 0)
			{
				fprintf(fp, "\n");
				mkdots(fp, nsets++);
			}
			int ihist = (p->scord[i][4]- p->scord[0][4])/dl;
			hist[ihist]++;
		}
		for(int i=0; i<nhist; i++)
			fprintf(fh, "%e %d\n", i*dl+dl/2., hist[i]);
		fclose(fp);
		fclose(fh);


	}

}

