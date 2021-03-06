#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

#include "pmla.hxx"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "Solver.hxx"
#include "tools.hxx"
#include "ObjectiveFunction.hxx"

#ifdef WIN32
#include <crtdbg.h>
#endif




pthread_mutex_t rtxMutex= PTHREAD_MUTEX_INITIALIZER;
inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }

const double clight=2.997824562e8;      // m/s
const double elcharge=1.6022e-19;       // Coulomb
const double elmass = 0.51104;

const double cut = 0.1;



double mu0 = 4.*M_PI*1.e-7;     // V*s / (A*m)
double eps0 = 1./(mu0*clight*clight); //  A*s / (V*m)
double I0= 17045.;      //      Alfven current ampere
// double I0= 4*M_PI*eps0*elmass*clight^3 / elcharge    // en.wikipedia.org/wiki/Perveance
// double elrad = elcharge^2 / (4*M_PI*eps0*elmass*clight^2) = elcharge*clight /I0

const double magnetSlice=0.01;
const int nBunchSlices = 20;

double uu[49] = { 1., 0., 0., 0., 0., 0., 0.,
	          0., 1., 0., 0., 0., 0., 0.,
		  0., 0., 1., 0., 0., 0., 0.,
		  0., 0., 0., 1., 0., 0., 0.,
		  0., 0., 0., 0., 1., 0., 0.,
		  0., 0., 0., 0., 0., 1., 0.,
		  0., 0., 0., 0., 0., 0., 1.};


class  Opttrack : public UnconstrainedObjectiveFunction
{
public:
	double  current[nBunchSlices], bunchlength[nBunchSlices],energy[nBunchSlices], IoverI0[nBunchSlices];
	double gammaL[nBunchSlices]; 
	double pmoment[nBunchSlices];
	double totMoment;
	double sigma[nBunchSlices][7][7];	// sigma matrix
	double sigmaSave[7][7];	// sigma matrix
	double mat[nBunchSlices][7][7];	// element matrix
	double matsp[7][7];	// space charge kick matriX
	double beami;
	double freq;
	double qbunch;
	double spcfac;
	int doSpaceCharge;
	int doSurv;
	int zeroDisp;
	int type[500];
	double  length[500], strength[500], sloc[500];
	int n;
	double locstart;
	int nfit, nparm,ngoal;
	int fitElem[100];
	int fitParm[100];
	double fitSave[100];
	double fitBest[100];
	double fitStep[100];
	double fitGoal[100];
	FILE * fs;
	FILE * fey;
	FILE * fex;
	FILE * fdy;
	FILE * fdx;
	FILE * fyp;
	FILE * fxp;
	FILE * fy;
	FILE * fx;
	double bestgoal;

	Opttrack(const char * file, double spcfac, int doSpaceCharge, int doSurv, int zeroDisp);
	virtual ~Opttrack();
	double track();
	void plot();
	void propagate(int i, double tmoment);
    	double eval(Vector v, int *nerror, int resource=0);
	void installBest();
	void saveBest();
};

Opttrack::~Opttrack()
{
	fclose(fx);
	fclose(fy);
	fclose(fs);
	fclose(fxp);
	fclose(fyp);
	fclose(fex);
	fclose(fey);
	fclose(fdx);
	fclose(fdy);
	fclose(fs);
}


Opttrack::Opttrack(const char * file, double spcfac, int doSpaceCharge, int doSurv, int zeroDisp)
{
	this->spcfac = spcfac;
	this->doSpaceCharge = doSpaceCharge;
	this->doSurv = doSurv;
	this->zeroDisp = zeroDisp;
	FILE * fi = fopen(file,"r");
	if(!fi)
	{
		perror(file);
		exit(-1);
	}





	Pmla *p = new Pmla("tape2.t2", "tape3.t3");

	beami = p->beami;
	freq  = p->freq*1e6;
	qbunch    = beami/freq;



	n=0;
// 	type = 0    	eof
// 	type = 1	drift
// 	type = 2	quad
//	type = 3	solenoid 
//	type = 4	bend
//	type = 5	edge
//	type = 6	rotate
//	length and fields are converted from cm and gauss to m and kG when read from file
//	energies are in MeV


	FILE * fm = fopen("slenv.madx","w");
	if(!fm)
	{
		perror("slenv.mad");
		exit(-1);
	}

	int driftnum=0;
	int quadnum=0;
	int bendnum=0;
	int solnum=0;
	int cavnum=0;
	char elname[1000][6];
	int elnum=0;

	nfit=0;
	nparm=0;
	int neo =0;
	int save=0;
	double loc =0;
	double len, app, pri, stre, angl, e1, e2, phas;
	char line[1255];
	char word[255];
	while (fgets(line, 255, fi))
	{
		for(char * l=line; *l; l++)
		{
			if( *l >= 'A' && *l <= 'Z') *l += 32;
			if( *l == ',') *l=' ';
		}
		int nw = sscanf(line,"%s", word);
		if(nw != 1) continue;


		if( !strcmp( word, "!@slenvstart"))
		{
#if 1
			printf("start at neo = %d\n", neo);
			if(neo >= p->nplots)
			{
				fprintf(stderr, "slenvstart at %d >= number of elemets in tape2 %d\n", neo,  p->nplots);
				exit(-1);
			}
			p->getElement(neo);
			int nbuf = p->nbuf;
			int nc=int(nbuf*cut/nBunchSlices);	// number of electrons to cut from each ed
			p->sortZ();
			totMoment=0.;
			for(int sl =0; sl < nBunchSlices; sl++)
			{
				double **scord = p->scord + sl*nc;
				int nbufsl = nc;	// number of electrons in this slice
				if(sl*nc+nbufsl > nbuf) nbufsl = nbuf - sl*nc;	// last interval may be shorter
				double totBuLen = p->scord[0][6] - p->scord[nbuf-1][6] ;
				p->calcSigmaMatrix(scord, nbufsl);
				memcpy(sigma[sl], p->sigmaMat2, 49*sizeof(double));
				sigma[sl][6][6] =1.;
				pmoment[sl] =  p->sigmaMat[6][5]; 
				totMoment += pmoment[sl] * nbufsl;
				bunchlength[sl] = (p->maxPart[4] - p->minPart[4]) /360. *clight/freq;
				IoverI0= qbunch*clight/(bunchlength[sl]*I0*nBunchSlices);
				if(zeroDisp)
				{
					for(int i=0; i<4; i++)
						sigma[sl][i][5] = sigma[sl][5][i] =0.;
					sigma[sl][5][5] =5e-3;
				}

				memcpy(sigmaSave[sl], sigma[sl], 49*sizeof(double));
			}
			totMoment /= nbuf;


// 			for testing: eliminates the dispersion
// 			for(int i=0; i<7; i++)
// 			{
// 				sigma[i][5]=0.;
// 				sigma[5][i]=0.;
// 			}
// 			sigma[5][5] =1.;
#else
			unit(sigma);
			double betx = 10;
			double alfx= 1.;
			double epsx=1e-5;
			double bety = 20;
			double alfy= -1.;
			double epsy=1e-5;
			sigma[0][0] = betx*epsx;
			sigma[1][1] = (1.+ alfx*alfx)/betx*epsx;
			sigma[0][1] = -alfx*epsx;
			sigma[1][0] = -alfx*epsx;
			sigma[2][2] = bety*epsy;
			sigma[3][3] = (1.+ alfy*alfy)/bety*epsy;
			sigma[2][3] = -alfy*epsy;
			sigma[3][2] = -alfy*epsy;
			pmoment = 2.;
#endif
// 			printf("energy = %e, gamma=%e, moment=%e\n", energy,gammaL,totMoment);
			locstart=loc;
			save=1;
			double epsx = sqrt(sigma[0][0]*sigma[1][1] - sigma[0][1]*sigma[1][0]);
			double epsy = sqrt(sigma[2][2]*sigma[3][3] - sigma[2][3]*sigma[3][2]);
			gammaL = sqrt( square(totMoment/elmass) +1.);
			double betagamma = totMoment/elmass;
			double gammaL1 = gammaL < 1.0001 ? 1.0001 : gammaL;
			printf("gamma %e pmoment %e elmass %e betagamma %e\n", gammaL,totMoment,elmass, betagamma);
			if(gammaL <=1.) gammaL=1.1;

#ifdef MADX
			fprintf(fm, "beam, mass = 0.00051104, charge = -1, gamma = %e, exn = %e, eyn = %e, sige = 0.001;\n", gammaL1, epsx*betagamma, epsy*betagamma);
			fprintf(fm, "bx := %e;\n",  sigma[0][0]/epsx);
			fprintf(fm, "ax := %e;\n", -sigma[0][1]/epsx);
			fprintf(fm, "by := %e;\n",  sigma[2][2]/epsy);
			fprintf(fm, "ay := %e;\n", -sigma[2][3]/epsy);
#else
			fprintf(fm, "beam, mass = 0.00051104, charge = -1, gamma = %e, exn = %e, eyn = %e, sige = 0.001\n", gammaL1, epsx*betagamma, epsy*betagamma);
			fprintf(fm, "bx := %e\n",  sigma[0][0]/epsx);
			fprintf(fm, "ax := %e\n", -sigma[0][1]/epsx);
			fprintf(fm, "by := %e\n",  sigma[2][2]/epsy);
			fprintf(fm, "ay := %e\n", -sigma[2][3]/epsy);
#endif
			continue;
		}

		if( !strcmp( word, "!@slenvstop"))
		{
			printf("stop  at neo = %d\n", neo);
			save=0;
			continue;
		}

		if( !strcmp( word, "!@slenvfit"))
		{
// 			printf("fit  at neo = %d\n", neo);
			fitElem[nfit]=n;
			int nw = sscanf(line,"%s%d%le", word, fitParm+nfit, fitStep+nfit);
			if(fitParm[nfit]+1 > nparm) nparm= fitParm[nfit]+1;
			if(nw != 3) 
			{
				fprintf(stderr, "missimg stepsize:%s", line);
				exit(-1);
			}
			printf("sfit n=%d nfit=%d nparm=%d fitparm=%d fitstep=%15.8e\n",
					n, nfit, nparm, fitParm[nfit],fitStep[nfit]);
			nfit++;
			continue;
		}


		if( !strcmp( word, "!@slenvgoal"))
		{
			int nw = sscanf(line,"%s%le", word, fitStep+ngoal);
			if(nw != 2) 
			{
				fprintf(stderr, "missimg stepsize:%s", line);
				exit(-1);
			}
			ngoal++;
			continue;
		}

		if( !strcmp( word, "cathode"))
		{
			neo++;
// 			printf("%-20s\t%d\n", word, neo);
			continue;
		}

		if( !strcmp( word, "drift"))
		{
			sscanf(line,"%s%lf", word, &len);
			len /= 100.;
			loc += len;
			if(save)
			{
				type[n]=1;
				sloc[n] = loc;
				length[n]=len;
				strength[n++]=0;
				sprintf(elname[elnum++], "d%d",driftnum);
#ifdef MADX
#else
#endif
				fprintf(fm, "d%d: drift, l=%f;\n", driftnum++, len);
			}
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "quad"))
		{
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri,&stre);
			len /= 100.;
			loc += len;
			if(save)
			{
				type[n]=2;
				sloc[n] = loc;
				length[n]=len;
				stre *= 100./1000.;		// convert from G/cm to kG/m
				strength[n++]=stre;
				double kstrength = stre / totMoment* clight*1e-7;
// 				printf("quad k1=%e, l=%f\n", kstrength, len);
				sprintf(elname[elnum++], "q%d",quadnum);
#ifdef MADX
				fprintf(fm, "kq%d:=%e;\n", quadnum, kstrength);
				fprintf(fm, "q%d: quadrupole, l=%f, k1:=kq%d;\n", quadnum, len, quadnum);
#else
				fprintf(fm, "kq%d:=%e\n", quadnum, kstrength);
				fprintf(fm, "q%d: quadrupole, l=%f, k1=kq%d\n", quadnum, len, quadnum);
#endif
				quadnum++;
			}
			neo++;
			continue;
		}

		if( !strcmp( word, "solenoid"))
		{
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri,&stre);
			len /= 100.;
			loc += len;
			if(save)
			{
				type[n]=3;
				sloc[n] = loc;
				length[n]=len;
				stre /= 1000.;		// convert from G to kG
				strength[n++]=stre;
				double kstrength = stre / totMoment* clight*1e-7;
// 				printf("solenoid ks=%e, l=%f\n", kstrength, len);
				sprintf(elname[elnum++], "s%d",solnum);
				double sol = stre / totMoment* clight*1e-7;
#ifdef MADX
				fprintf(fm, "ks%d:=%e;\n", solnum, sol);
				fprintf(fm, "s%d: solenoid, l=%f, ks:=ks%d;\n", solnum, len, solnum);
#else
				fprintf(fm, "ks%d:=%e\n", solnum, sol);
				fprintf(fm, "s%d: solenoid, l=%f, ks=ks%d\n", solnum, len, solnum);
#endif
				solnum++;
			}
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "bend"))
		{
			sscanf(line,"%s%lf%lf%lf%lf%lf%lf%lf", word, &len, &app, &pri, &stre, &angl, &e1, &e2);
			len /= 100.;
			angl *= M_PI/180;
			e1 *= M_PI/180;
			e2 *= M_PI/180;
// 			printf("%s", line);
// 			printf("bend len %e angl %e e1 %e e2 %e  enegr  %e\n", len,  angl, e1, e2, stre);
			if(save && e1 != 0.)
			{
				type[n]=5;		// edge
				sloc[n] = loc;
				length[n]=e1;
				strength[n++]=angl/len;   // 1/rho
			}
			loc += len;
			if(save)
			{
				type[n]=4;
				sloc[n] = loc;
				length[n]=len;
				strength[n]=angl;
// 				printf("bend length %e strength %e\n", length[n], strength[n]);
				n++;
				sprintf(elname[elnum++], "b%d",bendnum);
#ifdef MADX
				fprintf(fm, "b%d: sbend, l=%f, angle=%e,  e1=%e, e2=%e;\n", bendnum++, len,angl, e1,e2);
#else
				fprintf(fm, "b%d: sbend, l=%f, angle=%e, &\n e1=%e, e2=%e\n", bendnum++, len,angl, e1,e2);
#endif
			}
			if( save && e2 != 0.)
			{
				type[n]=5;
				sloc[n] = loc;
				length[n]=e2;
				strength[n++]=angl/len;
			}
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

// 	case 5: edge(mat, length[i], strength[i]);	// length is used for e1 and e2, strength is 1/rhoc
// void edge(double mat[7][7], double angle, double hx)
// {
// 	unit(mat);
// 	double foc = tan(angle)*hx;
// 	mat[1][0] =   foc;
// 	mat[3][2] =  -foc;
// }
		if( !strcmp( word, "cell"))
		{
			sscanf(line,"%s%lf%lf%lf%lf%lf", word, &len, &app, &pri,&phas, &stre);
			len /= 100.;
			neo++;
			loc += len;
			if(save)		// treat like a drift for now
			{
				type[n]=1;
				sloc[n] = loc;
				length[n]=len;
				strength[n++]=0;
				sprintf(elname[elnum++], "c%d",cavnum);
#ifdef MADX
				fprintf(fm, "c%d: quadrupole, l=%f, k1=0;\n", cavnum++, len);
#else
				fprintf(fm, "c%d: quadrupole, l=%f, k1=0\n", cavnum++, len);
#endif
			}
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "rotate"))
		{
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri,&stre);
			if(save)
			{
				type[n]=6;
				length[n]=0;
				sloc[n] = loc;
				strength[n++]=stre*M_PI/180.;
			}
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

	}
	fclose(fi);
#ifdef MADX
	fprintf(fm, "lin: line = ( \n");
	for(int i=0; i<elnum; i++)
	{
		fprintf(fm, "%s", elname[i]);
		if(i < elnum-1) fprintf(fm, ",");
		else 		fprintf(fm, ");");
		if( i%10 == 9 && i != elnum-1) fprintf(fm, "\n");
	}
	fprintf(fm, "\nuse, sequence=lin;\n");
	if(doSurv)
		fprintf(fm, "survey, file=\"slenv.survx\";\n");

	fprintf(fm, "\n\n!match,  betx=bx, alfx=ax, bety=by, alfy=ay;\n");
	fprintf(fm, "!vary, name=kq0, step=0.01;\n");
	fprintf(fm, "!constraint, range=m01, alfx=0, alfy=0;\n");
	fprintf(fm, "!lmdif, calls=5000, tolerance=1e-10;\n");
	fprintf(fm, "!endmatch;\n\n\n");

	fprintf(fm, "twiss, file=\"slenv.twissx\", betx:=bx, alfx:=ax, bety:=by, alfy:=ay;\n");
	fprintf(fm, "stop;\n");
#else
	fprintf(fm, "lin: line = ( &\n");
	for(int i=0; i<elnum; i++)
	{
		fprintf(fm, "%s", elname[i]);
		if(i < elnum-1) fprintf(fm, ",");
		else 		fprintf(fm, ")");
		if( i%10 == 9 && i != elnum-1) fprintf(fm, "&\n");
	}
	fprintf(fm, "\nuse, lin\nprint, full\n");
	if(doSurv)
		fprintf(fm, "survey, tape=\"slenv.survx\"\n");

	fprintf(fm, "\n\n!match,  betx=bx, alfx=ax, bety=by, alfy=ay\n");
	fprintf(fm, "!vary, kq6, step=0.01\n");
	fprintf(fm, "!constraint, range=m01, alfx=0, alfy=0\n");
	fprintf(fm, "!lmdif, calls=5000, tolerance=1e-10\n");
	fprintf(fm, "!endmatch\n\n\n");

	fprintf(fm, "twiss, file=\"slenv.twissx\", betx=bx, alfx=ax, bety=by, alfy=ay\n");
	fprintf(fm, "stop\n");
#endif
	fclose(fm);

	for(int i=0; i<nfit; i++)
	{
		int j=fitElem[i];
		if(type[j] ==  2) fitStep[i] /= 10.;	// quad g/cm -> kG/m
		if(type[j] ==  3) fitStep[i] /= 1000.;	// sol g -> kG
		fitBest[i] = fitSave[i] = strength[j];
		printf("tfit i=%d j=%d type[j]=%d fitstep=%15.8e\n", i, j, type[j], fitStep[i]);
	}
	bestgoal=1e33;
    	t=1;	// one quality function
	xOptimal.setSize(nparm); 	// sets the number of variables
	xStart.setSize(nparm);


	valueOptimal=0.0;
	
	for(int i=0; i<nparm; i++)  xStart[i]=0.;
	printf("end constructor\n");
}


void unit(double mat[7][7])
{
	memcpy(mat, uu, 49*sizeof(double));
}

void printMat(const char* name, double mat[7][7])
{
	printf("%s = \n", name);
	for(int i=0; i<7; i++)
	{
		for(int j=0; j<7; j++)
			printf(" %e  ", mat[i][j]);
		printf("\n");
	}
	printf("\n");
}
void transport(double mat[7][7], double sigma[7][7])
{
	double tep[7][7];
// 	printMat("mat", mat);
// 	printMat("sigma", sigma);

	for(int i=0; i<7; i++)
	{
		for(int j=0; j<7; j++)
		{
			tep[i][j] = 0.;
			for(int k=0; k<7; k++)
			{
				tep[i][j] += sigma[i][k] * mat[j][k];
			}
		}
	}

// 	printMat("tep", tep);
	for(int i=0; i<7; i++)
	{
		for(int j=0; j<7; j++)
		{
			sigma[i][j] = 0.;
			for(int k=0; k<7; k++)
				sigma[i][j] += mat[i][k] * tep[k][j];
		}
	}

// 	printMat("sigma", sigma);
// 	printf("\n\n\n");


}


void lintr3(double mat[7][7], double el, double qk, double hx, double hz)
{

	memcpy(mat, uu, 49*sizeof(double));	// set mat to unit matrix

//     x komponente

// 	printf("qk %e hx %e hz %e el %e \n", qk, hx, hz, el);

      	double kx=qk-hx*hx;
      	double sqk=sqrt(fabs(kx));
      	if(sqk > 1.e-10)
      	{
      		double phix=el*sqk;
		double s,c;
      		if(kx < 0. )
    		{	
   			s=sin(phix);
      			c=cos(phix);
      		}
		else
		{
   			s=sinh(phix);
      			c=cosh(phix);
		}
   		mat[0][0]=  c;
      		mat[0][1]=  s/sqk;
      		mat[0][5]= -hx*(1.-c)/(sqk*sqk);
      		mat[1][0]= -s*sqk;
      		mat[1][1]=  c;
      		mat[1][5]= -hx*s/sqk;
      		if(kx > 0. )
		{
   			mat[0][5]= -mat[0][5];
      			mat[1][0]= -mat[1][0];
      		}
   		mat[4][1]= -mat[0][5];
   		mat[4][0]= -mat[1][5];
	}
	else
	{
//     limes (k-1/rho**2)  gegen 0
      		mat[0][0]=  1.0;
      		mat[0][1]=  el;
      		mat[0][5]= -hx*el*el/2.;
      		mat[1][0]=  0.0;
      		mat[1][1]=  1.0;
      		mat[1][5]= -hx*el;
	}




//     z komponente


      	double kz= -qk-hz*hz;
      	sqk=sqrt(fabs(kz));
      	if(sqk > 1.e-10)
	{
      		double phiz=el*sqk;
		double s,c;
      		if(kz < 0. )
		{
      			s=sin(phiz);
      			c=cos(phiz);
		}
		else
		{
      			s=sinh(phiz);
      			c=cosh(phiz);
		}
      		mat[2][2]=  c;
      		mat[2][3]=  s/sqk;
      		mat[2][5]= -hz*(1.-c)/(sqk*sqk);
      		mat[3][2]= -s*sqk;
      		mat[3][3]=  c;
      		mat[3][5]= -hz*s/sqk;
      		if(kz > 0. )
		{
      			mat[2][5]= -mat[2][5];
      			mat[3][2]= -mat[3][2];
   		} 
   		mat[4][3]= -mat[2][5];
   		mat[4][2]= -mat[3][5];
	}
	else
	{
//     limes (k-1/rho**2)  gegen 0

    		mat[2][2]=  1.0;
      		mat[2][3]=  el;
      		mat[2][5]= -hz*el*el/2.;
      		mat[3][2]=  0.0;
      		mat[3][3]=  1.0;
      		mat[3][5]= -hz*el;
	}

// 	printMat("mat", mat);
// 	printf("lintra = len=%f, qk=%e, hx=%e\n", el, qk, hx);
// 	printf("det = %e  %e\n", mat[0][0]*mat[1][1] -  mat[0][1]*mat[1][0], mat[2][2]*mat[3][3] -  mat[2][3]*mat[3][2]);
// 	printf("\n\n\n");


}






void quad(double mat[7][7], double tlen, double gradient, double pmoment)
{
	double kstrength = gradient / pmoment* clight*1e-7;
// 	printf("quad  k=%e l = %f\n", kstrength, tlen);
	lintr3(mat, tlen, kstrength, 0., 0.);
}

void bend(double mat[7][7], double tlen, double angle)
{
	double drho = angle/tlen;
// 	printf("drho = %e\n", drho);
	lintr3(mat, tlen, 0., drho,  0.);
// 	printMat("bend mat", mat);
}
void edge(double mat[7][7], double angle, double hx)
{
	unit(mat);
	double foc = tan(angle)*hx;
// 	printf("edge angle %e hx %e foc %e\n", angle, hx, foc);
	mat[1][0] =   foc;
	mat[3][2] =  -foc;
}

		
void spacelense(double mat[7][7], double ksc, double tlen)
{
	unit(mat);
	double foc = ksc*tlen;
	mat[1][0] =   foc;
	mat[3][2] =   foc;
}

void solenoid(double mat[7][7], double tlen, double field, double pmoment)
{
	double sol = field / pmoment* clight*1e-7;
	memcpy(mat, uu, 49*sizeof(double));
      	if(sol  ==  0.)
      	{
	      mat[0][1]=mat[2][3] =tlen;
	      return;
      	}
//
//     CALCULATION OF SOLENOID ACCORDING TO :
//     "  G. RIPKEN : DESY INTERNAL REPORT  R1-70/5  1970  "
//     CROSSING ANGLE BETWEEN MAGNET AXIS AND CLOSED ORBIT NEGLECTED.


      	double ca=cos(sol*tlen)/2.;
      	double sa=sin(sol*tlen)/2.;
      	double cb=.5-ca; ca=.5+ca;


      	mat[0][0]=ca;
      	mat[1][1]=ca;
      	mat[2][2]=ca;
      	mat[3][3]=ca;

      	mat[2][0]=sa;
      	mat[3][1]=sa;
      	mat[0][2]=-sa;
      	mat[1][3]=-sa;


      	mat[0][1]=sa/sol*2.;
      	mat[2][3]=mat[0][1];
      	mat[1][0]=-sa*sol/2.;
      	mat[3][2]=mat[1][0];


      	mat[1][2]=cb*sol/2.;
      	mat[3][0]=-mat[1][2];
      	mat[2][1]=cb/sol*2.;
      	mat[0][3]=-mat[2][1];

// 	printMat("mat", mat);
// 	printf("sol = len=%f, field=%e, moment=%e\n", tlen, field, pmoment);
// 	printf("\n\n\n");

}

void Opttrack::plot()
{
// 	exit(0);
	fx = fopen("slenvx.agr","w");
	if(!fx)
	{
		perror("slenvx.agr");
		exit(-1);
	}

	fy = fopen("slenvy.agr","w");
	if(!fy)
	{
		perror("slenvy.agr");
		exit(-1);
	}

	fxp = fopen("slenvxp.agr","w");
	if(!fxp)
	{
		perror("slenvxp.agr");
		exit(-1);
	}

	fyp = fopen("slenvyp.agr","w");
	if(!fyp)
	{
		perror("slenvyp.agr");
		exit(-1);
	}

	FILE * fdx = fopen("slenvdx.agr","w");
	if(!fdx)
	{
		perror("slenvdx.agr");
		exit(-1);
	}

	FILE * fdy = fopen("slenvdy.agr","w");
	if(!fdy)
	{
		perror("slenvdy.agr");
		exit(-1);
	}

	fex = fopen("slenvex.agr","w");
	if(!fex)
	{
		perror("slenvex.agr");
		exit(-1);
	}

	fey = fopen("slenvey.agr","w");
	if(!fey)
	{
		perror("slenvey.agr");
		exit(-1);
	}

	fs = fopen("slenvs.agr","w");
	if(!fs)
	{
		perror("slenvs.agr");
		exit(-1);
	}


	memcpy(sigma, sigmaSave, 49*sizeof(double));
	double oldloc = locstart;
	double epsx= sqrt(sigma[0][0]*sigma[1][1] - sigma[1][0]*sigma[0][1]);
	double epsy= sqrt(sigma[2][2]*sigma[3][3] - sigma[3][2]*sigma[2][3]);
	fprintf(fex, "%e %e\n", oldloc, epsx);
	fprintf(fey, "%e %e\n", oldloc, epsy);
	fprintf(fdx, "%e %e\n", oldloc, sigma[0][5]/sigma[5][5]/1000.);
	fprintf(fdy, "%e %e\n", oldloc, sigma[2][5]/sigma[5][5]/1000.);
// 	fprintf(fx, "%e %e\n", oldloc, sigma[0][0]/epsx);
// 	fprintf(fy, "%e %e\n", oldloc, sigma[2][2]/epsy);
	if( sqrt(sigma[0][0]) < 0.1) fprintf(fx, "%e %e\n", oldloc, sqrt(sigma[0][0]));
	if( sqrt(sigma[2][2]) < 0.1) fprintf(fy, "%e %e\n", oldloc, sqrt(sigma[2][2]));
	fprintf(fxp, "%e %e\n", oldloc, -sigma[0][1]);
	fprintf(fyp, "%e %e\n", oldloc, -sigma[2][3]);
	fprintf(fs,  "%e  0\n", oldloc);
	for(int i=0; i<n; i++)
	{

		propagate(i, pmoment[sl]);
		double epsx = sqrt(sigma[0][0]*sigma[1][1] - sigma[0][1]*sigma[1][0]);
		double epsy = sqrt(sigma[2][2]*sigma[3][3] - sigma[2][3]*sigma[3][2]);
		fprintf(fex, "%e %e\n", sloc[i], epsx);
		fprintf(fey, "%e %e\n", sloc[i], epsy);
		fprintf(fdx, "%e %e\n", sloc[i], sigma[0][5]/sigma[5][5]/1000.);
		fprintf(fdy, "%e %e\n", sloc[i], sigma[2][5]/sigma[5][5]/1000.);
// 		if( sigma[0][0]/epsx < 10000) fprintf(fx, "%e %e\n", sloc[i], sigma[0][0]/epsx);
// 		if( sigma[2][2]/epsy < 10000) fprintf(fy, "%e %e\n", sloc[i], sigma[2][2]/epsy);
		if( sqrt(sigma[0][0]) < 0.1) fprintf(fx, "%e %e\n", sloc[i], sqrt(sigma[0][0]));
		if( sqrt(sigma[2][2]) < 0.1) fprintf(fy, "%e %e\n", sloc[i], sqrt(sigma[2][2]));
		fprintf(fxp, "%e %e\n", sloc[i], -sigma[0][1]);
		fprintf(fyp, "%e %e\n", sloc[i], -sigma[2][3]);
		fprintf(fs, "%e 0\n", oldloc);
		fprintf(fs, "%e %e\n", oldloc, type[i]*1e-4);
		fprintf(fs, "%e %e\n", sloc[i], type[i]*1e-4);
		fprintf(fs, "%e 0\n", sloc[i]);
		oldloc = sloc[i];
	}
	printMat("last sigma =", sigma);

}

void Opttrack::saveBest()
{
	for(int i=0; i<nfit; i++)
	{
		int j=fitElem[i];
		fitBest[i]= strength[j] ;
	}
}


void Opttrack::installBest()
{
	for(int i=0; i<nfit; i++)
	{
		int j=fitElem[i];
		strength[j] = fitBest[i];
		switch(type[j])
		{
			case 2:
				printf("!best quad %e\n",  strength[j]*10.);
				break;
			case 3:
				printf("!best sol %e\n",  strength[j]*1000.);
				break;
			default:
				break;
		}
	}
}

double Opttrack::eval(Vector X, int *nerror, int resource)
{
	if(nerror) *nerror=0;
    	double *p=X;
    	double *start=xStart;



	for(int i=0; i<nfit; i++)
	{
		double oldstrength = strength[fitElem[i]];
		strength[fitElem[i]]= fitSave[i] + p[fitParm[i]]*fitStep[i];
// 		printf("fit %3d %3d %3d %15.8e %15.8e %15.8e %15.8e %15.8e\n", i, fitElem[i], fitParm[i], fitSave[i],
// 				oldstrength, strength[fitElem[i]], fitStep[i], p[fitParm[i]]*fitStep[i]);
	}
	memcpy(sigma, sigmaSave, 49*sizeof(double));

	for(int i=0; i<n; i++)
	{
		propagate(i);
	}
	double goal = fabs(sigma[0][1]) +fabs(sigma[2][3]) + square(sigma[0][5]);
// 	double goal = square(sigma[0][1]) +square(sigma[2][3]);
// 	double goal = square(sigma[0][5]);
	char * nobest= strdup("   ");
	char * isbest= strdup("###");
	char * showbest=nobest;
        pthread_mutex_lock( &rtxMutex );
	if(goal < bestgoal)
	{
		showbest=isbest;
		bestgoal=goal;
		saveBest();
	}
	FILE * flog = stdout;
	fprintf(flog, "%sbest %2d %18.10e this %18.10e<- ", showbest, resource, bestgoal, goal);
	for(int i=0; i<nparm; i++)  fprintf(flog, "%18.10e ", p[i]);
	fprintf(flog, ": ");
	ngoal=1;
	for(int i=0; i<ngoal; i++)  fprintf(flog, "%18.10e ", goal);
	fprintf(flog, "\n");
	fflush(flog);
        pthread_mutex_unlock( &rtxMutex );

	return goal;
}


void Opttrack:: propagate(int i, double tmoment)
{
	int nn = int(length[i]/magnetSlice+0.5);
	if(nn <= 0) nn=1;
	if(type[i] == 5) nn = 1;	// dont spit edges
	double tlen=length[i]/nn;

	switch(type[i])
	{
		case 0: break;
		case 1: unit(mat);
			mat[0][1] = mat[2][3] = tlen;
			break;
		case 2: quad(mat, tlen, -strength[i], tmoment);
			break;
		case 3: solenoid(mat, tlen, strength[i], tmoment);
			break;
		case 4: bend(mat, tlen, strength[i]/nn);	// strength is the angle
			break;
		case 5: edge(mat, length[i], strength[i]);	// length is used for e1 and e2, strength is 1/rhoc
			break;
	}

	for(int j=0; j<nn; j++)
	{
		transport(mat, sigma);
		if(doSpaceCharge)
		{
			double ksc = spcfac * IoverI0 *elmass*elmass*elmass /( tmoment*tmoment*tmoment * (sigma[0][0]+sigma[2][2])   );
			spacelense(matsp, ksc, tlen);
			transport(matsp, sigma);
		}


// 		double tilt = - atan2(sigma[0][2], sigma[2][2] - sigma[0][0])/2.;
// 		double ctilt = cos(tilt);
// 		double stilt = sin(tilt);
// 		double sigxx= ctilt*ctilt*sigma[0][0] + 2*ctilt*stilt*sigma[0][2] + stilt*stilt*sigma[2][2];
// 		double sigyy= ctilt*ctilt*sigma[2][2] - 2*ctilt*stilt*sigma[0][2] + stilt*stilt*sigma[0][0];

	}
	
}


main(int argc, char ** argv)
{
	int doSpaceCharge=1;
	int doSurv=0;
	int doFit=0;
	int zeroDisp =0;
	while(argc > 2)
	{
		if( !strcmp(argv[1], "-f"))	// do fitting
		{
			doFit=1;
			argc--; argv++;
		}
		if( !strcmp(argv[1], "-l"))	//  ingnore space charge
		{
			doSpaceCharge=0;
			argc--; argv++;
		}
		if( !strcmp(argv[1], "-s"))	// include survey command in mad file
		{
			doSurv=1;
			argc--; argv++;
		}
		if( !strcmp(argv[1], "-d"))	// set the beginning dispersion to zero and the energy spread to 5.e-3
		{
			zeroDisp=1;
			argc--; argv++;
		}
	}

	double spcfac =1;
	if( argc > 2) spcfac = atof(argv[2]);

	Opttrack * op = new Opttrack(argv[1], spcfac, doSpaceCharge, doSurv, zeroDisp);

	if(doFit)
	{
    		double rhoStart=1e-0, rhoEnd=1e-5;
    		int niter=100000;
    		ObjectiveFunction *of = op;
		int maxthreads = of->setMaxThreads(1);
		CONDOR(rhoStart, rhoEnd, niter, of);
	}
	op->installBest();
	op->plot();

}

