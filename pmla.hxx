//	pmla class : object to read the parmew binary files
//	translation from the example.for file provided by
//	the parmela distribution
//	jorg kewisch 11.2002
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <math.h>
#include <sys/stat.h>
// #include "powell_class.hxx"

double dminv(double *a, int n);


#ifndef O_BINARY 
const int O_BINARY = 0;
#endif
class Pmla
{
	inline double square(double x) { return x*x;}
	static const double elmass = 0.51104;			// in MeV
	int fp2;		// file handle for the tape2 file
	int rdwr2;		// flag for tape2 file is open for rw
	char * elemfilename;
	int fp3;		// file handle for the tape3 file
	char * timefilename;
	int filehandle[10];	// up to 10 tape2 files
	int filesize[10];
	int ireclen;		// reclen of tape2 file
	char * buff;		// buffer to read data
	float * fbuff;		// same as buff, but float

	Pmla() {}		// never call the default constructor
public:
	int numrec;
	int npoints;
	int irun;
	char title[81];
	float pi;		// a famous number
	float twopi;		// the larger brother of a famous number
	float radian;		// 2*pi/360
	float clight;		// speed of light
	float erest[3];		// rest energy for the three species of particles
	float brhof;		
	float freq;		// frequency of base rf system
	float wavel;		// wave length
	float z0;		// initial position of ref particle
	float w0;		// initial energy of ref particle;
	float beami;		// beam current
	float dcon;		// a famous rat poison?
	float zlimit;		// particles with z > zlimit are eliminated
	int nel;		// number of elements in the beam line
	int nplots;		// number of elements already calculated
	double dwt;
	double wt0;
	double dum;
	double swt;
	double cwt;
	int nsteps;
	int nsc;
	int nout;
	int ngood;		// number of good particles
	double coolenergy;
	

// for tape2
// the distribution of particles as given by the input particles. number 0 is the reference particle.
//   COR[I][0]=x
//   COR[I][1]=beta gamma x
//   COR[I][2]=y
//   COR[I][3]=beta gamma y
//   COR[I][4]=Phase (degree)
//   COR[I][5]=beta gamma z
// for tape2
// distribution after each element
//   CORD[I][0]=x
//   CORD[I][1]=bgx  (it seems this is x' * p_total/pz) (right now it is used as x' and y', needs lookng at)
//   CORD[I][2]=y
//   CORD[I][3]=bgy
//   CORD[I][4]=Phase (degree)
//   CORD[I][5]=Energy (MeV)
//   CORD[I][6]=original index of the particle in the COR array

	float **icord;		// coordinates of the initial distribution [6] long
	float **cord;		// coordinates of the  distribution at a step see above [7] long
	double **dcord;		// coordinates of the  distribution at a step, converted to double, units in meters, dcord[][1] is x_prime [7]
	double **scord;		// pointers to sort dcord
	double **ford;		// temporary array, for conversions of the dcord array [7]
	char  *ic;
	char  *iicol;
	char  *ip;
	int   *icol;
	int   *map;			// maps saved elements to all elements
	float *zloc;			// s-position of the element
	float *el;
	char ** type;
	int * rotat;
	double * Bfield;		// magnetic field on axis from zout 

	int maxj;
	float rfsize;
	int numfile;

	int neo;		// current element read from file
	float wr;		// reference energy
	float pr;		// reference phase


// tape3 data
	int nstep;
	int nbuf;
	float wt;
	float ax,bx,xbx,ybx,exrms,xmaxrms;
	float ay,by,xby,yby,eyrms,ymaxrms;
	float az,bz,xbz,ybz,ezrms,zmaxrms;

//  results
	int *set;		// element has been calculated
        double *ener;		// energy avarage
        double *ener1;		// energy avarage (test)
        double *gammaE;		// gamma avarage 
        double *phas;		// phase avarage
        double *blen;		// buch length rms
        double *disl;		// horizontal dispersion of element;
        double *dislp;		// horizontal slope of dispersion at element
        double *dispx;		// horizontal dispersion of element;
        double *dispxp;		// horizontal slope of dispersion at element
        double *dispy;		// vertical dispersion of element;
        double *dispyp;		// vertical slope of dispersion at element
        double *orbx;		// avarage orbit
        double *orbxp;		// avarage orbit
        double *orby;		// avarage orbit
        double *orbyp;		// avarage orbit
        double *epsx;		// horizontal emittance at element
        double *epsx_n;		// normalized horizontal emittance at element
        double *betx;		// horizontal beta function at element
        double *alfx;		// horizontal alpha function at element
        double *epsy;		// vertical emittance function at element
        double *epsy_n;		// normalized vertical emittance at element
        double *bety;		// vertical beta function at element
        double *alfy;		// vertical alpha function at element
        double *epss;		// longitudinal emittance function at element
        double *magnetization;	// <beta*gamma*(xy'-yx')>
        double *emit4d;		// 4 d emittance
        double *emit4dsum;	// sum of 4 d slice emittance
	double *emitz;             // longitudinal emittance in  deg*keV
	double * temperature;	//  sqrt(gamma^2*x'^2 + gamma^2*y'^2 + dp/p^2)
	double * temperxp;	//  sqrt(gamma^2*x'^2)
	double * temperyp;	//  sqrt(gamma^2*y'^2)
        double *ttrans;		// transverse temperature
        double *ttranm;		// transverse temperature (magnetized)
	int slices;
	double *r_eps;		// radial emittance
	double *r_bet;		// radial beta
	double *r_alf;		// radial alfa
	double *r_rot;		// radial alfa
	double (*s_eps)[500];	// slice emittace
	double (*s_bet)[500];
	double (*s_alf)[500];
	double (*s_disp)[500];
	double (*s_dispp)[500];
	double (*s_rot)[500];	// slice vortex
        double *epsx_c;		// horizontal emittance at elementi after magnetization correction
        double *epsy_c;		// horizontal emittance at elementi after magnetization correction
        double *epsx_u;		// horizontal emittance at elementi after unrotation
        double *epsy_u;		// horizontal emittance at elementi after unrotation
        double *eps_dx;		// dispersion contribution to the emittance
        double *eps_lx;		// longitudinal position contribution to the emittance
        double *rl2;		// avarage square of larmour radii
        double *rot;		// rotation speed of the beam = (xy'-yx')/r2
        double *envx;		// x enveloppe
        double *envy;		// y enveloppe
        double *roundness;	// enveloppe roundness
        double *divx;		// x' enveloppe
        double *divy;		// y' enveloppe
        double *roundnessp;	// slope roundness
        double *slopex;		// x' enveloppe
        double *slopey;		// y' enveloppe
	double *xmax;
	double *ymax;
	double *zmax;
	double *tmax;
	double *xmin;
	double *ymin;
	double *zmin;
	double *tmin;
        double *sigma_p;	// energy spread
        double *coolsigma;	// sigma p relative to coolenergy
        double *dp_max;		// energy spread
        double *dp_min;		// energy spread
        double *sigma_l;	// bunch length
        double *sigma_t;	// bunch length in seconds
        double *dl_max;		// bunch length
        double *dl_min;		// bunch length




// matrix

	double sigmaMat[7][7];		// the raw sigma matrix: the columns are x,x',y,y',phase, momentum and 1 (one), colume 6 containes the averages
	double sigmaMat2[7][7];		// the sigma matrix with the avarage removed, colume 6 containes then zero
	double maxPart[7], minPart[7];
	double m11;
	double m12;
	double m13;
	double m14;
	double m21;
	double m22;
	double m23;
	double m24;
	double m31;
	double m32;
	double m33;
	double m34;
	double m41;
	double m42;
	double m43;
	double m44;
	double mmt[4][4];
	double mmf[4][4];
	double mdelta[4];

// sigmas
        double xm;
        double xpm;
        double xx;
        double xxp;
        double xpxp;
        double ym;
        double ypm;
        double yy;
        double yyp;
        double ypyp;
	double xy ;
	double xyp;
	double xpy;
	double xpyp;

        double dlm;
        double dldl;
        double dlxp;
        double dlxpp;

        double dpm ;
        double dpdp;
        double ddxp;
        double ddxpp;
        double ddyp;
        double ddypp;


	double dd;
	double ld;
	double ll;
	double xd;
	double xl;
	double xpd;
	double xpl;
	double yd;
	double yl;
	double ypd;
	double ypl;
	double flattness;


	double bestgoal;
	int gunEnd;		// zout stores end of gun index

	Pmla(const char * elemfilename, const char * timefilename);
	virtual ~Pmla() ;
	int zout(const char * zoutfilename, const char * linefilename, double size);
	int scatter(const char * path);
	int getElement(int ne); 	// fills the cord  and ic array for element ne
	int putElement(int ne, int what=0); 	// saves the cord  and ic array for element ne
						// save dcord when what == 0 else saves ford
	double getValue(const char * what);
	int    getCurve(const char * what);

	int getTime(int ne);		// fills the cord  and ic array for time step ne **does not work yet**
	void twiss(int print);		// calculates twiss parameters for last gotten element
	void twiss(int print, double ** pt, int np);		// calculates twiss parameters for last gotten element
	void twiss1(int print);		// calculates twiss parameters for last gotten element
	void twiss1(int print, double ** pt, int np);		// calculates twiss parameters for last gotten element (backup, not used)
	void calcSigmaMatrix(double ** pt, int np, int print);
	inline int mapint(float &x) { return *((int *) &x); } // interprets float as int
	void  freadbuff(int fp, float * buff, int nfloats);
// 	void  freadbuff(int fp, double * buff, int nfloats);
	int  freadbuffx(int fp, float * buff, int nfloats);
	void  fwritebuff(int fp, float * buff, int nfloats);
	void  readbuff(int fp, char * buff, int nbytes);
	void  writebuff(int fp, char * buff, int nbytes);
	double mag(double p[], int n, int print);
	double mag1(double p[], int n);
	void  testcord(int ne, double m11, double m12, double m21);
	void  particle( int & i, double R, double x, double xp, double y, double yp);
	int  sliceEmit(int print);
	void cutEmittance(double cut, int step);
	void cutHeadTail(int print, double headcut, double tailcut);
	void cutDP(int print, double headcut, double tailcut);
	void cutTemp(int print, double headcut, double tailcut);
	void cutTempT(int print, double headcut, double tailcut);
	void sortZ(char what);		// what is 'z' for longitudinal sorting, 'p' for momentum sorting
	void sumEmittance(double * a[], int good);
	void longFit(double ** pt, int np);


	void dump();
	void dumpLongit(const char * cwd, double cut);
	void impactTdump();
	void histogram();
	int  slicePhasePlot(char * file, int rotat);

	int  printCord(char * file);
	double  unrotate(int method, int remove_orbit, double rm);		// unrotate beam. method 0 is <xy'-yx'> == 0; method 1 is <rl^2> == minimum.

	void  saveDelta();
	double  linMatrix();
	void  readTimeFile();

};


