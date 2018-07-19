


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


// ===================================================================================================================
//
//       PARMELA FITTING PROGRAM
//
//       Jorg Kewisch, 2005 - 2014
//
//      syntax:	mainp.exe [-r] [-d] [-p] <parmela input file> 
//      flags:	-r	apply the recovery file. This allows to continue where the previous run ended
//      	-d	start the Xvfb program with display :60 or :61
//      	-p	run parmela only once
//
// 	The program uses all cores of the computer. It finds the number of cores from the /prog/cpuinfo file with the line:
// 	int maxthreads = of->setMaxThreads(0);
// 	On windows systems you have to chage thezero to the number of cores you have.
// 	The program creates a directory "calcdir" in which the calculations are done. Whenever a "better" solution is found
// 	the tape, savecor and output files are copied back into the current directory.
// 
//      the input file is a parmela input file that contains the fitting instructions in parmela comment lines
//	the fitting package condor by Ir. Frank Vanden Berghen, Univesite' Libre de Bruxelles is used for fitting
//
//	!@var <parameter number> <active> <initial step size> <lower bound> <upper bound>  	defines a fit parameter. parameter number starts with zero. 
//		<parameter number> must be consequetive integers starting from 0 and is used in the subs command 
//              <active> equal 1 means this variable is used, 0 means it is ignored
//		<initial step size> guess
//		<lower bound> <upper bound>  are boundaries for the step, not the final value 
//
//	!@goal <where> <what> <active> <goal> <weight>						defines the figure of merrit. any number of goals can be used
//		<where> is the parmela element number.
//		<what> describes what to optimize. see the getval function in pmla.cxx
//              <active> equal 1 means this goal is used, 0 means it is ignored
//		<goal> is the wanted value for <what>, normally zero.
//		<weight> is the weight of a goal when more than one goal is fitted. One is a good weight.
//	any number in the parmela input file can be optimized by inserting a subs statement above the line:
//
//	!@subs <number of lines> <word> <modifiers> [ <word> <modifiers> ... ]
//		<number of lines>	the modification applies to the next n lines
//		<word>			the modification applies to the nth word in those lines. A word is anything seperated by whitespace
//					start counting words with zero.
//
//  	<modifier> describes how the following  <number of lines> lines are modified. The word is read from those lines.
//  	if for example this word contains 1.5, then
//		 0       means that the parameter 0 is added to the initial value 1.5 
//		-0       means that the parameter 0 is subtracted to the initial value 1.5 
//		p0       means that the parameter 0  is a length and is converted into a phase shift using the values supplied by the @beta statement,
//				which is added to 1.5 
//		q0       means the same as p0 for a third harmonic cavity
//		b0       means that the parameter 0  is the change in radius of the cathode and 1.5 is the field on the cathode. The field will be adjusted 
//				so that the magnetization is constant
//		c0,s1    have to be used in pair to optimze A*cos(phi) and A*sin(phi) insteadt of A and Phi. The wordindex 4 and 5 is hardwired, this is quick and dirty and should be changed later
//
//	!@beta <Energy> <RF frequency>							defines the energy  and the RF frequency (same as on the first line)
//											to calculate  beta for p and q
//	!@cath <radius> <field>					`			defines the initial cathode radius and B-field for b
//
//	!@step		no optimizing, but stepping through all active @var values from minimim to maximum in step size
//
//	!@mad ***broken ***indicated the insertion of mad code. A second !@mad indicates the end of the mad code. All mad lines are parmela commens (!).
//	after that comes the equivalent parmela code which is ignored by this program. a third !@mad indicates the end of this section
//	(this is not implemented correctly at this time, left over from a previos version)
// 
// 
// 	!@cutl <fraction>
// 		For the evaluation uses only a fraction of the beam. "!@cutl 0.9" means that 5% particles are removed from the head and 5% from the tail
// 	!@cutp <fraction>
// 		As !@cutl, but the cut is made in momentum space.
// 
// 	!@switch  <name> <0|1|-1>
// 	!@if <name>
// 	!@endif <name>
// 		Conditional input. If the switch is set to 0 then the input between !@if and !@endif is not copied  into the parmela input file.
//		If the switch is set to -1 the input is included in the first evauation, and ignored in all following runs
// 		Can be nested. 
//
//	Other !@ codes are ignored by this program
// ===================================================================================================================
//	the program can be aborted in a orderly manner by creating a file named "stop" in the current directory , i.e.  "touch stop"
//	the stop file will be removed o exit
//	the program can be paused by creating a file named "pause" in the current directory. when the user deletes "pause" the optimization continues.
//	A file "plot" will be removed after running ~/RTX/main5. "plot will be removed.
// ===================================================================================================================
//
//	this program calls the condor program, which is licensed under GPL. I guess that makes this programm GPL too.
//
//	condor builds a map using (n+1)*(n+2)/2 calculations and then uses this map to optimize




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <dirent.h>
#include <time.h>



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

const double clight= 2.99792456210e8;	// in m/s
const double elmass = 0.51104;
inline double square(double x) { return x*x;}

const int maxFitVariables = 100;		// max number of fit parameters
const int maxFitGoals = 100;		// max number of fit goals
const int maxSwitches = 100;		// max number of fit goals
const int maxRunStore = 1000;	    	// max number of states in the recobvery buffer;

int whatsystem;
char * ParmelaDir;
char * Wine;
char * preParmela;
char * slenvProg;
char * main5;
pid_t pidXvfb;
const char * thislock ="aaaRTXmainpRuns";

int myexit(int err)
{
	if(pidXvfb > 0) kill(pidXvfb, SIGTERM);
	unlink(thislock);
	fflush(stdout);
	fflush(stderr);
	exit(err);
}


int splitline(char * line, char ** space,  char ** word)
{
	int n=0;
	char * p = line;
	char newSpace[1000], newWord[1000];
	while(1)
	{
		char * q = newSpace;
		while( *p > 0 && *p <= ' ') {*q++=*p++; }
		*q=0;
		if(*p == 0) return n;
		q = newWord;
		while(*p > ' ') {*q++=*p++;}
		*q=0;
		space[n]=strdup(newSpace);
		word [n]=strdup(newWord );
		n++;
	}	
}
int str2words(char * line, char * beg[], int length[], char * after[])
{
	int n=0;
	char * p = line;
	while(1)
	{
		while( *p > 0 && *p <= ' ') p++;
		if(*p == 0) return n;
		beg[n] = p;
		int l=0;
		while(*p > ' ') { l++; p++;}
		length[n]=l;
		after[n]=p;
		n++;
	}	
}


class M2condor : public UnconstrainedObjectiveFunction
{
		int nrun;
		double para[maxRunStore][maxFitVariables], resu[maxRunStore], subresu[maxRunStore][maxFitVariables];
		int var_defined[maxFitVariables];
		int var_used[maxFitVariables];
		double phFac;
		double cath_size, cath_field;
		double cut;
		int cutWhat;
		int has_beta, has_cath, has_slenv;
		char slenvLatt[300];
		char slenvCavt[300];
		char switchName [100][20];
		int  switchFlag[100];
		int nswitches;
		int onetime;
		void madfit(FILE * fo, FILE * fi, double pnow[]);
	public:
		char inputFileName[300];	// input file
		char cwd[300];			// name of current directory, is used in temp file names
		int fit_ne[maxFitGoals];		// number of element in parmela
		char fit_what[maxFitGoals][10];		// name of what to optimize, i.e. emit4d..
		double fit_val[maxFitGoals];		// desired value
		double fit_weight[maxFitGoals];		// weight relative to other goals
		double fit_this[maxFitGoals];		// result of the current evaluation

		char var_name[maxFitVariables][30];	// yupp, this is the name of the variable
		double var_mult[maxFitVariables];	// step size for the optimize
		double var_lower[maxFitVariables];	// boundary
		double var_upper[maxFitVariables];	// boundary
		double var_step[maxFitVariables];	// step size for the step function
		int var_active[maxFitVariables];	// is active?;
		int var2parm[maxFitVariables];		// translates active variables to parameters
		int parm2var[maxFitVariables];		// translates active parameters to variables 
		int nparm;	// number of active variables
		int nvar;	// total number of variabes
		int ngoals;	// number of goals
		int ngood;	// the number of good particles, which must not get lost 
		int ninput;	// number of particles specified on "input" lines
		int  stepping;
		double bestgoal;

		FILE * flog;
		M2condor(char * inpfilename, int recover, int onetime);
    		double eval(Vector v, int *nerror=NULL, int resource=0);
    		double eval();
		void step(int level);
};




M2condor::M2condor(char * inpfilename, int recover, int onetime) 
{
	this->onetime = onetime;
	time_t tim=time(0);
	char logfile[50];
	sprintf(logfile,"pmla_log.%ld.txt", tim);
	unlink("pmla_log.txt");
	symlink(logfile, "pmla_log.txt");
	flog = fopen(logfile, "w");
	if(!flog)
	{
		fprintf(stderr, "cant open pmla_log.txt\n");
		myexit(-1);
	}
	fprintf(flog, "Input file  %s,  pid %d \n", inpfilename, getpid());
       	FILE * fi = fopen(inpfilename, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", inpfilename);
		myexit(-1);
	}
	strcpy(inputFileName, inpfilename);

	int ls=strlen(inputFileName);
	if(strcmp(inputFileName+ls-4,".inp") )
	{
		fprintf(stderr, "input file name must end in .inp, is %s\n", inputFileName+ls-4);
		myexit(-1);
	}
	getcwd(cwd, 300);
	fprintf(flog, "dir is %s\n", cwd);
	char ctim[100];
	ctime_r(&tim, ctim);
	fprintf(flog, "time is %s\n", ctim);


	ngood = -1; 	// not set
	ninput =1;	// there is always the reference particle
	ngoals=0;
	nparm=0;
	nvar=0;
	stepping=0;
	nswitches = 0;
	phFac=0;
	cutWhat = 0;
	has_beta= has_cath=has_slenv=0;
	for(int i=0; i<maxFitVariables; i++) {var_defined[i]= 0; var_used[i]= 0; fit_val[i] =0.; fit_weight[i] =0.; fit_this[i] =0.;}

	char line[500];
	char * pl;
	char what[100];
	int sofar;
	int neo = 0;
	while( pl=fgets(line,499,fi))
	{
		for(char * l=line; *l; l++)
		{
			if( *l >= 'A' && *l <= 'Z') *l += 32;
			if( *l == ',') *l=' ';
		}
		if(line[0] == '!' && line[1] == '@')
		{
// 			printf("=%s\n", line);
			pl +=2;
			sscanf(pl,"%s%n", what, &sofar);
			pl += sofar;
//			!@subs <number of lines> <word> <modifiers> [ <word> <modifiers>  ... ]
			if( ! strcmp(what, "subs") ) 
			{
				int numlines, word;
				sscanf(pl, "%d%n", &numlines, &sofar);
				pl += sofar;
				while(sscanf(pl, "%d%s%n", &word, what, &sofar) == 2)
				{
// 					printf("word, what %d %s  %d\n",  word, what, sofar); 
					pl += sofar;
					char * w= what;
					if(*w != '+'  &&  *w != '-'  &&  *w != 'p'  &&  *w != 'q'  &&  *w != 'b'  &&  *w != 'c'  &&  *w == 's')
					{
						fprintf(stderr, "wrong operator in %s\n", w);
					}
					w++;
					int varno = -1;
					for(int i=0; i<nvar; i++)
					{
						if(  ! strcmp(w, var_name[i]) )
						{
							varno = i;
							break;
						}
					}
					if(varno == -1)
					{
						fprintf(stderr, "variable %s is not defined\n", w);
						myexit(-1);
					}

// 					printf("found parameter %d\n", parmno);
					continue;
				}
				continue;
			}

			if( ! strcmp(what, "var") ) 
			{
				int activ;
				double p_mult, p_lower, p_upper;
				char vname[100];
				int st=sscanf(pl, "%s%d%lf%lf%lf", vname, &activ, &p_mult, &p_lower, &p_upper);
				if(st != 5)  {fprintf(stderr, "error in var statement !!!\n"); myexit(-1);}

				if(nvar == maxFitVariables)  {fprintf(stderr, "too many variables  !!!\n"); myexit(-1);}
				for(int i=0; i<nvar; i++)
				{
					if(  ! strcmp(vname, var_name[i]) )
					{
						fprintf(stderr, "variable %s double defined\n", vname);
						myexit(-1);
					}
				}

				strcpy(var_name[nvar], vname);
				var_defined[nvar]=1;
				var_mult[nvar] = p_mult;
				var_lower[nvar] = p_lower;
				var_upper[nvar] = p_upper;
				var_active[nvar] = activ;
				if(activ) printf("parameter   %d %f %f %f\n", nvar, p_mult, p_lower, p_upper);
				continue;
			}

			if( ! strcmp(what, "step") ) 
			{
				stepping=1;
				continue;
			}

			if( ! strcmp(what, "goal") ) 
			{
				int active;
				sscanf(pl, "%d%d%s%lf%lf", fit_ne+ngoals, &active, fit_what[ngoals], fit_val+ngoals, fit_weight+ngoals); 
				if(active)
				{
					printf("condition %d %d %s %f %f\n",  fit_ne[ngoals], active, fit_what[ngoals], fit_val[ngoals], fit_weight[ngoals]); 
					ngoals++;
				}
				continue;
			}

			if( ! strcmp(what, "cutl") ) 
			{
				cutWhat = 2;
				sscanf(pl,"%lf", &cut);
				fprintf(flog, "cut longitudinal at %f %%\n", cut);
				continue;
			}

			if( ! strcmp(what, "cutp") ) 
			{
				cutWhat = 1;
				sscanf(pl,"%lf", &cut);
				fprintf(flog, "cut momentum at %f %%\n", cut);
				continue;
			}

			if( ! strcmp(what, "switch") ) 
			{
				char swn[20];
				int  nsw;
				sscanf(pl,"%s%d", swn, &nsw);
				printf(" found switch %s = %d\n", swn, nsw);
				for(int i=0; i< nswitches; i++)
				{
					if(  !strcmp(swn, switchName[i]) )
					{
						switchFlag[i]=nsw;
						goto foundswitch;
					}
				}
				if(nswitches == maxSwitches)
				{
					fprintf(stderr, "to many switches\n");
					myexit(-1);
				}
				switchFlag[nswitches]=nsw;
				strcpy(switchName[nswitches++], swn);

			foundswitch: 	continue;
			}

			if( ! strcmp(what, "beta") ) 
			{
				double ener_b, freq;
				sscanf(pl,"%lf%lf", &ener_b, &freq);
				double gamma=ener_b / elmass  +1.;
				double beta = sqrt(1.-1./(gamma*gamma));
				phFac= -360.*freq*1e6/(100.*clight*beta);	// 100 because the parmela input is in cm
				printf("beta  e %e gamma %e beta %e freq %e phFac %e\n", ener_b, gamma, beta, freq, phFac);
				has_beta=1;
				continue;
			}

			if( ! strcmp(what, "slenv") ) 
			{
				sscanf(pl,"%s%s", slenvLatt, slenvCavt);
				printf("slenvLatt is %s\n", slenvLatt);
				printf("slenvCavt is %s\n", slenvCavt);
				fprintf(flog, "slenvLatt is %s\n", slenvLatt);
				fprintf(flog, "slenvCavt is %s\n", slenvCavt);
				FILE * fp = fopen(slenvLatt, "r");
				if(!fp)
				{
					fprintf(stderr,"cant open %s for read\n", slenvLatt);
					myexit(-1);
				}
				fclose(fp);
				char file[300];
				sprintf(file, "%s01.T7", slenvCavt);
				fp = fopen(file, "r");
				if(!fp)
				{
					sprintf(file, "%s1.T7", slenvCavt);
					fp = fopen(file, "r");
					if(!fp)
					{
						fprintf(stderr,"cant open %s for read\n", file);
						myexit(-1);
					}
				}
				fclose(fp);
				has_slenv=1;
				continue;
			}

			if( ! strcmp(what, "cath") ) 
			{
				sscanf(pl,"%lf%lf", &cath_size, &cath_field);
				printf("cath size %f field %f\n", cath_size, cath_field);
				has_cath=1;
				continue;
			}
			fprintf(stderr, "ignored %s\n", line);
		}

		what[0]=0;
		sscanf(line, "%s", what);
		if( ! strcmp(what, "input")  ) 
		{
			int type, num;
			for(char * p = line; *p; p++) if(*p == ',') *p=' ';
			sscanf(line, "%s%d%d", what, &type, &num);
			ninput += num;
			continue;
		}
		if( !strcmp( what, "cathode"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "drift"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "quad"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "solenoid"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "bend"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "cell"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}

		if( !strcmp( what, "rotate"))
		{
			neo++;
// 			printf("%-20s\t%d\n", what, neo);
			continue;
		}


	}
	fclose(fi);
	printf("found %d particles on input lines\n", ninput);



	if(nvar <= 0 && (!onetime))
	{
		fprintf(stderr, "nvar = %d, no parameters specified or used\n", nvar);
		myexit(-1);
	}


	nparm =0;
	for(int i=0; i<nvar; i++)
	{
		parm2var[i] = -1;
		var2parm[i] = -1; 
	}
	for(int i=0; i<nvar; i++)
	{
		if(!var_used[i])    {fprintf(stderr, " parameter %d not used !!!\n", i); myexit(-1);}
		if(!var_defined[i]) {fprintf(stderr, " parameter %d not defined !!!\n", i); myexit(-1);}
		if(var_active[i]) 
		{
			parm2var[nparm]=i;
			var2parm[i] =nparm++; 
		}
	}
// 	for(int i=0; i<nvar; i++)
// 	{
// 		printf("var %d parm2var %d ,var2parm %d\n", i, parm2var[i],var2parm[i]);
// 	}
	for(int i=0; i< ngoals; i++) printf("cond %d %s %f %f\n", fit_ne[i], fit_what[i], fit_val[i], fit_weight[i]);










// ----------------------------------------
	nrun=0; 
	bestgoal=1e33;
	if(recover)
	{
		FILE * fp = fopen("recover.hex", "r");
		if(!fp)
		{
			fprintf(stderr,"cant open %s for read\n", "recover.hex");
			myexit(-1);
		}

		int nnparm, nngoals;
		char rfile[200];
		fscanf(fp, "%d%d%s",  &nnparm, &nngoals, rfile);
		if( strcmp(rfile, inputFileName)  )
		{
			fprintf(stderr,"Input file name %s is different from name in recover.hex %s\n", inputFileName, rfile);
			myexit(-1);
		}
		if(nnparm != nparm || nngoals != ngoals)
		{
			fprintf(stderr," parateter mismatch nparm=%d nnparm=%d ngoals=%d nngoals=%d\n", nparm, nnparm, ngoals, nngoals);
			myexit(-1);
		}
		fprintf(flog, "recovered nrun=%d, nparm=%d, ngoals=%d\n", nrun, nnparm, nngoals);

		char seperator[10];
		double val;
		nrun=0;
		while(1)
		{
			int it[2];
			double * tt=(double*)it;
			for(int j=0; j<nparm; j++)
			{
				if(fscanf(fp, "%x%x", it, it+1) != 2) goto ende;
				para[nrun][j]= *tt;
			}
			for(int j=0; j<ngoals; j++)
			{
				if(fscanf(fp, "%x%x", it, it+1) != 2) goto ende;
				subresu[nrun][j]= *tt;
			}
			if(fscanf(fp, "%x%x%s%le", it, it+1, seperator, &val) != 4) goto ende;
			if(seperator[0] != '=') goto ende;
			resu[nrun] = *tt;
			char * nobest= strdup("   ");
			char * isbest= strdup("###");
			char * showbest=nobest;
			if( resu[nrun] < bestgoal) 
			{
				showbest=isbest;
				bestgoal= resu[nrun];
			}
			fprintf(flog, "%sbest *old* %2d %18.10e this %18.10e<- ", showbest, -1, bestgoal, resu[nrun]);
			for(int j=0; j<nparm; j++)  fprintf(flog, "%18.10e ", para[nrun][j]);
			fprintf(flog, ": ");
			for(int j=0; j<ngoals; j++)  fprintf(flog, "%18.10e ", subresu[nrun][j]);
			fprintf(flog, "\n");
			nrun++;
		}
	ende:   fclose(fp);
		fprintf(flog, "recovered %d calculations\n", nrun);
		fflush(flog);
	}
// -----------------------------------------
    	t=1;	// one quality function
//     	strcpy(name,"mainp");


	if(nparm == 0)
	{
		fprintf(stderr, "No variables set, exiting\n");
		myexit(-1);
	}
	if(ngoals == 0)
	{
		fprintf(stderr, "No goals set, exiting\n");
		myexit(-1);
	}
	xOptimal.setSize(nparm); 	// sets the number of variables
	xStart.setSize(nparm);


	valueOptimal=0.0;
	
	for(int i=0; i<nparm; i++)  xStart[i]=0.;
	printf("end constructor\n");
}



void M2condor::step(int level)
{
	int steps=int(var_step[level]+0.5);
	for(int i=0; i<steps; i++)
	{
		double x= var_lower[level]+ i*(var_upper[level]- var_lower[level])/(steps-1);
		xStart[level]=x;
		if(level < nparm-1)
		{
			step(level+1);
		}
		else
		{
			int error;
			double r= eval(xStart, &error, 0);    // use parallel in future
			if(error) fprintf(stderr, "eval returned error %d\n", error);
			fprintf(flog, "step ");
			for(int j=0; j< level; j++) fprintf(flog, "%e ", double(xStart[level]));
			fprintf(flog, "=> ");
			for(int j=0; j< ngoals; j++) fprintf(flog, "%e ", fit_this[j]);
			fprintf(flog, "\n");
		}
			
	}
}





double M2condor::eval()
{
	for(int i=0; i<nparm; i++)  xStart[i]=0;
	int error;
	double r= eval(xStart, &error, 0);
	if(error) fprintf(stderr, "eval returned error %d\n", error);
	return r;
}



double M2condor::eval(Vector X, int *nerror, int resource)
{
	if(nerror) *nerror=0;
    	double *p=X;
    	double *start=xStart;







        pthread_mutex_lock( &rtxMutex );
	fprintf(flog, "run %d,%2d start =", nrun,resource);
	for(int i=0; i<nparm; i++)  fprintf(flog, " %20.10e ", start[i]);
	fprintf(flog, "\n");
	fprintf(flog, "run %d,%2d parms =", nrun,resource);
	for(int i=0; i<nparm; i++)  fprintf(flog, " %20.10e ", p[i]);
	fprintf(flog, "\n");
	double pnow[nparm];
	for(int i=0; i<nparm; i++) pnow[i]=p[i]*var_mult[parm2var[i]];
	fprintf(flog, "run %d,%2d trunc =", nrun,resource);
	for(int i=0; i<nparm; i++)  fprintf(flog, " %20.10e ", pnow[i]);
	fprintf(flog, "\n");

//	if outside limits return error
	int boundErr=0;
	for(int i=0; i<nparm; i++)
	{
		int j=parm2var[i];
		if(pnow[i] < var_lower[j]   || pnow[i] > var_upper[j])
		{
			fprintf(flog, "parameter %d exceedes boundary\n", i);
			 boundErr=1;
		}
	}
	fflush(flog);
        pthread_mutex_unlock( &rtxMutex );
	if(boundErr)
	{
		if(nerror) *nerror=1;
		return 2.e33;
	}



//	check if we have done that before
	for(int i=0; i<nrun; i++)
	{
		for(int j=0; j<nparm; j++)
			if(fabs(p[j] -para[i][j]) > 1.e-8) goto nextj;
		fprintf(flog, "done that in step %d, %20.10e\n", i, resu[i]);
		return resu[i];

nextj:	;
	}


	char resourcedir[20];
	sprintf(resourcedir, "calcdirs/%03d/", resource);

	// make temp file names. We have the original input "inputFileName" in the main directory. inputFileName is a class member
	// we make a copy with substitutions "tempInputFileName" in the resource directory, This file is named "calcdirs/resource/pmla_temp.inp".
	char tempInputFileName[300];
	sprintf(tempInputFileName, "pmla_temp.inp");
	// we need this file name with and without the calcdirs/resource part
	char dirTempInputFileName[300];
	sprintf(dirTempInputFileName, "calcdirs/%03d/%s", resource, tempInputFileName);
	// we make a copy without comment lines  "strippedInputFileName" in the resource directory,
	// because parmela croakes on to many comment lines
	// This file is named "pmla_dir1_dir2_resource.inp".
	// dir1 and dir2 are the last two levels of directories
	// this file name is argument to parmela and
	// it is visible in the ps -ef command, so it helps that it includes the dir name
	char strippedInputFileName[300];
	int le=strlen(cwd);
	char *p1=cwd, *p2=cwd;
	for(char *p3=cwd; *p3; p3++) if(*p3 == '/') { p1=p2; p2=p3;  *p3='_';}
	sprintf(strippedInputFileName,"pmla_%s_%03d.inp", p1+1, resource);
	// we need this file name with and without the calcdirs/resource part
	char dirStrippedInputFileName[300];
	sprintf(dirStrippedInputFileName,"calcdirs/%03d/%s", resource, strippedInputFileName);

// 	printf("%2d  inputFileName		%s\n", resource, inputFileName);
// 	printf("%2d  tempInputFileName		%s\n", resource, tempInputFileName);
// 	printf("%2d  dirTempInputFileName	%s\n", resource, dirTempInputFileName);
// 	printf("%2d  strippedInputFileName	%s\n", resource, strippedInputFileName);
// 	printf("%2d  dirStrippedInputFileName	%s\n", resource, dirStrippedInputFileName);
	








        FILE * fi = fopen(inputFileName, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", inputFileName);
		myexit(-1);
	}
	

        FILE * fo = fopen(dirTempInputFileName, "w");
	if(!fo)
	{
		fprintf(stderr, "cant open %s\n", dirTempInputFileName);
		myexit(-1);
	}

	char line[500];
	char * pl;
	char what[100];
	int sofar;
	while( pl=fgets(line,499,fi))
	{
		for(char * l=line; *l; l++) if(*l < ' ') *l=0;
		fprintf(fo, "%s\n", line);

//		!@subs <number of lines> <word> <modifiers> [ <word> <modifiers> ... ]
		if( ! strncmp(line, "!@subs",6) ) 
		{
// 			printf("line %s\n", line);
			pl +=6;
			int numlines;
			int numsubs=0;
			int words[maxFitVariables];
			int bs[maxFitVariables];
			double adds[maxFitVariables];
			int clike = -1, slike=-1;
			sscanf(pl, "%d%n", &numlines, &sofar);
			pl += sofar;
			while(sscanf(pl, "%d%s%n", words+numsubs, what, &sofar) == 2)
			{
				pl += sofar;
				double minus = 1.;
				int b=0;
				char * w= what;
				if(*w == '+')  { minus=  1.;}
				if(*w == '-')  { minus= -1.;}
				if(*w == 'p')  { minus= phFac;   if(!has_beta) { fprintf(stderr, "beta is not defined for p modifier\n"); myexit(-1); }}
				if(*w == 'q')  { minus= 3.*phFac;if(!has_beta) { fprintf(stderr, "beta is not defined for p modifier\n"); myexit(-1); }}
				if(*w == 'b')  { b= 1; if(!has_cath) { fprintf(stderr, "cath is not defined for p modifier\n"); myexit(-1); }}
				if(*w == 'c')  { clike=-2;}
				if(*w == 's')  { slike=-2;}
				w++;
				int varno = -1;
				for(int i=0; i<nvar; i++)
				{
					if(  ! strcmp(w, var_name[i]) )
					{
						varno = i;
						break;
					}
				}
				if( var_active[varno] )
				{
					adds[numsubs] = minus*pnow[var2parm[varno]];
					bs[numsubs] = b;
					numsubs++;
					if(clike == -2) clike=varno;
					if(slike == -2) slike=varno;
				}
			}
// 			for(int i=0; i<numsubs; i++)
// 				printf("subs %d lines %d word %d adds %e bs %d\n", numsubs, numlines, words[i], adds[i], bs[i]);



			for(int i=0; i<numlines; i++)
			{
				pl=fgets(line,499,fi);
				char * space[80];  char * word[80];
				int nwords = splitline(line,  space,  word);
// 				printf("nwords %d\n", nwords);
				if(clike >= 0 )
				{
					if( slike < 0)
					{
						fprintf(stderr,"c-parameter set but no s-parameter\n");
						myexit(-1);
					}
					double oldvalueP, oldvalueA;
					int st=sscanf(word[4],"%lf",&oldvalueP);
					if(st != 1) { fprintf(stderr,"not a number in substitution:\n%s\n",line ); myexit(-1);}
					st=sscanf(word[5],"%lf",&oldvalueA);
					if(st != 1) { fprintf(stderr,"not a number in substitution:\n%s\n",line ); myexit(-1);}

					double a0 = oldvalueA * cos(M_PI/180.*oldvalueP);
					double b0 = oldvalueA * sin(M_PI/180.*oldvalueP);
					double a1 = a0 + pnow[var2parm[clike]];
					double b1 = b0 + pnow[var2parm[slike]];
					double Anew = sqrt(a1*a1 +  b1*b1);
					double Pnew = atan2(b1,a1);
					for(int iw = 0; iw < 4; iw++)
					{
						fprintf(fo, "%s", space[iw]);
						fprintf(fo, "%s", word[iw]);
					}
					fprintf(fo, " %20.10e     %20.10e ", Pnew, Anew);
					for(int iw = 6; iw < nwords; iw++)
					{
						fprintf(fo, "%s", space[iw]);
						fprintf(fo, "%s", word[iw]);
					}
				}
				else
				for(int iw = 0; iw < nwords; iw++)
				{
					fprintf(fo, "%s", space[iw]);

					int iv;
					for(iv = 0; iv < numsubs; iv++)
					{
						if(words[iv] == iw) goto found;
					}
					fprintf(fo, "%s", word[iw]);
					continue;
			found:		
					double oldvalue;
					int st=sscanf(word[iw],"%lf",&oldvalue);
					if(st != 1) { fprintf(stderr,"not a number in substitution:\n%s\n",line ); myexit(-1);}
					if(bs[iv])
					{
						double cath_s1 = cath_size+adds[iv];
						double new_field= cath_field*cath_size*cath_size/(cath_s1*cath_s1);
// 						printf("cath_size=%f add=%f cath_field=%f cath_s1=%f new_field=%f\n",
// 								cath_size,adds[iv],cath_field,cath_s1,new_field);
						fprintf(fo, " %20.10e ", new_field);
					}
					else
					{
						fprintf(fo, " %20.10e ", adds[iv]+oldvalue);
					}
				}
				fprintf(fo,"\n");
				for(int k=0; k<nwords; k++)
				{
					free(word[k]);
					free(space[k]);
				}
			}
		}
	}
	fclose(fi);
	fclose(fo);

        fi = fopen(dirTempInputFileName, "r");
	if(!fi)
	{
		fprintf(stderr, "cant open %s\n", dirTempInputFileName);
		myexit(-1);
	}
        fo = fopen(dirStrippedInputFileName, "w");
	if(!fo)
	{
		fprintf(stderr, "cant open %s\n", dirStrippedInputFileName);
		myexit(-1);
	}
	while( fgets(line,499,fi))
	{
		if( !strncmp(line, "!@if", 4) )
		{
			char ifFlag1[20],  ifFlag2[20];
			sscanf(line+4, "%s", ifFlag1);
			printf("skip %s\n", ifFlag1);
			int onOff=-1;
			for(int i=0; i<nswitches; i++)
			{
				if( !strcmp(ifFlag1, switchName[i]) ) onOff=switchFlag[i];
			}
			if(onOff < 0)
			{
				fprintf(stderr, "switch '%s' not found \n", ifFlag1);
				myexit(-1);
			}
			if(!onOff)
			{
				while( fgets(line,499,fi))
				{
					if( !strncmp(line, "!@endif", 7) )
					{
						sscanf(line+7, "%s", ifFlag2);
						printf("skip %s \n", ifFlag2);
						if( !strcmp(ifFlag1, ifFlag2)) break;
					}
				}
			}
		}

		if(line[0] == '!') continue;
		for(char * l=line; *l; l++) if(*l > ' ') { fprintf(fo,"%s", line); break;}	// find the first non-blank non-control char and print the rest
	}
	fclose(fi);
	fclose(fo);

	char command[200];
#if 1
	printf("dir is %s, input is %s\n", cwd, inputFileName);
	sprintf(command,"rm -f calcdirs/%03d/SCGRID", resource);
// 	printf("execute %s\n", command);
	system(command);

	sprintf(command,"cd calcdirs/%03d ; %s  %s", resource, preParmela, tempInputFileName);
// 	printf("execute %s\n", command);
	system(command);


	sprintf(command,"cd calcdirs/%03d ; %s  %s/parmela.exe  %s", resource, Wine, ParmelaDir, strippedInputFileName); // here comments are stripped
// 	printf("execute %s\n", command);
	system(command);
	printf("done pm\n");
#else
	sprintf(command,"vi  %s", dirStrippedInputFileName); // here comments are stripped
	system(command);
	printf("done pm\n");
#endif




	char t2nam[300], t3nam[300];
	sprintf(t2nam, "calcdirs/%03d/tape2.t2", resource);
	sprintf(t3nam, "calcdirs/%03d/tape3.t3", resource);
	Pmla *pmla = new Pmla(t2nam, t3nam);

	// the very first time: get the number of particles at the start of the beam line.
	if(ngood < 0)
	{
		pmla->getElement(1);
		ngood = pmla->ngood;
	}




	pmla->slices=50;
	int oldne= -1;
	double goal=0.;
	if(!onetime) for(int i=0; i<ngoals; i++)
	{
		if(fit_ne[i] != oldne) 	// do the getElement onbly once for each element
		{
			if(fit_ne[i] > pmla->nplots)
			{
				fprintf(stderr,"goal at %d > nplots %d\n", fit_ne[i] , pmla->nplots);
				fprintf(flog  ,"goal at %d > nplots %d\n", fit_ne[i] , pmla->nplots);
				delete pmla;
				myexit(-1);
			}
			pmla->getElement(fit_ne[i]);
// 			if(pmla->ngood < ngood)	// if we have lost any particles dont count this run
			if(pmla->ngood < ninput)	// if we have lost any particles dont count this run
			{
				printf("particles lost, ngood = %d, should be %d\n", ngood, pmla->ngood);
				fprintf(flog, "particles lost, ngood = %d, should be %d\n", ngood, pmla->ngood);
				if(nerror) *nerror=1;
				return 1e33;
			}
                	switch(cutWhat)
                	{   
                        	case 0:
                                	pmla->twiss(7);
                                	break;
                        	case 1:
                                	pmla->cutDP(7, cut/2., cut/2.);
                                	break;
                        	case 2:
                                	pmla->cutHeadTail(7, cut/2., cut/2.);
                                	break;
                	}   

			pmla->sliceEmit(0);
			oldne=fit_ne[i];
		}

		if( strcmp(fit_what[i], "slenv") )
		{
			fit_this[i] = pmla->getValue(fit_what[i]);
		}
		else
		{
			// this needs work, later
			sprintf(command,"cd calcdirs/%03d ; %s -m -9 %s slice.list.%d  %s", resource, slenvProg, slenvCavt, fit_ne[i], slenvLatt);
			printf("execute %s\n", command);
			fflush(stdout);
			system(command);

			char latt1pp[300];
			sprintf(latt1pp, "calcdirs/%03d/lattice1.txt", resource);
        		FILE * fi = fopen(latt1pp, "r");
			if(!fi)
			{
				fprintf(stderr, "cant open %s\n", "slenv_result.txt");
				myexit(-1);
			}
			char word[100]; word[0]=0; int ret=0;
			while(  strcmp(word, "result") &&  ret != EOF  ) ret=fscanf(fi,"%s", word);
			if(ret == EOF)
			{
				fprintf(stderr, "no result in lattice1.txt\n");
				myexit(-1);
			}
			fscanf(fi,"%lf",fit_this+i);
			fclose(fi);
		}
		printf("value %s %e\n", fit_what[i], fit_this[i]);
		double diff = (fit_val[i]-fit_this[i])*fit_weight[i];
		goal += square(diff);
		fprintf(flog, "value %s at %d = %20.10e\n", fit_what[i], fit_ne[i], fit_this[i]);
	}
	delete pmla;





        pthread_mutex_lock( &rtxMutex );
	time_t tim = time(0);
	char * timc =ctime(&tim);
	timc[17]=0;
	char * nobest= strdup("   ");
	char * isbest= strdup("###");
	char * showbest=nobest;
	if(goal < bestgoal)
	{
		showbest=isbest;
		bestgoal=goal;
		if(!onetime)
		{
			sprintf(command,"mv  -f %s pmla_best.inp", dirTempInputFileName);
			printf("execute %s\n", command);
			system(command);
		}
		sprintf(command,"mv  -f calcdirs/%03d/lattice1.txt slenv_best.latt", resource);
		system(command);
		sprintf(command,"cp  -f calcdirs/%03d/tape* .", resource);
		printf("execute %s\n", command);
		system(command);
		sprintf(command,"cp  -f calcdirs/%03d/SAV* .", resource);
		printf("execute %s\n", command);
		system(command);
		sprintf(command,"mv  -f calcdirs/%03d/OUTPAR.TXT* .", resource);
		printf("execute %s\n", command);
		system(command);
// 		system(main5);
	}
	fprintf(flog, "%sbest %5s %2d %18.10e this %18.10e<- ", showbest, timc+12, resource, bestgoal, goal);
	for(int i=0; i<nparm; i++)  fprintf(flog, "%18.10e ", p[i]);
	fprintf(flog, ": ");
	for(int i=0; i<ngoals; i++)  fprintf(flog, "%18.10e ", fit_this[i]);
	fprintf(flog, "\n");
	fflush(flog);

	for(int j=0; j<nparm; j++) para[nrun][j] = p[j];
	for(int j=0; j<ngoals; j++) subresu[nrun][j] = fit_this[j];
	resu[nrun++] = goal;

	// update recover  -- needs work
	FILE * fpr = fopen("recover.hex", "w");
	if(!fpr)
	{
		fprintf(stderr,"cant open %s for write\n", "recover.hex");
		myexit(-1);
	}

	fprintf(fpr, "%d %d %s\n",  nparm, ngoals, inputFileName);
	for(int i=0; i<nrun; i++)
	{
		for(int j=0; j<nparm; j++) 
		{
			int * it=(int*) (para[i]+j);
			fprintf(fpr,"%8x %8x ", it[0], it[1]);
		}
		for(int j=0; j<ngoals; j++)
		{
			int * it=(int*) (subresu[i]+j);
			fprintf(fpr,"%8x %8x ", it[0], it[1]);
		}
		int * it=(int*) (resu+i);
		fprintf(fpr,"%8x %8x = %e\n", it[0], it[1],  resu[i]);
	}
	fclose(fpr);


	FILE * fs = fopen("plot", "r");
	if(fs)
	{
		fclose(fs);
		system("rm plot");
		system(main5);
		fprintf(stderr, "****made plot*****\n");
	}
        pthread_mutex_unlock( &rtxMutex );
	while(fs = fopen("pause", "r"))
	{
		fclose(fs);
		fprintf(stderr, "pause file found\n");
		sleep(30);
	}
        fs = fopen("stop", "r");
	if(fs)
	{
		fclose(fs);
		system("rm stop");
		fprintf(flog, "stop file found\n");
		myexit(0);
	}

	// turn of negative switches after the first run
	for(int i=0; i<nswitches; i++)
		if( switchFlag[i] < 0)  switchFlag[i] = 0;
   	updateCounter(goal,X);
	return goal;
}							




void M2condor::madfit(FILE * fo, FILE * fi, double pnow[])
{

	char line[512];
//  copy mad part to temp file with substitutions
	char tmpm[300];	// temp  file
	char tmp2[300];	// temp1 file
	int ls=strlen(inputFileName);
	strcpy(tmpm, inputFileName);
	strcpy(tmp2, inputFileName);
	strcpy(tmpm+ls-4, "_temp.mad");
	strcpy(tmp2+ls-4, "_temp2.mad");
	FILE * fp = fopen(tmp2,"w");
	if(!fp)
	{
		fprintf(stderr,"cant open %s for write\n", tmp2);
		myexit(-1);
	}

	while(1)
	{
		int c= fgetc(fi);
		if( c == EOF) break;
		if(c == '\r') continue;
		if(c == '\t') c = ' ';
		if( c == ']') break;

		if( c == '$')
		{
			double minus = 1.;
			double dummy;
			int b=0;
			int c= fgetc(fi);
			if(c == '-')  {c= fgetc(fi); minus= -1.;}
			if(c == 'p')  {c= fgetc(fi); minus= phFac;}
			if(c == 'q')  {c= fgetc(fi); minus= 3.*phFac;}
			if(c == 'b')  {c= fgetc(fi); b= 1;}
			c -=  '0';
			double add = minus*pnow[c];
			char num[50];
			fscanf(fi, "%s", num); 
			while(num[0] == '$')
			{
				double minus = 1.;
				int c= num[1];
				if(c == '-')  {c= num[2]; minus= -1.;}
				if(c == 'p')  {c= num[2]; minus= phFac;}
				if(c == 'q')  {c= num[2]; minus= 3.*phFac;}
				c -=  '0';
				add += minus*pnow[c];
				fscanf(fi, "%s", num); 
			}
			sscanf(num,"%lf", &dummy);

			if(b)
			{
				double cath_s1 = cath_size+add;
				double new_field= cath_field*cath_size*cath_size/(cath_s1*cath_s1);
				printf("cath_size=%f add=%f cath_field=%f cath_s1=%f new_field=%f\n", cath_size,add,cath_field,cath_s1,new_field);
				fprintf(fp, "%16.10e", new_field);
			}
			else
			{
				add += dummy;
				fprintf(fp, "%16.10e", add);
			}
// 			printf("add =  %20.10e ", add);
		}
		else

		fputc(c, fp);

	}
	fclose(fp);

// discard the parmela translation of the mad part
	while(1)
	{
		int c=fgetc(fi);
		if(c == EOF ) break;
		if(c == ']') break;
	}
	fgets(line,499,fi);

// copy (substituted) tmp2 file back into the to  the parmela input file (for saving pmal_best.inp)
// Also  stripp the !mad marker and copy into pmla_temp.mad 
	fp = fopen(tmp2,"r");
	if(!fp)
	{
		fprintf(stderr,"cant open %s for read\n", tmp2);
		myexit(-1);
	}
	FILE * fq = fopen(tmpm,"w");
	if(!fq)
	{
		fprintf(stderr,"cant open %s for write\n", tmpm);
		myexit(-1);
	}
	while( fgets(line,499,fp))
	{
		if( strncmp(line, "!mad ", 5)) break;
		fprintf(fq,"%s", line+5);
		fprintf(fo,"%s", line);
	}
	fprintf(fo,"!]\n");
	fclose(fp);
	fclose(fq);





	char command[200];
	sprintf(command,"~/MAD8windows/mad8 < %s > pmla_temp.echo", tmpm);
	printf("execute %s\n", command);
	system(command);

	printf("mv print pmla_temp.madout\n");
	system("mv print pmla_temp.madout");

	char tmp3[300]; strcpy(tmp3, tmpm); tmp3[strlen(tmp3)-4]=0;
	sprintf(command, "~/MAD8windows/madPict.exe -madout -pmla %s", tmp3);
	printf("execute %s\n", command);
	system(command);

//  now copy the translated lattice into the output file
	strcat(tmp3, ".pmla");
	fp = fopen(tmp3,"r");
	if(!fp)
	{
		fprintf(stderr,"cant open %s for read\n", tmp3);
		myexit(-1);
	}
	while(1)
	{
		int c= fgetc(fp);
		if(c == '\r') continue;
		if( c == EOF) break;
		fputc(c, fo);
	}
	fclose(fp);
	fprintf(fo,"!]\n");

}





main(int argc, char ** argv, char ** env)
{
	int recover =0;
	int onetime =0;
	int display =0;
        while(argc > 2)
        {
		if( !strcmp(argv[1], "-r") )
		{
			recover=1;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-d") )
		{
			display=1;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-o") )
		{
			onetime=0;
			argv++; argc--;
		}
		if( !strcmp(argv[1], "-p") )
		{
			onetime=1;
			argv++; argc--;
		}
	}

        if(argc != 2)
        {
                fprintf(stderr, "usage: %s <input file>\n", argv[0]);
                myexit(-1);
        }


	pid_t myPid = getpid();
	FILE * fpr=fopen(thislock, "r");
	if(fpr)
	{
		fprintf(stderr, "file %s exists, mainp may already be running\n", thislock);
		exit(-1);
	}
	fpr=fopen(thislock, "w");
	if(!fpr)
	{
		perror(thislock);
		exit(-1);
	}
	fprintf(fpr, "%d", myPid);
	fclose(fpr);

	unlink("stop");
	unlink("pause");
	unlink("plot");







	char * home = getenv("HOME");
	if(!home)
	{
		fprintf(stderr, "environment variable HOME not set\n");
		myexit(-1);
	}
	
	
	char xfile[300];
	char rfile[300];
	
	
	if(display)
	{
		for(display=60; display < 65; display++)
		{
	
			sprintf(rfile,"/tmp/.RTX%d-lock", display);
			sprintf(xfile,"/tmp/.X%d-lock", display);
			FILE * fpr=fopen(rfile, "r");
			if(fpr)
			{
				pid_t rpid = -1;
				pid_t xpid = -1;
				fscanf(fpr, "%d%d", &rpid, &xpid);
				fclose(fpr);
				if(rpid > 0)
				{
					int err = kill(rpid, 0);
					if(err ==  0) continue;		// display in use
				}
				unlink(rfile);
			}
			goto isset;
		}
		fprintf(stderr, "both displays are in use\n");
		myexit(-1);
				

isset:		printf("display is %d\n", display);
		FILE * fpx=fopen(xfile, "r");
		if(fpx) 
		{
			pid_t xpid = -1;
			fscanf(fpx, "%d", &xpid);
			fclose(fpx);
			if(xpid > 0)
			{
				int err = kill(xpid, SIGTERM);
				if(err ==  0) unlink(xfile);
				else sleep(4);
			}
		}
	
	
	
		char * arg[6];
		arg[0] = strdup("/usr/bin/Xvfb");
		arg[1] = strdup("-sp");
		arg[2] = new char[200]; strcpy(arg[2], home); strcat(arg[2], "/SecurityPolicy");
		arg[3] = new char[200]; sprintf(arg[3], ":%d", display);
		arg[4]=0;
		setenv("DISPLAY", arg[3],1);
	
    		pidXvfb= vfork();
    		if(pidXvfb == 0)
    		{
         		execve(arg[0], arg, env);
            		_exit(2);
     		}
		

	             
		FILE * fpr=fopen(rfile, "w");
		if(fpr)
		{
			fprintf(fpr, "%d %d\n", myPid, pidXvfb);
			fclose(fpr);
		}
	}





	char * OS = getenv("OS");
	if(OS && !strcmp(OS, "Windows_NT") )
	{
		Wine = strdup("");
		ParmelaDir= strdup("/cygdrive/i/LANL");
		preParmela= strdup("~/RTX/preParmela.exe");
		slenvProg= strdup("~/SLENV/mii.exe");
		main5 = strdup("~/RTX/main5.exe > main5.out");
		whatsystem = 0;
	}
	else
	{
		Wine =strdup("wine");
		ParmelaDir=getenv("ParmelaDir");
		preParmela=strdup("~/RTX/preParmela");
		slenvProg=strdup("~/SLENV/mii");
		main5 = strdup("~/RTX/main5 > main5.out");
		setenv("DISPLAY", "localhost:60.0", 1);
		whatsystem = 1;
	}

	char inpfilename[300];
        strcpy(inpfilename,argv[1]);


    	double rhoStart=1e-0, rhoEnd=1e-5;
    	int niter=100000;


    	M2condor *mof = new M2condor(inpfilename, recover, onetime);
    	ObjectiveFunction *of = mof;

	// set number of threads. argument = 0 means look at /proc/cpuinfo 
	// negative means /proc/cpuinfon but reserve -n threads for other things
	int maxthreads = of->setMaxThreads(0);
	int nparm=mof->nparm;
	int neededThreads = nparm*(nparm-1)/2;
	if(neededThreads < nparm) neededThreads= nparm;
	printf("available threads %d\n", maxthreads);
	printf("parameters %d needed threads %d\n",  mof->nparm, neededThreads);
	if(maxthreads > neededThreads) maxthreads = of->setMaxThreads(neededThreads);
	if(onetime)
	{
		maxthreads =1;
		maxthreads = of->setMaxThreads(maxthreads);
	}
	printf("will use %d threads\n", maxthreads);

	// create resources
	system("rm -rf calcdirs");
	printf("rm -rf calcdirs");
	system("mkdir calcdirs");
	for(int i=0; i<maxthreads; i++)
	{
		char cmd[200];
		sprintf(cmd, "mkdir calcdirs/%03d", i);
		printf("%s\n", cmd);
		system(cmd);
		sprintf(cmd, "cp tape* SAV* calcdirs/%03d", i);
		printf("%s\n", cmd);
		system(cmd);
	}       


	// if there are *.T7 files or LANL.INI in the directory, create symbolic links in each calc dir
        DIR * d = opendir(".");
        if(!d)
        {   
                fprintf(stderr, "Can't open '.'\n");
                exit(-1);
        }   

        while( struct dirent * dd = readdir(d))
        {   
		char * s = dd->d_name;
		int l = strlen(s);
		if( (!strcmp(s+l-3, ".T7"))  || (!strcmp(s, "LANL.INI"))  )
		{
                	printf(" link %s\n", s);
			char t[300]; strcpy(t,"../../"); strcat(t, s);		// symlink is tricky. "t" is the name of the file relative to the destination dir of the link
				
			for(int i=0; i<maxthreads; i++)
			{
				char r[300];
				if(whatsystem)
				{
					sprintf(r, "calcdirs/%03d/%s", i, s);
					symlink(t, r);
				}
				else
				{
					sprintf(r, "cp %s calcdirs/%03d", s, i);  	// fucking windows does not understang symlinks, need to make a copy
					system(r);
				}
			}
		}
        }   
        closedir(d);




	if(onetime)
	{
		mof->eval();
	}
	else
	if(mof->stepping)
	{
		mof->step(0);
	}
	else
	{
		CONDOR(rhoStart, rhoEnd, niter, of);
	}


    	delete of;
	unlink(thislock);

}
