#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }

const double clight=2.997824562e8;      // m/s
const double elcharge=1.6022e-19;       // Coulomb
const double elmass = 0.51104;


const char * types[] = {"eof", "drift", "quad", "solenoid", "bend", "cavity", "rotate", "zero"};

int main(int argc, char ** argv)
{
	int dodisp=0;
	int dogpt =0;
	FILE *  fgpt = 0;
	FILE *  fscr = 0;
	if( !strcmp(argv[1], "-d"))
	{
		dodisp=1;
		argv++; argc--;
	}
	if( !strcmp(argv[1], "-g"))
	{
		dogpt=1;
		argv++; argc--;
	}
	FILE * fi = fopen(argv[1],"r");
	if(!fi)
	{
		perror(argv[1]);
		exit(-1);
	}

//      char cwd[200];
//      getcwd(cwd, 200);
// 	for(char * p= cwd+28; *p; p++) if(*p == '/') *p ='_';

	char filename[200];
	strcpy(filename, "survey_");
	strcat(filename, argv[1]);
	int l = strlen(filename)-4;
	strcpy(filename+l, ".agr");
	fprintf(stderr, "outfile = %s\n", filename);

	FILE * fo = fopen(filename,"w");
	if(!fo)
	{
		perror(filename);
		exit(-1);
	}

	strcpy(filename, "center_");
	strcat(filename, argv[1]);
	l = strlen(filename)-4;
	strcpy(filename+l, ".agr");
	fprintf(stderr, "centerfile = %s\n", filename);

	FILE * fc = fopen(filename,"w");
	if(!fc)
	{
		perror(filename);
		exit(-1);
	}






// 	type = 0    	eof
// 	type = 1	drift
// 	type = 2	quad
//	type = 3	solenoid 
//	type = 4	bend
//	type = 5	cavity
//	type = 6	rotate
//	type = 7	!@surveyzero
//	length and fields are converted from cm and gauss to m and kG when read from file
//	energies are in MeV



	int n=0;
	int driftnum=0;
	int quadnum=0;
	int bendnum=0;
	int solnum=0;
	int cavnum=0;
	int elnum=0;
	double x=0.;
	double y=0.;
	double xdir=1.;
	double ydir=0.;

	int neo =0;
	int surveyzero=0;
	int save=0;
	double loc =0;
	double len, app, pri, stre, angl, e1, e2, phas;
	char line[1255];
	char word[255];
	int type[1000];
	double  length[1000], strength[1000], sloc[1000], edge1[1000], edge2[1000];
	char elname[1000][6];
	while (fgets(line, 255, fi))
	{
		for(char * l=line; *l; l++)
		{
			if( *l >= 'A' && *l <= 'Z') *l += 32;
			if( *l == ',') *l=' ';
		}
		int nw = sscanf(line,"%s", word);
		if(nw != 1) continue;


		if( !strcmp( word, "!@surveyzero"))
		{
			printf("found surveyzero n=%d\n", n);
			type[n++]=7;
			surveyzero=1;
			continue;
		}


		if( !strcmp( word, "!@survey"))
		{
			double angle;
			int nw = sscanf(line,"%s%lf%lf%lf", word, &x, &y, &angle);
			if(nw != 4)
			{
				fprintf(stderr, "error in %sShould be !@survey <x> <y> <angle>\n", line);
				exit(-1);
			}
			angle *= M_PI/180.;
			xdir = cos(angle);
			ydir = sin(angle);
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
			type[n]=1;
			sloc[n] = loc;
			length[n]=len;
			strength[n++]=0;
			sprintf(elname[elnum++], "d%d",driftnum++);
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "quad"))
		{
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri,&stre);
			len /= 100.;
			loc += len;
			type[n]=2;
			sloc[n] = loc;
			length[n]=len;
			stre *= 1e-2;		// convert from G/cm to T/m
			strength[n++]=stre;
			sprintf(elname[elnum++], "q%d",quadnum++);
			neo++;
			continue;
		}

		if( !strcmp( word, "solenoid"))
		{
			printf("%s", line);
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri, &stre);
			printf("%s   %lf   %lf   %lf   %lf\n\n", word, len, app, pri, stre);
			len /= 100.;
			loc += len;
			type[n]=3;
			sloc[n] = loc;
			length[n]=len;
			stre /= 10000.;		// convert from G to tesla
			strength[n++]=stre;
			sprintf(elname[elnum++], "s%d",solnum++);
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
			loc += len;
			type[n]=4;
			sloc[n] = loc;
			length[n]=len;
			strength[n]=angl;
			edge1[n]=e1;
			edge2[n]=e2;
// 			printf("bend length %e strength %e\n", length[n], strength[n]);
			n++;
			sprintf(elname[elnum++], "b%d",bendnum++);
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "cell"))
		{
			sscanf(line,"%s%lf%lf%lf%lf%lf", word, &len, &app, &pri,&phas, &stre);
			len /= 100.;
			neo++;
			loc += len;
			type[n]=5;
			sloc[n] = loc;
			length[n]=len;
			strength[n++]=0;
			sprintf(elname[elnum++], "c%d",cavnum);
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

		if( !strcmp( word, "rotate"))
		{
			sscanf(line,"%s%lf%lf%lf%lf", word, &len, &app, &pri,&stre);
			type[n]=6;
			length[n]=0;
			sloc[n] = loc;
			strength[n++]=stre*M_PI/180.;
			neo++;
// 			printf("%-20s\t%d\t%f\t%f\n", word, neo, len, loc);
			continue;
		}

// 		printf("%s not recognized\n", word);
// 		exit(-1);

	}
	fclose(fi);
	type[n]=0;

	printf("n = %d\n", n);

	if(surveyzero)
	{
		x=y=ydir=0.;xdir=1.;
		for(int i=0; i<n; i++)
		{
			if(type[i] == 4)
			{
				double x0=x;
				double y0=y;
				double rho=length[i]/strength[i];
				double c = cos(strength[i]);
				double s = sin(strength[i]);
				x = x0 - (c -1.) *rho*ydir + s*rho*xdir; 
				y = y0 + (c -1.) *rho*xdir + s*rho*ydir; 
				double xdir1 = xdir*c + ydir*s;
				double ydir1 = ydir*c - xdir*s;
				xdir = xdir1;
				ydir = ydir1;
			}
			else
			if(type[i] == 7)
			{
				double x1 =  x*xdir + y*ydir;
				double y1 = -x*ydir + y*xdir; 
				x = -x1;
				y = -y1;
				ydir = -ydir;
				printf("!@surveyzero at (%f, %f) (%f, %f)\n", x,y,xdir,ydir);
				break;
			}
			else
			{
				x += length[i]*xdir;
				y += length[i]*ydir;
			}
// 			printf("%e %e #%s\n", x, y,elname[i]);
		}
	}
	double xzero =x;
	double yzero =y;
	double xdirzero =xdir;
	double ydirzero =ydir;


// 	type = 0    	eof
// 	type = 1	drift
// 	type = 2	quad
//	type = 3	solenoid 
//	type = 4	bend
//	type = 5	cavity
//	type = 6	rotate
//	type = 7	!@surveyzero



	// hardwired parameters
	const double qwa[8] = { 0., 0., 0.05, 0.1, 0.1, 0.2, 0., 0.};	// width of a magnet





	fprintf(fo, "%e %e\n",      x, y);
	for(int i=0; i<n; i++)
	{
// 		printf("typ %d name %s\n", type[i], elname[i]);
		double qw = qwa[type[i]];
		switch(type[i])
		{
			case 1:		// drift
				x += length[i]*xdir;
				y += length[i]*ydir;
				fprintf(fo, "%e %e #%s\n", x, y,elname[i]);
				break;
			case 2:		//quad
			case 3:		//solenoid
			case 5:		//cavity
			{
				double x0=x;
				double y0=y;
				x += length[i]*xdir;
				y += length[i]*ydir;
				double p1x = x0+qw*ydir;
				double p1y = y0-qw*xdir;
				double p2x = x0-qw*ydir;
				double p2y = y0+qw*xdir;
				double p3x = x +qw*ydir;
				double p3y = y -qw*xdir;
				double p4x = x -qw*ydir;
				double p4y = y +qw*xdir;
				fprintf(fo, "%e %e #%s\n",  p1x, p1y,elname[i]);
				fprintf(fo, "%e %e\n",      p3x, p3y);
				fprintf(fo, "%e %e\n",      p4x, p4y);
				fprintf(fo, "%e %e\n",      p2x, p2y);
				fprintf(fo, "%e %e\n",      p1x, p1y);
				fprintf(fo, "%e %e\n",      p3x, p3y);
				fprintf(fo, "%e %e\n",        x,   y);
				double xc = x0+ length[i]*xdir*0.5;
				double yc = y0+ length[i]*ydir*0.5;
				fprintf(fc, "%-12s\t%f\t%f\n", types[type[i]],   x,   yc);
				break;
			}
			case 4:		//bend
			{
				double x0=x;
				double y0=y;
				double rho=length[i]/strength[i];
				double c = cos(strength[i]);
				double s = sin(strength[i]);
				x = x0 - (c -1.) *rho*ydir + s*rho*xdir; 
				y = y0 + (c -1.) *rho*xdir + s*rho*ydir; 

				double p1x = x0+qw*ydir;
				double p1y = y0-qw*xdir;
				double p2x = x0-qw*ydir;
				double p2y = y0+qw*xdir;
				double xdir1 = xdir*c + ydir*s;
				double ydir1 = ydir*c - xdir*s;
				xdir = xdir1;
				ydir = ydir1;
				double p3x = x+qw*ydir;
				double p3y = y-qw*xdir;
				double p4x = x-qw*ydir;
				double p4y = y+qw*xdir;
				fprintf(fo, "%e %e #%s\n",  p1x, p1y,elname[i]);
				fprintf(fo, "%e %e\n",      p3x, p3y);
				fprintf(fo, "%e %e\n",      p4x, p4y);
				fprintf(fo, "%e %e\n",      p2x, p2y);
				fprintf(fo, "%e %e\n",      p1x, p1y);
				fprintf(fo, "%e %e\n",      p3x, p3y);
				fprintf(fo, "%e %e\n",        x,   y);
				double dle = rho* tan(strength[i]/2.);
				double xc = x0+ dle*xdir;
				double yc = y0+ dle*ydir;
				fprintf(fc, "%-12s\t%f\t%f\n", types[type[i]],   x,   yc);
			}
		}
	}

	if(dodisp)
	{
		char command[1000];
		sprintf(command, "xmgrace -geometry 980x750+1000 -fixed 750 500 -noask  %s &", filename);
		system(command);
	}



	x = xzero;
	y = yzero;
	xdir = xdirzero;
	ydir = ydirzero;


	// hardwired parameters
	const double edge = 0.05 ;
	const double step =0.00005 ;
	const double dr   =0.001 ;
	const int rstep =20 ;

	if(dogpt)
	{
		char gptfile[200];
		strcpy(gptfile, argv[1]);
		strcpy(gptfile+strlen(gptfile)-3, "gpt");
		fgpt = fopen(gptfile,"w");
		if(!fgpt)
		{
			perror(gptfile);
			exit(-1);
		}

		strcpy(gptfile+strlen(gptfile)-3, "screen");
		fscr = fopen(gptfile,"w");
		if(!fscr)
		{
			perror(gptfile);
			exit(-1);
		}

		int nsol=0;
		int nquad =0;
		for(int i=0; i<n; i++)
		{
			printf("y %d %d %f %f\n", i, type[i], length[i], strength[i]);
		}
		for(int i=0; i<n; i++)
		{
			switch(type[i])
			{
				case 2:		//quad
				{
					fprintf(fgpt, "quad%02dField= %e ;\n", nquad++, strength[i]);
					break;
				}
				case 3:		//solenoid
				{
					fprintf(fgpt, "sol%02dField= %e ;\n", nsol++, strength[i] );
					printf("x %d %d %f %f\n", i, type[i], length[i], strength[i]);
					break;
				}
			}
		}



		nsol=0;
		nquad=0;
		double stot=0.;
		for(int i=0; i<n; i++)
		{
// 			printf("typ %d name %s\n", type[i], elname[i]);
			stot += length[i];
			switch(type[i])
			{

				case 1:		// drift
					x += length[i]*xdir;
					y += length[i]*ydir;
					fprintf(fscr, "screen(\"wcs\", %f, 0, %f,    %f, 0, %f,     0, 1, 0,       %f);\n",  y-ydir*stot, x-xdir*stot, xdir, ydir, stot);
					break;
				case 2:		//quad
				{
					double x0=x;
					double y0=y;
					x += length[i]*xdir;
					y += length[i]*ydir;
					fprintf(fgpt, "#quadrupole(\"wcs\", %f, 0, %f,    %f, 0, %f,  0, 1, 0, %f, quad%02dField);\n", (y+y0)/2., (x+x0)/2., xdir, ydir, length[i], nquad++);
					break;
				}
				case 5:		//cavity
				{
					double x0=x;
					double y0=y;
					x += length[i]*xdir;
					y += length[i]*ydir;
					fprintf(fgpt, "# %s(\"wcs\", %f, 0, %f,    %f, 0, %f,  0, 1, 0, l=%f s=%f);\n", types[type[i]], (y+y0)/2., (x+x0)/2., xdir, ydir, length[i], strength[i]);
					break;
				}
				case 3:		//solenoid
				{
					// make a box solenoid field map
					char filename1d[200];
					char filenamegdf[200];
					char filenamepa7[200];
					sprintf(filename1d, "%03dmmSol.1d", int(length[i]*1000.+0.5));
					sprintf(filenamegdf, "%03dmmSol.gdf", int(length[i]*1000.+0.5));
					sprintf(filenamepa7, "%03dmmSol.pa7", int(length[i]*1000.+0.5));
					FILE * fg = fopen(filename1d, "w");
					if(!fg)
					{
						perror(filename1d);
						exit(-1);
					}
					FILE * f7 = fopen(filenamepa7, "w");
					if(!f7)
					{
						perror(filenamepa7);
						exit(-1);
					}
 	 				fprintf(fg, "z Bz\n");
					double ltot =  length[i] + edge;
					int n = ltot/step;
					if( ! (n&1) ) n++;
					int m=n/2;
					fprintf(f7, "%e %e %d\n", 0., rstep*dr,  rstep);
					fprintf(f7, "%e %e %d\n", -ltot*100., ltot*100.,  2*m);
					double dz = (length[i]+edge)/(2.*m);

					for(int j= -m; j<= m; j++)
					{
						double z = j*dz; 
						double zz = fabs(z);
						int  in = ( (2.*zz) < (length[i]-edge));
	
						double Bz = in ?  1. : (1. - sin( M_PI*(zz-length[i]) /edge) )/2.;

						double Bzp  = in ? 0. : - cos( M_PI*(zz-length[i]) /edge) *M_PI/edge/2.;           if( z < 0) Bzp = - Bzp;
						double Bzpp = in ? 0. :   sin( M_PI*(zz-length[i]) /edge) *M_PI*M_PI/edge/edge/2.; if( z < 0) Bzpp = - Bzpp;

						if(!in) fprintf(fg, "%e %e\n", z, Bz);

						for(int k=0; k<=rstep; k++)
						{
							double r = k*dr;
							double Br1 = (- 0.5 *r * Bzp)     *1e4*strength[i] ;		// Br = -1/2 r Bz'(0) convert back to gauss
							double Bz1 = (Bz - 0.25 *r*r*Bzpp)*1e4*strength[i] ;		// Bz = Bz(0) - 1/4 r^2 Bz"(0)
							fprintf(f7, "%e %e\n", Br1, Bz1);
						}
					}
					fclose(fg);
					fclose(f7);
					char cmd[300];
					sprintf(cmd, "asci2gdf -o  %s %s", filenamegdf, filename1d);
					printf("%s\n", cmd);
					system(cmd);
	
	
					double x0=x;
					double y0=y;
					x += length[i]*xdir;
					y += length[i]*ydir;
					fprintf(fgpt, "map1D_B(\"wcs\", %f, 0., %f,    %f, 0., %f,  0., 1., 0., \"%s\", \"z\", \"Bz\", sol%02dField );\n", (y+y0)/2., (x+x0)/2., xdir, ydir, filenamegdf, nsol++);
					fprintf(fscr, "screen(\"wcs\", %f, 0, %f,    %f, 0, %f,     0, 1, 0,       %f);  #sol\n",  y, x-stot, xdir, ydir, stot);
					break;
				}
				case 4:		//bend
				{
					// proceed to the magnet middle
					{
					double x0=x;
					double y0=y;
					double rho=length[i]/strength[i];
					double c = cos(strength[i]/2.);
					double s = sin(strength[i]/2.);
					double x1 = x0 - (c -1.) *rho*ydir + s*rho*xdir; 
					double y1 = y0 + (c -1.) *rho*xdir + s*rho*ydir; 
					double xdir1 = xdir*c + ydir*s;
					double ydir1 = ydir*c - xdir*s;
					fprintf(fgpt, "#%s(\"wcs\", %f, 0, %f,    %f, 0, %f,  0, 1, 0, l=%f s=%f\n", "dipole", y1, x1, ydir1, xdir1, length[i], strength[i]);
					}
	
					// do the whole magnet
					double x0=x;
					double y0=y;
					double rho=length[i]/strength[i];
	// 				printf("rho %e str %e\n", rho, strength[i]);
					double c = cos(strength[i]);
					double s = sin(strength[i]);
					x = x0 - (c -1.) *rho*ydir + s*rho*xdir; 
					y = y0 + (c -1.) *rho*xdir + s*rho*ydir; 
	
					double xdir1 = xdir*c + ydir*s;
					double ydir1 = ydir*c - xdir*s;
					xdir = xdir1;
					ydir = ydir1;
				}
			}
		}
	}
}



