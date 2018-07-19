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




main(int argc, char ** argv)
{
	char path[200];
	int n1= -1, n2 = -1;

	char cwd[200];
	getcwd(cwd, 200);

	double coolenergy = 1.6;  	// MeV
	int cutWhat=0;
	double cut= 0.; 
	char * cutparm1=0;
	char * cutparm2=0;
	while(argc > 2 && argv[1][0] == '-') 
	{
		if( ! strcmp(argv[1], "-cool")) 
		{
			coolenergy = atof(argv[2]);
			argv += 2; argc -= 2; continue;
		}
		if( ! strcmp(argv[1], "-cp"))
		{
			cutWhat = 1;
			cut = atof(argv[2]);
			cutparm1= argv[1]; cutparm2 = argv[2];
			argv += 2; argc -= 2; continue;
		}
		if( ! strcmp(argv[1], "-cl"))
		{
			cutWhat = 2;
			cut = atof(argv[2]);
			cutparm1= argv[1]; cutparm2 = argv[2];
			argv += 2; argc -= 2; continue;
		}
		if( ! strcmp(argv[1], "-ct"))
		{
			cutWhat = 3;
			cut = atof(argv[2]);
			cutparm1= argv[1]; cutparm2 = argv[2];
			argv += 2; argc -= 2; continue;
		}
	}
	if( cutWhat)
	{
		if( cut < 0. || cut > 1.) 
		{
			fprintf(stderr, "cut fraction must be between 0 and 1, 0 means no cut\n");
			exit(-1);
		}
		strcat(cwd, " ");
		strcat(cwd, cutparm1);
		strcat(cwd, " ");
		strcat(cwd, cutparm2);
	}

	if(argc > 3)
	{
		fprintf(stderr, "usage: %s [ -c[plt] <fraction> ]  [<first>  [last]  ] \n", argv[0]);
		exit(-1);
	}
		
	if(argc >= 2) n1 = atoi(argv[1]);
	if(argc >= 3) n2 = atoi(argv[2]); else n2 = n1+1;



	char elemfilename[100];
	strcpy(elemfilename,"TAPE2.T2");
	char timefilename[100];
	strcpy(timefilename,"TAPE3.T3");
	Pmla *p = new Pmla(elemfilename, timefilename);
	char zoutfilename[100];
	strcpy(zoutfilename,"OUTPAR");
	strcat(zoutfilename,".TXT");
	p->zout(zoutfilename, "line.agr", 2e-7);
	p->zout(zoutfilename, "line100.agr", 1e-4);
	p->zout(zoutfilename, "lined.agr", 1e-3);
	p->zout(zoutfilename, "linee.agr", 1e-1);

	p->coolenergy = coolenergy;

	int nplots = p->nplots;

// 	p->saveDelta();
	p->slices=20;
// 	p->readTimeFile();

// 	p->twiss1();

	int print = (argc >= 2) ? 7 : 0;

	if(n1 < 0) { n1 = 1; n2 =  nplots; }
	if(n2 >=  nplots)  n2 =  nplots; 





	for(int icut =0; icut < 2; icut++)
	{

		if( (!icut) && (!cutWhat) ) continue; 
		for(int ne =n1; ne < n2; ne++)
		{
			if(ne >= nplots) break;
// 			printf("\n\n%d \n", ne);
			if(p->getElement(ne)) continue;
// 			p->testcord(ne, 1., lr, lr);
// 			p->scatter("scatter");

// 			det =p->linMatrix();
// 			if(argc >= 2) p->twiss1(7);

			switch(cutWhat)
			{
				case 0:
					p->twiss(print);
					break;
				case 1:
					p->cutDP(print, cut/2., cut/2.);
					break;
				case 2:
// 					p->cutHeadTail(7, cut/2., cut/2.);
					p->cutHeadTail(print, 0., cut);
					break;
				case 3:
					p->cutTempT(print, cut/2., cut/2.);
					break;
			}





//    			unrotate(int method, int remove_orbit, double rm)
// 			double rm1=p->unrotate(0,0, 0.);
//  			double rm2=p->unrotate(1,0, 0.);
// 			double rm1=p->unrotate(0,1, 0.); // remove orbit, changs method
// 			printf("rot= %e\n", p->rot[ne]);
// 			double rm2=p->unrotate(1,1, 0.);
// 			p->twiss(7);
// 			p->putElement(ne, 0);
// 
// 
// 

// 			if(argc == 2) p->impactTdump();
			if(argc == 2) p->sliceEmit(1);
			if(argc == 2) p->dump();
			if(argc >= 2) p->dumpLongit(cwd, cut);
// 			if(argc == 2)
// 			{
// 				double cut; int steps;
// 				printf("cut  (0.99) , steps? "); scanf("%lf%d", &cut, &steps);
// 				p->cutEmittance(cut, steps);
// 			}
// 			if(argc == 2) p->slicePhasePlot("slice_rrp.agr", 1);

// 			p->cutHeadTail(cut);
// 			p->cutDP(cut);
		}
		if(argc > 1) exit(0);
	
		Xmgr * xmgr;
		char title[200];
		char xfile[200];

// 		for(int i=2; i<nplots; i++) 
// 		{
// 			if( p->zloc[i] > 0. && (p->epsx_n[i] <= 0. || p->epsx_n[i] <= 0.)  )
// 			{
// 				nplots = i;
// 			}
// 		}	


// calculate invariant envelope

// 		double q=3.2e-9;
// 		double I0=17000;
// 		double sigma[nplots], siginv[nplots], sigina[nplots];
// 		for(int n=1; n<nplots-1; n++)
// 		{
// 			sigma[n]=sqrt(   (p->epsx[n])*(p->betx[n])+(p->epsy[n])*(p->bety[n]));
// 			double I=q*clight/(2.*sqrt(2.*M_PI)*p->sigma_l[n]/1.18);
// 			double gm = p->ener[n-1]/elmass+1.;
// 			double ga = p->ener[n  ]/elmass+1.;
// 			double gp = p->ener[n+1]/elmass+1.;
// 			double gammPrime=(gp - gm)/(p->zloc[n+1]-p->zloc[n-1]);
// 			if(gammPrime < 1.e-4) siginv[n] = 0.; else siginv[n] = 4.*sqrt(8.*I/(3.*I0*ga))/gammPrime;
// 			if(siginv[n] > 0.05) siginv[n]=0.;
// 			double gammPrimeA=(p->ener[p->gunEnd]/elmass)/(p->zloc[p->gunEnd]);	// gamma starts with 1
// 			if(gammPrimeA < 1.e-4) sigina[n] = 0.; else sigina[n] = 4.*sqrt(8.*I/(3.*I0*ga))/gammPrimeA;
// 			if(sigina[n] > 0.05) sigina[n]=0.;
// 			printf("z=%f i=%e g=%e %e %e l=%e gap=%e sig=%e\n", p->zloc[n],I, gm, ga, gp, p->sigma_l[n], gammPrime, siginv[n]);
// 		}
// 
// 
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "sigma.agr");
//      	xmgr = new Xmgr(xfile, "envelope", "Path length [m]","sigma [m]");
// 		xmgr->subTitle(cwd);
// 		xmgr->plot(p->zloc+1, sigma+1,  p->rotat+1, nplots-2);	//   actual envelope
// 		xmgr->plot(p->zloc+1, siginv+1, p->rotat+1, nplots-2);	//   invariant envelope 
// 		xmgr->plot(p->zloc+1, sigina+1, p->rotat+1, nplots-2);	//   invariante envelope with avarage gamma-prime
// 		delete xmgr;
// 


// 		strcpy(title, "Slice Betas"); 
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "beta.agr");
//      	xmgr = new Xmgr(xfile, title, "Path length [m]","Beta [m]");
// 		xmgr->subTitle(cwd);
// // 		xmgr->plot(p->zloc+2, p->r_bet+2, p->rotat+2, nplots-2);
// 		xmgr->plot(p->zloc+2, p->r_bet+2, p->rotat+2, 2);
// 		for(int s =0; s < p->slices; s++)
// 		{
// 			xmgr->plot(p->zloc+2, ((double *) p->s_bet)+s+20, p->rotat+2, nplots-2, 10);
// 		}
// 		delete xmgr;
// 
// 
// 		strcpy(title, "Slice Alphas/Betas");
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "alfa.agr");
//     		    xmgr = new Xmgr(xfile, title, "Path length [m]","Alpha");
// // 		xmgr->plot(p->zloc+2, p->r_alf+2, p->rotat+2, nplots-2);
// 		xmgr->plot(p->zloc+2, p->r_alf+2, p->rotat+2, 2);
// 		for(int s =0; s < p->slices; s++)
// 		{
// 			xmgr->plot(p->zloc+2, ((double *) p->s_alf)+s+20, p->rotat+2, nplots-2, 10);
// 		}
// 		delete xmgr;
// 
// 
// 		strcpy(title, "Slice Emittances");
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "emit.agr");
//     		    xmgr = new Xmgr(xfile, title, "Path length [m]","Emittance [m]");
// 		xmgr->frame(0., 250., 0., 2e-6);
// // 		xmgr->plot(p->zloc+2, p->r_eps+2, p->rotat+2, nplots-2);
// 		xmgr->plot(p->zloc+2, p->r_eps+2, p->rotat+2, 2);
// 		for(int s =0; s < p->slices; s++)
// 		{
// 			xmgr->plot(p->zloc+2, ((double *) p->s_eps)+s+20, p->rotat+2, nplots-2, 10);
// 		}
// 		delete xmgr;
// 
// 		strcpy(title, "Slice Vortex");
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "vort.agr");
// 	        xmgr = new Xmgr(xfile, title, "Path length [m]","<xyp-yxp> [m]");
// //		xmgr->frame(0., 250., 0., 2e-6);
// 		xmgr->plot(p->zloc+2, p->r_rot+2, p->rotat+2, nplots-2);
// 		for(int s =0; s < p->slices; s++)
// 		{
// 			xmgr->plot(p->zloc+2, ((double *) p->s_rot)+s+20, p->rotat+2, nplots-2, 10);
// 		}
// 		delete xmgr;
// 
// 
// 		strcpy(title, "Magnetic Field ");
// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "field.agr");
// 	        xmgr = new Xmgr(xfile, title, "Path length [m]","B [Gauss]");
// 		xmgr->plot(p->zloc+2, p->Bfield+2, p->rotat+2, nplots-2);
// 		delete xmgr;





		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "disp.agr");
        	xmgr = new Xmgr(xfile, "Dispersion", "Path length [m]","Dispersion[m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->dispx+3, p->rotat+3, nplots-2);
		xmgr->plot(p->zloc+3, p->dispy+3, p->rotat+3, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "dispp.agr");
        	xmgr = new Xmgr(xfile, "Dispersion", "Path length [m]","Dispersion slope");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->dispxp+3, p->rotat+3, nplots-2);
		xmgr->plot(p->zloc+3, p->dispyp+3, p->rotat+3, nplots-2);
		delete xmgr;



		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "magnet.agr");
        	xmgr = new Xmgr(xfile, "Magnetization", "Path length [m]","Normalized Magnetization [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->magnetization+1, p->rotat+1, nplots-2);
		delete xmgr;


		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "emit4d.agr");
        	xmgr = new Xmgr(xfile, "4d emittance sqrt", "Path length [m]","Normalized Emittance [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->emit4d+1, p->rotat+1, nplots-2);
		delete xmgr;




// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "emitu.agr");
//         	xmgr = new Xmgr(xfile, "Emittances", "Path length [m]","Normalized Emittance [m]");
// 		xmgr->subTitle(cwd);
// 		xmgr->plot(p->zloc+1, p->epsx_u+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->epsy_u+1, p->rotat+1, nplots-2);
// 		delete xmgr;

//
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "emitn.agr");
        	xmgr = new Xmgr(xfile, "Emittances", "Path length [m]","Normalized Emittance [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->epsx_n+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->epsy_n+1, p->rotat+1, nplots-2);
		delete xmgr;

//
//
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "sigt.agr");
        	xmgr = new Xmgr(xfile, "Bunch length", "Path length [m]","sigma_t[s]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->sigma_t+1, p->rotat+1, nplots-2);
		delete xmgr;
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "msigt.agr");
        	xmgr = new Xmgr(xfile, "Bunch length", "Path length [m]","sigma_t[s]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->tmax+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->tmin+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->dl_max+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->dl_min+1, p->rotat+1, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "sigl.agr");
        	xmgr = new Xmgr(xfile, "Bunch length", "Path length [m]","sigma_l[m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->sigma_l+1, p->rotat+1, nplots-2);
		delete xmgr;
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "msigl.agr");
        	xmgr = new Xmgr(xfile, "Bunch length", "Path length [m]","sigma_l[m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->zmax+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->zmin+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->dl_max+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->dl_min+1, p->rotat+1, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "emits.agr");
        	xmgr = new Xmgr(xfile, "Longitudinal Emittance", "Path length [m]","epss[m*MeV]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->epss+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		delete xmgr;

//
	
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "orbit.agr");
        	xmgr = new Xmgr(xfile, "Orbit", "Path length [m]","Offset [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->orbx+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->orby+1, p->rotat+1, nplots-2);
		delete xmgr;

//

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "temperature.agr");
        	xmgr = new Xmgr(xfile, "Temperature", "Path length [m]","");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->temperature+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		xmgr->plot(p->zloc+3, p->temperxp+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		xmgr->plot(p->zloc+3, p->temperyp+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		xmgr->plot(p->zloc+3, p->sigma_p +3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		delete xmgr;

//

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "cool.agr");
        	xmgr = new Xmgr(xfile, "RMS Energy deviation", "Path length [m]","sigma_E/E");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->coolsigma+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		delete xmgr;

//

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "sigp.agr");
        	xmgr = new Xmgr(xfile, "RMS Energy spread", "Path length [m]","sigma_E/E");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->sigma_p+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		delete xmgr;

//
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "msigp.agr");
        	xmgr = new Xmgr(xfile, "Max Energy spread", "Path length [m]","sigma_p/p");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+3, p->dp_max+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		xmgr->plot(p->zloc+3, p->dp_min+3, p->rotat+3, nplots-2);	// we skip the first 2 points because they screw up the scale
		delete xmgr;

//
		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "div.agr");
        	xmgr = new Xmgr(xfile, "Divergence", "Path length [m]","RMS divergence [rad]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->divx+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->divy+1, p->rotat+1, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "env.agr");
        	xmgr = new Xmgr(xfile, "RMS Envelopes", "Path length [m]","RMS size [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->envx+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->envy+1, p->rotat+1, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "menv.agr");
        	xmgr = new Xmgr(xfile, "Max Envelopes", "Path length [m]","RMS size [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->xmax+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->ymax+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->xmin+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->ymin+1, p->rotat+1, nplots-2);
		delete xmgr;

// 		for(int i=1; i<nplots-2; i++)
// 		{
// 				double rms = (p->envx[i] > p->envy[i]) ? p->envx[i] : p->envy[i];
// 				double mx  = (p->xmax[i] > p->ymax[i]) ? p->xmax[i] : p->ymax[i];
// 				mx  = (p->xmin[i] > mx) ? p->xmin[i] : mx;
// 				mx  = (p->ymin[i] > mx) ? p->ymin[i] : mx;
// 				printf("\t %12.5f \t %12.5f    %d %s %12.5f\n", rms, mx, i, p->type[i], p->zloc[i]);
// 		}


		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "ener.agr");
        	xmgr = new Xmgr(xfile, "Avarage energy", "Path length [m]","Energy [MeV]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->ener+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->ener1+1, p->rotat+1, nplots-2);
		delete xmgr;

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "gammaE.agr");
        	xmgr = new Xmgr(xfile, "Avarage gammaE", "Path length [m]","gamma ");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->gammaE+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->gammaE+1, p->rotat+1, nplots-2);
		delete xmgr;

//

// 		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
// 		strcat(xfile, "emitc.agr");
//      	xmgr = new Xmgr(xfile, "Thermal emittances", "Path length [m]","Normalized Emittance [m]");
// 		xmgr->subTitle(cwd);
// 		xmgr->plot(p->zloc+1, p->epsx_c+1, p->rotat+1, nplots-2);
// 		xmgr->plot(p->zloc+1, p->epsy_c+1, p->rotat+1, nplots-2);
// 		delete xmgr;

//

		if(icut) xfile[0]=0; else strcpy(xfile, "c-");
		strcat(xfile, "twiss.agr");
        	xmgr = new Xmgr(xfile, "Twiss Parameters", "Path length [m]","Betas, Disp, Dispp [m]");
		xmgr->subTitle(cwd);
		xmgr->plot(p->zloc+1, p->betx+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->bety+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->alfx+1, p->rotat+1, nplots-2);
		xmgr->plot(p->zloc+1, p->alfy+1, p->rotat+1, nplots-2);
		delete xmgr;


		getcwd(cwd, 200);
		cutWhat =0;
	}
}

