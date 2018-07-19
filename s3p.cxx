main(int argc, char ** argv)
{

		FILE * fpl;
		FILE * fpe;
		FILE * fpm;
		FILE * fpb;
		FILE * fpr;
		FILE * fps;
		FILE * fpg;
		FILE * fpc;
		FILE * fpd;
		FILE * fpp;
		FILE * fpde;
		FILE * fppp;
		FILE * fpp1;
		FILE * fpp2;
		FILE * fpp4;
// 		FILE * fppa;
		FILE * fpti;


		char file[100];
		sprintf(file,"slice.list.%d", neo);
		fpl = fopen(file,"w");
		if(!fpl)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}


	if(print)
	{
		sprintf(file,"slice.tips.%d", neo);
		fpti = fopen(file,"w");
		if(!fpti)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.magnet.%d", neo);
		fpm = fopen(file,"w");
		if(!fpm)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.emit4d.%d", neo);
		fpe = fopen(file,"w");
		if(!fpe)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.sigmar.%d", neo);
		fpr = fopen(file,"w");
		if(!fpr)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.sigmrp.%d", neo);
		fps = fopen(file,"w");
		if(!fps)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.gamma0.%d", neo);
		fpg = fopen(file,"w");
		if(!fpg)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.curren.%d", neo);
		fpc = fopen(file,"w");
		if(!fpc)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.disper.%d", neo);
		fpd = fopen(file,"w");
		if(!fpd)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.dispep.%d", neo);
		fpb = fopen(file,"w");
		if(!fpb)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.phase.%d", neo);
		fpp = fopen(file,"w");
		if(!fpp)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.density.%d", neo);
		fpde = fopen(file,"w");
		if(!fpde)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.phasep.%d", neo);
		fppp = fopen(file,"w");
		if(!fppp)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.p100.%d", neo);
		fpp1 = fopen(file,"w");
		if(!fpp1)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		fprintf(fpp1, "@    s0 symbol 9\n@    s0 symbol size 0.410000\n@    s0 line type 0\n");
		fprintf(fpp1, "@    s1 symbol 9\n@    s1 symbol size 0.410000\n@    s1 line type 0\n");
		fprintf(fpp1, "@    s2 symbol 9\n@    s2 symbol size 0.410000\n@    s2 line type 0\n");
		fprintf(fpp1, "@    s3 symbol 9\n@    s3 symbol size 0.410000\n@    s3 line type 0\n");
		fprintf(fpp1, "@    s4 symbol 9\n@    s4 symbol size 0.410000\n@    s4 line type 0\n");
		sprintf(file,"slice.p250.%d", neo);
		fpp2 = fopen(file,"w");
		if(!fpp2)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		fprintf(fpp2, "@    s0 symbol 9\n@    s0 symbol size 0.410000\n@    s0 line type 0\n");
		sprintf(file,"slice.p400.%d", neo);
		fpp4 = fopen(file,"w");
		if(!fpp4)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		fprintf(fpp4, "@    s0 symbol 9\n@    s0 symbol size 0.410000\n@    s0 line type 0\n");
// 		sprintf(file,"slice.pall.%d", neo);
// 		fppa = fopen(file,"w");
// 		if(!fppa)
// 		{
// 			fprintf(stderr,"cant open file %s\n", file);
// 			return -1;
// 		}
// 		fprintf(fppa, "@    s0 symbol 9\n@    s0 symbol size 0.410000\n@    s0 line type 0\n");
	}
	




	memcpy(ford, dcord, nbuf*7*sizeof(double));
	if(slices <= 0) slices=1;
	double lmin=1e30, lmax=-1e30;
	for(int i=0; i<nbuf; i++)
	{
// 		int j=int(ford[i][6])-1; // index
// 		if(j < 0 || j >= npoints)
// 		{
// 			fprintf(stderr, "wrong index in sliveemittance %d %d nbuf %d\n", i, j, nbuf);
// 			exit(-1);
// 		}
// 		ford[i][4] =  - double(icord[j][4]);

		if(ford[i][4] > lmax) lmax = ford[i][4];
		if(ford[i][4] < lmin) lmin = ford[i][4];
		
	}
	printf("lmin lmax %e %e nbuf %d ngood %d npoints %d\n", lmin, lmax, nbuf, ngood, npoints);

	double delta = (lmax - lmin)/slices;
	double * intrval = new double[slices+1];
	int * numpart = new int[slices];
	for(int s=0; s<slices; s++)
	{
		intrval[s]=lmin+s*delta;
		numpart[s] = 0;
	}
	intrval[slices]=lmax;


	for(int i=0; i< nbuf; i++)
	{
		int s = int((ford[i][4]-lmin)/delta);
		if(s < 0 || s >= slices) printf(" fort lmin lmax delta %e %e %e %e\n", ford[i][4], lmin , lmax, delta);
		numpart[s]++;
	}

	int maxpart=0; 
	int maxintv=0;
	for(int s=0; s<slices; s++)
		if(numpart[s] > maxpart)
		{
			maxpart = numpart[s];
			maxintv = s;
		}
	printf("maxpart %d maxintv %d\n", maxpart, maxintv);

	for(int s=0; s<maxintv; s++)
		if( numpart[s] < 2) lmin = intrval[s+1];

	for(int s=slices-1; s>maxintv; s--)
		if( numpart[s] < 2) lmax = intrval[s];

	delta = (lmax - lmin)/slices;
	for(int s=0; s<slices; s++) intrval[s]=lmin+s*delta;
	intrval[slices]=lmax;


	for(int s=0; s<slices; s++) s_rot[neo][s] = 0.;

	double dpm=0.;
        for(int i = 0; i < nbuf; i++) dpm += ford[i][5];
        dpm /= nbuf;


	emit4dsum[neo]=0.;
	double arr =0.;
	double arm =0.;
	double arrp =0.;
	double arprp =0.;
	double arot =0.;
	int nii=0;

	printf("lmin lmax %e %e nbuf %d ngood %d npoints %d\n", lmin, lmax, nbuf, ngood, npoints);
// 	if(print)
		fprintf(fpl, "nslices %d beami %e  freq  %e wavel %e lmin %e lmax %e zloc %e\n", slices, beami, freq, wavel, lmin, lmax, zloc[neo]);
	for(int s=0; s<slices; s++) 
	{
		int ni=0;
		double rrot  = 0.;
		double thetaprime = 0.;
        	double M     = 0.;
        	double rm    = 0.;
        	double rr    = 0.;
        	double rrp   = 0.;
        	double rpm   = 0.;
        	double rprp  = 0.;
		double ene   = 0.;
        	double xm    = 0.;
        	double xpm   = 0.;
        	double xx    = 0.;
        	double xxp   = 0.;
        	double xpxp  = 0.;
        	double ym    = 0.;
        	double ypm   = 0.;
        	double yy    = 0.;
        	double yyp   = 0.;
        	double ypyp  = 0.;
		double xy    = 0.;
		double xyp   = 0.;
		double xpy   = 0.;
		double xpyp  = 0.;
		double ddxp  = 0.;
		double ddxpp = 0.;
		double phase = 0.;
		double phasep= 0.;
		double energ = 0.;

		double s_time=(intrval[s] + intrval[s+1]-2.*lmin)/(2.*360.*freq*1e6);

		for(int i=0; i< nbuf; i++)
		{

			if(ford[i][4] >= intrval[s] && ford[i][4] <= intrval[s+1])
			{
         			double xt =ford[i][0];
	    			double xpt=ford[i][1];
	    			double yt =ford[i][2];
	    			double ypt=ford[i][3];
				double gamma= ford[i][5]/ erest[0]+1.;
				energ += ford[i][5];
				double betagamma =  sqrt(gamma*gamma-1);
				double p =  erest[0]*betagamma;
				double r2 = xt*xt+yt*yt;
				if(r2 > 1.e-20)
				{
					double r = sqrt(r2);
					rrot +=   r2*Bfield[neo]*clight*1.e-6;
					M += (xt*ypt-yt*xpt)*betagamma;
					thetaprime += (xt*ypt-yt*xpt)*betagamma /r2;

					double rpr= (xt*xpt+yt*ypt);
					double rp = rpr/r;
					double phas = atan(rpr/r2);
					phase += phas;
					rm += r;
					rr += r2;
					rrp += r*rp;
					rpm += rp;
					rprp += rp*rp;
					ene += ford[i][5];
            				xm   += xt;
            				xpm  += xpt;
            				xx   += xt*xt;
            				xxp  += xt*xpt;
            				xpxp += xpt*xpt;
            				ym   += yt;
            				ypm  += ypt;
            				yy   += yt*yt;
            				yyp  += yt*ypt;
            				ypyp += ypt*ypt;
	    				xy   += xt*yt;
	    				xyp  += xt*ypt;
	    				xpy  += xpt*yt;
	    				xpyp += xpt*ypt;
	    				double c5 = (ford[i][5]-dpm);
            				ddxp  += ford[i][0]*c5;
            				ddxpp += ford[i][1]*c5;
					phasep += phas*c5;
					ni++;
					if(print)
					{
						if(s == 100) fprintf(fpp1,"%e %e\n", r ,rp*betagamma);
						if(s == 250) fprintf(fpp2,"%e %e\n", r ,rp*betagamma);
						if(s == 400) fprintf(fpp4,"%e %e\n", r ,rp*betagamma);
// 						fprintf(fppa,"%e %e\n", r ,rp);
					}
				}
			}
		}
// 		printf("interval %d ni %d\n", s, ni);
		if(ni < 2)
		{
			fprintf(fpl, "%e, %e, %e, %e, %e, %e\n", 0., 0., 1., 0., 0., 0.);
			continue;
		}
		rrot /= ni;
		thetaprime /= ni;
		M /=ni;

		nii += ni;
        	rr /= ni;
        	rm /= ni;
        	rrp /= ni;
        	rprp /= ni;
        	ene /= ni;

        	xx /= ni;
        	xm /= ni;
        	xpm /= ni;
        	xxp /= ni;
        	xpxp /= ni;

        	yy /= ni;
        	ym /= ni;
        	ypm /= ni;
        	yyp /= ni;
        	ypyp /= ni;

		xy /= ni;
		xyp /= ni;
		xpy /= ni;
		xpyp /= ni;

		phase /= ni;
		phasep /= ni;
		energ /= ni;


		double gamma= ene/ erest[0]+1.;
		double betagamma = sqrt(gamma*gamma-1.);
		double emit4d2 = 
	    		xx*xpxp*yy*ypyp
	   		-xx*xpxp*yyp*yyp
	   		-xx*xpy*xpy*ypyp
	   		+2.*xx*xpy*xpyp*yyp
	   		-xx*xpyp*xpyp*yy
	   		-xxp*xxp*yy*ypyp
	   		+xxp*xxp*yyp*yyp
	   		+2.*xxp*xpy*xy*ypyp
	   		-2.*xxp*xpy*xyp*yyp
	   		-2.*xxp*xpyp*xy*yyp
	   		+2.*xxp*xpyp*xyp*yy
	   		-xpxp*xy*xy*ypyp
	   		+2.*xy*xpxp*xyp*yyp
	   		-2.*xy*xpyp*xyp*xpy
	   		-xpxp*xyp*xyp*yy
	   		+xyp*xyp*xpy*xpy
	   		+xy*xy*xpyp*xpyp
	   	;
		if(emit4d2 <= 0.)
		{
			fprintf(fpl, "%e, %e, %e, %e, %e, %e\n", 0., 0., 1., 0., 0., 0.);
			continue;
		}
		double emit4d = betagamma*sqrt(sqrt(emit4d2));

		double sigr=sqrt(rr);
		double sigrp= sqrt(rprp);
		double curr= double(slices*ni)/double(nbuf);
		double density=double(slices*ni)/rr;
        	s_rot[neo][0] = rrot;
        	s_rot[neo][1] = M;
        	s_rot[neo][2] = rr;
        	s_rot[neo][4] = Bfield[neo]/10000.;
		double re2 = rr*rprp -rrp*rrp;
        	s_eps[neo][s] = sqrt(re2);
		if(s_eps[neo][s] < 1.e-20) s_eps[neo][s]=1.e-20;
        	s_bet[neo][s] =   rr/s_eps[neo][s];
        	s_alf[neo][s] = -rrp/s_eps[neo][s];
        	s_alf[neo][s] /=  s_bet[neo][s];
        	s_disp[neo][s] = dpdp == 0 ? 0 : ddxp/dpdp;
        	s_dispp[neo][s] = dpdp == 0 ? 0 : ddxpp/dpdp;
        	phasep  = dpdp == 0 ? 0 : phasep/dpdp;
		emit4dsum[neo] += emit4d*emit4d*curr;
		arr +=rr*curr;
		arm +=rm*curr;
		arrp +=rrp*curr;
		arprp +=rprp*curr;
		arot += rrot+M;
 
		double sigrdc = ni ? sigr/ni : 0.;
		double sigrpdc = ni ? sigrp/ni : 0.;

		// wheta is the right definition?
// 		M = thetaprime*rr;

// 			fprintf(fpl, "%e, %e, %e, %e, %e, %e\n", sigr, sigrp, gamma, curr, -thetaprime*rr, emit4d);
			fprintf(fpl, "%e, %e, %e, %e, %e, %e\n", sigr, sigrp, gamma, curr, -M, emit4d);
		if(print)
		{

// 			printf("el=%d slice=%d e=%e b=%e a=%e  d=%e  dp=%e\n", neo, s, s_eps[neo][s], s_bet[neo][s], s_alf[neo][s],s_disp[neo][s],s_dispp[neo][s]);
// 			printf("ni=%d arr=%e arprp=%e \n", ni, arr, arprp);

// 			fprintf(fpe,"%d %e\n", s , emit4d);
// // 			fprintf(fpm,"%d %e\n", s ,  -thetaprime);
// 			fprintf(fpm,"%d %e\n", s ,  M);
// 			fprintf(fpr,"%d %e\n", s , sigr);
// 			fprintf(fps,"%d %e\n", s , sigrp);
// 			fprintf(fpg,"%d %e\n", s , gamma);
// 			fprintf(fpc,"%d %e\n", s , curr);
// 			fprintf(fpd,"%d %e\n", s , s_disp[neo][s]);
// 			fprintf(fpb,"%d %e\n", s , s_dispp[neo][s]);
// 			fprintf(fpp,"%d %e\n", s , phase);
// 			fprintf(fppp,"%d %e\n", s , phasep);
// 			fprintf(fpde,"%d %e\n", s , density);

			fprintf(fpe,"%e %e\n", s_time , emit4d);
// 			fprintf(fpm,"%e %e\n", s_time ,  -thetaprime);
			fprintf(fpm,"%e %e\n", s_time ,  M);
			fprintf(fpr,"%e %e\n", s_time , sigr);
			fprintf(fps,"%e %e\n", s_time , sigrp);
			fprintf(fpg,"%e %e\n", s_time , gamma);
			fprintf(fpc,"%e %e\n", s_time , curr);
			fprintf(fpd,"%e %e\n", s_time , s_disp[neo][s]);
			fprintf(fpb,"%e %e\n", s_time , s_dispp[neo][s]);
			fprintf(fpp,"%e %e\n", s_time , phase);
			fprintf(fppp,"%e %e\n", s_time , phasep);
			fprintf(fpde,"%e %e\n", s_time , density);

			fprintf(fpti,"%e %e\n", sigrdc, sigrpdc*betagamma);
		}
	}


	arr = sqrt(arr/slices);
	arprp = sqrt(arprp/slices);
	arot /= nbuf;
	printf("nbuf=%d arr=%e arprp=%e \n", nbuf, arr, arprp);

	emit4dsum[neo] /= slices;
	emit4dsum[neo] = sqrt(emit4dsum[neo]);
	double re2 = arr*arprp -arrp*arrp;
	r_eps[neo] = sqrt(re2);
	if(r_eps[neo] < 1.e-20) r_eps[neo]=1.e-20;
        r_bet[neo] =   arr/r_eps[neo];
        r_alf[neo] =  -arrp/r_eps[neo];
        r_rot[neo] =  arot;
	printf("r_bet[neo],r_alf[neo],r_eps[neo], emit4dsum = %f %f %e %e\n", r_bet[neo],r_alf[neo],r_eps[neo],emit4dsum[neo]);

	
		fclose(fpl);
	if(print)
	{
		fclose(fpe);
		fclose(fpm);
		fclose(fpr);
		fclose(fps);
		fclose(fpg);
		fclose(fpc);
		fclose(fpd);
		fclose(fpb);
		fclose(fpp);
		fclose(fppp);
		fclose(fpp1);
		fclose(fpp2);
		fclose(fpp4);
// 		fclose(fppa);
	}

	delete [] intrval;
	return 0;

}
