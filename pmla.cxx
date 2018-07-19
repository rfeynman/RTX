//	pmla class : object to read the parmew binary files
//	translation from the example.for file provided by
//	the parmela distribution
//	jorg kewisch 11.2002
//
#include "pmla.hxx"
double dminv(double * a, int n);


void  Pmla::readbuff(int fp, char * buff, int nbytes)
{
		int toread = nbytes;
// 		printf("toread %d \n", toread);
		char * bb = buff;
		while(toread > 0)
		{
// printf("toread %d \n", toread);
			int rdnow = toread > 2000 ? 2000 : toread;
			int rbytes = read(fp, bb, rdnow);
			if(rbytes < 0)
			{
				printf("read incomplete, errno = %d\n",
					errno);
				exit(-1);
			}
			toread -= rbytes;
			bb += rbytes;
// printf("toread %d nbytes %d errno %d\n", toread, rbytes, errno);
		}
}
	

// void  Pmla::freadbuff(int fp, double * buff, int nfloats)
// {
// }
	




void  Pmla::writebuff(int fp, char * buff, int nbytes)
{
		int toread = nbytes;
// 		printf("towrite %d \n", toread);
		char * bb = buff;
		while(toread > 0)
		{
// printf("toread %d \n", toread);
			int rdnow = toread > 2000 ? 2000 : toread;
			int rbytes = write(fp, bb, rdnow);
			if(rbytes < 0)
			{
				printf("write incomplete, errno = %d\n",
					errno);
				exit(-1);
			}
			toread -= rbytes;
			bb += rbytes;
// printf("toread %d nbytes %d errno %d\n", toread, rbytes, errno);
		}
}





void  Pmla::freadbuff(int fp, float * buff, int nfloats)
{
		int toread = 4*nfloats;
// 		printf("toread %d \n", toread);
		char * bb = (char *) buff;
		while(toread > 0)
		{
// printf("toread %d \n", toread);
			int rdnow = toread > 2000 ? 2000 : toread;
			int rbytes = read(fp, bb, rdnow);
			if(rbytes < 0)
			{
				printf("read incomplete, errno = %d\n",
					errno);
				exit(-1);
			}
			toread -= rbytes;
			bb += rbytes;
// printf("toread %d nbytes %d errno %d\n", toread, rbytes, errno);
		}
}



int  Pmla::freadbuffx(int fp, float * buff, int nfloats)
{
		int toread = 4*nfloats;
// 		printf("toread %d \n", toread);
		char * bb = (char *) buff;
		while(toread > 0)
		{
// printf("toread %d \n", toread);
			int rdnow = toread > 2000 ? 2000 : toread;
			int rbytes = read(fp, bb, rdnow);
			if(rbytes <= 0) return errno+100000;
			toread -= rbytes;
			bb += rbytes;
		}
		return 0;
}





void  Pmla::fwritebuff(int fp, float * buff, int nfloats)
{
		int toread = 4*nfloats;
// 		printf("towrite %d \n", toread);
		char * bb = (char *) buff;
		while(toread > 0)
		{
//		 printf("toread %d \n", toread);
			int rdnow = toread > 2000 ? 2000 : toread;
			int rbytes = write(fp, bb, rdnow);
			if(rbytes < 0)
			{
				printf("write incomplete, errno = %d\n",
					errno);
				exit(-1);
			}
			toread -= rbytes;
			bb += rbytes;
// printf("toread %d nbytes %d errno %d\n", toread, rbytes, errno);
		}
}

void  Pmla::readTimeFile()
{
//	the tape3 file is a binary sequential unformated fortran file. The Lahey compiler writes each record as 
//	the number of bytes of data, the data itself and the number of bytes again. Therfore the first and last four
//	bytes are ignored.
//	
//	  The IP array contains particle numbers (1, 2, or 3) defined on the
//	  CHARGE line. The format of the CHARGE line is:
//	  CHARGE Charge, ParticleType
//	  where Charge is an interger from -127 to +127 and ParticeleType is 1, 2, or 3
//	  The charge spcifies the particle type for the following INPUT line in Parmela
//	  input file. a different mass maybe associated with each ParticleType
//	
//	
//	  If NUM is equal to 5, then the file is of type TAPE3 and
//	  coordinates for time step number NE will appear in array CORD.
//	  The coordinates at each time step are stored in array CORD as
//	  follows:
//	
//	  CORD(1,I)=x
//	  CORD(2,I)=y
//	  CORD(3,I)=z
//	  CORD(4,I)=original index of the particle in the COR array
//	  CORD(5,I)=Kinetic Energy (keV)
//	
//	  Note that CORD(4,I) is always an integer number even though it is
//	  stored in a real array as a real number.
//	
//	 for out5 time steping
//	 ipt=-1 means process all particles
//	 ipt=0  means other particles i.e. cord(5*n-1)=0
//	 ipt= 1 means process only beam particles i.e. cord(5*n-1).ne.0
//	 cord (1-3) is x,y,z : cord(4) is orginal particle number :cord (5) is KE.
//	 use orginial particle number to identify particles type and charge.


	if( !fp3)
	{
		fprintf(stderr,"tape3 file is not open, change pmla.cxx\n");
		exit(-1);
	}
	FILE * fp1 = fopen("part1.dat", "w");
	if(!fp1)
	{
		fprintf(stderr," cant open %s\n", "part1.dat");
		return;
	}

	FILE * fpx = fopen("twx.dat", "w");
	if(!fpx)
	{
		fprintf(stderr," cant open %s\n", "twx.dat");
		return;
	}

	FILE * fpy = fopen("twy.dat", "w");
	if(!fpy)
	{
		fprintf(stderr," cant open %s\n", "twy.dat");
		return;
	}

	FILE * fpz = fopen("twz.dat", "w");
	if(!fpz)
	{
		fprintf(stderr," cant open %s\n", "twz.dat");
		return;
	}

	fprintf(fpx, "# ax,bx,xbx,ybx,exrms,xmaxrms\n");
	fprintf(fpy, "# ay,by,xby,yby,eyrms,ymaxrms\n");
	fprintf(fpz, "# az,bz,xbz,ybz,ezrms,zmaxrms\n");


	float line[20];
	while(!freadbuffx(fp3, line, 4))
	{


		int kk=0;
		int rechead= mapint( line[kk++]);	// record length
// 		printf("rechead = %d\n",rechead);
		int nstep= mapint( line[kk++]);		// steps
// 		printf("nstep = %d\n",nstep);
		int nbuf= mapint( line[kk++]);		// number of particles
		printf("nbuf = %d\n",nbuf);
		float wt=line[kk++];			// ?
		printf("wt = %f",wt);

		float (*cord)[5]  = new float[nbuf][5];
		freadbuff(fp3,(float *)  cord, 5*nbuf);

		char *ic= new char[nbuf];
		readbuff(fp3, ic, nbuf);
		char *ip= new char[nbuf];
		readbuff(fp3, ip, nbuf);
		fprintf(fp1, "%d %d %f %e %e %e %e %e %d %d \n", nsteps, nbuf, wt, cord[1][0], cord[1][1], cord[1][2], cord[1][3], cord[1][4], ic[1], ip[1] );

		freadbuff(fp3, line, 19);
		rechead= mapint( line[18]);
// 		printf("rechead = %d\n",rechead);
		rechead= nbuf*22+12+18*4;
// 		printf("bytes = %d\n",rechead);
	
		ax	= line[ 0];
		bx	= line[ 1];
		xbx	= line[ 2];
		ybx	= line[ 3];
		exrms	= line[ 4];
		xmaxrms	= line[ 5];

		ay	= line[ 6];
		by	= line[ 7];
		xby	= line[ 8];
		yby	= line[ 9];
		eyrms	= line[10];
		ymaxrms	= line[11];

		az	= line[12];
		bz	= line[13];
		xbz	= line[14];
		ybz	= line[15];
		ezrms	= line[16];
		zmaxrms	= line[17];

		fprintf(fpx, "%e %e %e %e %e %e\n", ax,bx,xbx,ybx,exrms,xmaxrms);
		fprintf(fpy, "%e %e %e %e %e %e\n", ay,by,xby,yby,eyrms,ymaxrms);
		fprintf(fpz, "%e %e %e %e %e %e\n", az,bz,xbz,ybz,ezrms,zmaxrms);
		printf("wt = %f\n",wt);
	}
}
	


Pmla::Pmla(const char * elemfilename, const char * timefilename)
{
	this->elemfilename = strdup(elemfilename);
	fp2 = open(elemfilename, O_RDONLY | O_BINARY);
	if( fp2 < 0)
	{
		fprintf(stderr, "Cant open file %s\n", elemfilename);
		exit(-1);
	}
	this->timefilename = strdup(timefilename);
// 	fp3 = open(timefilename, O_RDONLY | O_BINARY);
// 	if( fp3 < 0)
// 	{
// 		fprintf(stderr, "Cant open file %s\n", timefilename);
// 		exit(-1);
// 	}

	fp3 = 0;



	rdwr2=0;
	filehandle[0] = fp2;
	float line[57];
	freadbuff(fp2, line, 57);


	int kk=0;
	ireclen= mapint( line[kk++]);
// 	printf("ireclen = %d\n",ireclen);
	numrec = mapint( line[kk++]);
// 	printf("numrec  = %d\n",numrec );
	npoints= mapint( line[kk++]);
// 	printf("npoints = %d\n",npoints);
	irun   = mapint( line[kk++]);
// 	printf("irun    = %d\n",irun   );

	memcpy(title, line+kk, 80); title[80]=0; kk += 20;
// 	printf("\"%s\"\n", title);

	pi       =  line[kk  ];
// 	printf("pi      = %f\n",pi     );
	twopi    =  line[kk+1];
// 	printf("twopi   = %f\n",twopi  );
	radian   =  line[kk+2];
// 	printf("radian  = %f\n",radian );
	clight   =  line[kk+3];
// 	printf("clight  = %f\n",clight );
	erest[0]       =  line[kk+4]; 
	erest[1]       =  line[kk+5]; 
	erest[2]       =  line[kk+6]; 
// 	printf("erest   = %f  %f  %f\n", erest[0],erest[1],erest[2]     );
	brhof    =  line[kk+7];
// 	printf("brhof   = %f\n",brhof  );
	freq     =  line[kk+8];
// 	printf("freq    = %f\n",freq   );
	wavel    =  line[kk+9];
// 	printf("wavel   = %f\n",wavel  );
	z0       =  line[kk+10];
// 	printf("z0      = %f\n",z0     );
	w0       =  line[kk+11];
// 	printf("w0      = %f\n",w0     );
	beami    =  line[kk+12];
// 	printf("beami   = %f\n",beami  );
	dcon     =  line[kk+13];
// 	printf("dcon    = %f\n",dcon   );
	zlimit   =  line[kk+14]; kk += 16;
// 	printf("zlimit  = %f\n",zlimit );
	nel    = mapint( line[kk++]);
// 	printf("nel     = %d\n",nel    );



	dwt      = *( (double *) (line+kk));
// 	printf("dwt     = %f\n",dwt    );
	wt0      = *( (double *) (line+kk+ 2));
// 	printf("wt0     = %f\n",wt0    );
	dum      = *( (double *) (line+kk+ 4));
// 	printf("dum     = %f\n",dum    );
	swt      = *( (double *) (line+kk+ 6));
// 	printf("swt     = %f\n",swt    );
	cwt      = *( (double *) (line+kk+ 8));
// 	printf("cwt     = %f\n",cwt    );
	nsteps = mapint( line[kk+10]);
// 	printf("nsteps  = %d\n",nsteps );
	nsc    = mapint( line[kk+11]);
// 	printf("nsc     = %d\n",nsc    );
	nout   = mapint( line[kk+12]); kk += 13;
// 	printf("nout    = %d\n",nout   );


	ngood  = mapint( line[kk++]);
// 	printf("ngood   = %d\n",ngood  );

// 	printf("kk      = %d\n", kk );
	filesize[0] = lseek(fp2, 0, SEEK_END);
	nplots=((filesize[0])/ireclen)/numrec;
	printf("size    = %d nplots = %d\n", filesize[0], nplots);
	if(nplots == 0) return;

	lseek(fp2,0,SEEK_SET);

	icord          = new float*[npoints];
	cord           = new float*[npoints];
	dcord          = new double*[npoints];
	scord          = new double*[npoints];
	ford           = new double*[npoints];
	for(int i      =0; i< npoints; i++)
	{
		icord[i]       = new float[6];
		 cord[i]       = new float[7];
		dcord[i]       = new double[7];
		 ford[i]       = new double[7];
	}

	ic             = new char [npoints];
	iicol          = new char [npoints+1];
	ip             = new char [npoints+1];
	icol           = new int  [npoints+1];
	map            = new int  [nel+1];
	zloc           = new float[nel+1];
	el             = new float[nel+1];
	type           = new char*[nel+1];
	rotat          = new int  [nel+1];
        Bfield         = new double[nel+1];

        ener           = new double[nel+1];
        ener1          = new double[nel+1];
        gammaE         = new double[nel+1];
        phas           = new double[nel+1];
        blen           = new double[nel+1];
        disl           = new double[nel+1];
        dislp          = new double[nel+1];
        dispx          = new double[nel+1];
        dispxp         = new double[nel+1];
        dispy          = new double[nel+1];
        dispyp         = new double[nel+1];
        orbx           = new double[nel+1];
        orby           = new double[nel+1];
        orbxp          = new double[nel+1];
        orbyp          = new double[nel+1];
        epsx           = new double[nel+1];
        epsx_n         = new double[nel+1];
        betx           = new double[nel+1];
        alfx           = new double[nel+1];
        epsy           = new double[nel+1];
        epsy_n         = new double[nel+1];
        bety           = new double[nel+1];
        alfy           = new double[nel+1];
        epss           = new double[nel+1];
        magnetization  = new double[nel+1];
        emit4d         = new double[nel+1];
        emit4dsum      = new double[nel+1];
	emitz          = new double[nel+1];
	temperature    = new double[nel+1];
	temperxp       = new double[nel+1];
	temperyp       = new double[nel+1];
        ttrans         = new double[nel+1];
        ttranm         = new double[nel+1];
        eps_dx         = new double[nel+1];
        eps_lx         = new double[nel+1];
        epsx_c         = new double[nel+1];
        epsy_c         = new double[nel+1];
        epsx_u         = new double[nel+1];
        epsy_u         = new double[nel+1];
        rl2            = new double[nel+1];
        rot            = new double[nel+1];
        envx           = new double[nel+1];
        envy           = new double[nel+1];
        divx           = new double[nel+1];
        divy           = new double[nel+1];
        roundness      = new double[nel+1];
        roundnessp     = new double[nel+1];
        slopex         = new double[nel+1];
        slopey         = new double[nel+1];
	xmax           = new double[nel+1];
	ymax           = new double[nel+1];
	zmax           = new double[nel+1];
	tmax           = new double[nel+1];
	xmin           = new double[nel+1];
	ymin           = new double[nel+1];
	zmin           = new double[nel+1];
	tmin           = new double[nel+1];
        sigma_p        = new double[nel+1];
        coolsigma      = new double[nel+1];
        dp_max         = new double[nel+1];
        dp_min         = new double[nel+1];
        sigma_l        = new double[nel+1];
        sigma_t        = new double[nel+1];
        dl_max         = new double[nel+1];
        dl_min         = new double[nel+1];
	set            = new int   [nel+1];

	
	r_eps          = new double[nel+1];
	r_bet          = new double[nel+1];
	r_alf          = new double[nel+1];
	r_rot          = new double[nel+1];
	s_eps          = new double[nel+1][500];
	s_bet          = new double[nel+1][500];
	s_alf          = new double[nel+1][500];
	s_disp         = new double[nel+1][500];
	s_dispp        = new double[nel+1][500];
	s_rot          = new double[nel+1][500];

	buff = new char[ireclen+10];
	fbuff = (float*) buff;
	readbuff(fp2, buff, ireclen);

	for(int i=0; i <= nel; i++)
	{
		zloc[i]=fbuff[kk++]/100.;	// convert to metres
		el[i]=fbuff[kk++];
		map[i]= mapint(fbuff[kk++]);
		type[i]=NULL;
		rotat[i]=1;
	}

// 	for(int i=0; i<= nel; i++)
// 		printf("%6d %15f %15f %6d\n", i, zloc[i], el[i], map[i]);

	int left= numrec ==1 ? npoints : (ireclen/4-57-3*(nel+1))/7;

	for(int ii=0; ii < left; ii++)
	{
		memcpy(icord[ii], fbuff+kk, 6*sizeof(float));
		kk += 6;
	}
	int kk1 = 4*kk;
	memcpy( ip, buff+kk1, left); kk1 += left;
	memcpy( iicol, buff+kk1, left); kk1 += left;


	for(int irec=1; irec < numrec; irec++)
	{

// 		printf(" irec = %d\n", irec);
		lseek(fp2,ireclen*irec,SEEK_SET);
		readbuff(fp2, buff, ireclen);
		int ibigin = left;
		int more = (ireclen/4)/7;
		left += more;
		if(left > npoints)
		{
			left = npoints;
			more = left - ibigin;
		}
// 		printf("irec = %d ibigin = %d more=%d left=%d \n",irec, ibigin,more,left);
		kk1 = 0;
		for(int ii=0; ii < more; ii++)
		{
			memcpy(icord[ibigin+ii], fbuff+kk1, 6*sizeof(float));
			kk1 += 6;
		}
		kk1 = more*6*sizeof(float);
		memcpy((ip+ibigin), buff+kk1, more); kk1 += more;
		memcpy((iicol+ibigin), buff+kk1, more); kk1 += more;
	}


	maxj=1;
	for(int i=0; i<nel; i++) if(map[i] != 0) maxj++;
	rfsize = float(maxj)*float(numrec)*float(ireclen);
	numfile=int(rfsize/2.09e9+1);
// 	printf("maxj=%d, numfile=%d rfsize=%f numrec=%d ireclen=%d\n", maxj, numfile, rfsize,numrec,ireclen);
	maxj = int(rfsize/(float(numrec)*float(ireclen))/numfile)+1;
	if(numfile > 10) numfile = 10;
// 	printf("maxj=%d, numfile=%d rfsize=%f\n", maxj, numfile, rfsize);

	struct stat statbuf;

	for(int numf=1; numf < numfile; numf++)
	{
		char file[100];
		strcpy(file, elemfilename);
		int l = strlen(elemfilename);
		file[l]= numf+64; file[l+1]=0;
		filehandle[numf] = open(file, O_RDONLY);
		if( fp2 < 0)
		{
			printf(" cant open file %s\n", file);
			exit(-1);
		}
		fstat(filehandle[numf], &statbuf);
		filesize[numf]=statbuf.st_size;
		printf("open file %s, size %d\n", file, filesize[numf]);
	}
	gunEnd=-1;
	
}



Pmla::~Pmla()
{
	for(int i=0; i< npoints; i++)
	{
		delete [] icord[i];
		delete []  cord[i];
		delete [] dcord[i];
		delete []  ford[i];
	}
	delete [] icord   ;
	delete [] cord  ;
	delete [] dcord ;
	delete [] scord ;
	delete [] ford  ;

	delete [] ic    ;
	delete [] iicol ;
	delete [] ip    ;
	delete [] icol  ;
	delete [] map   ;
	delete [] zloc  ;
	delete [] el    ;
	for(int i=0; i <= nel; i++) free(type[i]);
	delete [] type  ;
	delete [] rotat ;
        delete [] Bfield;
        delete [] ener  ;
        delete [] ener1 ;
        delete [] gammaE;
        delete [] phas  ;
        delete [] blen  ;
        delete [] disl  ;
        delete [] dislp ;
        delete [] dispx ;
        delete [] dispxp;
        delete [] dispy ;
        delete [] dispyp;
        delete [] orbx  ;
        delete [] orby  ;
        delete [] orbxp ;
        delete [] orbyp ;
        delete [] epsx  ;
        delete [] epsx_n;
        delete [] betx  ;
        delete [] alfx  ;
        delete [] epsy  ;
        delete [] epsy_n;
        delete [] bety  ;
        delete [] alfy  ;
        delete [] epss  ;
        delete [] magnetization;
        delete [] emit4d;
        delete [] emit4dsum;
	delete [] emitz;
	delete [] temperature;
	delete [] temperxp;
	delete [] temperyp;
        delete [] ttrans;
        delete [] ttranm;
        delete [] eps_dx;
        delete [] eps_lx;
        delete [] epsx_c;
        delete [] epsy_c;
        delete [] epsx_u;
        delete [] epsy_u;
        delete [] rl2;
        delete [] rot;
        delete [] envx;
        delete [] envy;
        delete [] divx;
        delete [] divy;
        delete [] roundness;
        delete [] roundnessp;
        delete [] slopex;
        delete [] slopey;
	delete [] xmax;
	delete [] ymax;
	delete [] zmax;
	delete [] tmax;
	delete [] xmin;
	delete [] ymin;
	delete [] zmin;
	delete [] tmin;
        delete [] sigma_p;
        delete [] coolsigma;
        delete [] dp_max;
        delete [] dp_min;
        delete [] sigma_l;
        delete [] sigma_t;
        delete [] dl_max;
        delete [] dl_min;

	delete [] set   ;
	delete [] r_eps ;
	delete [] r_bet ;
	delete [] r_alf ;
	delete [] r_rot ;
	delete [] s_eps ;
	delete [] s_bet ;
	delete [] s_alf ;
	delete [] s_disp ;
	delete [] s_dispp ;
	delete [] s_rot ;
	delete [] buff ;
	close(fp2);
// 	close(fp3);
}




int Pmla::getElement(int ne)
{
	if(ne == 0)
	{
        	for(int i = 0; i < npoints; i++)
		{
        		dcord[i][0]=double(icord[i][0])/100.;
        		dcord[i][1]=double(icord[i][1]);
        		dcord[i][2]=double(icord[i][2])/100.;
        		dcord[i][3]=double(icord[i][3]);
        		dcord[i][4]=double(icord[i][4]);
        		dcord[i][5]=double(icord[i][5]);
			dcord[i][6]=double(i+1);
		}
		nbuf=npoints;
		neo=ne;

		return 0;
	}



	int nnx = -1;
	set[ne]=0;

	for(int i=0; i< nel; i++) if(map[i] == ne) {nnx=i; break;}
	if (nnx < 0)
	{
		printf("elemnet %d not found\n", ne);
		return(-1);
	}
	int numf=nnx/maxj;
	set[ne]=0;
	if(numf > 10) return(-1);
	
	int irecs=(nnx-numf*maxj)*numrec;
	int irece=irecs+numrec;
// 	printf("nnx=%d, irecs=%d irece=%d, nplots=%d, filezi=%d, ireclen=%d, numrec=%d\n",
// 			nnx, irecs, irece, nplots, filesize[0],ireclen,numrec);

	int errs = lseek(filehandle[numf],ireclen*irecs,SEEK_SET);
	if(errs < 0)
	{
			printf(" lseek error\n");
			return -1;
	}
	readbuff(filehandle[numf], buff, ireclen);
	nbuf= mapint(fbuff[0]);
	neo= mapint(fbuff[1]);
	pr = fbuff[2];
	wr = fbuff[3];
// 	printf("nbuf=%d neo=%d pr=%e wr=%e\n", nbuf, neo, pr, wr);
	if( neo != ne)
	{
		printf("****************************** neo=%d and ne=%d differ\n", neo, ne);
		return(-1);
	}
// 	printf("****************************** neo=%d and ne=%d same\n", neo, ne);

	int left= ( ireclen/4-4)/8; if(left > nbuf) left = nbuf; 
// 	printf("left=%d \n", left);

	for(int ii=0; ii < left; ii++)
		memcpy(cord[ii], fbuff+4+ii*7, 7*sizeof(float));
	int kk1 = (7*left+4)*sizeof(float);
	memcpy( ic, buff+kk1, left);



	for(int irec=irecs+1; irec < irece; irec++)
	{

// 		printf(" irec = %d\n", irec);
		lseek(fp2,ireclen*irec,SEEK_SET);
		int ibigin = left;
		int more = (ireclen/4)/8;
		left += more;
		if(left > nbuf)
		{
			left = nbuf;
			more = left - ibigin;
		}
// 		printf("ibigin = %d more=%d left=%d \n",ibigin,more,left);
		for(int ii=0; ii < more; ii++)
			freadbuff(fp2, cord[ibigin+ii], 7);
		readbuff(fp2, ic+ibigin, more);
// 		for(int i=ibigin; i<left; i++)
// 			printf("%6d %11f %11f %11f %11f %11f %11f %11f\n", i, 
// 				cord[i][0], cord[i][1], cord[i][2], 
// 				cord[i][3], cord[i][4], cord[i][5],
// 				cord[i][6]);
	}
	if(rotat[ne] < 0 )
	{
		for(int i=0; i<nbuf; i++) for(int j=0; j<4; j++) cord[i][j] = - cord[i][j];
	}
//	for(int i=0; i<nbuf; i++) icol[i] = iicol[int(cord[i][6])-1]-ic[i];

// 	printf("+++ %d +++ cord +++\n", ne);
// 
// 	for(int i=0; i<nbuf; i++)
// 	{
// 		double xp = cord[i][1]/sqrt( 1. + cord[i][1]*cord[i][1] + cord[i][3]*cord[i][3]);
// 		printf("%6d %11f %11f %11f %11f %11f %11f %11f : %e\n", i, 
// 			cord[i][0], cord[i][1], cord[i][2], 
// 			cord[i][3], cord[i][4], cord[i][5], cord[i][6], xp);
// 	}
		
	for(int i=0; i<nbuf; i++)
	{
		for(int j=0; j<7; j++) dcord[i][j] =  cord[i][j];
		dcord[i][0] *= 0.01; // change units from cm to m
		dcord[i][2] *= 0.01;
// 	    	double ppz=sqrt( 1. + dcord[i][1]*dcord[i][1] + dcord[i][3]*dcord[i][3]);
	    	double ppz=1.;
		dcord[i][1] /= ppz; // change  whatever to x' and y'
		dcord[i][3] /= ppz;
// 		if(dcord[i][6] < 7) printf("%11f   zero  %11f   %11f %11f %11f %11f %11f %11f\n", dcord[i][6], zloc[ne],
// 			dcord[i][0], dcord[i][1], dcord[i][2], 
// 			dcord[i][3], dcord[i][4], dcord[i][5]);
	}
	return 0;
}




int Pmla::putElement(int ne, int what)
{
	if(!rdwr2)
	{
		close(fp2);
		fp2 = open(elemfilename, O_RDWR   | O_BINARY);
		if( fp2 < 0)
		{
			fprintf(stderr, "Cant open file %s for write\n", elemfilename);
			perror("Why? ");
			exit(-1);
		}
		rdwr2=1;
	}

	int nnx = -1;
	set[ne]=0;

	for(int i=0; i< nel; i++) if(map[i] == ne) {nnx=i; break;}
	if (nnx < 0)
	{
		printf("elemnet %d not found\n", ne);
		return(-1);
	}
	int numf=nnx/maxj;
	set[ne]=0;
	if(numf > 10) return(-1);
	
	for(int i=0; i<nbuf; i++)
	{
		if(what) for(int j=0; j<7; j++) cord[i][j] =  ford[i][j];
		else     for(int j=0; j<7; j++) cord[i][j] =  dcord[i][j];
		cord[i][0] *= 100.; // change units from m to cm
		cord[i][2] *= 100.;
	    	double ppz=sqrt( 1. - cord[i][1]*cord[i][1] - cord[i][3]*cord[i][3]);
		cord[i][1] /= ppz; // change  whatever to x' and y'
		cord[i][3] /= ppz;
	}
	int irecs=(nnx-numf*maxj)*numrec;
	int irece=irecs+numrec;
// 	printf("nnx=%d, irecs=%d irece=%d, nplots=%d, filezi=%d, ireclen=%d, numrec=%d\n",
// 			nnx, irecs, irece, nplots, filesize[0],ireclen,numrec);

	int errs = lseek(filehandle[numf],ireclen*irecs,SEEK_SET);
	if(errs < 0)
	{
			printf(" lseek error\n");
			return -1;
	}
	readbuff(filehandle[numf], buff, ireclen);
	nbuf= mapint(fbuff[0]);
	neo= mapint(fbuff[1]);
	pr = fbuff[2];
	wr = fbuff[3];
// 	printf("nbuf=%d neo=%d pr=%e wr=%e\n", nbuf, neo, pr, wr);
	if( neo != ne)
	{
		printf("****************************** neo=%d and ne=%d differ\n", neo, ne);
		return(-1);
	}

	int left= ( ireclen/4-4)/8; if(left > nbuf) left = nbuf; 

	for(int ii=0; ii < left; ii++)
		memcpy(fbuff+4+ii*7, cord[ii], 7*sizeof(float));
	int kk1 = 4*(7*left+4);
	errs = lseek(filehandle[numf],ireclen*irecs,SEEK_SET);
	writebuff(filehandle[numf], buff, ireclen);



	for(int irec=irecs+1; irec < irece; irec++)
	{

// 		printf(" irec = %d\n", irec);
		lseek(fp2,ireclen*irec,SEEK_SET);
		int ibigin = left;
		int more = (ireclen/4)/8;
		left += more;
		if(left > nbuf)
		{
			left = nbuf;
			more = left - ibigin;
		}
// 		printf("ibigin = %d more=%d left=%d \n",ibigin,more,left);
		fwritebuff(fp2, cord[ibigin], 7*more);
		writebuff(fp2, ic+ibigin, more);
	}
		
	return 0;
}





double Pmla::getValue(const char * what)
{
	if(!strcmp(what,"ener")) return ener[neo];
	if(!strcmp(what,"gammaE")) return gammaE[neo];
	if(!strcmp(what,"phas")) return phas[neo];
	if(!strcmp(what,"blen")) return blen[neo];
	if(!strcmp(what,"disl")) return disl[neo];
	if(!strcmp(what,"dislp"))  return dislp[neo];
	if(!strcmp(what,"dispx"))   return dispx[neo];
	if(!strcmp(what,"dispxp"))  return dispxp[neo];
	if(!strcmp(what,"dispy"))   return dispy[neo];
	if(!strcmp(what,"dispyp"))  return dispyp[neo];
	if(!strcmp(what,"disp+"))   return s_disp[neo][slices-2];
	if(!strcmp(what,"dispp+"))  return s_dispp[neo][slices-2];
	if(!strcmp(what,"disp0"))   return s_disp[neo][(slices+1)/2];
	if(!strcmp(what,"dispp0"))  return s_dispp[neo][(slices+1)/2];
	if(!strcmp(what,"disp-"))   return s_disp[neo][1];
	if(!strcmp(what,"dispp-"))  return s_dispp[neo][1];
	if(!strcmp(what,"orbx"))   return orbx[neo];
	if(!strcmp(what,"orby"))   return orby[neo];
	if(!strcmp(what,"orbxp"))   return orbxp[neo];
	if(!strcmp(what,"orbyp"))   return orbyp[neo];
	if(!strcmp(what,"epsx"))   return epsx[neo];
	if(!strcmp(what,"epsx_n")) return epsx_n[neo];
	if(!strcmp(what,"betx")) return betx[neo];
	if(!strcmp(what,"alfx")) return alfx[neo];
	if(!strcmp(what,"epsy")) return epsy[neo];
	if(!strcmp(what,"epsy_n")) return epsy_n[neo];
	if(!strcmp(what,"bety")) return bety[neo];
	if(!strcmp(what,"alfy")) return alfy[neo];
	if(!strcmp(what,"epss"))   return epss[neo];
	if(!strcmp(what,"epssdl"))   return epss[neo]/sigma_l[neo]; 	// this is the enegy spread after removing the chirp
	if(!strcmp(what,"magnet")) return magnetization[neo];
	if(!strcmp(what,"emit4d")) return emit4d[neo];
	if(!strcmp(what,"emit4dsum")) return emit4dsum[neo];
	if(!strcmp(what,"emitz")) return emitz[neo];
	if(!strcmp(what,"temper")) return temperature[neo];
	if(!strcmp(what,"temperx")) return temperxp[neo];
	if(!strcmp(what,"tempery")) return temperyp[neo];
	if(!strcmp(what,"ttrans")) return ttrans[neo];
	if(!strcmp(what,"ttranm")) return ttranm[neo];
	if(!strcmp(what,"r_eps")) return r_eps[neo];
	if(!strcmp(what,"r_bet")) return r_bet[neo];
	if(!strcmp(what,"r_alf")) return r_alf[neo];
	if(!strcmp(what,"r_rot")) return r_rot[neo];
	if(!strcmp(what,"epsx_c")) return epsx_c[neo];
	if(!strcmp(what,"epsy_c")) return epsy_c[neo];
	if(!strcmp(what,"epsx_u")) return epsx_u[neo];
	if(!strcmp(what,"epsy_u")) return epsy_u[neo];
	if(!strcmp(what,"eps_dx")) return eps_dx[neo];
	if(!strcmp(what,"eps_lx")) return eps_lx[neo];
	if(!strcmp(what,"rl2")) return rl2[neo];
	if(!strcmp(what,"rot")) return rot[neo];
	if(!strcmp(what,"envx")) return envx[neo];
	if(!strcmp(what,"envy")) return envy[neo];
	if(!strcmp(what,"divx")) return divx[neo];
	if(!strcmp(what,"divy")) return divy[neo];
	if(!strcmp(what,"round")) return roundness[neo];
	if(!strcmp(what,"roundp")) return roundnessp[neo];
	if(!strcmp(what,"slopex")) return slopex[neo];
	if(!strcmp(what,"slopey")) return slopey[neo];
	if(!strcmp(what,"xmax")) return xmax[neo];
	if(!strcmp(what,"ymax")) return ymax[neo];
	if(!strcmp(what,"zmax")) return zmax[neo];
	if(!strcmp(what,"tmax")) return tmax[neo];
	if(!strcmp(what,"xmin")) return xmin[neo];
	if(!strcmp(what,"ymin")) return ymin[neo];
	if(!strcmp(what,"zmin")) return zmin[neo];
	if(!strcmp(what,"tmin")) return tmin[neo];
	if(!strcmp(what,"sigma_p")) return sigma_p[neo];
	if(!strcmp(what,"coolsigma")) return coolsigma[neo];
	if(!strcmp(what,"dp_max")) return dp_max[neo];
	if(!strcmp(what,"dp_min")) return dp_min[neo];
	if(!strcmp(what,"sigma_l")) return sigma_l[neo];
	if(!strcmp(what,"sigma_t")) return sigma_t[neo];
	if(!strcmp(what,"dl_max")) return dl_max[neo];
	if(!strcmp(what,"dl_min")) return dl_min[neo];
	if(!strcmp(what,"epss_dl")) return epss[neo]/sigma_l[neo];



	if(!strcmp(what,"xm")) return xm;
	if(!strcmp(what,"xpm")) return xpm;
	if(!strcmp(what,"xx")) return xx;
	if(!strcmp(what,"xxp")) return xxp;
	if(!strcmp(what,"xpxp")) return xpxp;
	if(!strcmp(what,"ym")) return ym;
	if(!strcmp(what,"ypm")) return ypm;
	if(!strcmp(what,"yy")) return yy;
	if(!strcmp(what,"yyp")) return yyp;
	if(!strcmp(what,"ypyp")) return ypyp;
	if(!strcmp(what,"xy")) return xy;
	if(!strcmp(what,"xyp")) return xyp;
	if(!strcmp(what,"xpy")) return xpy;
	if(!strcmp(what,"xpyp")) return xpyp;
	if(!strcmp(what,"coupl")) return xy*xy+xpyp*xpyp+xyp*xyp+xpy*xpy;

	if(!strcmp(what,"dlm")) return dlm;
	if(!strcmp(what,"dldl")) return dldl;
	if(!strcmp(what,"dlxp")) return dlxp;
	if(!strcmp(what,"dlxpp")) return dlxpp;

	if(!strcmp(what,"dpm")) return dpm;
	if(!strcmp(what,"dpdp")) return dpdp;
	if(!strcmp(what,"ddxp")) return ddxp;
	if(!strcmp(what,"ddxpp")) return ddxpp;



	if(!strcmp(what,"dd")) return dd;
	if(!strcmp(what,"ld")) return ld;
	if(!strcmp(what,"ll")) return ll;
	if(!strcmp(what,"xd")) return xd;
	if(!strcmp(what,"xl")) return xl;
	if(!strcmp(what,"xpd")) return xpd;
	if(!strcmp(what,"xpl")) return xpl;
	if(!strcmp(what,"yd")) return yd;
	if(!strcmp(what,"yl")) return yl;
	if(!strcmp(what,"ypd")) return ypd;
	if(!strcmp(what,"ypl")) return ypl;

	if(!strcmp(what,"flatt")) return flattness;

	fprintf(stderr, "parameter %s is not defined in Pmla::getValue\n", what);
	printf(         "parameter %s is not defined in Pmla::getValue\n", what);
	exit(-1);
}

int Pmla::getCurve(const char * what)
{
	if(!strcmp(what,"ener")) 	return 0;
	if(!strcmp(what,"gammaE")) 	return 0;
	if(!strcmp(what,"phas")) 	return 0;
	if(!strcmp(what,"blen")) 	return 0;
	if(!strcmp(what,"disl")) 	return 0;
	if(!strcmp(what,"dislp"))  	return 0;
	if(!strcmp(what,"dispx"))   	return 0;
	if(!strcmp(what,"dispxp"))  	return 0;
	if(!strcmp(what,"dispy"))   	return 0;
	if(!strcmp(what,"dispyp"))  	return 0;
	if(!strcmp(what,"disp+"))   	return 0;
	if(!strcmp(what,"dispp+"))  	return 0;
	if(!strcmp(what,"disp0"))   	return 0;
	if(!strcmp(what,"dispp0"))  	return 0;
	if(!strcmp(what,"disp-"))   	return 0;
	if(!strcmp(what,"dispp-"))  	return 0;
	if(!strcmp(what,"orbx"))   	return 0;
	if(!strcmp(what,"orby"))   	return 0;
	if(!strcmp(what,"orbxp"))   	return 0;
	if(!strcmp(what,"orbyp"))   	return 0;
	if(!strcmp(what,"epsx"))   	return 1;
	if(!strcmp(what,"epsx_n")) 	return 1;
	if(!strcmp(what,"betx")) 	return 0;
	if(!strcmp(what,"alfx")) 	return 0;
	if(!strcmp(what,"epsy")) 	return 1;
	if(!strcmp(what,"epsy_n")) 	return 1;
	if(!strcmp(what,"bety")) 	return 0;
	if(!strcmp(what,"alfy")) 	return 0;
	if(!strcmp(what,"epss"))   	return 1;
	if(!strcmp(what,"epssdl"))   	return 1; 	// this is the enegy spread after removing the chirp
	if(!strcmp(what,"magnet")) 	return 0;
	if(!strcmp(what,"emit4d")) 	return 1;
	if(!strcmp(what,"emit4dsum")) 	return 0;
	if(!strcmp(what,"emitz")) 	return 1;
	if(!strcmp(what,"temper")) 	return 1;
	if(!strcmp(what,"temperx")) 	return 1;
	if(!strcmp(what,"tempery")) 	return 1;
	if(!strcmp(what,"ttrans")) 	return 1;
	if(!strcmp(what,"ttranm")) 	return 0;
	if(!strcmp(what,"r_eps")) 	return 0;
	if(!strcmp(what,"r_bet")) 	return 0;
	if(!strcmp(what,"r_alf")) 	return 0;
	if(!strcmp(what,"r_rot")) 	return 0;
	if(!strcmp(what,"epsx_c")) 	return 0;
	if(!strcmp(what,"epsy_c")) 	return 0;
	if(!strcmp(what,"epsx_u")) 	return 0;
	if(!strcmp(what,"epsy_u")) 	return 0;
	if(!strcmp(what,"eps_dx")) 	return 0;
	if(!strcmp(what,"eps_lx")) 	return 0;
	if(!strcmp(what,"rl2")) 	return 0;
	if(!strcmp(what,"rot")) 	return 0;
	if(!strcmp(what,"envx")) 	return 0;
	if(!strcmp(what,"envy")) 	return 0;
	if(!strcmp(what,"divx")) 	return 0;
	if(!strcmp(what,"divy")) 	return 0;
	if(!strcmp(what,"round")) 	return 0;
	if(!strcmp(what,"roundp")) 	return 0;
	if(!strcmp(what,"slopex")) 	return 0;
	if(!strcmp(what,"slopey")) 	return 0;
	if(!strcmp(what,"xmax")) 	return 0;
	if(!strcmp(what,"ymax")) 	return 0;
	if(!strcmp(what,"zmax")) 	return 0;
	if(!strcmp(what,"tmax")) 	return 0;
	if(!strcmp(what,"xmin")) 	return 0;
	if(!strcmp(what,"ymin")) 	return 0;
	if(!strcmp(what,"zmin")) 	return 0;
	if(!strcmp(what,"tmin")) 	return 0;
	if(!strcmp(what,"sigma_p")) 	return 1;
	if(!strcmp(what,"coolsigma")) 	return 1;
	if(!strcmp(what,"dp_max")) 	return 0;
	if(!strcmp(what,"dp_min")) 	return 0;
	if(!strcmp(what,"sigma_l")) 	return 1;
	if(!strcmp(what,"sigma_t")) 	return 1;
	if(!strcmp(what,"dl_max")) 	return 0;
	if(!strcmp(what,"dl_min")) 	return 0;
	if(!strcmp(what,"epss_dl")) 	return 1;



	if(!strcmp(what,"xm")) 		return 0;
	if(!strcmp(what,"xpm")) 	return 0;
	if(!strcmp(what,"xx")) 		return 0;
	if(!strcmp(what,"xxp")) 	return 0;
	if(!strcmp(what,"xpxp")) 	return 0;
	if(!strcmp(what,"ym")) 		return 0;
	if(!strcmp(what,"ypm")) 	return 0;
	if(!strcmp(what,"yy")) 		return 0;
	if(!strcmp(what,"yyp")) 	return 0;
	if(!strcmp(what,"ypyp")) 	return 0;
	if(!strcmp(what,"xy")) 		return 0;
	if(!strcmp(what,"xyp")) 	return 0;
	if(!strcmp(what,"xpy")) 	return 0;
	if(!strcmp(what,"xpyp")) 	return 0;
	if(!strcmp(what,"coupl")) 	return 0;

	if(!strcmp(what,"dlm")) 	return 0;
	if(!strcmp(what,"dldl")) 	return 0;
	if(!strcmp(what,"dlxp")) 	return 0;
	if(!strcmp(what,"dlxpp")) 	return 0;

	if(!strcmp(what,"dpm")) 	return 0;
	if(!strcmp(what,"dpdp")) 	return 0;
	if(!strcmp(what,"ddxp")) 	return 0;
	if(!strcmp(what,"ddxpp")) 	return 0;



	if(!strcmp(what,"dd")) 		return 0;
	if(!strcmp(what,"ld")) 		return 0;
	if(!strcmp(what,"ll")) 		return 0;
	if(!strcmp(what,"xd")) 		return 0;
	if(!strcmp(what,"xl")) 		return 0;
	if(!strcmp(what,"xpd")) 	return 0;
	if(!strcmp(what,"xpl")) 	return 0;
	if(!strcmp(what,"yd")) 		return 0;
	if(!strcmp(what,"yl")) 		return 0;
	if(!strcmp(what,"ypd")) 	return 0;
	if(!strcmp(what,"ypl")) 	return 0;

	if(!strcmp(what,"flatt")) 	return 0;

	fprintf(stderr, "parameter %s is not defined in Pmla::getValue\n", what);
	exit(-1);
}


void Pmla::longFit(double ** pt, int np)
{
// 	np=1;
// 	pt++;
	const int nv = 2;
	double jac[nv][nv];
	double pp[nv];
	double cf[nv];

	memset(jac,0, nv*nv*sizeof(double));
	memset(pp, 0,   nv*sizeof(double));
	memset(cf, 0,   nv*sizeof(double));

	double dpm = 0.;
        for(int n = 0; n < np; n++) dpm += pt[n][5];
        dpm /= np;
        for(int n = 0; n < np; n++)
        {
		double jc[nv];
		jc[0]=sin(pt[n][4]*M_PI/180.);
// 		jc[0]=   (pt[n][4]*M_PI/180.);
// 		jc[0]=pt[n][4];
		jc[1]=cos(pt[n][4]*M_PI/180.);
// 		jc[2]=sin(pt[n][4]*3.*M_PI/180.);
// 		jc[3]=cos(pt[n][4]*3.*M_PI/180.);
// 		jc[nv-1]=1.;

		for(int i=0; i<nv; i++)
		{
			for(int j=0; j<nv; j++) jac[i][j] += jc[i]*jc[j];
			pp[i] += (pt[n][5]-dpm) * jc[i];
// 			pp[i] += (pt[n][5]) * jc[i];
		}
	}
	for(int i=0; i<nv; i++)
	{
		for(int j=0; j<nv; j++)
		{
// 			jac[i][j] /= np;
			printf("%e  ", jac[i][j]);
		}
// 		pp[i] /= np;
		printf(" = %e \n ", pp[i]);
	}

	jac[0][0] = 1./jac[0][0];	// inversion of a 1x1  matrix
// 	double det = dminv( (double *) jac, nv);
// 	printf("det = %e\n", det);
// 	if( fabs(det) < 1e-10)
// 	{
// 		fprintf(stderr, "in Pmla::longFit neo= %d, det = %e\n", neo, det);
// 		return;
// 	}
	for(int i=0; i<nv; i++)
	{
		for(int j=0; j<nv; j++)
		{
			cf[i] = jac[i][j]*pp[j];
		}
	}
	printf("lfit %e %e %e %e %e \n", cf[0], cf[1],cf[2],cf[3],cf[4]);

	
	FILE * fa = fopen("_a.agr", "w");
	if(!fa)
	{
		fprintf(stderr," cant open %s\n", "_a.agr");
		return;
	}
        fprintf(fa, "@with g0\n");
        fprintf(fa, "@    s%d symbol 9\n",0);
        fprintf(fa, "@    s%d symbol size 0.060000\n",0);
        fprintf(fa, "@    s%d line type 0\n",0);


        for(int n = 0; n < np; n++)
        {
		double jc[nv];
		jc[0]=sin(pt[n][4]*M_PI/180.);
// 		jc[0]=   (pt[n][4]*M_PI/180.);
// 		jc[0]=pt[n][4];
		jc[1]=cos(pt[n][4]*M_PI/180.);
// 		jc[2]=sin(pt[n][4]*3.*M_PI/180.);
// 		jc[3]=cos(pt[n][4]*3.*M_PI/180.);
// 		jc[nv-1]=1.;

		for(int i=0; i<nv; i++) pt[n][5] -= cf[i]*jc[i];
		fprintf(fa, "%e %e\n", pt[n][4],pt[n][5]);
	}
	fclose(fa);


}




void Pmla::twiss(int print)
{
	twiss(print, dcord, nbuf);
}




void Pmla::twiss(int print, double ** pt, int np)
{
        xm = 0;		// < x >
        xpm = 0;	// < x' >
        xx = 0;		// < x^2 >
        xxp = 0;	// < x * x' >
        xpxp = 0;	// < x'^2 >
        ym = 0;		// < y >
        ypm = 0;	// < y' >
        yy = 0;		// < y^2 >
        yyp = 0;	// < y * y' >
        ypyp=0;		// < y'^2 >
	xy  = 0;	// < x * y >
	xyp  = 0;	// < x * y' >
	xpy  = 0;	// < x' * y' >
	xpyp  = 0;

        dlm  = 0;	// < phi >
        dldl  = 0;      // < delta phi^2 >
        dlxp = 0;	// < x * delta phi >
        dlxpp= 0;	// < x' * delta phi>

	double ee = 0.;
	double de_cool = 0.;

        dpm  = 0;	// < E >
        dpdp = 0;	// < delta E^2 >
        ddxp = 0;	// < x*delta E >
        ddxpp= 0;	// < x' * delta E >
        ddyp = 0;	// < y * delta E >
        ddypp= 0;	// < y' * delta E >


	dd =0.;		//  < (delta p/p)^2 >
	ld =0.;		//  < delta phi * delta p/p >
	ll =0.;		//  < delta phi ^2 >
	xd =0.;		//  < x * delta p/p >
	xl =0.;		//  < x * delta phi >
	xpd =0.;	//  < x' * delta p/p >
	xpl =0.;	//  < x' * delta phi >
	yd =0.;		//  < y * delta p/p >
	yl =0.;    	//  < y * delta phi >
	ypd =0.;	//  < y' * delta p/p >
	ypl =0.;	//  < y' * delta phi >

	
	xmax[neo]= -1e33;
	ymax[neo]= -1e33;
	zmax[neo]= -1e33;
	xmin[neo]=  1e33;
	ymin[neo]=  1e33;
	zmin[neo]=  1e33;
	dp_max[neo]=  -1e33;
	dp_min[neo]=   1e33;
	dl_max[neo]=  -1e33;
	dl_min[neo]=   1e33;

        for(int i = 0; i < np; i++)
        {
            dlm += pt[i][4];
            dpm += pt[i][5];
	    ee += pt[i][5]*pt[i][5];
	    if(pt[i][0] > xmax[neo]) xmax[neo] = pt[i][0];
	    if(pt[i][0] < xmin[neo]) xmin[neo] = pt[i][0];

	    if(pt[i][2] > ymax[neo]) ymax[neo] = pt[i][2];
	    if(pt[i][2] < ymin[neo]) ymin[neo] = pt[i][2];

	    if(pt[i][4] > dl_max[neo]) dl_max[neo] = pt[i][4];
	    if(pt[i][4] < dl_min[neo]) dl_min[neo] = pt[i][4];

	    if(pt[i][5] > dp_max[neo]) dp_max[neo] = pt[i][5];
	    if(pt[i][5] < dp_min[neo]) dp_min[neo] = pt[i][5];
	}
        dpm /= np;
	ee /= np;
        dlm /= np;
// 	printf("dmp %e dpmax %e dpmin %e\n", dpm, dp_max[neo], dp_min[neo]);

	dp_max[neo] = fabs(dp_max[neo] - dpm)/dpm;
	dp_min[neo] = fabs(dp_min[neo] - dpm)/dpm;
	dl_max[neo] = fabs(dl_max[neo] - dlm);
	dl_min[neo] = fabs(dl_min[neo] - dlm);

	double gamma = dpm/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);
	double betaE = betagamma/gamma;
	double gg = sqrt(square(dpm/erest[0])+1.);
	double Eavg = (gg-1.)*erest[0];
	ener1[neo]=Eavg;

	zmax[neo] = dl_max[neo]*wavel*betaE/360.*0.01; // length in meters
	zmin[neo] = dl_min[neo]*wavel*betaE/360.*0.01; // length in meters
	tmax[neo] = dl_max[neo]/(freq*1e6 *360.); // length in seconds
	tmin[neo] = dl_min[neo]/(freq*1e6 *360.); // length in seconds

	sigma_t[neo]=sqrt(dldl)/(freq*1e6 *360.); // length in seconds


        for(int i = 0; i < np; i++)
        {
	    double c4 = (pt[i][4]-dlm);
            dldl  += c4*c4;
            dlxp  += pt[i][0]*c4;
            dlxpp += pt[i][1]*c4;

	    double c5 = (pt[i][5]-dpm);
            dpdp  += c5*c5;
            ddxp  += pt[i][0]*c5;
            ddxpp += pt[i][1]*c5;
            ddyp  += pt[i][2]*c5;
            ddypp += pt[i][3]*c5;

        }
        disl[neo]   = dldl == 0 ? 0 : dlxp/dldl;
        dislp[neo]  = dldl == 0 ? 0 : dlxpp/dldl;
        dispx[neo]  = dpdp == 0 ? 0 : ddxp/dpdp;
        dispxp[neo] = dpdp == 0 ? 0 : ddxpp/dpdp;
        dispy[neo]  = dpdp == 0 ? 0 : ddyp/dpdp;
        dispyp[neo] = dpdp == 0 ? 0 : ddypp/dpdp;
        dpdp /= np;
        dldl /= np;



        for(int i = 0; i < np; i++)
        {
// 		printf("%f %f %f %f %f %f %f\n", pt[i][0], pt[i][1], pt[i][2], pt[i][3], pt[i][4], pt[i][5], pt[i][6]);

// 	    	double c5 = (pt[i][5]-dpm)/dpm;
// 	    	double xt =pt[i][0]-disp[neo]*c5;
// 	    	double xpt=pt[i][1]-dispp[neo]*c5;

         	double xt =pt[i][0];
	    	double xpt=pt[i][1];
	    	double yt =pt[i][2];
	    	double ypt=pt[i][3];
	    	double l = pt[i][4] - dlm;
	    	double d =(pt[i][5] - dpm)/dpm;
	    	double dc =(pt[i][5] - coolenergy)/coolenergy;
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

	    	dd   += d*d;
	    	ld   += l*d;
	    	ll   += l*l;
	    	xd   += xt*d;
	    	xl   += xt*l;
	    	xpd  += xpt*d;
	    	xpl  += xpt*l;
	    	yd   += yt*d;
	    	yl   += yt*l;
	    	ypd  += ypt*d;
	    	ypl  += ypt*l;

		de_cool += dc*dc;
        }
        xx /= np;
        xm /= np;
        xpm /= np;
        xxp /= np;
        xpxp /= np;

        yy /= np;
        ym /= np;
        ypm /= np;
        yyp /= np;
        ypyp /= np;

	xy /= np;
	xyp /= np;
	xpy /= np;
	xpyp /= np;


	dd  /= np;
	ld  /= np;
	ll  /= np;
	xd  /= np;
	xl  /= np;
	xpd  /= np;
	xpl  /= np;
	yd  /= np;
	yl  /= np;
	ypd  /= np;
	ypl  /= np;

	de_cool /= np;

	sigma_p[neo]=sqrt(dd);
	coolsigma[neo]=sqrt(de_cool);
// 	printf("sigma pppppppp %e %e\n", sigma_p[neo], sqrt(ee - dpm*dpm)/dpm);
	sigma_l[neo]=sqrt(dldl)*wavel*betaE/360.*0.01; // length in meters
	sigma_t[neo]=sqrt(dldl)/(freq*1e6 *360.); // length in seconds
	temperature[neo] = sqrt(gamma * gamma * ( xpxp + ypyp)  + dpdp/(dpm*dpm) );
	temperature[neo] = sqrt(gamma * gamma * ( xpxp + ypyp)  + dpdp/(dpm*dpm) + 5e-6/40. + 25e-8);
	temperxp[neo] = gamma * sqrt(xpxp);
	temperyp[neo] = gamma * sqrt(ypyp);

	double m6[6][6] = { xx, xxp, xy, xyp, xl, xd, xxp, xpxp, xpy, xpyp, xpl, xpd, xy, xpy, yy, yyp, yl, yd,
		xyp, xpyp, yyp, ypyp, ypl, ypd, xl, xpl, yl, ypl, ll, ld, xd, xpd, yd, ypd, ld, dd};
// 	printf("sigmat %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", neo,
// 	 	xx, xxp, xy, xyp, xl, xd, xxp, xpxp, xpy, xpyp, xpl, xpd, xy, xpy, yy, yyp, yl, yd,
// 		xyp, xpyp, yyp, ypyp, ypl, ypd, xl, xpl, yl, ypl, ll, ld, xd, xpd, yd, ypd, ld, dd);

	double emit6d = dminv( (double *) m6, 6);

	magnetization[neo] = betagamma*(xyp-xpy);

	emit4d[neo] = betagamma*sqrt(sqrt(
	    xx*xpxp*yy*ypyp
	   -xx*xpxp*yyp*yyp
	   -xx*xpy*xpy*ypyp
	   +2*xx*xpy*xpyp*yyp
	   -xx*xpyp*xpyp*yy
	   -xxp*xxp*yy*ypyp
	   +xxp*xxp*yyp*yyp
	   +2*xxp*xpy*xy*ypyp
	   -2*xxp*xpy*xyp*yyp
	   -2*xxp*xpyp*xy*yyp
	   +2*xxp*xpyp*xyp*yy
	   -xpxp*xy*xy*ypyp
	   +2*xy*xpxp*xyp*yyp
	   -2*xy*xpyp*xyp*xpy
	   -xpxp*xyp*xyp*yy
	   +xyp*xyp*xpy*xpy
	   +xy*xy*xpyp*xpyp
	   ));
	
	orbx[neo]=xm;
	orby[neo]=ym;
	orbxp[neo]=xpm;
	orbyp[neo]=ypm;
	emitz[neo]= sqrt( ll*dd -ld*ld )* dpm*1000.;
	double tempt =  xpxp+ypyp;
	tempt *= dpm*dpm*1e6/0.511000;
	eps_dx[neo]= sqrt((ddxpp*ddxpp*xpxp - 2*xxp*ddxp*ddxpp + ddxp*ddxp*xx)/dpdp);
	eps_lx[neo]= sqrt((dlxpp*dlxpp*xpxp - 2*xxp*dlxp*dlxpp + dlxp*dlxp*xx)/dldl);

// 	double mmt[4][4];
// 	for(int i=0; i<4; i++)
// 	{
// 		for(int j=0; j<4; j++)
// 		{
// 			mmt[i][j] = (pt[i][j]-pt[i+4][j])/(2.*icord[j][j]);
// 			printf("%15e  ", icord[i][j]);
// 		}
// 		printf("\n");
// 	}



	ener[neo] = dpm;
	gammaE[neo] = gamma;
	phas[neo] = dlm;
	blen[neo] = sqrt(ll);
//         epsx[neo] = sqrt(xx*xpxp -xxp*xxp);
//         betx[neo] =   xx/epsx[neo];
//         alfx[neo] = -xxp/epsx[neo];
        epsx[neo] = sqrt(xx*xpxp -xxp*xxp - xx*xpm*xpm -xpxp*xm*xm +2.*xm *xpm*xxp);
        epsx_n[neo] = epsx[neo]*betagamma;
        betx[neo] =   (xx-xm*xm)/epsx[neo];
        alfx[neo] = -(xxp-xm*xpm)/epsx[neo];
	envx[neo] = sqrt(xx-xm*xm);
	divx[neo] = sqrt(xpxp-xpm*xpm);
	slopex[neo] = (xxp-xm*xpm)/(xx-xm*xm);
//         epsy[neo] = sqrt(yy*ypyp -yyp*yyp);
//         bety[neo] =   yy/epsy[neo];
//         alfy[neo] = -yyp/epsy[neo];
        epsy[neo] = sqrt(yy*ypyp -yyp*yyp - yy*ypm*ypm -ypyp*ym*ym +2.*ym *ypm*yyp);
        epsy_n[neo] = epsy[neo]*betagamma;
        bety[neo] =   (yy-ym*ym)/epsy[neo];
        alfy[neo] = -(yyp-ym*ypm)/epsy[neo];
	envy[neo] = sqrt(yy-ym*ym);
	divy[neo] = sqrt(ypyp-ypm*ypm);
	slopey[neo] = (yyp-ym*ypm)/(yy-ym*ym);

	roundness[neo] = fabs(envx[neo] -envy[neo]);
// 	roundnessp[neo] = fabs(slopex[neo] -slopey[neo]);
	roundnessp[neo] = fabs(xxp - yyp);

// 	sigma_p[neo]=sqrt(dpdp)/dpm;
// 	sigma_l[neo]=sqrt(dldl)*wavel/360.*0.01; // length in meters
        epss[neo] = sqrt(ll*dd -ld*ld )*wavel/360.*0.01; // length in meters;
	double eprod = 0.5*epsx[neo]*epsy[neo];
	if(print & 1) 
	{
// 	printf("%4d %7s %10.4f %11.4e %11.4e  %11.4e   %10.4f %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f %10.4f\n",
// 			neo,type[neo],zloc[neo],emit4d[neo],eprod, emit6d,betx[neo],alfx[neo],epsx[neo],dispx[neo], dispxp[neo], bety[neo],alfy[neo],epsy[neo]);
// 			neo,type[neo],zloc[neo],emit4d[neo],eprod, ener[neo],betx[neo],alfx[neo],epsx[neo],dispx[neo], dispxp[neo], bety[neo],alfy[neo],epsy[neo]);
//
	}
	if(print & 2)
	{
		printf("%d np, ngood, nbuf    %d     %d     %d     %d\n", neo, np, ngood, nbuf, npoints);
		printf("%d zloc               %f\n", neo, zloc[neo]);
		printf("%d orbit              %15.7e     %15.7e     %15.7e    %15.7e\n", neo, xm, xpm, ym, ypm);
		printf("%d dispersion         %15.7e     %15.7e     %15.7e    %15.7e\n", neo, dispx[neo], dispxp[neo], dispy[neo], dispyp[neo]);
		printf("%d emittance x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",  neo,  epsx[neo],epsy[neo],emit4d[neo]/betagamma,eprod);
		printf("%d normalizd x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",  neo,  epsx_n[neo],epsy_n[neo],emit4d[neo],eprod*betagamma);
		printf("%d sigma_xy, slope_xy %15.7e     %15.7e     %15.7e    %15.7e\n",  neo,  envx[neo], envy[neo], slopex[neo], slopey[neo]);
		printf("%d e longitudinal     %15.7e\n", neo, epss[neo]);
		printf("%d gamma,beta         %15.7e     %15.7e\n", neo, gamma, sqrt(1-1/(gamma*gamma)));
		printf("%d disl               %15.7e\n", neo, disl[neo]);
		printf("%d dislp              %15.7e\n", neo, dislp[neo]);
		printf("%d dldl               %15.7e\n", neo, dldl);
		printf("%d dpdp               %15.7e\n", neo, dpdp);
		printf("%d EAvg               %15.7e\n", neo, Eavg);
		printf("%d dpAvg              %15.7e\n", neo, dpm);
		printf("%d coolsigma          %15.7e\n", neo, coolsigma[neo]);
		printf("%d phasAvg            %15.7e\n", neo, dlm);
		printf("%d sigma_p            %15.7e\n", neo, sigma_p[neo]);
		printf("%d sigma_l            %15.7e\n", neo, sigma_l[neo]);
		printf("%d xSigMinMax         %f %f %f\n", neo, sqrt(xx), xmin[neo], xmax[neo]);
		printf("%d ySigMinMax         %f %f %f\n", neo, sqrt(yy), ymin[neo], ymax[neo]);
		printf("%d zSigMinMax         %f %f %f\n", neo, sigma_l[neo], zmin[neo], zmax[neo]);
		printf("%d eSigMinMax         %f %f %f\n", neo, sigma_p[neo], dp_min[neo], dp_max[neo]);
		printf("xx, xpxp, yy, ypyp, tempt = %e %e %e %e %e\n", xx, xpxp, yy, ypyp, tempt);
		printf("sqrt(xx), sqrt(yy) = %e %e gamma=%f beta=%f\n", sqrt(xx), sqrt(yy), gamma, sqrt(1-1/(gamma*gamma)));
		printf("<x>, <y> = %e %e gamma=%f\n", xm, ym, gamma);
		printf("%d eps_dx, eps_dl = %e, %e\n", neo, eps_dx[neo], eps_lx[neo]);


		printf("%d ener[neo], = %f %f %f %f\n", neo,ener[neo],ener[neo],ener[neo],ener[neo]);
		printf("%d betx[neo],alfx[neo],epsx[neo] = %f %f %e \n", neo,betx[neo],alfx[neo],epsx[neo]);
		printf("%d bety[neo],alfy[neo],epsy[neo] = %f %f %e \n", neo,bety[neo],alfy[neo],epsy[neo]);
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xx,   xxp,  xy,   xyp  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xxp,  xpxp, xpy,  xpyp );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xy,   xpy,  yy,   yyp  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xyp,  xpyp, yyp,  ypyp );
		printf("\n%d M=%e\n",  neo,(xyp-xpy)*betagamma);

		printf("%d ze         %15.8e %15.8e %15.8e     %15.8e %15.8e %15.8e    %15.8e %15.8e\n", neo, sigma_l[neo], zmin[neo], zmax[neo], sigma_p[neo], dp_min[neo], dp_max[neo], epsx_n[neo],epsy_n[neo]);
	}
}





// this is a backup from before mixing energy and momentum
// gonna use that to clean up
void Pmla::calcSigmaMatrix(double ** pt, int np, int print)
{
	for(int j=0; j<7; j++)
	{
		for(int k=0; k<7; k++) sigmaMat[j][k] =0.;
		maxPart[j]  = -1e33;
		minPart[j]  =  1e33;
	}



        for(int i = 0; i < np; i++)
        {
		double savept5 = pt[i][5];
		pt[i][6] = 1.;
                double gammaL= pt[i][5]/erest[0] +1.; 
                pt[i][5] = elmass*sqrt(gammaL*gammaL-1.);

		for(int j=0; j<7; j++)
		{
	    		if(pt[i][j] > maxPart[j]) maxPart[j] = pt[i][j];
	    		if(pt[i][j] < minPart[j]) minPart[j] = pt[i][j];

			for(int k=j; k<7; k++) sigmaMat[j][k] += pt[i][j] * pt[i][k];
		}
		pt[i][5]=savept5;
	}


	for(int j=0; j<7; j++)
	{
		for(int k=j; k<7; k++)
		{
			sigmaMat[j][k] /= np;
			if( j != k) sigmaMat[k][j] = sigmaMat[j][k];
		}
	}

	for(int j=0; j<7; j++)
	{
		for(int k=j; k<7; k++)
		{
			sigmaMat2[k][j] = sigmaMat2[j][k] = sigmaMat[j][k] - sigmaMat[j][6]*sigmaMat[k][6];
		}
	}

	for(int j=0; j<7; j++)
	{
		sigmaMat2[5][j] /= sigmaMat[5][6];
		sigmaMat2[j][5] /= sigmaMat[6][5];
	}




	if(print)
	{
		printf("sigmaMat = \n");
		for(int j=0; j<7; j++)
		{
			for(int k=0; k<7; k++) printf("  %e  ", sigmaMat[j][k]);
			printf("\n");
		}
		printf("\n");


		printf("sigmaMat2 = \n");
		for(int j=0; j<7; j++)
		{
			for(int k=0; k<7; k++) printf("  %e  ", sigmaMat2[j][k]);
			printf("\n");
		}
		printf("\n");
	}
}


void Pmla::twiss1(int print)
{
	twiss1(print, dcord, nbuf);
}
// this is a backup from before mixing energy and momentum
// gonna use that to clean up
void Pmla::twiss1(int print, double ** pt, int np)
{
	calcSigmaMatrix(pt,  np, 0);
#if 0
	xmax[neo] = maxPart[0];
	ymax[neo] = maxPart[2];
	zmax[neo] = maxPart[4] * wavel/360.*0.01; // length in meters;
	xmin[neo] = minPart[0];
	ymin[neo] = minPart[2];
	zmin[neo] = minPart[4] * wavel/360.*0.01; // length in meters;



        xm    = sigmaMat[6][0];
        xpm   = sigmaMat[6][1];
        ym    = sigmaMat[6][2];
        ypm   = sigmaMat[6][3];
        dlm   = sigmaMat[6][4];
        dpm   = sigmaMat[6][5];


        xx    = sigmaMat2[0][0];
        xxp   = sigmaMat2[0][1];
        xpxp  = sigmaMat2[1][1];
        yy    = sigmaMat2[2][2];
        yyp   = sigmaMat2[2][3];
        ypyp  = sigmaMat2[3][3];
	xy    = sigmaMat2[0][2];
	xyp   = sigmaMat2[0][3];
	xpy   = sigmaMat2[1][2];
	xpyp  = sigmaMat2[1][3];


	dd    = sigmaMat2[5][5];
	ld    = sigmaMat2[4][5];
	ll    = sigmaMat2[4][4];
	xd    = sigmaMat2[0][5];
	xl    = sigmaMat2[0][4];
	xpd   = sigmaMat2[1][5];
	xpl   = sigmaMat2[1][4];
	yd    = sigmaMat2[2][5];
	yl    = sigmaMat2[2][4];
	ypd   = sigmaMat2[3][5];
	ypl   = sigmaMat2[3][4];



        dpdp  = sigmaMat2[5][5];
        dldl  = sigmaMat2[4][4];

        dlxp  = sigmaMat2[0][4];
        dlxpp = sigmaMat2[1][4];
        ddxp  = sigmaMat2[0][5];
        ddxpp = sigmaMat2[1][5];
        ddyp  = sigmaMat2[2][5];
        ddypp = sigmaMat2[3][5];

	double gamma= dpm/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);

        disl[neo]   = dldl == 0 ? 0 : dlxp/dldl;
        dislp[neo]  = dldl == 0 ? 0 : dlxpp/dldl;
        dispx[neo]  = dpdp == 0 ? 0 : ddxp/dpdp;
        dispxp[neo] = dpdp == 0 ? 0 : ddxpp/dpdp;
        dispy[neo]  = dpdp == 0 ? 0 : ddyp/dpdp;
        dispyp[neo] = dpdp == 0 ? 0 : ddypp/dpdp;

	sigma_p[neo]=sqrt(dpdp)/dpm;
	sigma_l[neo]=sqrt(dldl)*wavel/360.*0.01; // length in meters




	double m6[6][6] = { xx, xxp, xy, xyp, xl, xd, xxp, xpxp, xpy, xpyp, xpl, xpd, xy, xpy, yy, yyp, yl, yd,
		xyp, xpyp, yyp, ypyp, ypl, ypd, xl, xpl, yl, ypl, ll, ld, xd, xpd, yd, ypd, ld, dd};

	double emit6d = dminv( (double *) m6, 6);

	magnetization[neo] = betagamma*(xyp-xpy);

	emit4d[neo] = betagamma*sqrt(sqrt(
	    xx*xpxp*yy*ypyp
	   -xx*xpxp*yyp*yyp
	   -xx*xpy*xpy*ypyp
	   +2*xx*xpy*xpyp*yyp
	   -xx*xpyp*xpyp*yy
	   -xxp*xxp*yy*ypyp
	   +xxp*xxp*yyp*yyp
	   +2*xxp*xpy*xy*ypyp
	   -2*xxp*xpy*xyp*yyp
	   -2*xxp*xpyp*xy*yyp
	   +2*xxp*xpyp*xyp*yy
	   -xpxp*xy*xy*ypyp
	   +2*xy*xpxp*xyp*yyp
	   -2*xy*xpyp*xyp*xpy
	   -xpxp*xyp*xyp*yy
	   +xyp*xyp*xpy*xpy
	   +xy*xy*xpyp*xpyp
	   ));
	
	emitz[neo]= sqrt( ll*dd -ld*ld )* dpm*1000.;
	double tempt =  xpxp+ypyp;
	tempt *= dpm*dpm*1e6/0.511000;
	eps_dx[neo]= sqrt((ddxpp*ddxpp*xpxp - 2*xxp*ddxp*ddxpp + ddxp*ddxp*xx)/dpdp);
	eps_lx[neo]= sqrt((dlxpp*dlxpp*xpxp - 2*xxp*dlxp*dlxpp + dlxp*dlxp*xx)/dldl);



	ener[neo] = dpm;
	gammaE[neo] = gamma;
	phas[neo] = sigmaMat[6][4];
	blen[neo] = sqrt(sigmaMat2[4][4]);
        epsx[neo] = sqrt(xx*xpxp -xxp*xxp - xx*xpm*xpm -xpxp*xm*xm +2.*xm *xpm*xxp);
        epsx_n[neo] = epsx[neo]*betagamma;
        betx[neo] =   (xx-xm*xm)/epsx[neo];
        alfx[neo] = -(xxp-xm*xpm)/epsx[neo];
	envx[neo] = sqrt(xx-xm*xm);
	divx[neo] = sqrt(xpxp-xpm*xpm);
        epsy[neo] = sqrt(yy*ypyp -yyp*yyp - yy*ypm*ypm -ypyp*ym*ym +2.*ym *ypm*yyp);
        epsy_n[neo] = epsy[neo]*betagamma;
        bety[neo] =   (yy-ym*ym)/epsy[neo];
        alfy[neo] = -(yyp-ym*ypm)/epsy[neo];
	envy[neo] = sqrt(yy-ym*ym);
	divy[neo] = sqrt(ypyp-ypm*ypm);

        epss[neo] = sqrt(ll*dd -ld*ld )*wavel/360.*0.01; // length in meters;
	double eprod = 0.5*epsx[neo]*epsy[neo];
	if(print & 1) 
	{
// 	printf("%4d %7s %10.4f %11.4e %11.4e  %11.4e   %10.4f %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f %10.4f\n",
// 			neo,type[neo],zloc[neo],emit4d[neo],eprod, emit6d,betx[neo],alfx[neo],epsx[neo],dispx[neo], dispxp[neo], bety[neo],alfy[neo],epsy[neo]);
// 			neo,type[neo],zloc[neo],emit4d[neo],eprod, ener[neo],betx[neo],alfx[neo],epsx[neo],dispx[neo], dispxp[neo], bety[neo],alfy[neo],epsy[neo]);
//
	}
	if(print & 2)
	{
		printf("%d np, ngood, nbuf    %d     %d     %d     %d\n", neo, np, ngood, nbuf, npoints);
		printf("%d zloc               %f\n", neo, zloc[neo]);
		printf("%d orbit              %15.7e     %15.7e     %15.7e    %15.7e\n", neo, xm, xpm, ym, ypm);
		printf("%d dispersion         %15.7e     %15.7e     %15.7e    %15.7e\n", neo, dispx[neo], dispxp[neo], dispy[neo], dispyp[neo]);
		printf("%d emittance x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",  neo,  epsx[neo]*1e6,epsy[neo]*1e6,emit4d[neo]/betagamma*1e6,eprod*1e6);
		printf("%d normalizd x y 4 xy %15.7e     %15.7e     %15.7e    %15.7e\n",  neo,  epsx_n[neo]*1e6,epsy_n[neo]*1e6,emit4d[neo]*1e6,eprod*betagamma*1e6);
		printf("%d e longitudinal     %15.7e\n", neo, epss[neo]);
		printf("%d gamma,beta         %15.7e     %15.7e\n", neo, gamma, sqrt(1-1/(gamma*gamma)));
		printf("%d disl               %15.7e\n", neo, disl[neo]);
		printf("%d dislp              %15.7e\n", neo, dislp[neo]);
		printf("%d dldl               %15.7e\n", neo, dldl);
		printf("%d dpdp               %15.7e\n", neo, dpdp);
		printf("%d dpAvg              %15.7e\n", neo, dpm);
		printf("%d phasAvg            %15.7e\n", neo, dlm);
		printf("%d sigma_p            %15.7e\n", neo, sigma_p[neo]);
		printf("%d sigma_l            %15.7e\n", neo, sigma_l[neo]);
		printf("%d xminmax            %f %f\n", neo, xmin[neo], xmax[neo]);
		printf("%d yminmax            %f %f\n", neo, ymin[neo], ymax[neo]);
		printf("%d zminmax            %f %f\n", neo, zmin[neo], zmax[neo]);
		printf("xx, xpxp, yy, ypyp, tempt = %f %f %f %f %f\n", xx, xpxp, yy, ypyp, tempt);
		printf("sqrt(xx), sqrt(yy) = %f %f gamma=%f beta=%f\n", sqrt(xx), sqrt(yy), gamma, sqrt(1-1/(gamma*gamma)));
		printf("<x>, <y> = %f %f gamma=%f\n", xm, ym, gamma);
		printf("%d eps_dx, eps_dl = %e, %e\n", neo, eps_dx[neo], eps_lx[neo]);


		printf("%d ener[neo], = %f %f %f %f\n", neo,ener[neo],ener[neo],ener[neo],ener[neo]);
		printf("%d betx[neo],alfx[neo],epsx[neo] = %f %f %e \n", neo,betx[neo],alfx[neo],epsx[neo]);
		printf("%d bety[neo],alfy[neo],epsy[neo] = %f %f %e \n", neo,bety[neo],alfy[neo],epsy[neo]);
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xx,   xxp,  xy,   xyp  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xxp,  xpxp, xpy,  xpyp );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xy,   xpy,  yy,   yyp  );
		printf("sigma\t %15.7e   %15.7e  %15.7e  %15.7e\n", xyp,  xpyp, yyp,  ypyp );
		printf("\n%d M=%e\n",  neo,(xyp-xpy)*betagamma);
	}
#endif
}






void quicksort (double* a[], int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
    int i=lo, j=hi;
    double  x=a[(lo+hi)/2][6];

    //  partition
    do
    {    
            while (a[i][6]<x) i++; 
            while (a[j][6]>x) j--;
            if (i<=j)
            {
		    // swap
	            double *  h=a[i]; a[i]=a[j]; a[j]=h;
	            i++; j--;
            }
     } while (i<=j);

    //  recursion
    if (lo<j) quicksort(a, lo, j);
    if (i<hi) quicksort(a, i, hi);
}




void Pmla::sumEmittance(double * a[], int good)
{

	printf("\n good=%d\n\n", good);


        xm = 0;
        xpm = 0;
        xx = 0;
        xxp = 0;
        xpxp = 0;
        ym = 0;
        ypm = 0;
        yy = 0;
        yyp = 0;
        ypyp=0;
	xy  = 0;
	xyp  = 0;
	xpy  = 0;
	xpyp  = 0;

        dpm  = 0;
        dpdp = 0;



        for(int i = 0; i < good; i++)
        {
// 		printf("%f %f %f %f %f %f %f\n", dcord[i][0], dcord[i][1], dcord[i][2], dcord[i][3], dcord[i][4], dcord[i][5], dcord[i][6]);


         	double xt =a[i][0];
	    	double xpt=a[i][1];
	    	double yt =a[i][2];
	    	double ypt=a[i][3];
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
            	dpm  += a[i][5];
        	dpdp += a[i][5]*a[i][5];
        }
        xx /= good;
        xm /= good;
        xpm /= good;
        xxp /= good;
        xpxp /= good;

        yy /= good;
        ym /= good;
        ypm /= good;
        yyp /= good;
        ypyp /= good;

	xy /= good;
	xyp /= good;
	xpy /= good;
	xpyp /= good;

	dpm /= good;
	dpdp /= good;
	sigma_p[neo]=sqrt(dpdp-dpm*dpm)/dpm;

	double gamma= dpm/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);
	emit4d[neo] = betagamma*sqrt(sqrt(
	    xx*xpxp*yy*ypyp
	   -xx*xpxp*yyp*yyp
	   -xx*xpy*xpy*ypyp
	   +2*xx*xpy*xpyp*yyp
	   -xx*xpyp*xpyp*yy
	   -xxp*xxp*yy*ypyp
	   +xxp*xxp*yyp*yyp
	   +2*xxp*xpy*xy*ypyp
	   -2*xxp*xpy*xyp*yyp
	   -2*xxp*xpyp*xy*yyp
	   +2*xxp*xpyp*xyp*yy
	   -xpxp*xy*xy*ypyp
	   +2*xy*xpxp*xyp*yyp
	   -2*xy*xpyp*xyp*xpy
	   -xpxp*xyp*xyp*yy
	   +xyp*xyp*xpy*xpy
	   +xy*xy*xpyp*xpyp
	   ));
        double epsxb = betagamma*sqrt(xx*xpxp -xxp*xxp - xx*xpm*xpm -xpxp*xm*xm +2.*xm *xpm*xxp);
        double epsyb = betagamma*sqrt(yy*ypyp -yyp*yyp - yy*ypm*ypm -ypyp*ym*ym +2.*ym *ypm*yyp);
	double perc = double(good)/double(nbuf);
	fprintf(stderr, "emittances: %e %e %e sigma %e at %f %%\n", epsxb, epsyb, emit4d[neo], sigma_p[neo], perc);
}


void Pmla::cutHeadTail(int print, double headcut, double tailcut)
{
	sortZ('z');
	int nhc=int(nbuf*headcut);
	int ntc=int(nbuf*tailcut);
	twiss(print, scord+nhc, nbuf-nhc-ntc);
}


void Pmla::cutDP(int print, double headcut, double tailcut)
{
//         	for(int i = 0; i < nbuf; i++)
//         	{
// 			dcord[i][6] = -dcord[i][4] ;
// 		}

// 	sortZ('p');
// 	int nhc=int(nbuf*headcut);
// 	int ntc=int(nbuf*tailcut);
// 	twiss(print, scord+nhc, nbuf-nhc-ntc);
	sortZ('a');
	int nhc=0;
	int ntc=int(nbuf*(tailcut+headcut));
	twiss(print, scord+nhc, nbuf-nhc-ntc);
}

void Pmla::cutTempT(int print, double headcut, double tailcut)
{

	sortZ('t');
	int nhc=int(nbuf*headcut);
	int ntc=int(nbuf*tailcut);
	twiss(print, scord+nhc, nbuf-nhc-ntc);
}


void Pmla::sortZ(char what)
{

	if(nbuf <= 0)
	{
		fprintf(stderr, "no good particles in Pmla::sortZ\n");
		return;
	}
	for(int i=0; i< nbuf; i++) scord[i] = dcord[i];
	if( what == 'x')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = dcord[i][0] ;
		}
	}
	else
	if( what == 'y')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = dcord[i][2] ;
		}
	}
	else
	if( what == 'r')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = dcord[i][0]*dcord[i][0] + dcord[i][2]*dcord[i][2] ;
		}
	}
	else
	if( what == 't')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			double xp2 = dcord[i][1]*dcord[i][1];
			double yp2 = dcord[i][3]*dcord[i][3];
			dcord[i][6] = xp2 > yp2 ? xp2 : yp2;
// 			dcord[i][6] = dcord[i][1]*dcord[i][1] + dcord[i][3]*dcord[i][3] ;
		}
	}
	else
	if( what == 'z')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = -dcord[i][4] ;
		}
	}
	else
	if( what == 'p')
	{
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = dcord[i][5] ;
		}
	}
	else
	if( what == 'a')
	{
		double pm = 0.;
        	for(int i = 0; i < nbuf; i++)
			pm += dcord[i][5];
		pm /= nbuf;
        	for(int i = 0; i < nbuf; i++)
        	{
			dcord[i][6] = square(dcord[i][5]-pm) ;
		}
	}
	quicksort(scord,0, nbuf-1);
}


void Pmla::cutTemp(int print, double headcut, double tailcut)
{

	if(nbuf <= 0)
	{
		fprintf(stderr, "no good particles in Pmla::cutHeadTail\n");
		return;
	}
	double dpm=0.;
	for(int i=0; i< nbuf; i++)
	{
		scord[i] = dcord[i];
        	dpm += dcord[i][5];
	}
        dpm /= nbuf;
	double gamma= dpm/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);
	int good=nbuf;
        for(int i = 0; i < good; i++)
        {
		dcord[i][6] = square(gamma*dcord[i][1]) + square(gamma*dcord[i][3]) + square( (dcord[i][5]-dpm) / dpm );
	}
	quicksort(scord,0, good-1);
	int nhc=int(nbuf*headcut);
	int ntc=int(nbuf*tailcut);
	twiss(print, scord+nhc, nbuf-nhc-ntc);
}



void Pmla::cutEmittance(double cut, int steps)
{
	char file[100];
	sprintf(file,"cut.%d.agr",  neo);
	FILE * fp = fopen(file, "w");
	if(!fp)
	{
		fprintf(stderr," cant open %s\n", file);
		return;
	}
	for(int s=0; s<= steps; s++)
		fprintf(fp, "@    s%d symbol 9\n@    s%d symbol size 0.210000\n@    s%d line type 0\n", s, s, s);


	double * a[nbuf];
	for(int i=0; i< nbuf; i++) a[i] = dcord[i];
	int good=nbuf;
	sumEmittance(a, good);


	for(int s=0; s< steps; s++)
	{
	   double fx, fxp, fx2, fxp2, fxxp;
	   int ix, ipx;
	   if(s & 1)
	   {
		fx  = 2.*(xm*xpxp - xpm*xxp);
		fxp = 2.*(xx*xpm - xm*xxp);
		fx2 = xpm*xpm - xpxp;
		fxp2 = xm*xm - xx;
		fxxp = 2.*(xxp - xm*xpm);
		ix=0; ipx=1;
	   }  else {
		fx  = 2.*(ym*ypyp - ypm*yyp);
		fxp = 2.*(yy*ypm - ym*yyp);
		fx2 = ypm*ypm - ypyp;
		fxp2 = ym*ym - yy;
		fxxp = 2.*(yyp - ym*ypm);
		ix=2; ipx=3;
	   } 

// 		fprintf(stderr, "fxfacts %e %e %e %e %e\n", fx, fxp, fx2, fxp2, fxxp);
// 		fprintf(stderr, "fyfacts %e %e %e %e %e\n", fy, fyp, fy2, fyp2, fyyp);

        	for(int i = 0; i < good; i++)
        	{
         		double xt =dcord[i][ix];
	    		double xpt=dcord[i][ipx];
// 	    		double yt =dcord[i][2];
// 	    		double ypt=dcord[i][3];
			double enx = fx *xt + fxp * xpt + fxxp * xt * xpt + fx2 * xt * xt + fxp2 * xpt * xpt;
// 			double eny = fy *yt + fyp * ypt + fyyp * yt * ypt + fy2 * yt * yt + fyp2 * ypt * ypt;
// 			dcord[i][6] = -enx - eny;
			dcord[i][6] = -enx ;
		}
	
		quicksort(a,0, good-1);
		int goody = int(good*cut);
        	for(int i = goody; i < good; i++)
        	{
         		double xt =dcord[i][0];
	    		double xpt=dcord[i][1];
			fprintf(fp, "%e %e\n", a[i][0], a[i][1]);
		}
		fprintf(fp, "&\n");
		good = goody;

		sumEmittance(a, good);
	}


        for(int i = 0; i < good; i++)
        {
         	double xt =a[i][0];
	    	double xpt=a[i][1];
		fprintf(fp, "%e %e\n", xt, xpt);
	}
	fclose(fp);

	
	sprintf(file,"act.%d.agr",  neo);
	fp = fopen(file, "w");
	if(!fp)
	{
		fprintf(stderr," cant open %s\n", file);
		return;
	}
	double gamma= dpm/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);

        double tepsx = sqrt(xx*xpxp -xxp*xxp - xx*xpm*xpm -xpxp*xm*xm +2.*xm *xpm*xxp);
        double tbetx =   (xx-xm*xm)/tepsx;
        double talfx = -(xxp-xm*xpm)/tepsx;
	double tgamx = (1. + talfx*talfx)/tbetx;

        double tepsy = sqrt(yy*ypyp -yyp*yyp - yy*ypm*ypm -ypyp*ym*ym +2.*ym *ypm*yyp);
        double tbety =   (yy-ym*ym)/tepsy;
        double talfy = -(yyp-ym*ypm)/tepsy;
	double tgamy = (1. + talfy*talfy)/tbety;

	double actmax=0.;
        for(int i = 0; i < nbuf; i++)
        {
         	double xt =a[i][0];
	    	double xpt=a[i][1];
		double act = betagamma*(tgamx * xt*xt + 2.*talfx* xt*xpt + tbetx*xpt*xpt);
		if(act > actmax) actmax = act;
	}
	int bins[100]; 
        for(int i = 0; i < 100; i++) bins[i]=0;
        for(int i = 0; i < nbuf; i++)
        {
         	double xt =a[i][0];
	    	double xpt=a[i][1];
		double act = betagamma*(tgamx * xt*xt + 2.*talfx* xt*xpt + tbetx*xpt*xpt);
		int ix = int(act/actmax*99.);
		if(ix < 0 || ix >= 100)
		printf("ix: %d act %e actmax %e\n", ix, act, actmax);
		for(int ixx=0; ixx<=ix; ixx++)  bins[ixx]++;
	}
	printf("actmax= %e\n", actmax);
	// this plots the particle density in action space
        for(int i = 0; i < 100; i++) fprintf(fp, "%e %e\n", i*actmax/99., double(bins[i])/nbuf);
	printf("actmax= %e\n", actmax);
	fclose(fp);


	double wbx0 = sqrt(tbetx);
	double wbx1 = sqrt(4e-6/tepsx);
	double ax0  = talfx;
	double ax1  = 0.;;

	double wby0 = sqrt(tbety);
	double wby1 = sqrt(4e-6/tepsy);
	double ay0  = talfy;
	double ay1  = 0.;;

	double phix=.0;
	double phiy=0.;
	double cox  = cos(phix);
	double six  = sin(phix);
	double coy  = cos(phiy);
	double siy  = sin(phiy);
// 	printf("c s c s: %e %e %e %e\n", cox, six, coy, siy);


	double mmf00 = wbx1/wbx0*(cox+ax0*six);
	double mmf01 = wbx0*wbx1*six;
	double mmf10 = -(cox*ax1-cox*ax0+six+six*ax0*ax1)/(wbx0*wbx1);
	double mmf11 = wbx0/wbx1*(cox-ax1*six);

	double mmf22 = wby1/wby0*(coy+ay0*siy);
	double mmf23 = wby0*wby1*siy;
	double mmf32 = -(coy*ay1-coy*ay0+siy+siy*ay0*ay1)/(wby0*wby1);
	double mmf33 = wby0*(coy-ay1*siy)/wby1;
// 	printf("mmmmmmm: %e %e %e %e\n", mmf10, mmf11, mmf32, mmf33);

	double rmax = 0.;
	double lmax = -1e30;
	double lmin =  1e30;
        for(int i = 0; i < nbuf; i++)
        {
		int f;
        	ford[i][0] =mmf00*a[i][0] + mmf01*a[i][1];
         	ford[i][1] =mmf10*a[i][0] + mmf11*a[i][1];
         	ford[i][2] =mmf22*a[i][2] + mmf23*a[i][3];
         	ford[i][3] =mmf32*a[i][2] + mmf33*a[i][3];
         	ford[i][4] =      a[i][4];
         	ford[i][5] =      a[i][5];
         	ford[i][6] =      a[i][6];
		double r2 = square(ford[i][0]) + square(ford[i][2]);
		if(r2 > rmax) rmax=r2;
		if(ford[i][4] < lmin) lmin=ford[i][4];
		if(ford[i][4] > lmax) lmax=ford[i][4];
	}
	rmax=sqrt(rmax)+1e-6; // make slightly bigger to include all particles in the bins
	double ldif= lmax-lmin+1e-6; // lmax is now the difference
	const int rbins=50;
	const int lbins=50;
	double ttrans[lbins][rbins];
	memset(ttrans, 0, rbins*lbins*sizeof(double));
	double ntrans[lbins][rbins];
	memset(ttrans, 0, rbins*lbins*sizeof(int));
        for(int i = 0; i < nbuf; i++)
        {
		double r = sqrt(square(ford[i][0]) + square(ford[i][2]));
		int ir=int(rbins*r/rmax);
		int il=int(lbins*(ford[i][4]-lmin)/ldif);
		double t = square(ford[i][1]) + square(ford[i][3]);
		ttrans[il][ir] += t;
		ntrans[il][ir] ++;
	}
	sprintf(file,"ttrans.%d.agr",  neo);
	fp = fopen(file, "w");
	if(!fp)
	{
		fprintf(stderr," cant open %s\n", file);
		return;
	}
        for(int ir = 0; ir < rbins; ir++)
        for(int il = 0; il < lbins; il++)
	{
		double t = (ntrans[il][ir] > 0) ? ttrans[il][ir]/ntrans[il][ir] : 0.;
		double r = (ir+0.5)/rbins*rmax;
		double l = (il+0.5)/lbins*ldif;
		fprintf(fp,"%e %e %e\n", l, r, t);
	}
	fclose(fp);
}













int Pmla::scatter(const char * path)
{
	char file[100];
	sprintf(file,"%s.xxp.%d.agr", path, neo);
	FILE * fpx = fopen(file, "w");
	if(!fpx)
	{
		fprintf(stderr," cant open %s\n", file);
		return(-1);
	}
	sprintf(file,"%s.yyp.%d.agr", path, neo);
	FILE * fpy = fopen(file, "w");
	if(!fpy)
	{
		fprintf(stderr," cant open %s\n", file);
		return(-1);
	}
        for(int i = 0; i < nbuf; i++)
	{
		fprintf(fpx,"%e %e\n", dcord[i][0], dcord[i][1]);
		fprintf(fpy,"%e %e\n", dcord[i][2], dcord[i][3]);
	}
	fclose(fpx);
	fclose(fpy);
	return 0;
}


int Pmla::zout(const char * zoutfilename, const char * linefilename, double size)
{
	char line[300];
	FILE * fp = fopen(zoutfilename, "r");
	if(!fp)
	{
		fprintf(stderr," cant open OUTPAR.TXT\n");
		return(-1);
	}
	FILE * fo = fopen(linefilename, "w");
	if(!fo)
	{
		fprintf(stderr," cant open line.agr\n");
		return(-1);
	}
	while( fgets(line, 299, fp))
	{
		if(!strncmp(line, " zout", 5)) goto found;
	}
	fprintf(stderr," ZOUT not found\n");
	return(-1);
found:	fgets(line, 299, fp);
	fprintf(fo, "0 0\n");

	int startGun=0;
	int r=1;
	while( fgets(line, 299, fp))
	{
		int i; float z1,z2, dz, B; char typ[20];
		char * li=line; while(*li) { if(*li >= 'A' && *li <= 'Z') *li += 32; li++;}
		int rea= sscanf(line,"%d%f%s%f%f%f", &i, &z1, typ, &z2, &dz, &B);
		if(i > nplots) break;
		if(rea == 5) B=0.;
		else
		if(rea != 6) break;
		z1 /= 100.;
		z2 /= 100.;
		dz /= 100.;
// 		printf("****%d   %f %f %s\n", i , z1, z2, typ);
		if( i > nel) break;

		free(type[i]); type[i] = strdup(typ);
		if(!strcmp(typ, "rotate")) r= -r;
		rotat[i] = r;
		Bfield[i]=B;
		if(!strcmp(typ, "sbload")) rotat[i]= 0;
	//	printf("rotat %d %d\n", i, rotat[i]);
		if( !strcmp(typ, "drift"))
		{
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		else
		if( !strcmp(typ, "bend"))
		{
			fprintf(fo, "%e %e\n", z1, 5*size);
			fprintf(fo, "%e %e\n", z2, 5*size);
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		else
		if( !strcmp(typ, "quad"))
		{
			fprintf(fo, "%e %e\n", z1, 7*size);
			fprintf(fo, "%e %e\n", z2, 7*size);
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		else
		if( !strcmp(typ, "solenoid"))
		{
			fprintf(fo, "%e %e\n", z1, 9*size);
			fprintf(fo, "%e %e\n", z2, 9*size);
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		if( !strcmp(typ, "cell"))
		{
			if(startGun == 0)  startGun=1;
			fprintf(fo, "%e %e\n", z1, 12*size);
			fprintf(fo, "%e %e\n", z2, 12*size);
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		else
		{
			fprintf(fo, "%e %e\n", z2, 0.);
		}
		if( startGun == 1 && strcmp(typ, "cell"))
		{
			gunEnd = i-1;
			startGun= 2;
		}
	}
// 	printf("end of gun index is %d\n", gunEnd);
// 	while( 1)
// 	{
// 		char line2[200];
// 		sscanf(line, "%s", line2);
// 		if(!strcmp(line2, "start")) goto found2;
// 		if(!fgets(line, 299, fp)) break;
// 	}
// 	fprintf(stderr," start keyword not found\n");
// 	return(-1);
// found2:	fgets(line, 299, fp);
	fclose(fp);
	fclose(fo);
	return 0;
}






double Pmla::mag(double p[], int NPARM, int print)
{

	double m4[4][4] = { xx, xxp, xy, xyp, xxp, xpxp, xpy, xpyp, xy, xpy, yy, yyp, xyp, xpyp, yyp, ypyp};

// 	printf("p[%d] = %e %e %e %e %e %e r=%f\n", NPARM, p[0], p[1], p[2], p[3], p[4], p[5], R);
	double wbx0 = sqrt( betx[neo]);
	double wby0 = sqrt( bety[neo]);
	double ax0  = alfx[neo];
	double ay0  = alfy[neo];
	

	
	double cox  = cos(p[0]);
	double six  = -sin(p[0]);
// 	double coy  = cos(p[0]);
// 	double siy  = -sin(p[0]);
	double coy  = 1.;
	double siy  = 0.;


	double wbx1 = 0.05+fabs(p[1]);
	double wby1 = 0.05+fabs(p[2]);
	
	double ax1=0.;
	double ay1=0.;
	if(NPARM == 5)
	{
		ax1  = p[3];
		ay1  = p[4];
	}

	if(NPARM == 7)
	{
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[0][0],  m4[0][1], m4[0][2],  m4[0][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[1][0],  m4[1][1], m4[1][2],  m4[1][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[2][0],  m4[2][1], m4[2][2],  m4[2][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n\n", m4[3][0],  m4[3][1], m4[3][2],  m4[3][3] );
		double a[4][4];
		a[0][0] = (m4[0][0] - m4[2][0] - m4[0][2] + m4[2][2])/2.;
		a[1][0] = (m4[1][0] - m4[3][0] - m4[1][2] + m4[3][2])/2.;
		a[2][0] = (m4[0][0] + m4[2][0] - m4[0][2] - m4[2][2])/2.;
		a[3][0] = (m4[1][0] + m4[3][0] - m4[1][2] - m4[3][2])/2.;

		a[0][1] = (m4[0][1] - m4[2][1] - m4[0][3] + m4[2][3])/2.;
		a[1][1] = (m4[1][1] - m4[3][1] - m4[1][3] + m4[3][3])/2.;
		a[2][1] = (m4[0][1] + m4[2][1] - m4[0][3] - m4[2][3])/2.;
		a[3][1] = (m4[1][1] + m4[3][1] - m4[1][3] - m4[3][3])/2.;

		a[0][2] = (m4[0][0] - m4[2][0] + m4[0][2] - m4[2][2])/2.;
		a[1][2] = (m4[1][0] - m4[3][0] + m4[1][2] - m4[3][2])/2.;
		a[2][2] = (m4[0][0] + m4[2][0] + m4[0][2] + m4[2][2])/2.;
		a[3][2] = (m4[1][0] + m4[3][0] + m4[1][2] + m4[3][2])/2.;

		a[0][3] = (m4[0][1] - m4[2][1] + m4[0][3] - m4[2][3])/2.;
		a[1][3] = (m4[1][1] - m4[3][1] + m4[1][3] - m4[3][3])/2.;
		a[2][3] = (m4[0][1] + m4[2][1] + m4[0][3] + m4[2][3])/2.;
		a[3][3] = (m4[1][1] + m4[3][1] + m4[1][3] + m4[3][3])/2.;



		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", a[0][0],  a[0][1], a[0][2],  a[0][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", a[1][0],  a[1][1], a[1][2],  a[1][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", a[2][0],  a[2][1], a[2][2],  a[2][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n\n", a[3][0],  a[3][1], a[3][2],  a[3][3] );



		m4[0][0] = ( a[0][0] + a[2][0] + a[0][2] + a[2][2])/2.;
		m4[1][0] = ( a[1][0] + a[3][0] + a[1][2] + a[3][2])/2.;
		m4[2][0] = (-a[0][0] + a[2][0] - a[0][2] + a[2][2])/2.;
		m4[3][0] = (-a[1][0] + a[3][0] - a[1][2] + a[3][2])/2.;

		m4[0][1] = ( a[0][1] + a[2][1] + a[0][3] + a[2][3])/2.;
		m4[1][1] = ( a[1][1] + a[3][1] + a[1][3] + a[3][3])/2.;
		m4[2][1] = (-a[0][1] + a[2][1] - a[0][3] + a[2][3])/2.;
		m4[3][1] = (-a[1][1] + a[3][1] - a[1][3] + a[3][3])/2.;

		m4[0][2] = (-a[0][0] - a[2][0] + a[0][2] + a[2][2])/2.;
		m4[1][2] = (-a[1][0] - a[3][0] + a[1][2] + a[3][2])/2.;
		m4[2][2] = ( a[0][0] - a[2][0] - a[0][2] + a[2][2])/2.;
		m4[3][2] = ( a[1][0] - a[3][0] - a[1][2] + a[3][2])/2.;

		m4[0][3] = (-a[0][1] - a[2][1] + a[0][3] + a[2][3])/2.;
		m4[1][3] = (-a[1][1] - a[3][1] + a[1][3] + a[3][3])/2.;
		m4[2][3] = ( a[0][1] - a[2][1] - a[0][3] + a[2][3])/2.;
		m4[3][3] = ( a[1][1] - a[3][1] - a[1][3] + a[3][3])/2.;




		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[0][0],  m4[0][1], m4[0][2],  m4[0][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[1][0],  m4[1][1], m4[1][2],  m4[1][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[2][0],  m4[2][1], m4[2][2],  m4[2][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", m4[3][0],  m4[3][1], m4[3][2],  m4[3][3] );




// 		ax1  = p[3];
// 		ay1  = p[4];
	}

// 	printf("abc %e %e %e %e %e %e %e %e\n", cox, six, coy, siy, wbx1, wby1, ax1, ay1);

	double mmf06 = p[1];
	double mmf16 = p[2];
	double mmf15 = p[3];
	double mmf05 = p[4];
	
	mmf[0][0] = wbx1/wbx0*(cox+ax0*six);
	mmf[0][1] = wbx0*wbx1*six;
	mmf[1][0] = -(cox*ax1-cox*ax0+six+six*ax0*ax1)/(wbx0*wbx1);
	mmf[1][1] = wbx0/wbx1*(cox-ax1*six);

	mmf[2][2] = wby1/wby0*(coy+ay0*siy);
	mmf[2][3] = wby0*wby1*siy;
	mmf[3][2] = -(coy*ay1-coy*ay0+siy+siy*ay0*ay1)/(wby0*wby1);
	mmf[3][3] = wby0*(coy-ay1*siy)/wby1;

	
	

	double aa[4][4];
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	{
		aa[i][j] =0.;
		for(int k=0; k<4; k++) aa[i][j] += mmf[i][k]*m4[k][j];
	}

	double bb[4][4];
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	{
		bb[i][j] =0.;
		for(int k=0; k<4; k++) bb[i][j] += aa[i][k]*mmf[j][k];
	}

	// solenoid strength
	double r=  -15.;
	for(int i=0; i<4; i++)
	{
		bb[1][i] += r*bb[2][i];
		bb[3][i] -= r*bb[0][i];
	}
	for(int i=0; i<4; i++)
	{
		bb[i][1] += r*bb[i][2];
		bb[i][3] -= r*bb[i][0];
	}




// 	double energ=ener[neo];
// 	double gamma= ener[neo]/ erest[0]+1.;
// 	double betagamma = sqrt(gamma*gamma-1.);

	double tt = bb[1][1] + bb[3][3]; 


// 	if(tt < bestgoal)
// 	{
// 		bestgoal = tt;
// 		printf("f  v2 = %e  %e phas/2pi=%e,   betx=%f bety=%f, alfx=%f, alfy=%f\n", tt, bestgoal, p[0]/2./M_PI, p[1], p[2], p[3], p[4]);
// 	}
	if(print)
	{
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", bb[0][0],  bb[0][1], bb[0][2],  bb[0][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", bb[1][0],  bb[1][1], bb[1][2],  bb[1][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", bb[2][0],  bb[2][1], bb[2][2],  bb[2][3] );
		printf("sigma b\t %15.7e   %15.7e  %15.7e  %15.7e\n", bb[3][0],  bb[3][1], bb[3][2],  bb[3][3] );
		printf("mmf    \t %15.7e   %15.7e  %15.7e  %15.7e\n", mmf[0][0],  mmf[0][1], mmf[0][2],  mmf[0][3] );
		printf("mmf    \t %15.7e   %15.7e  %15.7e  %15.7e\n", mmf[1][0],  mmf[1][1], mmf[1][2],  mmf[1][3] );
		printf("mmf    \t %15.7e   %15.7e  %15.7e  %15.7e\n", mmf[2][0],  mmf[2][1], mmf[2][2],  mmf[2][3] );
		printf("mmf    \t %15.7e   %15.7e  %15.7e  %15.7e\n", mmf[3][0],  mmf[3][1], mmf[3][2],  mmf[3][3] );
		printf("here   betx=%f, alfx=%f, bety=%f, alfy=%f\n", betx[neo],alfx[neo],bety[neo],alfy[neo]);
		printf("there  betx=%f, alfx=%f, bety=%f, alfy=%f\n", wbx1*wbx1,ax1,wby1*wby1,ay1);
		printf("delta nu = %f\n", p[0]/2./M_PI);
	}



	return tt;

}






double Pmla::mag1(double p[], int NPARM)
{

	double gamma= ener[neo]/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);
	printf("gamma = %f\n", gamma);
	double t= mag( p, NPARM, 1);
	t = (t)*betagamma*1e6;
	
	printf("mmf=\n");
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
			printf("%f  ", mmf[i][j]);
		printf("\n");
	}
	
	printf("p[%d] = %e %e %e %e %e %e t=%f\n", NPARM, p[0], p[1], p[2], p[3], p[4], p[5], t);

	return t;

}





void Pmla::histogram()
{
	const int nbins=100;
	int hist[nbins];
	for(int j=0; j<6; j++)
	{
		double lmin=0, lmax=0;
		for(int i=0; i<nbuf; i++)
		{
			if(ford[i][j] > lmax) lmax = ford[i][j];
			if(ford[i][j] < lmin) lmin = ford[i][j];
			
		}

		double delta = (lmax - lmin)/nbins;
		for(int s=0; s<=nbins; s++) hist[s]=0;
		for(int i=0; i<nbuf; i++)
		{
			int n=int((nbins-1)*(ford[i][j] - lmin)/(lmax - lmin));
			hist[n]++;
		}
		char file[10]; strcpy(file,"histoA"); file[5] += j;
		FILE * fp=fopen(file,"w");
		if(!fp)
		{
			fprintf(stderr, "cant open %s\n", file);
			return;
		}
		for(int s=0; s<nbins; s++) fprintf(fp, "%f %f\n", lmin + (lmax-lmin)*double(s)/double(nbins), double(hist[s])/nbuf);
		fclose(fp);



	}
}





void Pmla::dump()
{
	char file[100];
	sprintf(file, "dump.%d", neo);
	FILE * fp=fopen(file,"w");
	if(!fp)
	{
		fprintf(stderr, "cant open %s\n", file);
		return;
	}
	sortZ('z');
	fprintf(fp, "# %6d %6d\n", nbuf, neo);
	for(int i=0; i<nbuf; i++)
	{
		fprintf(fp, "%11e %11e %11e %11e %11e %11f %d\n",  
			scord[i][0], scord[i][1], scord[i][2], 
			scord[i][3], scord[i][4], scord[i][5], int(scord[i][6]));
// 		fprintf(fp, "%11e %11e\n",  dcord[i][4], dcord[i][5]);
// 		if(int(dcord[i][6]) == 0)
// 		{
// 			fprintf(stderr, "dump ref part %d\n", i);
// 			fprintf(stderr, "%6d %11e %11e %11e %11e %11e %11f %d\n", i, 
// 				dcord[i][0], dcord[i][1], dcord[i][2], 
// 				dcord[i][3], dcord[i][4], dcord[i][5], int(dcord[i][6]));
// 		}
	}
	fclose(fp);
}




void Pmla::dumpLongit(const char * cwd, double cut)
{
	char file[100];
	sprintf(file, "long.%d.agr", neo);
	FILE * fp=fopen(file,"w");
	if(!fp)
	{
		perror(file);
		return;
	}
	fprintf(fp, "# %6d %6d\n", nbuf, neo);
	fprintf(fp, "@with g0\n");
	fprintf(fp, "@ title \"Longitudinal phase space at %d\"\n", neo);
	fprintf(fp, "@ subtitle \"%s\"\n", cwd);
	fprintf(fp, "@ xaxis label \"Phase [deg]\"\n");
	fprintf(fp, "@ yaxis label \"Energy [MeV]\"\n");
	fprintf(fp, "@  s0 symbol 9\n");
	fprintf(fp, "@  s0 line type 0\n");
	fprintf(fp, "@  s0 symbol size 0.01\n"); 
	fprintf(fp, "@hardcopy device \"PNG\"\n");
	fprintf(fp, "@device \"PNG\" DPI 1200\n");
	sortZ('p');
	int nhc=int(nbuf*cut/2.);
	int ntc=int(nbuf*cut/2.);
// 	sortZ('a');
// 	int nhc=0;
// 	int ntc=int(nbuf*cut);
	for(int i=nhc; i<nbuf-nhc-ntc; i++)
	{
		fprintf(fp, "%11e %11e\n", scord[i][4], scord[i][5]);
// 		fprintf(fp, "%11e %11e\n", scord[i][1], scord[i][3]);
	}
	fclose(fp);
}





void Pmla::impactTdump()
{
	char file[100];
	sprintf(file, "impactTdump.%d", neo);
	FILE * fp=fopen(file,"w");
	if(!fp)
	{
		fprintf(stderr, "cant open %s\n", file);
		return;
	}
	fprintf(fp, "%6d\n", nbuf);
	double lmin=-1e33;
        for(int i = 0; i < nbuf; i++)
        {
            if( dcord[i][4] > lmin) lmin =  dcord[i][4];
	}
	dlm /= nbuf;
	for(int i=0; i<nbuf; i++)
	{
		double gamma= dcord[i][5]/ erest[0]+1.;
		double betagamma = sqrt(gamma*gamma-1.);
		fprintf(fp, "%11e %11e %11e %11e %11e %11f\n",  
			dcord[i][0], dcord[i][1]*betagamma,
			dcord[i][2], dcord[i][3]*betagamma,
			(lmin-dcord[i][4])*wavel/360.*0.01,
			betagamma );
	}
	fclose(fp);
}











int  Pmla::sliceEmit(int print)
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
		sprintf(file,"slice.tips.%04d", neo);
		fpti = fopen(file,"w");
		if(!fpti)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.magnet.%04d", neo);
		fpm = fopen(file,"w");
		if(!fpm)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.emit4d.%04d", neo);
		fpe = fopen(file,"w");
		if(!fpe)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.sigmar.%04d", neo);
		fpr = fopen(file,"w");
		if(!fpr)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.sigmrp.%04d", neo);
		fps = fopen(file,"w");
		if(!fps)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.gamma0.%04d", neo);
		fpg = fopen(file,"w");
		if(!fpg)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.curren.%04d", neo);
		fpc = fopen(file,"w");
		if(!fpc)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.disper.%04d", neo);
		fpd = fopen(file,"w");
		if(!fpd)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.dispep.%04d", neo);
		fpb = fopen(file,"w");
		if(!fpb)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.phase.%04d", neo);
		fpp = fopen(file,"w");
		if(!fpp)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.density.%04d", neo);
		fpde = fopen(file,"w");
		if(!fpde)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.phasep.%04d", neo);
		fppp = fopen(file,"w");
		if(!fppp)
		{
			fprintf(stderr,"cant open file %s\n", file);
			return -1;
		}
		sprintf(file,"slice.p100.%04d", neo);
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
	




	for(int ii=0; ii < nbuf; ii++)
		memcpy(ford[ii], dcord[ii], 7*sizeof(double));
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
	if(delta < 1e-13) delta = 360.;
	double * intrval = new double[slices+1];
	int * numpart = new int[slices+1];
	for(int s=0; s<slices; s++)
	{
		intrval[s]=lmin+s*delta;
		numpart[s] = 0;
	}
	intrval[slices]=lmax;
	numpart[slices]=0;


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
			if(print)fprintf(fpde,"%e %d\n", intrval[s] , numpart[s]);
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


	flattness =0.;
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
		flattness += square(ni - double(nbuf)/slices);
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
// 			fprintf(fpde,"%e %e\n", s_time , density);

			fprintf(fpti,"%e %e\n", sigrdc, sigrpdc*betagamma);
		}
	}


	arr = sqrt(arr/slices);
	arprp = sqrt(arprp/slices);
	arot /= nbuf;
// 	printf("nbuf=%d arr=%e arprp=%e \n", nbuf, arr, arprp);

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



int  Pmla::slicePhasePlot(char * file, int rotat)
{

	FILE * fp = fopen(file,"w");
	if(!fp)
	{
		fprintf(stderr,"cant open file %s\n", file);
		return -1;
	}
	for(int s=0; s<slices; s++) 
	{
		fprintf(fp, "@    s%d symbol %d\n",s, s+1);
		fprintf(fp, "@    s%d symbol size 0.160000\n", s);
		fprintf(fp, "@    s%d line type 0\n", s);
	}

	for(int ii=0; ii < nbuf; ii++)
		memcpy(ford[ii], dcord[ii], 7*sizeof(double));
	if(slices <= 0) slices=1;
	double lmin=0, lmax=0;
	for(int i=0; i<nbuf; i++)
	{
		if(ford[i][4] > lmax) lmax = ford[i][5];
		if(ford[i][4] < lmin) lmin = ford[i][5];
		
	}
	double delta = (lmax - lmin)/slices;
	double * intrval = new double[slices+1];
	for(int s=0; s<=slices; s++) intrval[s]=lmin+s*delta;


	double rotation =0.;
	if(rotat)
	{
		for(int i=0; i< nbuf; i++)
		{
			double r2 = ford[i][0]*ford[i][0]+ford[i][2]*ford[i][2];
			if(r2 > 1.e-8)
			{
				rotation += (ford[i][0]*ford[i][3]-ford[i][2]*ford[i][1])/r2;
			}
		}
		rotation /= nbuf;
		for(int i=0; i< nbuf; i++)
		{
			ford[i][1] += rotation*ford[i][2];
			ford[i][3] -= rotation*ford[i][0];
		}
	}




	for(int s=0; s<slices; s++) 
	{
		if(s) fprintf(fp, "&\n");
		for(int i=0; i<nbuf; i++)
		{
			if(ford[i][5] >= intrval[s] && ford[i][5] <= intrval[s+1])
			{
				double r2 = ford[i][0]*ford[i][0]+ford[i][2]*ford[i][2];
				double r = sqrt(r2);
				double rp = (ford[i][0]*ford[i][1]+ford[i][2]*ford[i][3])/r;
// 				fprintf(fp, "%e %e\n", r, rp);
				fprintf(fp, "%e %e\n", ford[i][0], ford[i][1]);
			}
		}
	}
	
	delete [] intrval;
	fclose(fp);
	return 0;

}

int  Pmla::printCord(char * file)
{
	FILE * fp = fopen(file, "w");
	if(!fp) return errno;
	fprintf(fp,"%d %f %e %e\n", nbuf, zloc[neo], pr, wr);
	for(int i=0; i< nbuf; i++)
	{
		for(int j=0; j< 7; j++) fprintf(fp,"%e  ", dcord[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 0;
}



double  Pmla::unrotate(int method, int remove_orbit, double rm)
{
//	method = 0 -> rm = (x*y' - y*x')/r2    (changs method)
//	method = 1 -> rm = (x'2+y'2)/(x*y' - y*x')	(vladimirs method)
//	method = 2 -> rm = argument

	if(remove_orbit)  for(int i=0; i< nbuf; i++)
	{
		ford[i][0] = dcord[i][0] -xm;
		ford[i][1] = dcord[i][1] -xpm;
		ford[i][2] = dcord[i][2] -ym;
		ford[i][3] = dcord[i][3] -ypm;
	}
	else  for(int i=0; i< nbuf; i++)
	{
		ford[i][0] = dcord[i][0];
		ford[i][1] = dcord[i][1];
		ford[i][2] = dcord[i][2];
		ford[i][3] = dcord[i][3];
	}






	double vt = 0.;
	double rotation =0.;
	if(method == 0)
	{
        	for(int i=0; i< nbuf; i++)
        	{
                	double r2 = ford[i][0]*ford[i][0]+ford[i][2]*ford[i][2];
                	if(r2 > 1.e-8)
                	{
				//  this is  dTheta/ds in units 1/m
                        	rotation += (ford[i][0]*ford[i][3]-ford[i][2]*ford[i][1])/r2;
                	}
//                 	vt +=        ford[i][0]*ford[i][0]+ford[i][2]*ford[i][2];
//                        	rotation += ford[i][0]*ford[i][3]-ford[i][2]*ford[i][1];
        	}
        	rm = rotation / nbuf;
//         	rm = rotation / vt;
	}
	else if(method == 1)
	{
		for(int i=0; i< nbuf; i++)
		{
			double x  = ford[i][0];
			double xp = ford[i][1];
			double y  = ford[i][2];
			double yp = ford[i][3];
			vt += xp*xp+yp*yp;
			rotation += x*yp - y*xp;
		}
		rm = vt/rotation;
	}

	rot[neo]=1./rm;
	double rm2 = rm*rm;

	for(int i=0; i< nbuf; i++)
	{
		ford[i][1] += rm*ford[i][2];
		ford[i][3] -= rm*ford[i][0];
	}
	printf("unrotate(%d,%d):  1/thetaPrime = %f\n", method, remove_orbit, rot[neo]);
	double xx=0;
	double yy=0;
	double xpxp=0;
	double ypyp=0;
	double xxp=0;
	double yyp=0;
        for(int i = 1; i < nbuf; i++)
        {
		double x  = ford[i][0];
		double xp = ford[i][1];
		double y  = ford[i][2];
		double yp = ford[i][3];
		xx   += x*x;
		yy   += y*y;
		xxp  += x*xp;
		yyp  += y*yp;
		xpxp += xp*xp;
		ypyp += yp*yp;
		double t = sqrt((xp*xp+yp*yp)/rm2);
		double r = sqrt(x*x+y*y);
// 		if(method ==2) printf("%e %e\n", r, t);

	}
	xx /= nbuf;
	yy /= nbuf;
	xxp /= nbuf;
	yyp /= nbuf;
	xpxp /= nbuf;
	ypyp /= nbuf;
        epsx_u[neo] = sqrt(xx*xpxp - xxp*xxp);
        epsy_u[neo] = sqrt(yy*ypyp - yyp*yyp);
	double t = xpxp+ypyp;
	double tt = t/(rm*rm);
        double betx =   xx/epsx_u[neo];
        double alfx = -xxp/epsx_u[neo];
        double bety =   yy/epsy_c[neo];
        double alfy = -yyp/epsy_c[neo];
	printf("unrot  betx %f, alfx %f, xx %f\n", betx, alfx, sqrt(xx));
	printf("unrot  bety %f, alfy %f, yy %f\n", bety, alfy, sqrt(yy));
	double gamma= ener[neo]/ erest[0]+1.;
	double betagamma = sqrt(gamma*gamma-1.);
        epsx_u[neo] *= betagamma;
        epsy_u[neo] *= betagamma;
	printf("unrot rl2 =%f %e\n", rm, tt);
	printf("unrot ex  =%f  %e\n", rm, epsx_u[neo]);
	printf("unrot ey  =%f  %e\n", rm, epsy_u[neo]);
	printf("unrot sx  =%f  %e\n", rm, sqrt(xx));
	printf("unrot sy  =%f  %e\n", rm, sqrt(yy));
	return rm;
}



double  Pmla::linMatrix()
{
       	for(int i=0; i<4; i++)
       	{
               	for(int j=0; j<4; j++)
               	{
                       	mmt[i][j] = ( dcord[i][j]+dcord[i+4][j])/mdelta[i]/2.;
                       	printf("%15e  ", mmt[i][j]);
               	}
               	printf("\n");
       	}
       	printf("\n");
	double mmd[4][4];
	memcpy(mmd,mmt,16*sizeof(double));
	double det= dminv((double *) mmd, 4);
	return det;
}



void  Pmla::saveDelta()
{
	if(!getElement(1)) twiss(0);
        for(int j=0; j<4; j++)
        {
		mdelta[j] = dcord[j][j];
        }
	printf("deltas = %f %f %f %f\n", mdelta[0], mdelta[1], mdelta[2], mdelta[3]);

}




void  Pmla::particle( int & i, double R, double x, double xn, double y, double yn)
{
	double xp = -0.5*R*y + xn;
	double yp =  0.5*R*x + yn;
	
	dcord[i][0] =   mmt[0][0] * x + mmt[0][1] * xp;
	dcord[i][1] =   mmt[1][0] * x + mmt[1][1] * xp;
	dcord[i][2] =   mmt[2][2] * y + mmt[2][3] * yp;
	dcord[i][3] =   mmt[3][2] * y + mmt[3][3] * yp;
	dcord[i][4] =  0.;
	dcord[i][5] =  2.;
	dcord[i][6] =  0.;
// 	printf("%f %f %f %f %f %f %f\n", dcord[i][0], dcord[i][1], dcord[i][2], dcord[i][3], dcord[i][4], dcord[i][5], dcord[i][6]);
			
	ford[i][0] =    x;
	ford[i][1] =    xn;
	ford[i][2] =    y;
	ford[i][3] =    yn;
	ford[i][4] =  0.;
	ford[i][5] =  2.;
	ford[i][6] =  0.;
	i++;
}
			

void  Pmla::testcord(int ne, double m11, double m12, double m21)
{
	srandom(76873215);
	neo = ne;
	double R = 0.5;

	double m22 =  ( 1.+ m21*m12)/m11;
	for(int i = 0; i<4; i++)
	for(int j = 0; j<4; j++) mmt[i][j]=0.;
	mmt[0][0] = mmt[2][2] = m11;
	mmt[0][1] = mmt[2][3] = m12;
	mmt[1][0] = mmt[3][2] = m21;
	mmt[1][1] = mmt[3][3] = m22;


// 	nbuf=npoints=4;
// 	for(int i=0; i<nbuf; i++)
// 	{
// 		double cc[6];
// 		for(int j=0; j <6; j++) cc[j]=0.;
// 		cc[i]=1.;
// 		cc[1] += -0.5*R*cc[2];
// 		cc[3] +=  0.5*R*cc[0];
// 
// 		dcord[i][0] =   m11 * cc[0] + m12 * cc[1];
// 		dcord[i][1] =   m21 * cc[0] + m22 * cc[1];
// 		dcord[i][2] =   m11 * cc[2] + m12 * cc[3];
// 		dcord[i][3] =   m21 * cc[2] + m22 * cc[3];
// 		dcord[i][4] =  0.;
// 		dcord[i][5] =  2.;
// 		dcord[i][6] =  0.;
// 
// 	}
// 	return;



	nbuf=npoints=4096;
	for(int i=0; i<nbuf; )
	{
    		double v    =  2.*M_PI*double(random())/(double(RAND_MAX)+1.);
        	double w    =  -2.0*log(sqrt(0.25*double(random())/(double(RAND_MAX)+1.)));
    		double vp   =  2.*M_PI*double(random())/(double(RAND_MAX)+1.);
        	double wp   =  -2.0*log(sqrt(0.25*double(random())/(double(RAND_MAX)+1.)));

		double x =  w*cos(v);
		double y =  w*sin(v);
		

		double noise=0.01;
		double xn =  noise*wp*cos(vp);
		double yn =  noise*wp*sin(vp);


		particle( i, R,  x,  xn,  y,  yn);
		particle( i, R, -x,  xn,  y,  yn);
		particle( i, R,  x,  xn, -y,  yn);
		particle( i, R, -x,  xn, -y,  yn);

		particle( i, R,  y,  xn,  x,  yn);
		particle( i, R, -y,  xn,  x,  yn);
		particle( i, R,  y,  xn, -x,  yn);
		particle( i, R, -y,  xn, -x,  yn);

		particle( i, R,  x,  -xn,  y,  -yn);
		particle( i, R, -x,  -xn,  y,  -yn);
		particle( i, R,  x,  -xn, -y,  -yn);
		particle( i, R, -x,  -xn, -y,  -yn);

		particle( i, R,  y,  -xn,  x,  -yn);
		particle( i, R, -y,  -xn,  x,  -yn);
		particle( i, R,  y,  -xn, -x,  -yn);
		particle( i, R, -y,  -xn, -x,  -yn);



			
	}

	double xx=0.;
	double yy=0.;
	double rot =0.;
	for(int i=0; i<nbuf; i++)
	{
		rot += dcord[i][0]*dcord[i][3] -  dcord[i][2]*dcord[i][1];
		yy += dcord[i][0]*dcord[i][0];
		xx += dcord[i][2]*dcord[i][2];
	}
	rot /= nbuf;
	printf("rot = %e\n", rot);
	xx /= nbuf;
	yy /= nbuf;
	rot /= (xx+yy);
	double flow = xx*yy*R;
	printf("made  xx = %e  yy = %e flow = %e R=%f rot = %e\n",  xx, yy, flow,R, rot);
}
