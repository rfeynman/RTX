#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


//	read a parmela save file and  set the avarage orbit and dispersion to zero.
//	The file is a fortran unformated file, i.e. each record is preceeded and followed by the length of the record in bytes
//	the file containes only 1 record with the length 16+54*npart, where npart is the number of particles.
//	the 16 bytes are the reference phase (double), the number of good (not lost) particles "ngood" and the number of particles at the start "npart".

main(int argc, char ** argv)
{
	int fi = open(argv[1],O_RDONLY);
	if(fi < 0)
	{
		perror(argv[1]);
		exit(0);
	}
	off_t size = lseek(fi, 0, SEEK_END);
	lseek(fi, 0, SEEK_SET);
	printf("size is %ld\n", size);
	char * buff = new char[size];
	ssize_t rsize = read(fi, buff, size);
	printf("read %ld\n", rsize);
	close(fi);

	double refPhase; memcpy(&refPhase, buff+4, 8);
	int ngood; memcpy(&ngood, buff+12, 4);
	int npart; memcpy(&npart, buff+16, 4);
	printf("refPhase = %e np = %d, ngood=%d\n", refPhase, npart, ngood);
	char buff2[54];
	char * pp;

	double xm, xpm, ym, ypm, dm;
	double xd, xpd, yd, ypd, dd;
	double dx, dpx, dy, dpy;


	pp = buff+20;
	double * p = (double *) pp;
	p[0]=0.;
	p[1]=0.;
	p[2]=0.;
	p[3]=0.;

	int fo = open(argv[1],O_WRONLY);
	if(fo < 0)
	{
		perror("argv[1]");
		exit(0);
	}
	ssize_t wsize = write(fo, buff, size);
	close(fo);


	exit(0);



}
