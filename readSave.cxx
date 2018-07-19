#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


//	read a parmela save file. The file is a fortran unformated file, i.e. each record is preceeded and followed by the length of the record in bytes
//	the file containes only 1 record with the length 16+54*npart, where npart is the number of particles.
//	the 16 bytes are the reference phase (double), the 

main(int argc, char ** argv)
{
	int fi = open("SAVECORa",O_RDONLY);
	if(fi < 0)
	{
		perror("SAVECORa");
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
	char * p = buff+20;
	char buff2[54];
	for(int ip = 0; ip < ngood; ip++)
	{
		memcpy(buff2, p, 54);
		p += 54;
		double cord[6];
		memcpy(cord, buff2, 6*8);
		printf("%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e  ", cord[0], cord[1], cord[2], cord[3], cord[4], cord[5]);
		int n1; memcpy(&n1, buff2+48, 4); printf("%10d ", n1);
		char a = buff2[52]; printf("%2d ", a);
		char b = buff2[52]; printf("%2d ", b);
		printf("\n");

	}
	exit(0);



}
