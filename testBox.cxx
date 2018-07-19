#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


const int million=1000000000;
double rfac = 2./RAND_MAX;
main(int argc, char ** argv)
{
	double s=0;
	for(int i=0; i< million; i++)
	{
		 double r = rfac*random()-1.;
		 s += r*r;
	}
	s /= million;
	printf("%e  %e\n", s, sqrt(s));
}
