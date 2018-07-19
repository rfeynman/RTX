// calculate phase shilft of cavities when they are moved longitudinally
// jorg kewisch 2006


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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

const double clight= 2.99792456210e10; /// in cm/sec

main(int argc, char ** argv)
{
        double ener_b, freq, distance, phas;
        printf("energy [MeV] ? "); scanf("%lf", &ener_b);
        printf("frequency [MHz] ? "); scanf("%lf", &freq);
        printf("distance [cm] ? "); scanf("%lf", &distance);
        double gamma=ener_b / 0.51104  +1.;
        double beta = sqrt(1.-1./(gamma*gamma));
        double phFac= -360.*freq*1e6/(clight*beta);
	double add= phFac*distance;
        printf("beta  e %e gamma %e beta %e freq %e phFac %e\n", ener_b, gamma, beta, freq, phFac);
	while(1)
	{
		printf("old phase ? "); 
		int err=scanf("%lf", &phas);
		if(err != 1)
		{
			char line[1000];
			scanf("%s", line);
			continue;
		}
		double phas1 = phas + add;
		double phas3 = phas + 3.*add;
		while(phas1 <    0.) phas1 += 360.;
		while(phas1 >= 360.) phas1 -= 360.;
		while(phas3 <    0.) phas3 += 360.;
		while(phas3 >= 360.) phas3 -= 360.;
		printf(" %20.10e  %20.10e\n", phas1, phas3);
	}

}
