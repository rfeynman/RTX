


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


#include "xmgr_class.hxx"



Xmgr::Xmgr(const char * file, const char * title, const char * xaxis, const char * yaxis)
{
	fp = fopen(file, "w");
	if(!fp)
	{
		fprintf(stderr, "Cant open %s!", file);
		exit(-1);
	}
	fprintf(fp,"@title \"%s\"\n", title);
	fprintf(fp,"@frame linewidth 2\n");
	fprintf(fp,"@xaxis label \"%s\"\n", xaxis);
	fprintf(fp,"@yaxis label \"%s\"\n", yaxis);
	fprintf(fp,"@map color 15 to (255, 255, 0), \"yellow\"\n");
	fprintf(fp,"@map color 5 to (0, 139, 0), \"green4\"\n");

	first=1;
}


Xmgr::~Xmgr()
{
	fprintf(fp, "@hardcopy device \"PNG\"\n");
	fprintf(fp, "@device \"PNG\" DPI 1200\n");

	fclose(fp);
}


void Xmgr::subTitle(const char * s)
{
	fprintf(fp, "@with g0\n");
	fprintf(fp, "@subtitle \"%s\"\n", s);
}

void Xmgr::frame(double xmin, double xmax, double ymin, double ymax)
{
	fprintf(fp,"@g0 on\n");
	fprintf(fp,"@with g0\n");
	fprintf(fp,"@    world xmin %f\n", xmin);
	fprintf(fp,"@    world xmax %f\n", xmax);
	fprintf(fp,"@    world ymin %f\n", ymin);
	fprintf(fp,"@    world ymax %f\n", ymax);
}

void Xmgr::plot(double * x, double * y, int * skip, int nelem, int yinc)
{
	if(!first) fprintf(fp,"&\n");
	first=0;
	int * sk = skip;
	int ski;

	for(int i=0; i< nelem; i++)
	{
		ski=0;
		if(sk) 
		{
			int ski = *skip;
		}
		if(!ski) fprintf(fp,"%e %e\n", *x, *y);
		x++;
		y += yinc;
		skip++;
	}
}

void Xmgr::plot(float * x, double * y, int * skip, int nelem, int yinc)
{
	if(!first) fprintf(fp,"&\n");
	first=0;
	int * sk = skip;
	int ski;

	for(int i=0; i< nelem; i++)
	{
		ski=0;
		if(sk) 
		{
			int ski = *skip;
		}
		if(!ski) fprintf(fp,"%e %e\n", *x, *y);
		x++;
		y += yinc;
		skip++;
	}
}

