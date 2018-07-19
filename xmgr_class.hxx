#include <stdio.h>
#include <stdlib.h>

class Xmgr
{
	FILE * fp;
	int first;
public:
	Xmgr(const char * file, const char * title, const char * xaxis, const char * yaxis);
	~Xmgr();
	void plot( double * x, double * y, int nelem, int yinc=1);
	void plot(double * x, double * y, int * skip, int nelem, int yinc=1);
	void plot(float * x, double * y, int * skip, int nelem, int yinc=1);
	void frame(double xmin, double xmax, double ymin, double ymax);
	void subTitle(const char * s);

};

