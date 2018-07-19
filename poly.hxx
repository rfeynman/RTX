//	fit of a polynom of order n to m data points (x[],y[])
//	coeff is a m long array of double pointers pointing to 
//	m+1 long arrays of doubles.
//	a is a n long array of doubles for the result.

int polyInit(double * x, int m, double * coef[]);

double polyFit(double * x, double * y, int m, double * coef[], double * a, int n);

int polyTrack(double tolerance, int maxRows, double * x, double * y, double * yp,
	double * a, int n, double x0, double xf, int nrows = 0, int first=0);

int polyTrackFill(double ymax, int maxRows, double * xr, double * yr, double * yp,
	double * a, int n, double x0, double xf, int nrows = 0, int first=0);

double polyValue(double x, double * a, int n);
double polyPrime(double x, double * a, int n);
double polyDoublePrime(double x, double * a, int n);

