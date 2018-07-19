#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <error.h>


main(int argc, char ** argv)
{
	printf("*%20.10e|\n", 1e-20);
	printf("*%20.10e|\n", -1e-20);
}
