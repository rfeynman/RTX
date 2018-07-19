#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <error.h>


main(int argc, char ** argv)
{
	FILE * fp = fopen(argv[1], "r");
	if(!fp)
	{
		perror(argv[1]);
		exit(-1);
	}

	FILE * fo = fopen("_traj.agr", "w");
	if(!fo)
	{
		perror("_traj.agr");
		exit(-1);
	}

	char line[1000];
	while( fgets(line, 1000, fp)  )
	{
		fgets(line, 1000, fp);
		fgets(line, 1000, fp);

		while( line[3] != ' ')
		{
			double x, y, z, t;
			sscanf(line, "%lf%lf%lf%lf", &x, &y, &z, &t);
			fgets(line, 1000, fp);
			fprintf(fo, "%e %ei\n", z, x);
		}
		fprintf(fo, "&\n");
	}




}
