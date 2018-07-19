#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


main(int argc, char ** argv)
{
	FILE * fp = fopen("OUTPAR.TXT", "r");
	if(!fp)
	{
		perror("OUTPAR.TXT");
		exit(-1);
	}

	char line[1000];

	while( fgets(line, 1000, fp) )
	{
		char start[20];
		sscanf(line, "%19s", start);
		if(!strcmp(start, "start")) goto found_start;
	}
	fprintf(stderr, "no start found in OUTPAR\n");
	exit(-1);

found_start:
	while( fgets(line, 1000, fp) )
	{
		int element;
		int num = sscanf(line, "%d", &element);
		if(num && (element > 0)) printf("%d ", element);
	}
	printf("\n");
}


