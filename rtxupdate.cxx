#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


main(int argc, char ** argv)
{
	for(int ff = 1; ff < argc; ff++)
	{
		char filename[200];
		strcpy(filename, argv[ff]);
		int lf = strlen(filename);
		if( strcmp(filename+lf-4, ".inp") ) continue;
		char filename2[200];
		strcpy(filename2, argv[ff]);
		strcpy(filename2+lf-4, ".rtxUpdate");
		char cmd[200], line[1000];
		sprintf(cmd, "mv %s %s", filename, filename2);
		system(cmd);

		FILE * fi = fopen(filename2, "r");
		if(!fi)
		{
			perror(filename2);
			exit(-1);
		}
		FILE * fo = fopen(filename, "w");
		if(!fo)
		{
			perror(filename);
			exit(-1);
		}

		while( fgets(line, 1000, fi) )
		{
			if( ! strncmp(line, "!@var", 5)  )
			{
				int sofar, active, num;
				char var[10];
				sscanf(line, "%s%d%d%n", var, &num, &active, &sofar);
				fprintf(fo, "!@var %1d  %-3d   %s", active, num, line+sofar);
			}
			else
			{
				fprintf(fo, "%s", line);
			}
		}
		fclose(fi);
		fclose(fo);
		sprintf(cmd, "touch -r  %s %s", filename2, filename);
		system(cmd);
	}




}
