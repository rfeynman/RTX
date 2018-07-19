#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>



main(int argc, char ** argv)
{
	DIR * d = opendir(".");
	if(!d)
	{
		fprintf(stderr, "Can't open '.'\n");
		exit(-1);
	}

	while( struct dirent * dd = readdir(d))
	{
		printf("%s\n", dd->d_name);
	}


	closedir(d);
}
