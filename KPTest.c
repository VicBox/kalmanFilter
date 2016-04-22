#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kalmanPosition.h"


FILE *inputfile, *measurementfile, *outputfile;


int main()
{
	inputfile = fopen("./u_in", "r");
	measurementfile = fopen("./Ground_y", "r");

	if (inputfile == NULL || measurementfile == NULL)
	{
		printf("*** data could not be opened!\n");
		return -1;
	}
	else
	{
		printf("\n\nfile open successful\n");	
	}

	outputfile = fopen("./output", "w+");
	if (outputfile != NULL)
	{
		printf("outputfile created!\n");
	}
	else
	{
		printf("create file error\n");
		return -1;	
	}

	double cam_pos, u_in = 0.0, x;
	double GQGt  = 0.01, Phi = 1, H = 1, R = 16, P = 0.25;

	fscanf( measurementfile, "%lf", &cam_pos);

	x = cam_pos;

	kalman_init( &GQGt, &Phi, &H, &R, &P, &x);

	while( 1 )
	{
		if(fscanf(measurementfile,"%lf", &cam_pos ) == -1 || fscanf( inputfile, "%lf", &u_in ) == -1)
		{
			break;	
		}

		kalman_step( &u_in, &cam_pos);
		fprintf( outputfile, "%f\n", *kalman_get_state() );
	}
	printf("\nCalculattion success\n");
	getchar();
	return 0;
}
