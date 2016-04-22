/*  kalmanPosition.h
    This is the include file for the kalman filtering functions defined
    in kalman.c.

    J. Watlington, 11/16/95
	
	modified by 
	ZY. Zhang 21/04/16
*/

/***********    Linear Kalman Filtering   *************/

extern void kalman_init( double *GQGt, double *Phi, double *H, double *R,
		double *P, double *x );

extern void kalman_step( double *u_in, double *z_in );

extern double *kalman_get_state( void ); 