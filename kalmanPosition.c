/*  kalmanPosition.c

    This file contains the code for a kalman filter.

    J. Watlington, 11/15/95 
	
	modified to 1_D space
	ZY. Zhang 21/04/16 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kalmanPosition.h"

//#define PRINT_DEBUG (0)


/*  The following are the global variables of the Kalman filters,
    used to point to data structures used throughout.     */

static double  *state_pre;         /* ptr to apriori state, x(-)     */
static double  *state_post;        /* ptr to aposteriori state , x(+) */

static double  *cov_pre;          /* ptr to apriori covariance , P(-) */
static double  *cov_post;         /* ptr to apriori covariance , P(-) */
static double  *sys_noise_cov;    /* system noise covariance  (GQGt)  */
static double  *mea_noise_cov;    /* measurement noise variance vector (R)  */

static double  *sys_transfer;     /* system transfer function (Phi)    */
static double  *mea_transfer;     /* measurement transfer function (H) */

static double  *kalman_gain;      /* The Kalman Gain  (K) */

int            global_step = 0;    /* the current step number (k) */
//int            measurement_size;   /* number of elems in measurement */
//int            state_size;         /* number of elements in state    */

/*  Temporary variables, declared statically to avoid lots of run-time
    memory allocation.      */

static double  *z_estimate;        /* a measurement_size x 1 vector */
//static double  *temp_state_state; /* a state_size x state_size  */
//static double  *temp_meas_state;  /* a measurement_size x state_size  */
//static double  *temp_meas_meas;   /* a measurement_size squared  */
//static double  *temp_meas_2;      /* another one ! */

/*  Prototypes of internal functions  */

static void update_system( double *z, double *x_minus,
			 double *kalman_gain, double *x_plus );
static void estimate_prob( double *P_post, double *Phi, double *GQGt,
			  double *P_pre );
static void update_prob(double *P_pre, double *R, double *H,
	double *P_post, double *K);

/****************  System Model Manifestations   *******************/
static void        apply_system(double *old_state, double *u, double *new_state);

/****************  Measurement Model Manifestations   **************/
static void        apply_measurement( double *new_state, double *est_measurement );



/******************************************************************

  Linear Kalman Filtering

  kalman_init()
  This function initializes the kalman filter.  Note that for a
  straight-forward (linear) Kalman filter, this is the only place that
  K and P are computed...      */

void kalman_init( double *GQGt, double *Phi, double *H, double *R,
		double *P, double *x)
{
  /*  Init the global variables using the arguments.  */

  //vec_copy( x, state_post, state_size );
  state_post = x;
  //mat_copy( P, cov_post, state_size, state_size );
  cov_post = P;

  sys_noise_cov = GQGt;
  mea_noise_cov = R;

  sys_transfer = Phi;
  mea_transfer = H;

  cov_pre = malloc( sizeof( double ) );
  state_pre = malloc( sizeof( double ) );
  kalman_gain = malloc( sizeof( double ) );
  z_estimate = malloc( sizeof( double ) );
  
  
}


/*  kalman_step()
    This function takes a set of measurements, and performs a single
    recursion of the straight-forward kalman filter.
*/

void kalman_step( double *u_in, double *z_in )
{
  /**************  Estimation Loop  ***************/

  estimate_prob( cov_post, sys_transfer, sys_noise_cov, cov_pre );
  apply_system( state_post, u_in, state_pre );
  
  update_prob( cov_pre, mea_noise_cov, mea_transfer, cov_post, kalman_gain );
  update_system( z_in, state_pre, kalman_gain, state_post );

  global_step++;

#ifdef PRINT_DEBUG
  printf("u = %f, z = %f, x = %f, step = %d\n", u_in, z_in, state_post, global_step);
#endif
}

/*  kalman_get_state
    This function returns a pointer to the current estimate (a posteriori)
    of the system state.         */

double *kalman_get_state( void )
{
  return( state_post );
}

/************************************************************

   Internal Functions, defined in order of appearance

*/

/* update_system()
   This function generates an updated version of the state estimate,
   based on what we know about the measurement system.       */
   
static void update_system( double *z, double *x_pre,
			 double *K, double *x_post )
{
#ifdef PRINT_DEBUG
  printf( "ekf: updating system\n" );
#endif

  apply_measurement( x_pre, z_estimate );
  //vec_sub( z, z_estimate, z_estimate, measurement_size );
  *z_estimate = ( *z ) - ( *z_estimate );
  //mat_mult_vector( K, z_estimate, x_post, state_size, measurement_size );
  *x_post = ( *K ) * ( *z_estimate );
  //vec_add( x_post, x_pre, x_post, state_size );
  *x_post = ( *x_post ) + ( *x_pre );
}

/* estimate_prob()
   This function estimates the change in the variance of the state
   variables, given the system transfer function.   */

static void estimate_prob( double *P_post, double *Phi, double *GQGt,
			  double *P_pre )
{
#ifdef PRINT_DEBUG
  printf( "ekf: estimating prob\n" );
#endif
  
  *P_pre = ( *Phi ) * ( *P_post ) * ( *Phi ) + ( *GQGt );
}


/*  update_prob()
    This function updates the state variable variances.
    Inputs:
    P_pre - the apriori probability matrix ( state x state )
    R     - measurement noise covariance   ( meas x meas )
    H     - the measurement transfer matrix ( meas x state )
    Outputs:
    P_post - the aposteriori probability matrix (state x state )
    K     - the Kalman gain matrix ( state x meas )
*/

static void update_prob( double *P_pre, double *R, double *H,
			double *P_post, double *K )
{
#ifdef PRINT_DEBUG
  printf( "ekf: updating prob\n" );
#endif
#ifdef DIV_DEBUG
  printf( "P = %f\n", P_pre );
#endif
  /*mat_mult( H, P_pre, temp_meas_state,
	   measurement_size, state_size, state_size );
  mat_mult_transpose( H, temp_meas_state, temp_meas_meas,
		     measurement_size, state_size, measurement_size );
  mat_add( temp_meas_meas, R, temp_meas_meas,
	  measurement_size, measurement_size );

  take_inverse( temp_meas_meas, temp_meas_2, measurement_size );
  */
/*  
#ifdef DIV_DEBUG
  printf( "1 / (HPH + R)", temp_meas_2,
	       measurement_size, measurement_size );
#endif
  mat_transpose_mult( temp_meas_state, temp_meas_2, K,
		     state_size, measurement_size, measurement_size );
*/
			 
  *K = (*P_pre) * (*H) / ( (*H) * (*P_pre) * (*H) + ( *R ) );

/*
  print_matrix( "Kalman Gain", K, state_size, measurement_size );

  mat_mult( K, temp_meas_state, temp_state_state,
	   state_size, measurement_size, state_size );
#ifdef PRINT_DEBUG
  printf( "ekf: updating prob 3\n" );
#endif
  mat_add( temp_state_state, P_pre, P_post, state_size, state_size );
*/
  
  *P_post = ( 1.0 - ( *K ) * ( *H )) * ( *P_pre );
}

/****************  System Model Manifestations   **************
  
  apply_system()
  This function is called to predict the next state given a
  current state.  It is an evaluation of the system's transfer
  function.  If used with a Linear Kalman filter, this should be a
  linear function, otherwise it may be non-linear.     */
  
void apply_system( double *old_state, double *u, double *new_state )
{
  *new_state = ( *sys_transfer ) * ( *old_state ) + ( *u ) ;
}

  /****************  Measurement Model Manifestations   **************
  
  apply_measurement()
  This function is called to predict the next measurement given a
  predicted state.  It is an evaluation of the measurement's transfer
  function.*/
void apply_measurement( double *new_state, double *est_measurement )
{  
  *est_measurement = ( *mea_transfer ) * ( *new_state );
}
