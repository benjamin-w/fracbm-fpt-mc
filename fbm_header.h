/* fracbm-fpt-mc (2019)
 * Authors: Benjamin Walter, Kay Wiese
 * Header File */

// LIBRARIES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <fftw3.h>
#include <lapacke.h>
#include <cblas.h>

// MACROS
#define IJ2K(a,b) (a+b*(b+1)/2) // Converts matrix indices
#define ARRAY_REALLOC_FACTOR 2.0 // Factor for realloc
#define MIN(a,b) ( (a < b) ? (a) : (b))
#define MAX(a,b) ( (a > b) ? (a) : (b))
#define ABS(a) ((a > 0) ? (a): (-a) )
#define ALLOC(p,n)  (p)=malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define FFT_ALLOC(p,n)  (p)=fftw_malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define REALLOC(p,n)  (p)=realloc( (p) , (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 

// STRUCT
// Complex numbers
typedef struct complex_z
{
	 double r;
	 double i;
} complex_z;

typedef struct bridge_process // This is the bridge process spanned by two midpoints. In its center, a midpoint will be generated
{
	double right_time;	// Time of right endpoint
	double right_value;	// Height of right endpoint
	double left_time;	// Time of left endpoint
	double left_value;	// Height of left endpoint
	int generation;		// Generation
	double threshold;	// They know where the upper threshold is
	double critical_strip;	// And the critical strip;
	int centre_critical;	// 1 - Center is in critical strip -> Sub-divide; 0 - Off the critical strip -> Stop
	int crossing_threshold; // 1 - The connection of both endpoints crosses the threshold
	struct bridge_process * left_sub_bridge;
	struct bridge_process * right_sub_bridge;
	struct bridge_process * parental_bridge;
	struct bridge_process * root_bridge;
	struct bridge_process * previous_root;
	struct bridge_process ** pointer_stack; // Only root node has this. List of all bridges attached to the root. Free every pointer in this list.
	int stack_top;
} bridge_process;

typedef struct triag_matrix
{
	/* This struct stores a sequence of points X_1, X_2, ... and their inverse correlation matrix. It can be dynamically managed as points are added. Convention for triagonal matrix: column major form and 'upper' triagonal form.*/
	long size; // This is the size of the matrix itself, that is 0, 1, ..., size - 1 are array indices for both vectors 'trajectory_x' and 'trajectory_t', and 0, 1, ..., size*(size+1)/2 - 1 are indices for upper triagonal matrix. Note that size=N, meanwhile fracbm has N+1 entries. So size is the number of 'free' points, the first one being fixed.
	long array_length; // This is the length of the array. If size gets to large, realloc. If size == array_length, this array is full.
	double * inv_corr_matrix; // ( size * ( size + 1) / 2) entries in symmetric inverse correlation matrix of all points already known
	double * trajectory_x; // All points of the trajectory already known. Length = size + 1 (X_0 = 0 doesn't count, and is neglected in inverse correlation matrix (null mode)). 
	double * trajectory_t; // And the corresponding time points. Length = size + 1
	double hurst_parameter;
} triag_matrix;

// FUNCTIONS
void initialise( fftw_complex** ,  fftw_complex** , fftw_complex** ,  fftw_complex** ,  double** , complex_z** , long N, gsl_rng**, const gsl_rng_type**, int);
void initialise_trajectory ( double**, long);
void initialise_inverse_correlation_matrix(double***, long );
double erfcinv(double);
void write_correlation_exponents(double*, long, double, double);
void write_correlation(fftw_complex*, double*, long);
void write_inverse_correlation_matrix(double **, long, double);
void generate_random_vector(complex_z*, fftw_complex*, fftw_complex*,long, double);
void set_to_zero(double*, long);
void integrate_noise(double*, fftw_complex*, double, double, long, double, int*, double);
void find_fpt(double*, double*, double, long, double, double**, double, int);
void fpt_to_zvar(double, double, double);
void initialise_critical_bridge(bridge_process**, double, double, double, double, double, double, bridge_process*);
void split_and_search_bridge(bridge_process*, int*, double*,  double);

bridge_process* check_this_bridge(bridge_process*, double*, double);
bridge_process* split_bridge(bridge_process*);
double generate_random_conditional_midpoint(double, double, double, double);
double crossing_time_of_bridge(bridge_process*);
double time_time_correlation(double, double, double);

void enlarge_QI(void);
void print_QI(void);
void print_bridge(bridge_process*);
void copy_QI(double**, int, double, double*, double);
void free_tree(bridge_process**);
void free_bridge(bridge_process**);

// GLOBAL VARIABLES
extern int max_generation;
extern double *gamma_N_vec;
extern double *g_vec;

// GSL RNG
extern const gsl_rng_type *T;
extern gsl_rng *r;
extern int seed;
extern triag_matrix *QI; //malloc somewhere
