/* This code simulates a fractional Brownian Motion on the interval [0,1] and finds the first passage time.
 * 
 * Authors: Benjamin Walter (Imperial College) , Kay Wiese (ENS Paris)
 * 
 * For any questions/feedback/comments --> b.walter16@imperial.ac.uk
 * 
 */


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

#include "fbm_header.h"

double time_GEN_MIDPOINT = 0.0;
double time_STEP_1 = 0.0, time_STEP_2 = 0.0, time_STEP_3 = 0.0, time_STEP_4 = 0.0, time_STEP_5= 0.0,time_STEP_6 = 0.0, time_STEP_7;

int main(int argc, char *argv[])
{
	setlinebuf(stdout);
	clock_t t_start = clock();
	// VARIABLES
	
	// physics
	// Default values, overwritten by getopt
	double frac_drift = 0.0;
	double lin_drift = 0.0;
	double hurst = 0.5; 	// Hurst parameter 0 < h < 1
	int g =8;		// 2^g is number of points
	
	// experiment
	int iteration = 10000;	// Size of ensemble
	double delta_passage_heigths = .1; // This is the gap between lines at which first passages will be read
	double epsilon = 1e-9; // Tolerance probability with which a midpoint might be falsenegative, i.e. traverse the boundary
	int iter;
	

	
	// observables
	double *passage_heights;
	double *first_passage_times;
	double *zvar; // z = m /(sqrt(2) * t^H)

	// GSL RNG
	const gsl_rng_type *T;
	gsl_rng *r;
	int seed = -1;

	// ARGUMENTS
	opterr = 0;
	int c = 0;
        while( (c = getopt (argc, argv, "h:g:S:I:m:n:E:") ) != -1)
	{                switch(c)
                        {
				case 'm':
					lin_drift = (( double) atof(optarg));
					break;
                        	case 'n':
					frac_drift = (( double) atof(optarg));
					break;
                                case 'h':
                                        hurst = atof(optarg);
                                        break;
				case 'g':
					g = atoi(optarg);
					break;
				case 'S':
					seed = atoi(optarg);
					break;
				case 'I':
					iteration = atoi(optarg);
					break;
				case 'E':
					epsilon = atof(optarg);
					break;
                       		default:
                                exit(EXIT_FAILURE);
                        }
	}// getopt ends

	//int i; // Index
	long N = ((long) pow(2,g));
	double invN = 1/(( double) N);
	double sigma = 1.0;//sqrt(2.0);

	printf("# FBM Simulator, Hurst: %g, N: %ld \n", hurst, N);	
	
	clock_t t_begin_alloc_pointers = clock();

	fftw_complex *correlation, *circulant_eigenvalues, *rndW, *fracGN;
       	double *correlation_exponents, *fracbm, **QInverseCorrelationCatalogue;
	complex_z *randomComplexGaussian;
	fftw_plan p1, p2 ;


	// Initialise observables
	initialise(&correlation, &circulant_eigenvalues, &rndW, &fracGN, &correlation_exponents, &randomComplexGaussian, N, &r, &T, seed);
	initialise_observables(&passage_heights, &first_passage_times, &zvar, delta_passage_heigths);
	initialise_trajectory(&fracbm, N);
	initialise_inverse_correlation_matrix(&QInverseCorrelationCatalogue, N);

	// Initialise FFT plans
	p1 = fftw_plan_dft_1d(2*N , correlation, circulant_eigenvalues, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(2*N, rndW, fracGN, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	clock_t t_end_alloc_pointers = clock();
	
	// Print out header
	printf("# FBM Simulator\n# Bridge with Matrix inversion\n# Hurst parameter: %g\n# Linear drift (mu): %g\n#Fractional drift (nu): %g\n# N: %ld\n# MAX generation: %i \n", hurst, lin_drift, frac_drift, N, MAX_GENERATION);	
	
	// Write correlation of noise
	write_correlation_exponents(&correlation_exponents, N, invN, hurst);
	write_correlation(&correlation, correlation_exponents, N);

	// Write Inverse of correlation matrix of FBM ('Q'(N)-matrix)
	write_inverse_correlation_matrix(&QInverseCorrelationCatalogue, N, hurst);	
		
	// FFT into circulant eigenvalues
	clock_t t_begin_fft_p1 = clock();
	fftw_execute(p1); 
	clock_t t_end_fft_p1 = clock();
	
	// Some time measurements
	double time_RNG = 0.0;
	double time_FFT = ( ((double) t_end_fft_p1) - ((double) t_begin_fft_p1));
	double time_INTEGRATE = 0.0;
	double time_FPT = 0.0;
	double time_ITERATION = 0.0;
	clock_t t_begin_random, t_end_random, t_begin_FFT_p2, t_end_FFT_p2,  t_begin_integrate, t_end_integrate, t_begin_find_fpt, t_end_find_fpt ;

	clock_t t_begin_iterate, t_end_iterate;

	int last_point_index; // Index of the first point to cross the barrier (=N, if this doesn't happen)

	for(iter = 0; iter < iteration; iter++)
	{
		t_begin_iterate = clock();

		/*corr_write_time = 0.0;
		corr_write_index = 0;*/

		t_begin_random = clock();
		generate_random_vector(&randomComplexGaussian, &rndW, circulant_eigenvalues, N, sigma, r);
		
		t_end_random = clock();

		/*if( (iter % 10000) == 0)
		{
			// Hard coded, clean up memory leak behind fftw
			fftw_destroy_plan(p2);
			fftw_cleanup();	
			p2 = fftw_plan_dft_1d(2*N, rndW, fracGN, FFTW_BACKWARD, FFTW_ESTIMATE);
		}*/
		t_begin_FFT_p2 = clock();
		fftw_execute(p2);	
		t_end_FFT_p2 = clock();
		
		// 
		// Reset first passage times
		
		set_to_zero(&first_passage_times, ( (int) (MSTEP + 1)) );

		// Integrate fractional Gaussain noise to fBM
		t_begin_integrate = clock();	
		integrate_noise(&fracbm, fracGN, lin_drift, frac_drift, N, hurst, &last_point_index, passage_heights);
		t_end_integrate = clock();	
		
		// Find maximum to recursive depth RECURSION_DEPTH
		t_begin_find_fpt = clock();
		find_fpt(fracbm, &first_passage_times, passage_heights, N, epsilon, r, QInverseCorrelationCatalogue, hurst, MIN((last_point_index), N));

		t_end_find_fpt = clock();

		// Convert first passage times into Laplace variables
		fpt_to_zvar(&zvar, passage_heights, first_passage_times, hurst);
		
		t_end_iterate = clock();

		time_FPT += (((double) t_end_find_fpt) - ((double) t_begin_find_fpt));
		time_RNG += (((double) t_end_random) - ((double) t_begin_random));
		time_FFT += (((double) t_end_FFT_p2) - ((double) t_begin_FFT_p2));
		time_INTEGRATE += (((double) t_end_integrate) - ((double) t_begin_integrate));		
		time_ITERATION += (((double) t_end_iterate) - ((double) t_begin_iterate));		
	}// End iteration	

	double time_total = (((double) t_end_iterate) - ((double) t_start))/CLOCKS_PER_SEC;

	printf("# ***************\n# TIME STATISTICS\n# Time total: %g\
			\n# Time spent on iteration: %g\n# Time spent on one realisation of 2^%i (average): %g\n# Time spent on FFTs: %g\n# Time spent on drawing random numbers: %g\n\
# Time spent on integrating noise: %g\n\
# Time spent on allocating memory: %g\n\
# Time spent on finding FPT: %g\n\
# Time spend generating midpoints %g\n\
# \n\
# Time Step 1:\t%g (%g %%)\n\
# Time Step 2:\t%g (%g %%)\n\
# Time Step 3:\t%g (%g %%)\n\
# Time Step 4:\t%g (%g %%)\n\
# Time Step 5:\t%g (%g %%)\n\
# Time Step 6:\t%g (%g %%)\n\
# Time Step 7:\t%g (%g %%)\n",\
		       	time_total ,\
		       	time_ITERATION/CLOCKS_PER_SEC,\
		       	g,\
		       	time_ITERATION/(CLOCKS_PER_SEC* (double) iteration),\
		       	time_FFT/CLOCKS_PER_SEC, time_RNG/CLOCKS_PER_SEC,\
		       	time_INTEGRATE/CLOCKS_PER_SEC,\
		       	(((double) t_end_alloc_pointers) - ((double) t_begin_alloc_pointers))/CLOCKS_PER_SEC,\
		       	time_FPT/CLOCKS_PER_SEC,\
			(time_GEN_MIDPOINT/CLOCKS_PER_SEC),\
													    (time_STEP_1/CLOCKS_PER_SEC), (100*time_STEP_1/time_GEN_MIDPOINT),\
							(time_STEP_2/CLOCKS_PER_SEC),(100*time_STEP_2/time_GEN_MIDPOINT),\
							(time_STEP_3/CLOCKS_PER_SEC),(100*time_STEP_3/time_GEN_MIDPOINT),\
							(time_STEP_4/CLOCKS_PER_SEC),(100*time_STEP_4/time_GEN_MIDPOINT),\
							 (time_STEP_5/CLOCKS_PER_SEC),(100*time_STEP_5/time_GEN_MIDPOINT),\
													     (time_STEP_6/CLOCKS_PER_SEC),(100*time_STEP_6/time_GEN_MIDPOINT),\
													      (time_STEP_7/CLOCKS_PER_SEC),(100*time_STEP_7/time_GEN_MIDPOINT));	       


	return 0;
}


