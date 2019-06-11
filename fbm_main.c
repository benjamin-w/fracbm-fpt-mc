/* This code simulates a fractional Brownian Motion on the interval [0,1] and finds the first passage time.
 * 
 * Authors: Benjamin Walter (Imperial College) , Kay Wiese (ENS Paris)
 * 
 * For any questions/feedback/comments --> b.walter16@imperial.ac.uk
 * 
 */

#include "fbm_header.h"

int main(int argc, char *argv[])
{
	setlinebuf(stdout);
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
					lin_drift = ( - ( double) atof(optarg)); // Sign chosen in accordance with paper. 
					break;
                        	case 'n':
					frac_drift = (- ( double) atof(optarg));
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
	
	
	// Print out header
	printf("# FBM Simulator\n# Bridge with Matrix inversion\n# Hurst parameter: %g\n# Linear drift (mu): %g\n#Fractional drift (nu): %g\n# N: %ld\n# MAX generation: %i \n", hurst, lin_drift, frac_drift, N, MAX_GENERATION);	
	
	// Write correlation of noise
	write_correlation_exponents(&correlation_exponents, N, invN, hurst);
	write_correlation(&correlation, correlation_exponents, N);

	// Write Inverse of correlation matrix of FBM ('Q'(N)-matrix)
	write_inverse_correlation_matrix(&QInverseCorrelationCatalogue, N, hurst);	
		
	// FFT into circulant eigenvalues
	fftw_execute(p1); 
	
	int last_point_index; // Index of the first point to cross the barrier (=N, if this doesn't happen)

	for(iter = 0; iter < iteration; iter++)
	{
		generate_random_vector(&randomComplexGaussian, &rndW, circulant_eigenvalues, N, sigma, r);
		
		fftw_execute(p2);	
		
		// Reset first passage times
		set_to_zero(&first_passage_times, ( (int) (MSTEP + 1)) );

		// Integrate fractional Gaussain noise to fBM
		integrate_noise(&fracbm, fracGN, lin_drift, frac_drift, N, hurst, &last_point_index, passage_heights);
		
		// Find maximum to recursive depth RECURSION_DEPTH
		find_fpt(fracbm, &first_passage_times, passage_heights, N, epsilon, r, QInverseCorrelationCatalogue, hurst, MIN((last_point_index), N));


		// Convert first passage times into Laplace variables
		fpt_to_zvar(&zvar, passage_heights, first_passage_times, hurst);
	}// End iteration	
	return 0;
}


