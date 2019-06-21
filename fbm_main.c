/* This code simulates a fractional Brownian Motion on the interval [0,1] and finds the first passage time.
 * 
 * Authors: Benjamin Walter (Imperial College) , Kay Wiese (ENS Paris)
 * 
 * For any questions/feedback/comments --> b.walter16@imperial.ac.uk
 * 
 */

#include "fbm_header.h"

int max_generation;
triag_matrix *QI;
// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;
int seed;



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
	max_generation = 8;
	
	// experiment
	int iteration = 10000;	// Size of ensemble
	double epsilon = 1e-9; // Tolerance probability with which a midpoint might be falsenegative, i.e. traverse the boundary
	int iter;
	

	
	// observables
	double passage_heights = 0.1;
	double first_passage_times = 0.0;
	int seed = -1;

	// ARGUMENTS
	opterr = 0;
	int c = 0;
        while( (c = getopt (argc, argv, "h:g:G:S:I:m:n:E:") ) != -1)
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
				case 'G':
					max_generation = atoi(optarg);
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
	initialise_trajectory(&fracbm, N);
	initialise_inverse_correlation_matrix(&QInverseCorrelationCatalogue, N);

	//Initialise global observables
	ALLOC(gamma_N_vec, pow(2, (g+max_generation))); // Maximal number of points possible
	ALLOC(g_vec, pow(2,(g+max_generation)));
	long array_length = 2*N;
	ALLOC(QI, 1);
	ALLOC(QI->inv_corr_matrix, (array_length*(array_length + 1) / 2));
	ALLOC(QI->trajectory_x, (array_length + 1));
	ALLOC(QI->trajectory_t, (array_length + 1));
	QI->array_length = array_length;
	QI->hurst_parameter = hurst;

	// Initialise FFT plans
	p1 = fftw_plan_dft_1d(2*N , correlation, circulant_eigenvalues, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(2*N, rndW, fracGN, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	
	// Print out header
	printf("# FBM Simulator\n# Bridge with Matrix inversion\n# Hurst parameter: %g\n# Linear drift (mu): %g\n#Fractional drift (nu): %g\n# Barrier height at %g\n# N: %ld\n# MAX generation: %i \n", hurst, lin_drift, frac_drift, passage_heights, N, max_generation);	
	
	// Write correlation of noise
	write_correlation_exponents(correlation_exponents, N, invN, hurst);
	write_correlation(correlation, correlation_exponents, N);
	
	// Write Inverse of correlation matrix of FBM ('Q'(N)-matrix)
	write_inverse_correlation_matrix(QInverseCorrelationCatalogue, N, hurst);	
		
	// FFT into circulant eigenvalues
	fftw_execute(p1); 
	
	int last_point_index; // Index of the first point to cross the barrier (=N, if this doesn't happen)

	for(iter = 0; iter < iteration; iter++)
	{
		
		generate_random_vector(randomComplexGaussian, rndW, circulant_eigenvalues, N, sigma);
		
		fftw_execute(p2);	
		
		// Reset first passage times
		first_passage_times = 0.0;
		// Integrate fractional Gaussain noise to fBM
		integrate_noise(fracbm, fracGN, lin_drift, frac_drift, N, hurst, &last_point_index, passage_heights);
		
		// Find maximum to recursive depth RECURSION_DEPTH
		find_fpt(fracbm, &first_passage_times, passage_heights, N, epsilon,  QInverseCorrelationCatalogue, hurst, last_point_index);
		

		// Convert first passage times into Laplace variables
		fpt_to_zvar(passage_heights, first_passage_times, hurst);
		
	}// End iteration

	return 0;
}


