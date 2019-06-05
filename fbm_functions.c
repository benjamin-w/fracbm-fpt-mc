#include "fbm_header.h"

void initialise( fftw_complex** correlation,  fftw_complex** circulant_eigenvalues,  fftw_complex** rndW,  fftw_complex** fracGN,  double** correlation_exponents, complex_z** randomComplexGaussian, long N, gsl_rng** r,const gsl_rng_type** T, int seed)
{
	// These are the objects that are N long (the increments)
	*correlation =  fftw_malloc( (2 * N ) * sizeof( fftw_complex ) );
	if(!(*correlation)){printf("correlation malloc failed.\n"); exit(1) ;};
	*circulant_eigenvalues = fftw_malloc( (2 * N ) * sizeof( fftw_complex) );
	if(!(*circulant_eigenvalues)){printf("circulant eigenvalues malloc failed.\n"); exit(1) ;};
	*rndW =fftw_malloc( 2 * N * sizeof( fftw_complex ));
	if(!(*rndW)){printf("rndW malloc failed.\n"); exit(1) ;};
	*fracGN = fftw_malloc( 2 * N * sizeof( fftw_complex) );
	if(!(*fracGN)){printf("fracGN malloc failed.\n"); exit(1) ;};
	*correlation_exponents =  calloc(( N + 1) , sizeof( double));
	if(!(*correlation_exponents)){printf("correlation_exponents malloc failed.\n"); exit(1) ;};
	*randomComplexGaussian = (complex_z *) malloc(N * sizeof(complex_z));
	if(!(*randomComplexGaussian)){printf("randomComplexGaussian malloc failed.\n"); exit(1) ;};
	
	/* Initialises GSL Random Generator */
        gsl_rng_env_setup();
        *T = gsl_rng_default;
        *r = gsl_rng_alloc (*T);
        if(seed==-1) seed = ((int) (((int) clock() ) % 100000));
        gsl_rng_set(*r, seed);
        printf("# SEED %i\n",seed);
}

void initialise_trajectory(  double ** fracbm, long N)
{
	// N increments give N+1 points
	*fracbm = malloc( (N + 1) * sizeof(  double));
	if(!(*fracbm)){printf("fracbm malloc failed.\n"); exit(1);}
}

void initialise_inverse_correlation_matrix(double*** QInverseCorrelation, long N )
{
	// This is inverse correlation matrix of X_1...X_N with X_0 = 0 fixed.
	
	long i;
	*QInverseCorrelation = malloc( N * sizeof(double*));
	for(i = 0; i < N; i++)
	{
		(*QInverseCorrelation)[i] = (double*) malloc( ((N * (N+1)) / 2) * sizeof(double));	
	}
}

void initialise_observables(double** passage_heights, double** first_passage_times, double** zvar, double delta_passage_heights)
{
	*passage_heights = malloc( (MSTEP + 1) *  sizeof(  double ));
	assert(*passage_heights != NULL);
	*first_passage_times = malloc((MSTEP + 1) * sizeof(  double ));
	assert(*first_passage_times !=NULL);
	*zvar = malloc((MSTEP + 1) * sizeof(  double ));
	assert(*zvar != NULL);

	int i;
	for(i=0; i <= MSTEP; i++)
	{
		(*passage_heights)[i] = (i * delta_passage_heights);
		printf("# Barrier %i at %g\n",i,(*passage_heights)[i]);

	}	
}

double erfcinv(double x)
{
	/* Approximation of the inverse erfc function followin Blair, Edwards, Johnson 1976 */
	double eta = - log(sqrt(M_PI) * x);
	double erfcinv_x = sqrt( eta - 0.5* log(eta) + (1/eta)*(0.25*log(eta) - 0.5));
	return erfcinv_x;
}

void write_correlation_exponents(double** correlation_exponents, long N, double invN, double hurst)
{
	int i;
	*correlation_exponents[0] = 0.0; // 0^0 = 0.
	for(i=1; i<=N; i++)
	{
		(*correlation_exponents)[i] =pow( (( double) i) * invN, ( (double) 2*hurst ));
	}
}

void write_correlation(fftw_complex** correlation, double* correlation_exponents, long N)
{
	(*correlation)[0][0] = 2*(correlation_exponents[1] );
	(*correlation)[0][1] = 0.0;
	(*correlation)[N][0] = 0.0;
	(*correlation)[N][1] = 0.0;
	
	int i;
	for(i=1; i<=N-1; i++)
	{
		(*correlation)[2*N - i][0] = (*correlation)[i][0] = (correlation_exponents[i+1] + correlation_exponents[i-1] - 2*correlation_exponents[i]);
		(*correlation)[2*N - i][1] = (*correlation)[i][1] = 0;
	}
}

void write_inverse_correlation_matrix(double*** Q, long N, double hurst)
{
	// This is a catalogue of N different correlation matrices where the n.th correlation matrix corresponds to the one of the first N points fixed.

	int NN = ((int)N); // ! Long is being casted int ! but just because LAPACKE needs to deal with it.
	double delta_t = (1/ ((double) N));
	
	int size; // This is the size of the matrix going to go from 1 to N
	
	int i,j;
	lapack_int info;
        
	for(size = 1; size <= NN; size++)
	{
		for(i = 1; i <= size; i++)
		{
			for(j=i; j <= size; j++)
			{
				((*Q)[size-1])[IJ2K((i-1),(j-1))] =  time_time_correlation( i * delta_t, j * delta_t, hurst);
			}
		}
		info=LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'U', size, (*Q)[size-1]);
		if(info != 0){printf("Lapack Cholesky decomposition of correlation matrix failed.\n"); exit(1);}
		info=LAPACKE_dpptri(LAPACK_COL_MAJOR, 'U', size, (*Q)[size-1]);
		if(info != 0){printf("Lapack inversion of correlation matrix failed.\n"); exit(1);}
	}
}

void generate_random_vector(complex_z** randomvector, fftw_complex** rndW, fftw_complex *circulant_eigenvalues, long N, double sigma, gsl_rng* r)
{
	// This routine creates the random vector used in Davies Harte generation of subgrid
	int i;
	double invN = 1/((double) N);

	for(i = 0; i < N; i++)
		{
			(*randomvector)[i].r= gsl_ran_gaussian_ziggurat(r,sigma);
			(*randomvector)[i].i = gsl_ran_gaussian_ziggurat(r,sigma);
		}

	(*rndW)[0][0] = sqrt(0.5* circulant_eigenvalues[0][0]*invN)*gsl_ran_gaussian_ziggurat(r,sigma);
	(*rndW)[0][1] = 0.0;

	(*rndW)[N][0] = sqrt(0.5*circulant_eigenvalues[N][0]*invN)*gsl_ran_gaussian_ziggurat(r,sigma);
	(*rndW)[N][1] = 0.0;
	
	for(i=1; i < N; i++)
	{
		(*rndW)[i][0] = (sqrt(0.25*circulant_eigenvalues[i][0]*invN) * (*randomvector)[i].r);
		(*rndW)[i][1] = sqrt(0.25*circulant_eigenvalues[i][0]*invN) * (*randomvector)[i].i;
		(*rndW)[2*N-i][0] =(*rndW)[i][0]; 
		(*rndW)[2*N-i][1] = -(*rndW)[i][1]; 
	}
}

void set_to_zero(double** pointer, int length)
{
	// generic function to "re-calloc" pointer
	int i;
	for(i = 0; i < length; i++)
	{
		(*pointer)[i] = 0.0;
	}
}

void integrate_noise(double** fracbm, fftw_complex* fracGN,double lin_drift, double frac_drift, long N, double hurst, int* last_point_index, double* passage_heights)
{	
	// Find the first point to jump over the barrier (if exists). Then throw away all points behind. Take the appropiate inverse matrix and pass it on.

	// Assuming passage_heights is only 1 big.
	
	(*fracbm)[0] = 0.0;
	double delta_t = (1/((double) N));
	*last_point_index = ((int) N);

	int i;
	for(i = 1; i <= N; i++)
	{
		(*fracbm)[i] = ((*fracbm)[i-1] + fracGN[i-1][0] + (lin_drift + frac_drift*(pow((((double)i)-0.5)*delta_t, (2*hurst - 1.0))))*delta_t); // Ito integration, fractional drift is evaluated at midpoint --> Check for consequences
		if( (*fracbm)[i] > passage_heights[1])
		{
			if(i < *last_point_index){*last_point_index = i;}
		}

		if( i == (*last_point_index)) break; 
	}
}

void find_fpt(double* fracbm, double** first_passage_times, double* passage_heights, long N, double epsilon, gsl_rng* r, double** QCatalogue, double hurst, int last_point_index)
{
	// Find FPT knowing that first passage happens in [0, last_point_index * delta_t]
	int i;
	double delta_t = (1/((double) N));
	int height_counter = 1;
	
	triag_matrix *QI;
	copy_QI(&QI, QCatalogue, last_point_index, hurst, &fracbm, delta_t); // Here a local copy of QI is created that is conditioned on last_point_index
	int fpt_found = 0;
	double critical_strip = (erfcinv(2*epsilon)*(sqrt( ( (4.0/pow(2.0,2*hurst)) - 1)))*pow(delta_t, hurst)); // See currentnotes.pdf for calculation

	bridge_process *critical_bridge;
	critical_bridge = NULL;

	// MISSING: for loop over all height_counters
	/* The tree algorithm */
	i = 0;
	while(fpt_found == 0)
	{
		i++; // Go to next bridge
		if(i > last_point_index) break; // if i>=N, then last box is always wrong.  

		// The criterion to split is whether either endpoint lies in the critical zone
		if ((MAX(fracbm[i],fracbm[i-1])) > (passage_heights[height_counter] - critical_strip))
		{
			// Bridge is critical
			initialise_critical_bridge(&critical_bridge, ((double) i*delta_t), fracbm[i], ((double)(i-1)*delta_t), fracbm[i-1],  passage_heights[height_counter], critical_strip, critical_bridge); 			
			// Start splitting bridge
			split_and_search_bridge(&critical_bridge, &fpt_found, &((*first_passage_times)[height_counter]), r, delta_t, &QI);
		}
	}

	if(fpt_found == 0)
	{
		((*first_passage_times)[height_counter]) = 1.0; 
	}

	bridge_process* point_bridge_old = critical_bridge;
	bridge_process* point_bridge_new = critical_bridge->previous_root;
	
	if(point_bridge_new != NULL){
		while(point_bridge_new != NULL)
		{
			//print_bridge(&point_bridge_old);
			free_tree(&(point_bridge_old->root_bridge));
			free_bridge(&(point_bridge_old->root_bridge)); // At the end you need to kill the root.
			point_bridge_old = point_bridge_new;
			point_bridge_new = point_bridge_old->previous_root;
		}
	}

	if(point_bridge_new == NULL){free_tree(&(point_bridge_old->root_bridge));free_bridge(&(point_bridge_old->root_bridge));} 
	
	free_QI(&QI); // Free entire structure
}

void fpt_to_zvar(double** zvar, double* passage_heights, double* first_passage_times, double hurst)
{
	int i;
	for(i = 1; i <= MSTEP; i++)
	{
		if(first_passage_times[i] > 0)
		{
			(*zvar)[i] = (passage_heights[i] / (sqrt(2.0) * pow(first_passage_times[i], hurst) ));
	
		}
		printf("%.12f\n", (*zvar)[i]);
	}
}

void initialise_critical_bridge(bridge_process** root_bridge, double rtime, double rvalue, double ltime, double lvalue, double threshold, double critical_strip, bridge_process* old_bridge)
{
	(*root_bridge) = (bridge_process*) malloc(sizeof(bridge_process));
	(*root_bridge)->right_time = rtime;
	(*root_bridge)->right_value = rvalue;
	(*root_bridge)->left_time = ltime;
	(*root_bridge)->left_value = lvalue;
	(*root_bridge)->generation = 0;
	(*root_bridge)->threshold = threshold;
	if( (rvalue > threshold) && (lvalue < threshold)){ (*root_bridge)->crossing_threshold = 1;} else {(*root_bridge)->crossing_threshold = 0;}
	(*root_bridge)->critical_strip = critical_strip; 
	(*root_bridge)->centre_critical = 1; // Because otherwise this function wouldn't have been invoked. 
	(*root_bridge)->left_sub_bridge = NULL;
	(*root_bridge)->right_sub_bridge = NULL;
	(*root_bridge)->parental_bridge = NULL;
	(*root_bridge)->root_bridge = (*root_bridge);
	if(old_bridge == NULL){(*root_bridge)->previous_root = NULL;}
	else{
		(*root_bridge)->previous_root = old_bridge->root_bridge;

	}; // Last Tree needs to get appended to then deleta all of them.
	(*root_bridge)->pointer_stack = (bridge_process **) malloc( (pow(2, MAX_GENERATION + 1) - 1) * sizeof( bridge_process* )); // There are up to 2^G - 1 nodes in a tree. 
	(*root_bridge)->pointer_stack[0] = (*root_bridge);
	(*root_bridge)->stack_top = 1;
}

void split_and_search_bridge(bridge_process** initial_bridge, int* fpt_found, double* fpt, gsl_rng* r, double delta_t, triag_matrix** QI)
{
	bridge_process *which_bridge_shall_i_check;
        which_bridge_shall_i_check= (*initial_bridge);
	(*fpt_found) = 0;
	double fpt_best_guess = 1.0; // In case that after the algorithm this value still has this value, we know that it hasn't been found yet and will keep fpt_found at 0, otherwise it will be overwritten by the best guess.
	while(which_bridge_shall_i_check)
	{
		which_bridge_shall_i_check = check_this_bridge(&which_bridge_shall_i_check, &fpt_best_guess, r, delta_t, QI);
	}


	if(fpt_best_guess != 1.0){*fpt_found = 1; *fpt = fpt_best_guess;}else{*fpt_found = 0;}
}

bridge_process* check_this_bridge(bridge_process** incoming_bridge, double* fpt_best_guess, gsl_rng* r, double delta_t, triag_matrix** QI)
{
	/* CHECK_THIS_BRIDGE
	 * COMMENTARY
	 * This function reads a bridge process (i.e. a fractional brownian bridge of which only the endpoints are known) and decides which bridge process to look at next, i.e. it returns a pointer to the next bridge process. If this bridge process is a subprocess (i.e. a bridge between an end- and a new mid-point), this bridge will be generated and evaluated before returning the pointer. This is realised in another external function.
	 * Each bridge process can be categorised into one out of eight states that are characterised by five properties.
	 * The five properties are binary
	 * A) Whether the process has less than two children (0) or not (1).
	 * B) Whether the bridge process' left endpoint is after the best estimate of the first passage time, no (0), yes (1)
	* C) Whether it crosses the threshold, i.e. contains a passage event, (0), or not (1).
	 * D) Whether the generation of the bridge process is maximal, i.e. the degree of subdivisions that were necessary to construct that bridge process have reached the maximal desired resolution, (1), or not (0).
	 * E) Whether it center is critical, i.e. whether with a reasonable probability ( > \epsilon) its midpoint can be expected to trigger a passage event, (0), or not (1).
	 * 
	 * The corresponding values of A)-E) determine the reply of the system. Let's denote the state by a word of up to five bits with '0' or '1'. If the word is shorter, it means that the subsequent letters do not affect the outcome any longer
	 * The process is terminated, whenever word=1, 01, 0001, 0011, or 00101
	 * The process leads to a split of the bridge process whenever word=00000 or 00100
	 * The process leads to a split of the right subprocess only whenever word=00001 (this is a pathological case with probability < epsilon
	 *
	 */
	
	if( (*incoming_bridge) == NULL)
	{
		return NULL;
	}

	bridge_process *outgoing_bridge;
	
	if( (*incoming_bridge)->right_sub_bridge != NULL)
	{
		// 1*
		// If this is the case, this means the process is already going back to its root. Both children have already been created
		outgoing_bridge = (*incoming_bridge)->parental_bridge;
	}
	else
	{
		// 0*
		if( ((*fpt_best_guess) - (*incoming_bridge)->left_time) > -delta_t )// This is the maximum distance where still theoretically the FPT could be improved
		{
			// 00*
			if(((*incoming_bridge)->crossing_threshold) == 1)
			{
				// 000*
				// Update FPT if a bridge with a crossing is passed.
				*fpt_best_guess = crossing_time_of_bridge(*incoming_bridge); // If there is a crossing, update FPT. It can be that the FPT of a finer bridge is *after* the coarser one, but that is fine, latter is always a better estimate!
				
				if( ((*incoming_bridge)->generation) < MAX_GENERATION)
				{
					// 0000*
					return split_bridge(incoming_bridge, r, QI);
				}
				else
				{
					// 0001*
					//outgoing_bridge = (*incoming_bridge)->parental_bridge;
					return NULL; // Yes, because if a FPT was found in the smallest bridge possible, it is not going to be improved. If you delete this line, in fact it will get worse if a parent had no right child and its left time is good enough. Long time undetected "bug".
				}
			}
			else
			{
				// 001*
				if( ((*incoming_bridge)->generation) < MAX_GENERATION)
				{
					// 0010
					if( ((*incoming_bridge)->centre_critical ))
					{
						// 00100
						return split_bridge(incoming_bridge, r, QI);
					}
					else
					{
						// 00101
						outgoing_bridge = (*incoming_bridge)->parental_bridge;
					}
				}
				else
				{
					// 0011
					outgoing_bridge = (*incoming_bridge)->parental_bridge;
				}
			}
		}
		else
		{
			//01*
			outgoing_bridge = (*incoming_bridge)->parental_bridge;
		}
	}
	return outgoing_bridge;
}

bridge_process* split_bridge(bridge_process** parent_bridge, gsl_rng* r, triag_matrix** QI)
{
	/* SPLIT_BRIDGE
	* COMMENTARY
	* This function recevies a bridge process from check_this_bridge with the task to return the left or right subbridge spanning the midpoint with either of the endpoints and with full information (as required by the struct bridge_process in fbm_header.h)
	* If the left child has not been generated yet, a midpoint has to be drawn according to a Gaussian distribution known from fBM bridges and conditioned on all previously known points
	(* If the left child already exists, and a right child is to be generated, it simply inherits the right endpoint of the left child as its left endpoint, now drawing has to be done.
	*/

	bridge_process *sub_process;
	sub_process = NULL;
	
	double hurst = ((*QI)->hurst_parameter);
	
	// As bridge arrives, it has either zero or one child. Check.
	if( (*parent_bridge)->left_sub_bridge == NULL)
	{
		// Make left child
		sub_process = (bridge_process*) malloc(sizeof(bridge_process));
		if(sub_process == NULL){printf("Allocation of sub process failed\n");exit(1);}
		sub_process->right_time = 0.5*( (*parent_bridge)->left_time + (*parent_bridge)->right_time);
		sub_process->right_value = generate_random_conditional_midpoint( ((*parent_bridge)->right_time), ((*parent_bridge)->right_value), ((*parent_bridge)->left_time), ((*parent_bridge)->left_value), r, QI);
		sub_process->left_time = (*parent_bridge)->left_time;
		sub_process->left_value = (*parent_bridge)->left_value;
		sub_process->generation = (((*parent_bridge)->generation)+1);
		sub_process->threshold = ((*parent_bridge)->threshold);
		sub_process->critical_strip = (pow(2,-hurst) * (*parent_bridge)->critical_strip);
		if( MAX( (sub_process->right_value), (sub_process->left_value)) > ((sub_process->threshold) - (sub_process->critical_strip))){ sub_process->centre_critical = 1;}else{sub_process->centre_critical = 0;}
		if( ((sub_process->left_value) < (sub_process->threshold)) && ((sub_process->right_value) > (sub_process->threshold))){sub_process->crossing_threshold = 1;}else{sub_process->crossing_threshold = 0;}
		
		sub_process->left_sub_bridge = NULL;
		sub_process->right_sub_bridge = NULL;
		// Link it to parent
		sub_process->parental_bridge = (*parent_bridge);
		(*parent_bridge)->left_sub_bridge = sub_process;
		
		// Link this child to its roots
		sub_process->root_bridge = (*parent_bridge)->root_bridge;
		((sub_process->root_bridge)->pointer_stack)[(sub_process->root_bridge)->stack_top] = sub_process;
		(sub_process->root_bridge)->stack_top++;
	
		// Blank all out that only belongs to root
		sub_process->previous_root = NULL;
		sub_process->pointer_stack = NULL;	
		sub_process->stack_top = 0;
	}
	else
	{
		// Make right child
		sub_process = (bridge_process*) malloc(sizeof(bridge_process));
		if(sub_process == NULL){printf("Allocation of sub process failed\n");exit(1);}
		sub_process->right_time = (*parent_bridge)->right_time;
		sub_process->right_value = (*parent_bridge)->right_value;
		sub_process->left_time = 0.5*( (*parent_bridge)->left_time + (*parent_bridge)->right_time);
		sub_process->left_value = (((*parent_bridge)->left_sub_bridge)->right_value);
		sub_process->generation = (((*parent_bridge)->generation)+1);
		sub_process->threshold = ((*parent_bridge)->threshold);
		sub_process->critical_strip = (pow(2,-hurst) * (*parent_bridge)->critical_strip); 
		if( MAX( (sub_process->right_value), (sub_process->left_value)) > ((sub_process->threshold) - (sub_process->critical_strip))){ sub_process->centre_critical = 1;}else{sub_process->centre_critical = 0;}
		if( ((sub_process->left_value) < (sub_process->threshold)) && ((sub_process->right_value) > (sub_process->threshold))){sub_process->crossing_threshold = 1;}else{sub_process->crossing_threshold = 0;}
		sub_process->left_sub_bridge = NULL;
		sub_process->right_sub_bridge = NULL;

		// Link it to parent
		sub_process->parental_bridge = (*parent_bridge);
		(*parent_bridge)->right_sub_bridge = sub_process;
		
		// Link this child to its roots
		sub_process->root_bridge = (*parent_bridge)->root_bridge;
		((sub_process->root_bridge)->pointer_stack)[(sub_process->root_bridge)->stack_top] = sub_process;
		(sub_process->root_bridge)->stack_top++;
		
		// Blank all out that only belongs to root
		sub_process->previous_root = NULL;
		sub_process->pointer_stack = NULL;	
		sub_process->stack_top = 0;
	}

	return sub_process;
}


double generate_random_conditional_midpoint( double right_time, double right_value, double left_time, double left_value, gsl_rng* r, triag_matrix** QI)
{
	double mean, sigma /*should be "\sigma^2" ! */, midpoint;
	/* What's the random midpoint conditioned on being a (fractional) Brownian motion, conditioned on all points known so far?
	The inverse correlation matrix of all pointsis stored in QI. The new point needs to be added according to mean and variance as found of that matrix. The matrix then needs to be enlarged to also include the new point. 
	Algorithm (with identical variable names) in 'fracbm/doc/currentnotes.pdf'
	*/
		
	/* Matrix inversion core */
	long i;
	double midtime = (0.50*(right_time + left_time));	
	long number_of_points = (*QI)->size;
	double *gamma_N_vec = (double*) malloc(number_of_points * sizeof(double));
	double *g_vec = (double*) malloc(number_of_points * sizeof(double));
	double hurst = (*QI)->hurst_parameter;
	
	// Step 1, Gamma vector
	for(i = 0; i < number_of_points; i++)
	{
		gamma_N_vec[i] = time_time_correlation(midtime, (*QI)->trajectory_t[i+1], hurst); // No cross correlation with t_0 = 0. \gamma_i = <t_i \tilde{t}> for t_i > 0
	}

	// Step 2, g-vector. g = Q * \gamma
	cblas_dspmv(CblasColMajor, CblasUpper, ((int) number_of_points), 1.0, ((*QI)->inv_corr_matrix),  gamma_N_vec, 1, 0, g_vec, 1);

	// Step 3, Compute mean; \mu = g * X
	mean = cblas_ddot(number_of_points, g_vec, 1, &((*QI)->trajectory_x[1]), 1); // Observe offset by one.
	
	// Step 4, Compute variance; \sigma^2 = 2*t - g*x => LOSS OF SIGNIFICANCE !
	sigma = (2.0 * pow(midtime,(2*hurst)));
	for(i = 0; i < number_of_points; i++)
	{
		sigma -= ( g_vec[i] * gamma_N_vec[i]);
	}

	// Step 5, Draw normal distributed midpoint
	midpoint = (mean + gsl_ran_gaussian_ziggurat(r, (sqrt(sigma)))); // the absolute value is only there for loss of significance errors

	// Step 6, Save new points, enlarge matrix and save new inverse correlation matrix
	double inv_sigma = (1/sigma);
	// First add sigma^{-2}*g*g^T on top of Q(N)
	cblas_dspr(CblasColMajor,CblasUpper, ((int) number_of_points), inv_sigma, g_vec, 1, (*QI)->inv_corr_matrix); 
	// Check if matrix needs to be enlarged.
	if( ((*QI)->size) >= ((*QI)->array_length) ){enlarge_array(QI);} 
	int new_index = ((*QI)->size); // The largest index so far was (*QI)->size , so now it\s one more
 	//Add new row with X_new (midpoint)
	(*QI)->trajectory_t[new_index+1] = midtime;
	(*QI)->trajectory_x[new_index+1] = midpoint;

	// The new column has coordinates new_index*(new_index+1)/2 , ..., (new_index + 1) * (new_index + 2)/2 - 1, that is the new column has new_index + 1 entries.
	for (i = (new_index*(new_index + 1)/2); i < (((new_index + 1)*(new_index + 2)/2) - 1); i ++)
	{
		// That's the column padded on the new row
		(*QI)->inv_corr_matrix[i] = -inv_sigma*g_vec[ i - (new_index*(new_index + 1)/2)];
	}
	(*QI)->inv_corr_matrix[i] = inv_sigma; // That's the new diagonal entry
	(*QI)->size = new_index + 1; // Enlarge size.
	
	/* Free malloc's */
	free(gamma_N_vec);
	free(g_vec);

	return midpoint;
}

double crossing_time_of_bridge(bridge_process* process)
{	
	// For a bridge that crosses the threshold, take the interesection of its linear interpolation with the threshold as first passage time.
	double fpt = 0.0; // Can never happen, as T = 1
	// Sanity check
	if( process->crossing_threshold == 1)
	{	
		fpt = ( (process->left_time) + (process->right_time - process->left_time)/(process->right_value - process->left_value)*(process->threshold - process->left_value));
	}
	return fpt;
}

double time_time_correlation(double ti, double tj, double hurst)
{
	if(hurst == 0.5){return (2*MIN(ti,tj));}
	else{return ( pow(fabs(ti),2*hurst) + pow(fabs(tj),2*hurst) - pow(fabs(ti-tj),2*hurst));}
}

void enlarge_array(triag_matrix** QI)
{
	long old_size = (*QI)->array_length;
	long new_size = ((long)(ARRAY_REALLOC_FACTOR * old_size));
	(*QI)->inv_corr_matrix = (double*) realloc( (*QI)->inv_corr_matrix, ((new_size*(new_size + 1)/2)*sizeof(double)) );
	(*QI)->trajectory_x = (double*) realloc( (*QI)->trajectory_x, ((new_size + 1) * sizeof(double)));
	(*QI)->trajectory_t = (double*) realloc( (*QI)->trajectory_t, ((new_size + 1) * sizeof(double)));
	(*QI)->array_length = new_size;
}

// Some useful functions for debugging. Not normally used

void print_QI(triag_matrix** QI)
{
	
	printf("===============QI============\nSIZE: %ld, ARRAY_LENGTH: %ld, ADDRESS: %p\n", (*QI)->size, (*QI)->array_length, *QI);
	long L = (*QI)->size;
	int i,j;
	for(i = 0; i <= L; i++)
	{
		printf("(%g,%g)\t",(*QI)->trajectory_t[i], (*QI)->trajectory_x[i]);	
	}
	printf("\n");
	printf("*********** Matrix values **************\n");
	for(i = 0; i < L; i++)
	{
		for(j = 0; j < i; j++)
		{
			printf("\t");
		}
		for(j = i; j < L; j++)
		{
			printf("%g\t",(*QI)->inv_corr_matrix[IJ2K(i,j)]);
		}
		printf("\n");
	}
	printf("***** Print QI End******\n");
}

void print_bridge(bridge_process** bridge)
{
	printf(" +++ Bridge +++ \n");
	printf("+Address: %p \n", *bridge);
	printf("+Root address: %p \n", (*bridge)->root_bridge);
	printf("+Previous root: %p \n", (*bridge)->previous_root);
	printf("+Node count : %i \n", ((*bridge)->root_bridge)->stack_top);
	int i, nodecount = ((*bridge)->root_bridge)->stack_top;
	printf("+Relatives: \n");
	for(i =0; i < nodecount; i++)
	{
		printf("->  %p [%g:%g] (%i)\n", ((*bridge)->root_bridge)->pointer_stack[i],  (((*bridge)->root_bridge)->pointer_stack[i])->left_time,  (((*bridge)->root_bridge)->pointer_stack[i])->right_time,  (((*bridge)->root_bridge)->pointer_stack[i])->generation);
	}
	printf(" +++++++++++ \n");
}

void copy_QI(triag_matrix** QI, double **Q, int last_point_index, double hurst, double** fracbm, double delta_t)
{
	int i,j;
	long array_length = 2*last_point_index; // Start off by making it twice as big, that should suffice for 90% of cases.
	*QI = malloc( sizeof(triag_matrix));;	
	(*QI)->size = last_point_index;
	(*QI)->array_length = array_length;
	(*QI)->inv_corr_matrix = (double*) malloc( (array_length * ( array_length+1) / 2) * sizeof(double));
	(*QI)->trajectory_x = (double*) malloc( (array_length+1) * sizeof(double));
	(*QI)->trajectory_t = (double*) malloc( (array_length+1) * sizeof(double));
	(*QI)->hurst_parameter = hurst;
	(*QI)->trajectory_x[0] = 0.0;
	(*QI)->trajectory_t[0] = 0.0;
	for(i = 0; i < last_point_index ; i++)
        {
		((*QI)->trajectory_x)[i+1] = (*fracbm)[i+1];
		((*QI)->trajectory_t)[i+1] = ((i+1)*delta_t);
		for(j = i; j < last_point_index ; j++)
		{
			((*QI)->inv_corr_matrix)[IJ2K(i,j)] = Q[last_point_index-1][IJ2K(i,j)]; // Here we are reading off the inverse correlation matrix from the previously found inverse correlation matrix.
		}
        }
}

void free_QI(triag_matrix** QI)
{
	free((*QI)->inv_corr_matrix);
	free((*QI)->trajectory_t);
	free((*QI)->trajectory_x);
	free(*QI);
}

void free_tree(bridge_process** bridge)
{
	// Frees a "tree", so a nested sequence of generated midpoints
	int nodecount = ((*bridge)->root_bridge)->stack_top;
	bridge_process *root=(*bridge)->root_bridge;
	int i;
	for(i = 1; i < nodecount; i++)
	{
		free_bridge(&(root->pointer_stack[i]));
	}
}

void free_bridge(bridge_process** bridge)
{
	free((*bridge)->pointer_stack);
	free(*bridge);
}
