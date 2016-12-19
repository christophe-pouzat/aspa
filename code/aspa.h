/** @file aspa.h
 *  @brief Function prototypes for the aspa (analyse des sequences 
 *         de potentiels d'action) library.
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_heapsort.h>

gsl_vector * read_spike_train(FILE * STREAM, double sampling_frequency);
