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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

/** @brief Structure holding spike trains related data.
 *
 *  Since we work most of the time in a multi-trial setting,
 *  the structure is designed to hold all the information
 *  allowing to go back and forth between actual spike time
 *  (the data stored in member spike_train) and within trial
 *  time.
*/
typedef struct
{
  size_t n_spikes; //<! Number of spikes
  double inter_trial_interval; //<! The inter trial interval (s)
  double * spike_train; //<! The spike train
  size_t aggregate; //<! Aggregation indicator (0 is no aggregation, number of aggregated trials otherwise)
} aspa_spike_train_data;

aspa_spike_train_data * aspa_get_spike_train_data(double sampling_frequency, double inter_trial_interval);

int aspa_spike_train_data_free(aspa_spike_train_data * st);
  
/** @brief Structure holding inter spike interval (isi)
 *         related data
 *
 *  Each isi is characterized by its value held in member isi, 
 *  the time of the spike forming the left boundary of the isi,
 *  held in member spike_time, the trial of origin, held in 
 *  member trial and its rank, held in member rank.
 *  The members spike_time, isi, trial and rank are all pointers
 *  to arrays of length n_isi.
 *  Member n_trials hold the number of trials  
*/
typedef struct
{
  size_t n_isi; //!< Number of isi
  size_t n_trials; //!< Number of trials
  double * spike_time; //!< Time of spike making the left boundary of an isi
  double * isi; //!< isi
  size_t * trial; //!< trial of origin
  size_t * rank; //!< isi rank
} aspa_isi_data;

aspa_isi_data * aspa_isi_data_alloc(size_t n_isi);

int aspa_isi_data_free(aspa_isi_data * isid);

gsl_vector * read_spike_train(FILE * STREAM, double sampling_frequency);

gsl_vector * get_isi(gsl_vector *st, double inter_trial_interval, size_t *n_trials);

gsl_matrix * get_isi_rank(gsl_vector *st, double inter_trial_interval, size_t *n_trials);

aspa_isi_data * aspa_get_isi_data(aspa_spike_train_data *st);
