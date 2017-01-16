/** @file aspa.h
 *  @brief Function prototypes for the aspa (analyse des sequences 
 *         de potentiels d'action) library.
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

gsl_vector * aspa_raw_fscanf(FILE * STREAM, double sampling_frequency);

/** @brief Structure holding arrays of gsl_vectors each vector containing
 *         a single trial spike train.
 *
 *  Since we work most of the time in a multi-trial setting,
 *  the structure is designed to hold arrays of gsl_vectors.
 *  sta stands for: spike train array.
*/
typedef struct
{
  size_t n_trials; //<! Number of trials
  double onset; //<! Stimulus onset time (s)
  double offset; //<! Stimulus offset time (s)
  double trial_duration; //<! Single trial duration (s)
  double * trial_start_time; //<! Vector holding the actual start time of each trial
  gsl_vector ** st; //<! The spike trains
} aspa_sta;

aspa_sta * aspa_sta_alloc(size_t n_trials, double onset, double offset, double trial_duration);

int aspa_sta_free(aspa_sta * sta);

gsl_vector * aspa_sta_get_st(const aspa_sta * sta, size_t st_index);

double aspa_sta_get_st_start(const aspa_sta * sta, size_t st_index);

int aspa_sta_set_st_start(aspa_sta * sta, size_t st_index, double time);

aspa_sta * aspa_sta_from_raw(gsl_vector * raw, double inter_trial_interval, double onset, double offset, double trial_duration);

int aspa_sta_fprintf(FILE * stream, const aspa_sta * sta, bool flat);

aspa_sta * aspa_sta_fscanf(FILE * STREAM);

int aspa_sta_fwrite(FILE * stream, const aspa_sta * sta, bool flat);

aspa_sta * aspa_sta_fread(FILE * STREAM);

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

/** @brief Structure holding arrays of aspa_spike_train_data.
 *
 *  Since we work most of the time in a multi-trial setting,
 *  the structure is designed to hold arrays of aspa_spike_train_data.
*/
typedef struct
{
  size_t n_trials; //<! Number of trials
  double * trial_start_time; //<! Vector holding the actual start time of each trial
  aspa_spike_train_data ** spike_train_data; //<! The spike train data
} aspa_spike_train_data_array;

aspa_spike_train_data_array * aspa_get_spike_train_data_array(aspa_spike_train_data * std);

int aspa_spike_train_data_array_free(aspa_spike_train_data_array * sta);

aspa_spike_train_data * aspa_aggregate_spike_train(const aspa_spike_train_data_array * sta);

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
