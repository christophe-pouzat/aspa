/** @file aspa.h
 *  @brief Function prototypes for the aspa (analyse des sequences 
 *         de potentiels d'action) library.
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include <stdlib.h>
#include <string.h>
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_histogram.h>

gsl_vector * aspa_raw_fscanf(FILE * STREAM, double sampling_frequency);

/** @brief Structure holding arrays of gsl_vectors each vector containing
 *         a single trial spike train.
 *
 *  Since we work most of the time in a multi-trial setting,
 *  the structure is designed to hold arrays of gsl_vectors.
 *  Some of the analysis carried out on these objects will 
 *  require an "aggregation" of several trials (after aligning
 *  them on a "reference time" like the stimulus onset). We will
 *  keep track of this aggregation with the n_aggregated member.
 *  The latter will be 1 if no aggregation has been performed and
 *  will contain the number of aggregated trials otherwise.
 *  sta stands for: spike train array.
*/
typedef struct
{
  size_t n_trials; //!< Number of trials
  size_t n_aggregated; //!< Number of "real trials" aggregated per "trial" 
  double onset; //!< Stimulus onset time (s)
  double offset; //!< Stimulus offset time (s)
  double trial_duration; //!< Single trial duration (s)
  double * trial_start_time; //!< Vector holding the actual start time of each trial
  gsl_vector ** st; //!< The spike trains
} aspa_sta;

aspa_sta * aspa_sta_alloc(size_t n_trials, size_t n_aggregated, double onset, double offset, double trial_duration);

int aspa_sta_free(aspa_sta * sta);

gsl_vector * aspa_sta_get_st(const aspa_sta * sta, size_t st_index);

double aspa_sta_get_st_start(const aspa_sta * sta, size_t st_index);

int aspa_sta_set_st_start(aspa_sta * sta, size_t st_index, double time);

aspa_sta * aspa_sta_from_raw(gsl_vector * raw, double inter_trial_interval, double onset, double offset, double trial_duration);

int aspa_sta_fprintf(FILE * stream, const aspa_sta * sta, bool flat);

aspa_sta * aspa_sta_fscanf(FILE * STREAM);

size_t aspa_sta_n_spikes(const aspa_sta * sta);

size_t aspa_sta_n_spikes_max(const aspa_sta * sta);

double aspa_sta_rate(const aspa_sta * sta);

gsl_vector * aspa_sta_isi(const aspa_sta * sta);

int aspa_sta_fwrite(FILE * stream, const aspa_sta * sta, bool flat);

aspa_sta * aspa_sta_fread(FILE * STREAM);

aspa_sta * aspa_sta_aggregate(const aspa_sta * sta);

void aspa_cp_plot_i(const aspa_sta * sta, bool flat, bool normalized);

int aspa_cp_plot_g(FILE * STREAM, const aspa_sta * sta, bool flat, bool normalized);

void aspa_raster_plot_i(const aspa_sta * sta);

int aspa_raster_plot_g(FILE * STREAM, const aspa_sta * sta);

/** Structure holding basic sample summary statistics
*/
typedef struct
{
  size_t n; //!< sample size
  double min; //!< sample minimum
  double max; //!< sample maximum
  double upperq; //!< sample upper quartile
  double lowerq; //!< sample lower quartile
  double mean; //!< sample mean
  double median; //!< sample median
  double mad; //!< sample mad
  double var; //!< sample variance
} aspa_fns;

aspa_fns aspa_fns_get(const gsl_vector * data);

int aspa_fns_fprintf(FILE * STREAM, aspa_fns * fns);

double aspa_lagged_spearman(const gsl_vector * data, size_t lag);

double aspa_cdf_K(int n,double d);

double aspa_cdf_Kplus(int n,double d);

double aspa_Kolmogorov_D(gsl_vector * data, bool sorted, char * what);

double aspa_cdf_norm_P(double x);

double aspa_cdf_norm_Q(double x);

double aspa_AndersonDarling_W2(gsl_vector * data, bool sorted);

double aspa_cdf_ADinf_P(double z);

double aspa_cdf_AD_P(int n,double z);

double aspa_errfix(int n, double x);

double aspa_adinf(double z);

int aspa_durbin_modification(const gsl_vector * seq, gsl_vector * res);

void aspa_lagged_rank_plot_i(const aspa_sta * sta, size_t lag);

int aspa_lagged_rank_plot_g(FILE * STREAM, const aspa_sta * sta, size_t lag);
