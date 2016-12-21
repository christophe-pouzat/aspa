#include "aspa.h"
#define default_length 1000

/** @brief Reads data from STREAM and allocates a gsl_vector
 *
 *  The data pointed to by STREAM are assumed to be organized
 *  in a single column (one spike time per line).
 *
 *  @param[in/out] STREAM pointer to an opened file or stdin
 *  @param[in] sampling_frequency as its name says (in Hz)
 *  @returns a pointer ot an initialized gsl_vector
*/
gsl_vector * read_spike_train(FILE * STREAM, double sampling_frequency)
{
  size_t buffer_length = default_length;
  double *buffer=calloc(buffer_length, sizeof(double));
  size_t counter=0;
  char line[BUFSIZ];
  while (fgets (line, BUFSIZ, STREAM))
  {
    buffer[counter] = atof(line)/sampling_frequency;
    counter++;
    if (counter>=buffer_length)
    {
      buffer_length*=2;
      buffer=realloc(buffer,buffer_length*sizeof(double));
    }
  }
  if (!feof(STREAM))
  {
    fprintf (stderr, "Reading problem\n");
    exit (EXIT_FAILURE);
  }
  gsl_vector * res = gsl_vector_alloc(counter);
  for (size_t i=0; i<counter; i++)
    gsl_vector_set(res,i,buffer[i]);
  free(buffer);
  return res;
}


/** @brief Allocates memory for an isi_data structure
 * 
 *  @param n_isi The number of isi.
 *  @return A pointer to an allocated isi_data structure. 
*/
aspa_isi_data * aspa_isi_data_alloc(size_t n_isi)
{
  aspa_isi_data * res = malloc(sizeof(aspa_isi_data));
  res->n_isi = n_isi;
  res->spike_time = malloc(n_isi*sizeof(double));
  res->isi = malloc(n_isi*sizeof(double));
  res->trial = malloc(n_isi*sizeof(size_t));
  res->rank = malloc(n_isi*sizeof(size_t));
  return res;
}

/** @brief Frees memory of isi_data structure
 * 
 *  @param isid The isi_data to free
 *  @return 0 if everyhing goes fine.
*/
int aspa_isi_data_free(aspa_isi_data * isid)
{
  free(isid->spike_time);
  free(isid->isi);
  free(isid->trial);
  free(isid->rank);
  free(isid);
  return 0;
}

/** @brief Get a vector of inter spike intervals (isi) from a spike train
 *
 *  The only subtlety is the potential presence of many trials assumed to 
 *  be regularly spaced every inter_trial_interval seconds. If such is the
 *  case the last spike of trial i and the first of trial i+1 should not
 *  be used to construct an isi.
 *
 *  @param[in] st a pointer to a gsl_vector containing the spike train
 *  @param[in] inter_trial_interval the elapsed time between trials
 *  @param[out] n_trials a pointer to a size_t, the number of trials
 *  @return A pointer to an allocated gsl_vector 
*/
gsl_vector * get_isi(gsl_vector *st, double inter_trial_interval, size_t *n_trials)
{
  size_t n_spikes=st->size; // Number of spikes in spike train
  double isi[n_spikes]; // Array containing isi
  *n_trials=0; // The number of trials
  double last_st = gsl_vector_get(st,0); // Time of last spike
  size_t trial_idx=floor(last_st/inter_trial_interval); // In which trial is the current spike
  (*n_trials)++; 
  size_t n_isi=0; // The number of isi
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = gsl_vector_get(st,i); // Time of spike i
    size_t current_idx = floor(current_st/inter_trial_interval);
    if (current_idx > trial_idx)
    { // If the trial of spike i is not the one of the previous spike
      // We change trial_idx
      trial_idx = current_idx;
      // We increase the number of "observed" trials
      (*n_trials)++;
    } else {
      // The current and previous spike give an isi value
      isi[n_isi] = current_st-last_st;
      n_isi++;
    }
    last_st=current_st;
  }
  gsl_vector *res=gsl_vector_alloc(n_isi);
  for (size_t i=0; i<n_isi; i++)
    gsl_vector_set(res,i,isi[i]);
  return res;
}

/** @brief Get a matrix whose first row contains the spike time of the left
 *         boundary of an isi, the second row contains the isi and the third
 *         row contains the rank of the isi.
 *
 *  The only subtlety is the potential presence of many trials assumed to 
 *  be regularly spaced every inter_trial_interval seconds. If such is the
 *  case the last spike of trial i and the first of trial i+1 should not
 *  be used to construct an isi.
 *
 *  @param[in] st a pointer to a gsl_vector containing the spike train
 *  @param[in] inter_trial_interval the elapsed time between trials
 *  @param[out] n_trials a pointer to a size_t, the number of trials
 *  @return A pointer to an allocated gsl_matrix 
*/
gsl_matrix * get_isi_rank(gsl_vector *st, double inter_trial_interval, size_t *n_trials)
{
  size_t n_spikes=st->size; // Number of spikes in spike train
  double isi[n_spikes]; // Array containing isi
  double left[n_spikes]; // Array containing spike time of the "left" spike of an isi
  *n_trials=0; // The number of trials
  double last_st = gsl_vector_get(st,0); // Time of last spike
  size_t trial_idx=floor(last_st/inter_trial_interval); // In which trial is the current spike
  (*n_trials)++; 
  size_t n_isi=0; // The number of isi
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = gsl_vector_get(st,i); // Time of spike i
    size_t current_idx = floor(current_st/inter_trial_interval);
    if (current_idx > trial_idx)
    { // If the trial of spike i is not the one of the previous spike
      // We change trial_idx
      trial_idx = current_idx;
      // We increase the number of "observed" trials
      (*n_trials)++;
    } else {
      // The current and previous spike give an isi value
      isi[n_isi] = current_st-last_st;
      left[n_isi] = last_st;
      n_isi++;
    }
    last_st=current_st;
  }
  gsl_matrix *res=gsl_matrix_alloc(3,n_isi);
  for (size_t i=0; i<n_isi; i++)
  {
    // Column 0 contains the time of the spike making the
    // left boundary of the isi
    gsl_matrix_set(res,0,i,left[i]);
    // Column 1 contains the isi
    gsl_matrix_set(res,1,i,isi[i]);
  }
  // Column 3 of res will contain the rank of the isi
  // We obtain the rank following Sec. 12.4 of the GSL manual 
  gsl_vector_const_view c1 = gsl_matrix_const_row(res,1);
  gsl_permutation * perm = gsl_permutation_alloc(n_isi);
  gsl_permutation * rank = gsl_permutation_alloc(n_isi);
  gsl_sort_vector_index(perm,&c1.vector);
  gsl_permutation_inverse(rank,perm);
  for (size_t i=0; i<n_isi; i++)
    gsl_matrix_set(res,2,i,(double) gsl_permutation_get(rank,i));
  gsl_permutation_free(perm);
  gsl_permutation_free(rank);
  return res;
}

/** @brief Allocates and initializes an aspa_isi_data structure.
 *
 *  The only subtlety is the potential presence of many trials assumed to 
 *  be regularly spaced every inter_trial_interval seconds. If such is the
 *  case the last spike of trial i and the first of trial i+1 should not
 *  be used to construct an isi.
 *
 *  @param[in] st a pointer to a gsl_vector containing the spike train
 *  @param[in] inter_trial_interval the elapsed time between trials
 *  @return A pointer to an allocated and initialized aspa_isi_data structure 
*/
aspa_isi_data * aspa_get_isi_data(gsl_vector *st, double inter_trial_interval)
{
  size_t n_spikes=st->size; // Number of spikes in spike train
  size_t trial[n_spikes]; // Trial
  double isi[n_spikes]; // Array containing isi
  double left[n_spikes]; // Array containing spike time of the "left" spike of an isi
  size_t n_trials=0; // The number of trials
  double last_st = gsl_vector_get(st,0); // Time of last spike
  size_t trial_idx=floor(last_st/inter_trial_interval); // In which trial is the current spike
  n_trials++; 
  size_t n_isi=0; // The number of isi
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = gsl_vector_get(st,i); // Time of spike i
    size_t current_idx = floor(current_st/inter_trial_interval);
    if (current_idx > trial_idx)
    { // If the trial of spike i is not the one of the previous spike
      // We change trial_idx
      trial_idx = current_idx;
      // We increase the number of "observed" trials
      n_trials++;
    } else {
      // The current and previous spike give an isi value
      isi[n_isi] = current_st-last_st;
      left[n_isi] = last_st;
      trial[n_isi] = trial_idx;
      n_isi++;
    }
    last_st=current_st;
  }
  aspa_isi_data * res = aspa_isi_data_alloc(n_isi);
  res->n_trials = n_trials;
  for (size_t i=0; i<n_isi; i++)
  {
    res->spike_time[i] = left[i];
    res->isi[i] = isi[i];
    res->trial[i] = trial[i];
  }
  // We obtain the rank following Sec. 12.4 of the GSL manual
  gsl_vector_const_view isi_gv = gsl_vector_const_view_array(res->isi,n_isi);
  gsl_permutation * perm = gsl_permutation_alloc(n_isi);
  gsl_permutation * rank = gsl_permutation_alloc(n_isi);
  gsl_sort_vector_index(perm,&isi_gv.vector);
  gsl_permutation_inverse(rank,perm);
  for (size_t i=0; i<n_isi; i++)
    res->rank[i] = gsl_permutation_get(rank,i);
  gsl_permutation_free(perm);
  gsl_permutation_free(rank);
  return res;
}
