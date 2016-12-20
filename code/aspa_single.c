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

