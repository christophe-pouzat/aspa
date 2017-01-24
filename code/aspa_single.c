/** @file aspa_single.c
 *  @brief Function definitions for single neuron spike trains
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include "aspa.h"
#define default_length 1000

/** @brief Reads data from stdin, allocates and intializes a
 *         gsl_vector
 *
 *  The data entered in stdin are assumed to be organized
 *  in a single column (one spike time per line).
 *
 *  @param[in] STREAM a pointer to an opened text file
 *  @param[in] sampling_frequency as its name says (in Hz)
 *  @returns a pointer to an initialized gsl_vector
*/
gsl_vector * aspa_raw_fscanf(FILE * STREAM,
			     double sampling_frequency)
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

/** @brief Allocates an aspa_sta
 *
 *  The function splits the times of its input into
 *  as many gsl_vector as there are trials where each new gsl_vector 
 *  contains the data from a single trial aligned on the stimulus 
 *  onset time.
 *
 *  @param[in] n_trials the number of trials 
 *  @param[in] onset stimulus onset time (in s)
 *  @param[in] offset stimulus offset time (in s)
 *  @param[in] trial_duration as its name says (in s)
 *  @returns a pointer to an allocated aspa_sta
*/
aspa_sta * aspa_sta_alloc(size_t n_trials, size_t n_aggregated, double onset, double offset, double trial_duration)
{
  aspa_sta * res = malloc(sizeof(aspa_sta));
  res->n_trials = n_trials;
  res->n_aggregated = n_aggregated;
  res->onset = onset;
  res->offset = offset;
  res->trial_duration = trial_duration;
  res->trial_start_time = malloc(n_trials*sizeof(double));
  res->st = malloc(n_trials*sizeof(gsl_vector *));
  return res;
}

/** @brief Frees an aspa_sta
 *
 *  @param[in/out] A pointer to an allocated aspa_sta structure
 *  @returns 0 if everything goes fine 
*/
int aspa_sta_free(aspa_sta * sta)
{
  free(sta->trial_start_time);
  for (size_t i=0; i < sta->n_trials; i++)
    gsl_vector_free(sta->st[i]);
  free(sta->st);
  free(sta);
  return 0;
}

/** @brief Returns a pointer to a given spike train of an 
 *         aspa_sta structure.
 *
 *  Range checking is performed, if argument st_index is too 
 *  large an error is generated.
 *  
 *  @param[in] sta A pointer to an allocated aspa_sta structure
 *  @param[in] st_index the index of the requested spike train
 *  returns a pointer to a gsl_vector
*/
gsl_vector * aspa_sta_get_st(const aspa_sta * sta, size_t st_index)
{
  assert (st_index < sta->n_trials);
  return sta->st[st_index];
}

/** @brief Returns "true" start time of a given spike train of an 
 *         aspa_sta structure.
 *
 *  Range checking is performed, if argument st_index is too 
 *  large an error is generated.
 *  
 *  @param[in] sta A pointer to an allocated aspa_sta structure
 *  @param[in] st_index the index of the requested spike train
 *  returns the trial start time
*/
double aspa_sta_get_st_start(const aspa_sta * sta, size_t st_index)
{
  assert (st_index < sta->n_trials);
  return sta->trial_start_time[st_index];
}

/** @brief Sets "true" start time of a given spike train of an 
 *         aspa_sta structure.
 *
 *  Range checking is performed, if argument st_index is too 
 *  large an error is generated.
 *  
 *  @param[in] sta A pointer to an allocated aspa_sta structure
 *  @param[in] st_index the index of the requested trial
 *  @param[in] time the onset of the trial
 *  returns 0 if everything goes fine
*/
int aspa_sta_set_st_start(aspa_sta * sta, size_t st_index, double time)
{
  assert (st_index < sta->n_trials);
  sta->trial_start_time[st_index]=time;
  return 0;
}

/** @brief Allocates and initializes an aspa_sta structure from
 *         a "flat" representation of the data
 *
 *  The first argument, raw, is assumed to contain the actual
 *  spike times of all the trials one after the other.
 *  The function splits the times of raw into
 *  as many gsl_vector as there are trials where each new gsl_vector
 *  contains the data from a single trial alligned on the stimulus 
 *  onset time.
 *
 *  @param[in] raw pointer to a gsl_vector containing the "flat" data
 *  @param[in] inter_trial_interval as its name says (in s)
 *  @param[in] onset stimulus onset time (in s)
 *  @param[in] offset stimulus offset time (in s)
 *  @param[in] trial_duration as its name says (in s)
 *  @returns a pointer to an initialized aspa_sta
*/
aspa_sta * aspa_sta_from_raw(gsl_vector * raw, double inter_trial_interval, double onset, double offset, double trial_duration)
{
  size_t n_spikes=raw->size; // Number of spikes in raw
  double iti = inter_trial_interval;
  // Find out the number of trials
  double last_st = gsl_vector_get(raw,0); // Time of last considered spike
  size_t trial_idx=floor(last_st/iti); // In which trial is the current spike
  size_t n_trials=1; // The number of trials
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = gsl_vector_get(raw,i); // Time of spike i
    size_t current_idx = floor(current_st/iti);
    if (current_idx > trial_idx)
    {
      trial_idx = current_idx;
      n_trials++;
    }
  }
  aspa_sta * res = aspa_sta_alloc(n_trials, 1, onset, offset, trial_duration);
  double current_train[n_spikes];
  size_t s_idx=0;
  for (size_t t_idx=0; t_idx<n_trials; t_idx++)
  {
    double current_st = gsl_vector_get(raw,s_idx);
    size_t current_idx = floor(current_st/iti);
    aspa_sta_set_st_start(res,t_idx,current_idx*iti);
    size_t within_index=0;
    while (floor(current_st/iti) == current_idx)
    {
      current_train[within_index] = current_st-current_idx*iti;
      within_index++;
      s_idx++;
      if (s_idx == n_spikes)
	break;
      current_st = gsl_vector_get(raw,s_idx);
    }
    res->st[t_idx] = gsl_vector_alloc(within_index);
    gsl_vector * st = aspa_sta_get_st(res,t_idx);
    for (size_t i=0; i<within_index; i++)
      gsl_vector_set(st,i,current_train[i]);
    if (s_idx == n_spikes)
	break;
  }
  return res;
}

/** @brief Prints to stream the content of an aspa_sta structure
 *
 *  The printing "format" is selected throught the boolean variable
 *  flat. If the latter is set to true, the spike times are written
 *  one after the other (one per line) and the times are the actual
 *  ones, not the "within trial" times. If flat is set to false, a 
 *  header whose first five lines start with:
 *  \# Number of trials:
 *  \# Number of aggregated trials:
 *  \# Stimulus onset:
 *  \# Stimulus offset:
 *  \# Single trial duration:
 *  is printed first. The spike times (within trial times) of each
 *  trial are printed next with two blank lines separating the
 *  different trials (suitable for gnuplot). Each trial starts with 
 *  the following three comments lines:
 *  \# Start of trial:
 *  \# Trial start time:
 *  \# Number of spikes:
 *  and ends with:
 *  \# End of trial:
 *
 *  @param[in/out] stream a pointer to an opened text file
 *  @param[in] sta pointer to the aspa_sta structure to be written
 *  @param[in] flat boolean indicator controlling what is written
 *  @returns0 if successful  
*/
int aspa_sta_fprintf(FILE * stream, const aspa_sta * sta, bool flat)
{
  if (flat == false) {
    fprintf(stream,"# Number of trials: %d\n",(int) sta->n_trials);
    fprintf(stream,"# Number of aggregated trials: %d\n",(int) sta->n_aggregated);
    fprintf(stream,"# Stimulus onset: %g (s)\n",sta->onset);
    fprintf(stream,"# Stimulus offset: %g (s)\n",sta->offset);
    fprintf(stream,"# Single trial duration: %g (s)\n",sta->trial_duration);
    fprintf(stream,"\n\n");
  }
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    if (flat == false) {
      fprintf(stream,"# Start of trial: %d\n",(int) t_idx);
      fprintf(stream,"# Trial start time: %g (s)\n", aspa_sta_get_st_start(sta,t_idx));
      fprintf(stream,"# Number of spikes: %d\n",(int) st->size);
      for (size_t s_idx=0; s_idx < st->size; s_idx++)
	fprintf(stream,"%g\n", gsl_vector_get(st,s_idx));
      fprintf(stream,"# End of trial: %d\n",(int) t_idx);
      fprintf(stream,"\n\n");
    } else {
      double t_start = aspa_sta_get_st_start(sta,t_idx);
      for (size_t s_idx=0; s_idx < st->size; s_idx++)
	fprintf(stream,"%g\n", gsl_vector_get(st,s_idx)+t_start);
    }
  }
  return 0;
}

/** @brief Scans multiple trials from a text file and return the
 *         result in an allocated pointer to an aspa_sta structure
 *
 *  The expected format of the input is assumed to follow the following
 *  layout: a header whose first five lines start with:
 *  \# Number of trials:
 *  \# Number of aggregated trials:
 *  \# Stimulus onset:
 *  \# Stimulus offset:
 *  \# Single trial duration:
 *  The spike times (within trial times) of each
 *  trial should be found next with two blank lines separating the
 *  different trials (suitable for gnuplot). Each trial starts with 
 *  the following three comments lines:
 *  \# Start of trial:
 *  \# Trial start time:
 *  \# Number of spikes:
 *  and ends with:
 *  \# End of trial:  
 *
 *  @param[in/out] stream a pointer to an opened text file
 *  @returnsa pointer to an allocated aspa_sta structure
*/
aspa_sta * aspa_sta_fscanf(FILE * STREAM)
{
  char buffer[256];
  char value[128];
  // Read line per line
  fgets(buffer, sizeof(buffer), STREAM);
  sscanf(buffer, "# Number of trials:  %127s", value);
  size_t n_trials = atoi(value);
  fgets(buffer, sizeof(buffer), STREAM);
  sscanf(buffer, "# Number of aggregated trials:  %127s", value);
  size_t n_aggregated = atoi(value);
  fgets(buffer, sizeof(buffer), STREAM);
  sscanf(buffer, "# Stimulus onset:  %127s", value);
  double onset = atof(value);
  fgets(buffer, sizeof(buffer), STREAM);
  sscanf(buffer, "# Stimulus offset:  %127s", value);
  double offset = atof(value);
  fgets(buffer, sizeof(buffer), STREAM);
  sscanf(buffer, "# Single trial duration:  %127s", value);
  double trial_duration = atof(value);
  aspa_sta * res = aspa_sta_alloc(n_trials, n_aggregated, onset, offset, trial_duration);
  for (size_t t_idx=0; t_idx < n_trials; t_idx++)
  {
    // Read two blank lines
    fgets(buffer, sizeof(buffer), STREAM);
    fgets(buffer, sizeof(buffer), STREAM);
    // Read line with trial number
    fgets(buffer, sizeof(buffer), STREAM);
    // Read line with trial start time
    fgets(buffer, sizeof(buffer), STREAM);
    sscanf(buffer, "# Trial start time:  %127s", value);
    aspa_sta_set_st_start(res,t_idx,(double) atof(value));
    // Read line with the number of spikes
    fgets(buffer, sizeof(buffer), STREAM);
    sscanf(buffer, "# Number of spikes:  %127s", value);
    size_t n_spikes = atoi(value);
    // Allocate spike times vector
    res->st[t_idx] = gsl_vector_alloc(n_spikes);
    // Loop over the spike times
    for (size_t s_idx=0; s_idx < n_spikes; s_idx++)
    {
      float spike_time;
      fgets(buffer, sizeof(buffer), STREAM);
      sscanf(buffer,"%f",&spike_time);
      //fscanf(STREAM,"%f",&spike_time);
      gsl_vector_set(res->st[t_idx],s_idx,(double) spike_time);
    }
    // Read line with trial number
    fgets(buffer, sizeof(buffer), STREAM);
  }
  return res;
}

/** @brief Returns the total number of spikes contained in 
 *         an aspa_sta structure
 *
 *  @param[in] sta a pointer to an aspa_sta structure
 *  @result the total number of spikes
*/
size_t aspa_sta_n_spikes(const aspa_sta * sta)
{
  size_t n_total = 0;
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    n_total += st->size;
  }
  return n_total;
}

/** @brief Returns means number of spikes per second
 *
 *  @aparam[in] sta a pointer to an aspa_sta structure
 *  @returnsthe mean spike rate
*/
double aspa_sta_rate(const aspa_sta * sta)
{
  double total_obs_time = sta->n_trials*sta->trial_duration;
  return aspa_sta_n_spikes(sta)/total_obs_time/sta->n_aggregated;
}

/** @brief Return a gsl_vector containing the inter spike intervals
 *         (ISI) of an aspa_sta structure.
 *
 *  The ISI from each trial are obtained and put together, one after
 *  the other.
 *
 *  @param[in] sta a pointer to an apsa_sta structure
 *  @returnsa pointer to a gsl_vector with the ISI
*/
gsl_vector * aspa_sta_isi(const aspa_sta * sta)
{
  gsl_vector * isi = gsl_vector_alloc(aspa_sta_n_spikes(sta)-sta->n_trials);
  size_t isi_idx=0;
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    double present = gsl_vector_get(st,0);
    for (size_t i=0; i < (st->size-1); i++)
    {
      double next = gsl_vector_get(st,i+1);
      gsl_vector_set(isi,isi_idx,next-present);
      present=next;
      isi_idx++;
    }
  }
  return isi;
}

/** @brief Prints in binary to stream the content of an aspa_sta structure
 *
 *  What is printed is selected throught the boolean variable
 *  flat. If the latter is set to true, the number of spikes is written first
 *  as a size_t followed by the spike times 
 *  one after the other--the times are the actual
 *  ones, not the "within trial" times--. 
 *  If flat is set to false, a size_t with the number of trials is written first
 *  followed by the number of aggregated trials (size_t), the stimulus onset, 
 *  offset and the single trial duration as doubles.
 *  Then, for each trial, a trial start time (double), the number of spikes in the 
 *  trial (size_t) followed by the within trials spike times.
 *
 *  @param[in/out] stream a pointer to an opened text file
 *  @param[in] sta pointer to the aspa_sta structure to be written
 *  @param[in] flat boolean indicator controlling what is written
 *  @returns0 if successful  
*/
int aspa_sta_fwrite(FILE * stream, const aspa_sta * sta, bool flat)
{
  if (flat == true) {
    // find out the total number of spikes
    size_t n_total = aspa_sta_n_spikes(sta);
    fwrite(&n_total,sizeof(size_t),1,stream);
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      double t_start = aspa_sta_get_st_start(sta,t_idx);
      double spike_time;
      for (size_t s_idx=0; s_idx < st->size; s_idx++)
      {
	spike_time = gsl_vector_get(st,s_idx)+t_start; 
	fwrite(&spike_time,sizeof(double),1,stream);
      }
    }
  } else {
    fwrite(&(sta->n_trials),sizeof(size_t),1,stream);
    fwrite(&(sta->n_aggregated),sizeof(size_t),1,stream);
    fwrite(&(sta->onset),sizeof(double),1,stream);
    fwrite(&(sta->offset),sizeof(double),1,stream);
    fwrite(&(sta->trial_duration),sizeof(double),1,stream);
  }
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    double t_start = aspa_sta_get_st_start(sta,t_idx);
    fwrite(&t_start,sizeof(double),1,stream);
    fwrite(&(st->size),sizeof(size_t),1,stream);
    gsl_vector_fwrite(stream,st);
  }
  return 0;
}

/** @brief Reads multiple trials from a binary file and return the
 *         result in an allocated pointer to an aspa_sta structure
 *
 *  The expected format of the input is assumed to follow the following
 *  layout: 
 *  a size_t with the number of trials followed a size_t with the number
 *  of aggregated trials, by the stimulus onset, 
 *  offset and the single trial duration as doubles.
 *  Then, for each trial, a trial start time (double), the number of spikes in the 
 *  trial (size_t) followed by the within trials spike times.
 *
 *  @param[in/out] stream a pointer to an opened text file
 *  @returnsa pointer to an allocated aspa_sta structure
*/
aspa_sta * aspa_sta_fread(FILE * STREAM)
{
  size_t n_trials;
  fread(&n_trials, sizeof(size_t),1,STREAM);
  size_t n_aggregated;
  fread(&n_aggregated, sizeof(size_t),1,STREAM);
  double onset;
  fread(&onset, sizeof(double),1,STREAM);
  double offset;
  fread(&offset, sizeof(double),1,STREAM);
  double trial_duration;
  fread(&trial_duration, sizeof(double),1,STREAM);
  aspa_sta * res = aspa_sta_alloc(n_trials, n_aggregated, onset, offset, trial_duration);
  for (size_t t_idx=0; t_idx < n_trials; t_idx++)
  {
    double start_time;
    fread(&start_time, sizeof(double),1,STREAM);
    aspa_sta_set_st_start(res,t_idx,start_time);
    size_t n_spikes;
    fread(&n_spikes, sizeof(size_t),1,STREAM);
    // Allocate spike times vector
    res->st[t_idx] = gsl_vector_alloc(n_spikes);
    gsl_vector * st = aspa_sta_get_st(res,t_idx);
    gsl_vector_fread(STREAM,st);
  }
  return res;
}

/** @brief Aggregates many trials of a spike train
 *
 *  @param[in] sta pointer to the aspa_sta to aggregate
 *  @returnsa pointer to new "aggregated" aspa_sta if everyhing goes fine.
*/
aspa_sta * aspa_sta_aggregate(const aspa_sta * sta)
{
  aspa_sta * res = aspa_sta_alloc(1, sta->n_trials, sta->onset, sta->offset, sta->trial_duration);
  aspa_sta_set_st_start(res,0,aspa_sta_get_st_start(sta,0));
  size_t n_trials = sta->n_trials;
  res->st[0] = gsl_vector_alloc(aspa_sta_n_spikes(sta));
  gsl_vector * rst = aspa_sta_get_st(res,0);
  size_t s_idx=0;
  for (size_t t_idx=0; t_idx<n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    for (size_t i=0; i < st->size; i++)
    {
      gsl_vector_set(rst,s_idx,gsl_vector_get(st,i));
      s_idx++;
    }
  }
  gsl_sort_vector(rst);
  return res;
}

/** @brief Plots the observed counting process assiociated with
 *         an aspa_sta structure
 *
 *  The "observed counting process" (OCP) of an aspa_sta structure can
 *  be built in several ways: 
 *   - If the process is not aggregated the counts vs real time can
 *   be displayed (a single OCP is then shown). (param flat -> true)
 *   - If the process is not aggregated the counts vs within trial
 *   time can be displayed (as many OCP as trials are then shown).
 *   - If the process as been aggregated the mean OCP where the 
 *   count increase due to each spike is 1 / number_of_trials can
 *   be displayed (param normalized -> true).
 *
 *  @param[in] sta a pointer to the sta structure
 *  @param[in] flat boolean controlling if actual or within trial
 *             time is used
 *  @param[in] normalized boolean controlling if the mean OCP is
 *             is displayed 
 *  @returnsnothing the function is only used for its side effect
*/
void aspa_cp_plot_i(const aspa_sta * sta, bool flat, bool normalized)
{
  FILE *gp=NULL;
  if (!gp)
    gp = popen("gnuplot -persist","w");
  if (!gp) {
    printf("Couldn't open Gnuplot.\n");
    return;
  }
  fprintf(gp,"set term qt; set grid; unset key\n");
  fprintf(gp,"set xlabel 'Time (s)'\n");
  if (normalized == true || flat == true) {
    if (flat == true)
    {
      fprintf(gp,"set title 'Observed counting process'\n");
      fprintf(gp,"set ylabel 'Events count'\n");
    } else {
      fprintf(gp,"set title 'Observed mean counting process'\n");
      fprintf(gp,"set ylabel 'Mean events count'\n");
    }
    fprintf(gp," plot '-' u 1:2 with steps\n");
    double step = 1.0/sta->n_aggregated;
    double s_idx = step;
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      double start_time = aspa_sta_get_st_start(sta,t_idx);
      for (size_t i=0; i < st->size; i++)
      {
	fprintf(gp,"%g %g\n", gsl_vector_get(st,i)+start_time, s_idx);
	s_idx+=step;
      }
    }
  } else {
    fprintf(gp,"set title 'Observed counting processes'\n");
    fprintf(gp,"set ylabel 'Events count'\n");
    fprintf(gp," plot '-' u 1:2 with steps\n");
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      int s_idx = 1;
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      for (size_t i=0; i < st->size; i++)
      {
	fprintf(gp,"%g %d\n", gsl_vector_get(st,i), (int) s_idx);
	s_idx++;
      }
      fprintf(gp,"\n\n");
    }
  }
  fprintf(gp,"e\n");
  fflush(gp);
  pclose(gp);
}

/** @brief Writes the observed counting process assiociated with
 *         an aspa_sta structure to a file in gnuplot friendly format
 *
 *  The "observed counting process" (OCP) of an aspa_sta structure can
 *  be built in several ways: 
 *   - If the process is not aggregated the counts vs real time can
 *   be displayed (a single OCP is then shown). (param flat -> true)
 *   - If the process is not aggregated the counts vs within trial
 *   time can be displayed (as many OCP as trials are then shown).
 *   - If the process as been aggregated the mean OCP where the 
 *   count increase due to each spike is 1 / number_of_trials can
 *   be displayed (param normalized -> true).
 *
 *  @param[in/out] STREAM an open file
 *  @param[in] sta a pointer to the sta structure
 *  @param[in] flat boolean controlling if actual or within trial
 *             time is used
 *  @param[in] normalized boolean controlling if the mean OCP is
 *             is displayed 
 *  @returns0 if everything goes fine
*/
int aspa_cp_plot_g(FILE * STREAM, const aspa_sta * sta, bool flat, bool normalized)
{
  if (normalized == true || flat == true) {
    double step = 1.0/sta->n_aggregated;
    double s_idx = step;
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      double start_time = aspa_sta_get_st_start(sta,t_idx);
      for (size_t i=0; i < st->size; i++)
      {
	fprintf(STREAM,"%g %g\n", gsl_vector_get(st,i)+start_time, s_idx);
	s_idx+=step;
      }
    }
  } else {
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      int s_idx = 1;
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      for (size_t i=0; i < st->size; i++)
      {
	fprintf(STREAM,"%g %d\n", gsl_vector_get(st,i), (int) s_idx);
	s_idx++;
      }
      fprintf(STREAM,"\n\n");
    }
  }
  return 0;
}

/** @brief Generates a raster plot from an aspa_sta structure
 *
 *  Gnuplot is called. An interactive window pops up.
 *
 *  @param[in] sta a pointer to an aspa_sta structure
 *  @returns nothing the function is only used for its side effect
*/
void aspa_raster_plot_i(const aspa_sta * sta)
{
  FILE *gp=NULL;
  if (!gp)
    gp = popen("gnuplot -persist","w");
  if (!gp) {
    printf("Couldn't open Gnuplot.\n");
    return;
  }
  fprintf(gp,"set term qt; set grid; unset key\n");
  fprintf(gp,"set xlabel 'Time (s)'\n");
  fprintf(gp,"set ylabel 'Trial'\n");
  fprintf(gp," plot '-' u 1:2 with dots\n");
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    for (size_t i=0; i < st->size; i++)
      fprintf(gp,"%g %d\n", gsl_vector_get(st,i), (int) t_idx+1);
    fprintf(gp,"\n\n");
  }
  fprintf(gp,"e\n");
  fflush(gp);
  pclose(gp);
}

/** @brief Writes a raster plot from an aspa_sta structure
 *         to a file in a gnuplot friendly format
 *
 *  @param[in/out] STREAM an open file
 *  @param[in] sta a pointer to an aspa_sta structure
 *  @returns 0 if everything goes fine
*/
int aspa_raster_plot_g(FILE * STREAM, const aspa_sta * sta)
{
  for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
  {
    gsl_vector * st = aspa_sta_get_st(sta,t_idx);
    for (size_t i=0; i < st->size; i++)
      fprintf(STREAM,"%g %d\n", gsl_vector_get(st,i), (int) t_idx+1);
    fprintf(STREAM,"\n\n");
  }
  return 0;
}

/** @brief Compute five number summary of a gsl_vector
 *         as well as some other basic statistics
 *
 *  The [five number summary](https://en.wikipedia.org/wiki/Five-number_summary)
 *  as well as the [median absolute deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation)
 *  (MAD), the mean and variance are computed.
 *
 *  @param[in] data a pointer to a `gsl_vector`
 *  @returnsan `aspa_fns` structure
*/
aspa_fns aspa_fns_get(const gsl_vector * data)
{
  size_t n = data->size;
  gsl_vector * tmp = gsl_vector_alloc(n);
  gsl_vector_memcpy(tmp,data);
  gsl_sort_vector(tmp);
  double median = gsl_stats_median_from_sorted_data(tmp->data,1,n);
  double upperq = gsl_stats_quantile_from_sorted_data(tmp->data,1,n,0.75);
  double lowerq = gsl_stats_quantile_from_sorted_data(tmp->data,1,n,0.25);
  double min = gsl_vector_get(tmp,0);
  double max = gsl_vector_get(tmp,n-1);
  double mean = gsl_stats_mean(tmp->data,1,n);
  double var = gsl_stats_variance(tmp->data,1,n);
  gsl_vector_add_constant(tmp,-median);
  for (size_t i=0; i<n; i++)
    gsl_vector_set(tmp,i,fabs(gsl_vector_get(tmp,i)));
  gsl_sort_vector(tmp);
  double mad = 1.4826*gsl_stats_median_from_sorted_data(tmp->data,1,n);
  gsl_vector_free(tmp);
  return (aspa_fns) {.n=n,.mean=mean,.min=min,
      .max=max,.upperq=upperq,.lowerq=lowerq,
      .median=median,.mad=mad,.var=var};
}

/** @brief Prints to stream the content of an `aspa_fns` structure
 *
 *  @param[in/out] `STREAM` a pointer to an open file
 *  @param[in] `fns` an `aspa_fns` structure
 *  @returns0 if everything goes fine
*/
int aspa_fns_fprintf(FILE * STREAM, aspa_fns * fns)
{
  fprintf(STREAM,"The sample contains %d elements.\n",
	  (int) fns->n);
  fprintf(STREAM,"The mean and SD are   : %5.4f and %5.4f.\n",
	  fns->mean, sqrt(fns->var));
  fprintf(STREAM,"The median and MAD are: %5.4f and %5.4f.\n",
	  fns->median, fns->mad);
  fprintf(STREAM,"The five number summary:\n");
  fprintf(STREAM,"Min.   1st qrt Median 3rd qrt Max. \n");
  fprintf(STREAM,"%6.4f %6.4f  %6.4f %6.4f  %6.4f\n",
	 fns->min,fns->lowerq,fns->median,fns->upperq,fns->max);
  return 0;
}

/** @brief Computes lagged spearman correlation from a gsl_vector
 *
 *  The [Spearman rank correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
 *  is just the Pearson correlation on the ranks.
 *
 *  @param[in] data a pointer to a gsl_vector
 *  @param[in] lag the lag at which the correlation is computed
 *  @returnsthe correlation coefficient
*/
double aspa_lagged_spearman(const gsl_vector * data, size_t lag)
{
  size_t n = data->size;
  assert (lag < n); // make sure the lag is small enough
  gsl_vector_const_view lagged = gsl_vector_const_subvector(data,lag,n-lag);
  double work[2*(n-lag)];
  return gsl_stats_spearman(data->data,1,(&lagged.vector)->data,1,n-lag,work);
}

/** @brief Generates a lagged rank plot
 *
 *  The data are ranked and the rank of the (i+l)th
 *  isi is plotted agains the rank of the ith isi,
 *  where l is the lag.
 *  Gnuplot is called. An interactive window pops up.
 *
 *  @param[in] sta a pointer to an aspa_sta structure
 *  @param[in] lag the lag 
 *  @returns nothing the function is only used for its side effect
*/
void aspa_lagged_rank_plot_i(const aspa_sta * sta, size_t lag)
{
  gsl_vector * isi = aspa_sta_isi(sta);
  size_t n = isi->size;
  assert (lag < n); // make sure the lag is small enough
  gsl_permutation * perm = gsl_permutation_alloc(n);
  gsl_permutation * rank = gsl_permutation_alloc(n);
  gsl_sort_vector_index(perm,isi);
  gsl_permutation_inverse(rank,perm);
  FILE *gp=NULL;
  if (!gp)
    gp = popen("gnuplot -persist","w");
  if (!gp) {
    printf("Couldn't open Gnuplot.\n");
    return;
  }
  fprintf(gp,"set term qt; set grid; unset key\n");
  fprintf(gp,"set xlabel 'Rank i'\n");
  fprintf(gp,"set ylabel 'Rank i + %d'\n",(int) lag);
  fprintf(gp," plot '-' u 1:2 with dots\n");
  for (size_t i=0; i < n-lag-1; i++)
    fprintf(gp,"%d %d\n", (int) rank->data[i], (int) rank->data[i+lag]);
  fprintf(gp,"e\n");
  fflush(gp);
  pclose(gp);
  gsl_permutation_free(perm);
  gsl_permutation_free(rank);
  gsl_vector_free(isi);
}
