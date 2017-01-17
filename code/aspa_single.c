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
 *  @return 0 if successful  
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
 *  @return a pointer to an allocated aspa_sta structure
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
 *  @return 0 if successful  
*/
int aspa_sta_fwrite(FILE * stream, const aspa_sta * sta, bool flat)
{
  if (flat == true) {
    // find out the total number of spikes
    size_t n_total = 0;
    for (size_t t_idx=0; t_idx < sta->n_trials; t_idx++)
    {
      gsl_vector * st = aspa_sta_get_st(sta,t_idx);
      n_total += st->size;
    }
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
 *  @return a pointer to an allocated aspa_sta structure
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

/** @brief Reads data from stdin, allocates and intializes an
 *         aspa_spike_train_data structure
 *
 *  The data entered in stdin are assumed to be organized
 *  in a single column (one spike time per line).
 *
 *  @param[in] sampling_frequency as its name says (in Hz)
 *  @param[in] inter_trial_interval as its name says (in s)
 *  @returns a pointer to an initialized aspa_spike_train_data
*/
aspa_spike_train_data * aspa_get_spike_train_data(double sampling_frequency,
						  double inter_trial_interval)
{
  size_t buffer_length = default_length;
  double *buffer=calloc(buffer_length, sizeof(double));
  size_t counter=0;
  char line[BUFSIZ];
  while (fgets (line, BUFSIZ, stdin))
  {
    buffer[counter] = atof(line)/sampling_frequency;
    counter++;
    if (counter>=buffer_length)
    {
      buffer_length*=2;
      buffer=realloc(buffer,buffer_length*sizeof(double));
    }
  }
  if (!feof(stdin))
  {
    fprintf (stderr, "Reading problem\n");
    exit (EXIT_FAILURE);
  }
  aspa_spike_train_data * res = malloc(sizeof(aspa_spike_train_data));
  res->inter_trial_interval = inter_trial_interval;
  res->n_spikes = counter;
  res->aggregate = 0;
  res->spike_train = malloc(counter*sizeof(double));
  for (size_t i=0; i<counter; i++)
    res->spike_train[i]=buffer[i];
  free(buffer);
  return res;
}

/** @brief Frees memory of aspa_spike_train_data structure
 * 
 *  @param st The aspa_spike_train_data to free
 *  @return 0 if everyhing goes fine.
*/
int aspa_spike_train_data_free(aspa_spike_train_data * st)
{
  free(st->spike_train);
  free(st);
  return 0;
}

/** @brief Creates an aspa_spike_train_data_array from an aspa_spike_train_data
 *
 *  The function splits the times of the spike_train member of its input into
 *  as many aspa_spike_train_data as there are trials where each new aspa_spike_train_data
 *  contains the offset data from a single trial
 *
 *  @param[in] st a pointer to an aspa_spike_train_data structure
 *  @returns a pointer to an initialized aspa_spike_train_data_array
*/
aspa_spike_train_data_array * aspa_get_spike_train_data_array(aspa_spike_train_data * st)
{
  aspa_spike_train_data_array * res = malloc(sizeof(aspa_spike_train_data_array));
  size_t n_spikes=st->n_spikes; // Number of spikes in spike train
  double iti = st->inter_trial_interval;
  // Find ot the number of trials
  double last_st = st->spike_train[0]; // Time of last spike
  size_t trial_idx=floor(last_st/iti); // In which trial is the current spike
  size_t n_trials=1; // The number of trials
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = st->spike_train[i]; // Time of spike i
    size_t current_idx = floor(current_st/iti);
    if (current_idx > trial_idx)
      n_trials++;
  }
  res->n_trials = n_trials;
  res->trial_start_time = malloc(n_trials*sizeof(double));
  res->spike_train_data = malloc(n_trials*sizeof(aspa_spike_train_data *));
  double current_train[n_spikes];
  size_t s_idx=0;
  for (size_t t_idx=0; t_idx<n_trials; t_idx++)
  {
    double current_st = st->spike_train[s_idx];
    size_t current_idx = floor(current_st/iti);
    res->trial_start_time[t_idx] = current_idx*iti;
    size_t within_index=0;
    while (floor(current_st/iti) == current_idx && s_idx<n_spikes)
    {
      current_train[within_index] = current_st-current_idx*iti;
      within_index++;
      s_idx++;
      current_st = st->spike_train[s_idx];
    }
    res->spike_train_data[t_idx] = malloc(sizeof(aspa_spike_train_data));
    res->spike_train_data[t_idx]->n_spikes=within_index;
    res->spike_train_data[t_idx]->inter_trial_interval=iti;
    res->spike_train_data[t_idx]->spike_train=malloc(within_index*sizeof(double));
    res->spike_train_data[t_idx]->aggregate=0;
    for (size_t i=0; i<within_index; i++)
      res->spike_train_data[t_idx]->spike_train[i] = current_train[i];
  }
  return res;
}

/** @brief Frees memory of aspa_spike_train_data_array structure
 * 
 *  @param sta The aspa_spike_train_data_array to free
 *  @return 0 if everyhing goes fine.
*/
int aspa_spike_train_data_array_free(aspa_spike_train_data_array * sta)
{
  free(sta->trial_start_time);
  for (size_t i=0; i<sta->n_trials; i++)
    free(sta->spike_train_data[i]);
  free(sta->spike_train_data);
  free(sta);
  return 0;
}

/** @brief Aggregates many trials of a spike train
 *
 *  @param[in] st The aspa_spike_train_data to aggregate
 *  @return a new "aggregated" aspa_spike_train_data if everyhing goes fine.
*/
//aspa_spike_train_data * aspa_aggregate_spike_train(const aspa_spike_train_data * st)
aspa_spike_train_data * aspa_aggregate_spike_train(const aspa_spike_train_data_array * sta)
{
  aspa_spike_train_data * res = malloc(sizeof(aspa_spike_train_data));
  size_t n_trials = sta->n_trials;
  size_t n_spikes=0; // Number of spikes in spike train
  for (size_t i=0; i<n_trials; i++)
    n_spikes+=(sta->spike_train_data[i])->n_spikes;
  //size_t n_spikes=st->n_spikes; // Number of spikes in spike train
  res->n_spikes = n_spikes;
  //double iti = st->inter_trial_interval;
  double iti = (sta->spike_train_data[0])->inter_trial_interval;
  res->inter_trial_interval = iti;
  res->spike_train = malloc(n_spikes*sizeof(double));
  size_t s_idx=0;
  for (size_t i=0; i<n_trials; i++)
  {
    for (size_t j=0; j<(sta->spike_train_data[i])->n_spikes; j++)
    {
      res->spike_train[s_idx] = (sta->spike_train_data[i])->spike_train[j];
      s_idx++;
    }
  }
  //size_t n_trials=0; // The number of trials
  //double last_st = st->spike_train[0]; // Time of last spike
  //size_t trial_idx=floor(last_st/iti); // In which trial is the current spike
  //res->spike_train[0] = st->spike_train[0]-iti*trial_idx;
  //n_trials++; 
  //for (size_t i=1; i<n_spikes; i++)
  //{
  //  double current_st = st->spike_train[i]; // Time of spike i
  //  size_t current_idx = floor(current_st/iti);
  //  if (current_idx > trial_idx)
  //  { // If the trial of spike i is not the one of the previous spike
  //    // We change trial_idx
  //    trial_idx = current_idx;
  //    // We increase the number of "observed" trials
  //    n_trials++;
  //  }
  //  res->spike_train[i] = current_st-iti*trial_idx;
  //  last_st=current_st;
  //}
  res->aggregate = n_trials;
  gsl_sort(res->spike_train,1,n_spikes);
  return res;
}

/** @brief Reads data from STREAM and allocates a gsl_vector
 *
 *  The data pointed to by STREAM are assumed to be organized
 *  in a single column (one spike time per line).
 *
 *  @param[in/out] STREAM pointer to an opened file or stdin
 *  @param[in] sampling_frequency as its name says (in Hz)
 *  @returns a pointer to an initialized gsl_vector
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
 *  @param[in] st a pointer to a aspa_spike_train_data containing the spike train
 *  @return A pointer to an allocated and initialized aspa_isi_data structure 
*/
aspa_isi_data * aspa_get_isi_data(aspa_spike_train_data *st)
{
  size_t n_spikes=st->n_spikes; // Number of spikes in spike train
  double inter_trial_interval=st->inter_trial_interval;
  size_t trial[n_spikes]; // Trial
  double isi[n_spikes]; // Array containing isi
  double left[n_spikes]; // Array containing spike time of the "left" spike of an isi
  size_t n_trials=0; // The number of trials
  double last_st = st->spike_train[0]; // Time of last spike
  size_t trial_idx=floor(last_st/inter_trial_interval); // In which trial is the current spike
  n_trials++; 
  size_t n_isi=0; // The number of isi
  for (size_t i=1; i<n_spikes; i++)
  {
    double current_st = st->spike_train[i]; // Time of spike i
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
