/** @file aspa_read_spike_train.c
 *  @brief User program for reading spike train data
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include "aspa.h"

#include <string.h>
#include <getopt.h>

int read_args(int argc, char ** argv,
	      double * inter_trial_interval,
	      double * stim_onset,
	      double * stim_offset,
	      double * trial_duration,
	      double * sample2second,
	      size_t * in_bin,
	      size_t * out_bin);

void print_usage();

int main(int argc, char ** argv)
{
  size_t in_bin,out_bin;
  double inter_trial_interval=0;
  double stim_onset,stim_offset,sample2second;
  double trial_duration=0;
  int status = read_args(argc,argv,&inter_trial_interval,
			 &stim_onset,&stim_offset,&trial_duration,
			 &sample2second,&in_bin,&out_bin);
  if (status == -1) exit (EXIT_FAILURE);
  aspa_sta * sta;
  if (trial_duration > 0)
  { // Read flat test file with spike times one after the other
    gsl_vector * st_flat = aspa_raw_fscanf(stdin,sample2second);
    sta = aspa_sta_from_raw(st_flat, inter_trial_interval,
			    stim_onset, stim_offset,
			    trial_duration);
    gsl_vector_free(st_flat);
  }
  else
  {
    if (in_bin == 0)
      sta = aspa_sta_fscanf(stdin);
    else
      sta = aspa_sta_fread(stdin);
  }
  
  if (out_bin == 0)
    aspa_sta_fprintf(stdout,sta,false);
  else
    aspa_sta_fwrite(stdout,sta,false);
  
  aspa_sta_free(sta);
  return 0;
}

/** @brief Reads command line arguments.
 *  
 *  @param[in] argc argument of main
 *  @param[in] argv argument of main
 *  @param[out] inter_trial_interval the inter trial interval (in s)
 *  @param[out] stim_onset stimulus onset (in s), can do without, defaut 0
 *  @param[out] stim_offset stimulus offset (in s), can do without, defaut 0
 *  @param[out] trial_duration the duration of individual trials (in s)
 *  @param[out] sample2second sample to second correction 
 *              (sampling rate if data not already in s, default 1)
 *  @param[out] in_bin input format, "txt" or "bin" (default "txt")
 *  @param[out] out_bin output format, "txt" or "bin" (default "txt")
 *  @return 0 when everything goes fine
*/
int read_args(int argc, char ** argv,
	      double * inter_trial_interval,
	      double * stim_onset,
	      double * stim_offset,
	      double * trial_duration,
	      double * sample2second,
	      size_t * in_bin,
	      size_t * out_bin)
{
  // Define default values
  *stim_onset=0;
  *stim_offset=0;
  *sample2second=15000;
  *in_bin=0;
  *out_bin=0;
  {int opt;
    static struct option long_options[] = {
      {"in_bin",no_argument,NULL,'i'},
      {"out_bin",no_argument,NULL,'o'},
      {"sample2second",optional_argument,NULL,'s'},
      {"trial_duration",optional_argument,NULL,'d'},
      {"inter_trial_interval",optional_argument,NULL,'t'},
      {"stim_onset",optional_argument,NULL,'m'},
      {"stim_offset",optional_argument,NULL,'f'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"hios:d:t:m:f:",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 's':
      {
	float s=atof(optarg);
	if (s <= 0)
	{
	  fprintf(stderr,"The sample to seconds conversion factor should be > 0.\n");
	  return -1;
	}
	*sample2second = (double) s; 
      }
      break;
      case 'd':
      {
	float d=atof(optarg);
	if (d <= 0)
	{
	  fprintf(stderr,"Trial duration should be > 0.\n");
	  return -1;
	}
	*trial_duration = (double) d; 
      }
      break;
      case 't':
      {
	float t=atof(optarg);
	if (t <= 0)
	{
	  fprintf(stderr,"The inter trial interval should be > 0.\n");
	  return -1;
	}
	*inter_trial_interval = (double) t; 
      }
      break;
      case 'm':
      {
	float m=atof(optarg);
	*stim_onset = (double) m; 
      }
      break;
      case 'f':
      {
	float f=atof(optarg);
	*stim_offset = (double) f; 
      }
      break;
      case 'i': *in_bin=1;
	break;
      case 'o': *out_bin=1;
	break;
      case 'h': print_usage();
	return -1;
      default : print_usage();
	return -1;
      }
    }
  }
  if (*trial_duration > *inter_trial_interval)
  {
    fprintf(stderr,"Trial duration cannot be larger than inter trial interval.\n");
    return -1;
  }
  if (*stim_offset < *stim_onset)
  {
    fprintf(stderr,"Stim offset must be larger than stim onset.\n");
    return -1;
  }
  return 0;
}

/** @brief Prints usage to command line.
 * 
*/
void print_usage()
{
  printf("Usage: \n"
	 "  --in_bin: specify binary data input\n"
	 "  --out_bin: specify binary data output\n"
	 "  --sample2second <positive real>: the factor by which times\n"
	 "  in input data are divided in order get spike times in seconds\n"
	 "  used only when reading 'raw' data (default 15000)\n"
	 "  --inter_trial_interval <positive real>: the inter trial\n"
	 "  interval (in s) used only when reading 'raw' data\n"
	 "  --trial_duration <positive real>: the recorded duration\n"
	 "  (in s) of each trial used only when reading 'raw' data\n"
	 "  --stim_onset <real>: the stimulus onset time\n"
	 "  (in s) if that makes sense, used only when reading 'raw' data\n"
	 "  --stim_offset <real>: the stimulus offset time\n"
	 "  (in s) if that makes sense, used only when reading 'raw' data\n"
	 "\n");
}
