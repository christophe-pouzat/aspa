/** @file aspa_mst_plot.c
 *  @brief User program for plotting multi trials spike trains
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

#include <string.h>
#include <getopt.h>

#define NUM_WHAT 4
// Define allowed values for what
char * good_what[] = {"raster",
		      "cp_rt",
		      "cp_wt",
		      "cp_norm"};

int read_args(int argc, char ** argv,
	      size_t * in_bin,
	      char * what,
	      size_t * text);

void print_usage();

int main(int argc, char ** argv)
{
  size_t in_bin, text;
  char what[8];
  int status = read_args(argc,argv,&in_bin,what,&text);
  if (status == -1) exit (EXIT_FAILURE);
  aspa_sta * sta;
  if (in_bin == 0)
      sta = aspa_sta_fscanf(stdin);
  else
    sta = aspa_sta_fread(stdin);

  if (text == 0)
  { // Interactive use of gnuplot
    if (strcmp(what,good_what[0])==0)
      aspa_raster_plot_i(sta); // raster plot
    if (strcmp(what,good_what[1])==0)
      aspa_cp_plot_i(sta, true, false); // cp_rt
    if (strcmp(what,good_what[2])==0)
      aspa_cp_plot_i(sta, false, false); // cp_wt
    if (strcmp(what,good_what[3])==0)
    { // normalized counting process, cp_norm
      if (sta->n_aggregated == 1)
      { // must aggregate first
	aspa_sta * asta = aspa_sta_aggregate(sta);
	aspa_cp_plot_i(asta,false,true);
	aspa_sta_free(asta);
      }
      else
      {
	aspa_cp_plot_i(sta,true,true);
      }
    }
  }
  else
  { // Print to stdout
    if (strcmp(what,good_what[0])==0)
      aspa_raster_plot_g(stdout,sta); // raster plot
    if (strcmp(what,good_what[1])==0)
      aspa_cp_plot_g(stdout,sta, true, false); // cp_rt
    if (strcmp(what,good_what[2])==0)
      aspa_cp_plot_g(stdout,sta, false, false); // cp_wt
    if (strcmp(what,good_what[3])==0)
    { // normalized counting process, cp_norm
      if (sta->n_aggregated == 1)
      { // must aggregate first
	aspa_sta * asta = aspa_sta_aggregate(sta);
	aspa_cp_plot_g(stdout,asta,false,true);
	aspa_sta_free(asta);
      }
      else
      {
	aspa_cp_plot_g(stdout,sta,true,true);
      }
    }
  }
  aspa_sta_free(sta);
  return 0;
}

/** @brief Reads command line arguments.
 *  
 *  @param[in] argc argument of main
 *  @param[in] argv argument of main
 *  @param[out] in_bin input format, O for "txt" 1 for "bin" (default 0)
 *  @param[out] what the type of plot, one of:
 *              "raster", "cp_rt", "cp_wt", "cp_norm"
 *  @param[out] text output, O for interactive window 1 for "text" (default 0)
 *  @returns 0 when everything goes fine
*/
int read_args(int argc, char ** argv,
	      size_t * in_bin,
	      char * what,
	      size_t * text)
{
  // Define default values
  *in_bin=0;
  *text=0;
  strcpy(what,"cp_rt");
  {int opt;
    static struct option long_options[] = {
      {"in_bin",no_argument,NULL,'i'},
      {"text",no_argument,NULL,'t'},
      {"what",optional_argument,NULL,'w'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"hitw:",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 'w':
      {
	strcpy(what,optarg);
	size_t i;
	for (i=0; i<NUM_WHAT; i++)
	{
	  if (strcmp(what,good_what[i]) == 0) 
	    break;
	}
	if (i == NUM_WHAT)
	{
	  fprintf(stderr,"Unknown type of plot: %s\n",what);
	  return -1;
	}
      }
      break;
      case 'i': *in_bin=1;
	break;
      case 't': *text=1;
	break;
      case 'h': print_usage();
	return -1;
      default : print_usage();
	return -1;
      }
    }
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
	 "  --text: specify text output\n"
	 "  --what <string>: one of 'raster', 'cp_rt', 'cp_wt',\n"
	 "  'cp_norm', the type of plot (see bellow)\n"
	 "\n"
	 "An interactive lot is generated.\n"
	 "If what is set to 'raster' a raster plot is generated.\n"
	 "If what is set to 'cp_rt' the observed counting process\n"
	 "in 'real' time is generated, that is trial appear one after\n"
	 "the other.\n"
	 "If what is set to 'cp_wt' the observed counting processes\n"
	 "corresponding to each trial are shown on the 'within trial time.\n"
	 "If what is set to 'cp_norm' the normalized aggregated counting\n"
	 "process is displayed (normalization means here that the step size\n"
	 "due to each spike in each trial is 1/number of trials; in a sense\n"
	 "the 'mean' counting process is displayed).\n");
}
