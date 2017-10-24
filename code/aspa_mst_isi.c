/** @file aspa_mst_isi.c
 *  @brief User program for getting ISI from multi trials spike trains
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

#include <getopt.h>

int read_args(int argc, char ** argv,
	      size_t * in_bin);

int main(int argc, char ** argv)
{
  size_t in_bin;
  int status = read_args(argc,argv,&in_bin);
  if (status == -1) exit (EXIT_FAILURE);
  aspa_sta * sta;
  if (in_bin == 0)
      sta = aspa_sta_fscanf(stdin);
  else
    sta = aspa_sta_fread(stdin);

  gsl_vector * isi = aspa_sta_isi(sta);
  for (size_t i=0; i<isi->size; i++)
    fprintf(stdout,"%g\n",gsl_vector_get(isi,i));
  
  gsl_vector_free(isi);
  aspa_sta_free(sta);
  return 0;
}

int read_args(int argc, char ** argv,
	      size_t * in_bin)
{
  static char usage[] = \
    "usage: %s [-b --bin] [-h --help]\n\n"
    "  -b --bin: the data read from the 'stdin' are in binary format.\n"
    "  -h --help: prints this message.\n"
    " The program reads data from the 'stdin' (default in text format)\n"
    " most likely resulting from a call to 'aspa_read_spike_train',\n"
    " gets the inter spike intervals and writes them in text format\n"
    " to the stdout\n\n";
  // Define default values
  *in_bin=0;
  {int opt;
    static struct option long_options[] = {
      {"bin",no_argument,NULL,'b'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"hb",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 'b': *in_bin=1;
	break;
      case 'h': printf(usage,argv[0]);
	return -1;
      default : fprintf(stderr,usage,argv[0]);
	return -1;
      }
    }
  }
  return 0;
}
