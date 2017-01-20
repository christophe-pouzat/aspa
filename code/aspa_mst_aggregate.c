/** @file aspa_mst_aggregate.c
 *  @brief User program for aggregating multi trial spike trains
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include "aspa.h"

#include <string.h>
#include <getopt.h>

int read_args(int argc, char ** argv,
	      size_t * in_bin,
	      size_t * out_bin);

void print_usage();

int main(int argc, char ** argv)
{
  size_t in_bin,out_bin;
  int status = read_args(argc,argv,&in_bin,&out_bin);
  if (status == -1) exit (EXIT_FAILURE);
  aspa_sta * sta;
  if (in_bin == 0)
      sta = aspa_sta_fscanf(stdin);
  else
    sta = aspa_sta_fread(stdin);
  aspa_sta * asta = aspa_sta_aggregate(sta);
  aspa_sta_free(sta);
  if (out_bin == 0)
    aspa_sta_fprintf(stdout,asta,false);
  else
    aspa_sta_fwrite(stdout,asta,false);
  aspa_sta_free(asta);
  return 0;
}

/** @brief Reads command line arguments.
 *  
 *  @param[in] argc argument of main
 *  @param[in] argv argument of main
 *  @param[out] in_bin input format, O for "txt" 1 for "bin" (default 0)
 *  @param[out] out_bin output format, O for "txt" 1 for "bin" (default 0)
 *  @return 0 when everything goes fine
*/
int read_args(int argc, char ** argv,
	      size_t * in_bin,
	      size_t * out_bin)
{
  // Define default values
  *in_bin=0;
  *out_bin=0;
  {int opt;
    static struct option long_options[] = {
      {"in_bin",no_argument,NULL,'i'},
      {"out_bin",no_argument,NULL,'o'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"hio",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
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
	 "\n"
	 "Aggregates several trials into a single one.\n");
}
