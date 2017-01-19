#include "aspa.h"

#include <string.h>
#include <getopt.h>

int read_args(int argc, char ** argv,
	      size_t * in_bin);

void print_usage();

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
  aspa_fns isi_fns = aspa_fns_get(isi);
  if (sta->n_aggregated == 1)
    fprintf(stdout,"Data from %d trials.\n", (int) sta->n_trials);
  else
    fprintf(stdout,"Data from %d aggregated trials.\n", (int) sta->n_aggregated);
  fprintf(stdout,"The mean rate is: %4g Hz.\n", aspa_sta_rate(sta));
  fprintf(stdout,"The inter spike interval statistics are:\n"),
  aspa_fns_fprintf(stdout,&isi_fns);
  double src = aspa_lagged_spearman(isi, 1);
  fprintf(stdout,"A 95%% confidence interval for the lag 1 Spearman rank correlation is: [%g,%g].\n",
	 src-1.96*0.6325/sqrt(isi->size-1),src+1.96*0.6325/sqrt(isi->size-1));
  gsl_vector_free(isi);
  aspa_sta_free(sta);
  return 0;
}

/** @brief Reads command line arguments.
 *  
 *  @param[in] argc argument of main
 *  @param[in] argv argument of main
 *  @param[out] in_bin input format, O for "txt" 1 for "bin" (default 0)
 *  @return 0 when everything goes fine
*/
int read_args(int argc, char ** argv,
	      size_t * in_bin)
{
  // Define default values
  *in_bin=0;
  {int opt;
    static struct option long_options[] = {
      {"in_bin",no_argument,NULL,'i'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"hi",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 'i': *in_bin=1;
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
	 "\n"
	 "Returns five number summary and additional stats.\n");
}
