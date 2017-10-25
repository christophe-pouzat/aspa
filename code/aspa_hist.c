/** @file aspa_hist.c
 *  @brief User program for creating an histogram
 *
 *  The data are read from the `stdin` in text format. A log transformation can be applied before
 *  contructing the histogram. The resulting histogram is written to the `stdout`.
 *  The log transformation should be used when dealing with positive random variables since
 *  their densities are often skewed to the right.
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

#include <getopt.h>

int read_args(int argc, char ** argv,
	      size_t * use_log,
	      size_t * n_bins,
	      size_t * prob);

int main(int argc, char ** argv)
{
  size_t use_log, n_bins, prob;
  int status = read_args(argc,argv,&use_log,&n_bins,&prob);
  if (status == -1) exit (EXIT_FAILURE);

  int nb;
  fscanf(stdin,"%d",&nb); // Get the sample size
  fprintf(stderr,"Sample size: %d\n", nb);
  size_t n = (size_t) nb;
  if (use_log)
    fprintf(stderr,"Using a log transformation of the data.\n");
  
  double x[n];
  double y;
  fscanf(stdin,"%lg",&y); // Read first element
  x[0] = y;
  double xmin=x[0],xmax=x[0];
  size_t x_idx=1;  
  while (fscanf (stdin, "%lg", &y) == 1) { // Read the other elements
    x[x_idx] = y;
    if (x[x_idx] < xmin)
      xmin = x[x_idx];
    if (x[x_idx] > xmax)
      xmax = x[x_idx];
    x_idx+=1;
  }
  if (x_idx != n) {
    fprintf(stderr,"The number of elements read %d is wrong.\n", (int) x_idx);
    return -1;
  }
  xmin = xmin-DBL_EPSILON;
  xmax = xmax+DBL_EPSILON;
  gsl_histogram * hist= gsl_histogram_alloc(n_bins);
  double range[n_bins+1];
  if (use_log) {
    if (xmin <= 0) {
      fprintf(stderr,"Negative values, cannot use log transform!\n");
      return -1;
    }
    double delta = (log(xmax)-log(xmin))/n_bins;
    range[0] = log(xmin);
    for (size_t i=1; i<n_bins+1; i++)
      range[i] = range[i-1]+delta;
    for (size_t i=0; i<n_bins+1; i++)
      range[i] = exp(range[i]);
  } else {
    double delta = (xmax-xmin)/n_bins;
    range[0] = xmin;
    for (size_t i=1; i<n_bins+1; i++)
      range[i] = range[i-1]+delta;
  }
  gsl_histogram_set_ranges(hist, range, n_bins+1);
  for (size_t i=0; i<n; i++)
    gsl_histogram_increment(hist, x[i]);

  if (prob) { // Normalize the histogram
    for (size_t bin_idx=0; bin_idx<n_bins; bin_idx++)
      hist->bin[bin_idx] /= n*(hist->range[bin_idx+1]-hist->range[bin_idx]);
  }
  gsl_histogram_fprintf (stdout, hist, "%g", "%g");
  gsl_histogram_free(hist);
  return 0;
}

int read_args(int argc, char ** argv,
	      size_t * use_log,
	      size_t * n_bins,
	      size_t * prob)
{
  static char usage[] = \
    "usage: %s -n --n_bins=integer [-l --log] \n"
    "          [-p --prob] [-h --help]\n\n"
    "  -n --n_bins <positive integer>: the number of bins to use.\n"
    "  -l --log: should the log of the observations be used?\n"
    "  -p --prob: should the result be normalized so that the\n"
    "     the histogram integral is one?\n"
    "  -h --help: prints this message.\n"
    " The program reads data from the 'stdin' (in text format),\n"
    " the first line should contain the number of observations (integer)\n"
    " the following lines should contain the observations, one per line\n"
    " in decimal notation. If a log transformed of the data is requested\n"
    " it is applied first. Then a histogram with 'n_bins' bins is contructed.\n"
    " If a log transformation was used, the bins boundaries are transformed\n"
    " back to the original scale. If argument 'prob' is specified, the histogram\n"
    " is normalized such that its integral is one (it becomes of proper PDF\n"
    " estimator), otherwise the bin counts are kept.\n"
    " The program prints to the 'stdout' on 3 columns: the left bin boundary;\n"
    " the right bin boundary; the bin count or bin frequency (depending on the\n"
    " specification of argument 'prob').\n\n";
  // Define default values
  *use_log=0;
  *prob=0;
  {int opt;
    static struct option long_options[] = {
      {"log",no_argument,NULL,'l'},
      {"n_bins",required_argument,NULL,'n'},
      {"prob",optional_argument,NULL,'p'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"lhn:p",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 'l': *use_log=1;
	break;
      case 'n': *n_bins = (size_t) atoi(optarg);
	break;
      case 'p': *prob = 1;
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
