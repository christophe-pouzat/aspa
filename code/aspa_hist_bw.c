/** @file aspa_hist_bw.c
 *  @brief User program for selecting an histogram binwidth by cross-validation
 *
 *  The cross-validation score of Mats Rudemo "Empirical Choice of Histograms and Kernel Density Estimators"
 *  _Scandinavian Journal of Statistics_ **9**:65-78, 1982.
 *  The program computes the cross-validation score (that approximates the integrated mean squared error)
 *  for a range of bin number and prints to the `stdout` on two columns the binwidth and score.
 *  The data are read from the `stdin` in text format. A log transformation can be applied before
 *  contructing the histograms.
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

#include <getopt.h>

int read_args(int argc, char ** argv,
	      size_t * use_log,
	      size_t * from,
	      size_t * to);

int main(int argc, char ** argv)
{
  size_t use_log, from, to;
  int status = read_args(argc,argv,&use_log,&from,&to);
  if (status == -1) exit (EXIT_FAILURE);

  int nb;
  fscanf(stdin,"%d",&nb); // Get the sample size
  fprintf(stderr,"Sample size: %d\n", nb);
  size_t n = (size_t) nb;
  // Check values of to and from
  if (to < from) {
    from = 2;
    to = n / 10;
  }
  fprintf(stderr,"Exploring now number of bins between %d and %d.\n",
	  (int) from, (int) to);
  if (use_log)
    fprintf(stderr,"Using a log transformation of the data.\n");
  
  double x[n];
  double y;
  fscanf(stdin,"%lg",&y); // Read first element
  if (use_log) { // If log transformation applied
    if (y < 0) { // Make sure it makes sense
      fprintf(stderr,"Negative number, cannot take the log!\n");
      return -1;
    }
    x[0] = log(y); // Do it
  } else {
    x[0] = y;
  }
  double xmin=x[0],xmax=x[0];
  size_t x_idx=1;  
  while (fscanf (stdin, "%lg", &y) == 1) {
    if (use_log) { // If log transformation applied
      if (y < 0) { // Make sure it makes sense
	fprintf(stderr,"Negative number, cannot take the log!\n");
	return -1;
      }
      x[x_idx] = log(y); // Do it
    } else {
      x[x_idx] = y;
    }
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
  fprintf(stderr,"Histograms will be built between %g and %g.\n",
	  xmin,xmax);
  gsl_histogram * hist;
  size_t mbest;
  double jbest;
  double norm_factor = 1.0/(double) n;
  for (size_t m=from; m<=to; m++) {
    hist = gsl_histogram_alloc(m);
    gsl_histogram_set_ranges_uniform(hist, xmin, xmax);
    for (size_t i=0; i<n; i++)
      gsl_histogram_increment(hist, x[i]);
    gsl_histogram_scale(hist,norm_factor);
    gsl_histogram_mul(hist,hist);
    double p_hat_square = gsl_histogram_sum(hist);
    double j_score = (2.0-(n+1.0)*p_hat_square)*m/(n-1.0);
    fprintf(stdout,"%d %g\n",(int) m, j_score);
    if (m == from) {
      mbest=m;
      jbest=j_score;
    } else {
      if (j_score < jbest) {
	mbest=m;
	jbest=j_score;
      }
    }
    gsl_histogram_free(hist);
  }
  fprintf(stdout,"\n");
  fprintf(stderr,"The best number of bins is: %d giving a score of %g\n",(int) mbest,jbest);
  return 0;
}

int read_args(int argc, char ** argv,
	      size_t * use_log,
	      size_t * from,
	      size_t * to)
{
  static char usage[] = \
    "usage: %s [-l --log] [-f --from=integer]\n"
    "          [-t --to=integer] [-h --help]\n\n"
    "  -l --log: should the log of the observations be used?\n"
    "  -f --from <positive integer>: the smallest number of bins to explore\n"
    "     (default set to 2).\n"
    "  -t --to <positive integer>: the largest number of bins to explore\n"
    "     (default set to the number of observations divided by 5).\n"
    "  -h --help: prints this message.\n"
    " The program reads data from the 'stdin' (in text format),\n"
    " the first line should contain the number of observations (integer)\n"
    " the following lines should contain the observations, one per line\n"
    " in decimal notation. If a log transformed of the data is requested\n"
    " it is applied first. Then a sequence of histograms with number of\n"
    " bins variying between from and to (inclusive) is contructed. For each\n"
    " histogram, the empirical bin probability (bin count divided by sample\n"
    " size), p_hat_i for each bin i is obtained and the score is computed as\n"
    " follows: (2-(n+1) * sum(p_hat_i^2))*m/(n-1) where m is the number of\n"
    " bins and where the summation is done over the bins. The best number\n"
    " of bins is the one giving the smallest score.\n"
    " The program prints to the 'stdout' the number of bins and the\n"
    " corresponding cross-validation scores on two columns.\n\n";
  // Define default values
  *use_log=0;
  *from=2;
  *to=1;
  {int opt;
    static struct option long_options[] = {
      {"log",no_argument,NULL,'l'},
      {"from",optional_argument,NULL,'f'},
      {"to",optional_argument,NULL,'t'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,"lhf:t:",long_options,\
			      &long_index)) != -1) {
      switch(opt) {
      case 'l': *use_log=1;
	break;
      case 'f': *from = (size_t) atoi(optarg);
	break;
      case 't': *to = (size_t) atoi(optarg);
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
