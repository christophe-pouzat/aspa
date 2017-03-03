/** @file aspa_cdf_norm_test.c
 *  @brief User program for testing functions aspa_cdf_norm_P and aspa_cdf_norm_Q
 *
 *  Compares the outputs of G. Marsaglia [J.Stat.Software. 11(4): 1-11](https://www.jstatsoft.org/article/view/v011i04),
 *  `Phi` and `cPhi` with the ones of `gsl_cdf_gaussian_P` and `gsl_cdf_gaussian_Q` 
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

int main()
{
  char gsl_P[] = "gsl_cdf_gaussian_P";
  char gsl_Q[] = "gsl_cdf_gaussian_Q";
  char Marsaglia_P[] = "Marsaglia (2004) Phi";
  char Marsaglia_Q[] = "Marsaglia (2004) cPhi";
  char X[] = "x";
  char diff_P[] = "Diff. gaussian_P-Phi";
  char diff_Q[] = "Diff. gaussian_Q-cPhi";
  printf("----------------------------------------------------------------------------------------------------------\n"
	 "Comparison between the Gaussian CDF and complementary CDF given by\n"
	 "functions gsl_cdf_gaussian_P and gsl_cdf_gaussian_Q of the GSL and\n"
	 "functions Phi and cPhi of G. Marsaglia [J. Stat. Software 11(4): 1-11]:\n"
	 "----------------------------------------------------------------------------------------------------------\n");
  printf("%5s %21s %21s %21s %21s %21s %21s\n", X, gsl_P, Marsaglia_P, diff_P, gsl_Q, Marsaglia_Q, diff_Q);
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 0.1, gsl_cdf_gaussian_P(0.1,1.), aspa_cdf_norm_P(0.1),
	 gsl_cdf_gaussian_P(0.1,1.)-aspa_cdf_norm_P(0.1),
	 gsl_cdf_gaussian_Q(0.1,1.), aspa_cdf_norm_Q(0.1),gsl_cdf_gaussian_Q(0.1,1.)-aspa_cdf_norm_Q(0.1));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 0.25, gsl_cdf_gaussian_P(0.25,1.), aspa_cdf_norm_P(0.25),
	 gsl_cdf_gaussian_P(0.25,1.)-aspa_cdf_norm_P(0.25),
	 gsl_cdf_gaussian_Q(0.25,1.), aspa_cdf_norm_Q(0.25),gsl_cdf_gaussian_Q(0.25,1.)-aspa_cdf_norm_Q(0.25));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 0.5, gsl_cdf_gaussian_P(0.5,1.), aspa_cdf_norm_P(0.5),
	 gsl_cdf_gaussian_P(0.5,1.)-aspa_cdf_norm_P(0.5),
	 gsl_cdf_gaussian_Q(0.5,1.), aspa_cdf_norm_Q(0.5),gsl_cdf_gaussian_Q(0.5,1.)-aspa_cdf_norm_Q(0.5));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 1., gsl_cdf_gaussian_P(1.,1.), aspa_cdf_norm_P(1.),
	 gsl_cdf_gaussian_P(1.,1.)-aspa_cdf_norm_P(1.),
	 gsl_cdf_gaussian_Q(1.,1.), aspa_cdf_norm_Q(1.),gsl_cdf_gaussian_Q(1.,1.)-aspa_cdf_norm_Q(1.));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 1.5, gsl_cdf_gaussian_P(1.5,1.), aspa_cdf_norm_P(1.5),
	 gsl_cdf_gaussian_P(1.5,1.)-aspa_cdf_norm_P(1.5),
	 gsl_cdf_gaussian_Q(1.5,1.), aspa_cdf_norm_Q(1.5),gsl_cdf_gaussian_Q(1.5,1.)-aspa_cdf_norm_Q(1.5));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 1.96, gsl_cdf_gaussian_P(1.96,1.), aspa_cdf_norm_P(1.96),
	 gsl_cdf_gaussian_P(1.96,1.)-aspa_cdf_norm_P(1.96),
	 gsl_cdf_gaussian_Q(1.96,1.), aspa_cdf_norm_Q(1.96),gsl_cdf_gaussian_Q(1.96,1.)-aspa_cdf_norm_Q(1.96));
  printf("%5g %21g %21g %21g %21g %21g %21g\n", 1.96, gsl_cdf_gaussian_P(3.,1.), aspa_cdf_norm_P(3.),
	 gsl_cdf_gaussian_P(3.,1.)-aspa_cdf_norm_P(3.),
	 gsl_cdf_gaussian_Q(3.,1.), aspa_cdf_norm_Q(3.),gsl_cdf_gaussian_Q(3.,1.)-aspa_cdf_norm_Q(3.));
  printf("----------------------------------------------------------------------------------------------------------\n");
  return 0;
}
