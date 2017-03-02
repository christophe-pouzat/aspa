/** @file aspa_cdf_K_test.c
 *  @brief User program for testing function aspa_cdf_K
 *
 *  Compares the output of G. Marsaglia, Wai Wan 
 *  Tsang and Jingbo Wong, [J.Stat.Software. 8(18): 1-4](https://www.jstatsoft.org/article/view/v008i18),
 *  with the some elements of the table of Z. Birnbaum (1952) JASA 47(229): 425-441
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

int main()
{
  char Birnbaum[] = "Birnbaum (1952)";
  char Marsaglia[] = "Marsaglia et al (2003)";
  char N[] = "N";
  char e[] = "e";
  printf("----------------------------------------------------------------------------------------------------------\n"
	 "Comparison between some table entries of\n"
	 "Birnbaum (1952) Numerical Tabulation of the Distribution of Kolmogorov's Statistic for Finite Sample Size,\n"
	 "JASA 47(229): 425-441, with the output of the code of G. Marsaglia, Wai Wan Tsang and Jingbo Wong (2003)\n"
	 "J. Stat. Software 8(18): 1-4\n"
	 "----------------------------------------------------------------------------------------------------------\n");
  printf("%5s %5s %23s %23s\n", N, e, Birnbaum,Marsaglia);
  printf("%5d %5d %23g %23g\n", 6, 4, 0.99623, aspa_cdf_K(6,4./6));
  printf("%5d %5d %23g %23g\n", 17, 3, 0.39630, aspa_cdf_K(17,3./17));
  printf("%5d %5d %23g %23g\n", 29, 9, 0.99441, aspa_cdf_K(29,9./29));
  printf("%5d %5d %23g %23g\n", 34, 6, 0.78663, aspa_cdf_K(34,6./34));
  printf("%5d %5d %23g %23g\n", 45, 13, 0.99919, aspa_cdf_K(45,13./45));
  printf("%5d %5d %23g %23g\n", 56, 8, 0.81552, aspa_cdf_K(56,8./56));
  printf("%5d %5d %23g %23g\n", 67, 3, 0.00154, aspa_cdf_K(67,3./67));
  printf("%5d %5d %23g %23g\n", 77, 14, 0.98936, aspa_cdf_K(77,14./77));
  printf("%5d %5d %23g %23g\n", 84, 6, 0.24247, aspa_cdf_K(84,6./84));
  printf("%5d %5d %23g %23g\n", 98, 10, 0.75771, aspa_cdf_K(98,10./98));
  printf("----------------------------------------------------------------------------------------------------------\n");
  return 0;
}
