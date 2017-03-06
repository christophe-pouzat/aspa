/** @file aspa_cdf_K_test.c
 *  @brief User program for testing function aspa_cdf_K
 *
 *  Compares the output of G. Marsaglia, Wai Wan 
 *  Tsang and Jingbo Wong, [J.Stat.Software. 8(18): 1-4](https://www.jstatsoft.org/article/view/v008i18),
 *  with the some elements of the table of Z. Birnbaum (1952) JASA 47(229): 425-441 and reproduce table 1
 *  of Birnbaum and Tingey (1951) One-sided confidence contours for probability distribution functions
 *  _The Annals of Mathematical Statistics_ __22__: 592-596.
 * 
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
  char diff[] = "Difference B-M";
  printf("----------------------------------------------------------------------------------------------------------\n"
	 "Comparison between some table entries of\n"
	 "Birnbaum (1952) Numerical Tabulation of the Distribution of Kolmogorov's Statistic for Finite Sample Size,\n"
	 "JASA 47(229): 425-441, with the output of the code of G. Marsaglia, Wai Wan Tsang and Jingbo Wong (2003)\n"
	 "J. Stat. Software 8(18): 1-4\n"
	 "----------------------------------------------------------------------------------------------------------\n");
  printf("%5s %5s %23s %23s %23s\n", N, e, Birnbaum, Marsaglia, diff);
  printf("%5d %5d %23g %23g %23g\n", 6, 4, 0.99623, aspa_cdf_K(6,4./6), 0.99623-aspa_cdf_K(6,4./6));
  printf("%5d %5d %23g %23g %23g\n", 17, 3, 0.39630, aspa_cdf_K(17,3./17), 0.39630-aspa_cdf_K(17,3./17));
  printf("%5d %5d %23g %23g %23g\n", 29, 9, 0.99441, aspa_cdf_K(29,9./29), 0.99441-aspa_cdf_K(29,9./29));
  printf("%5d %5d %23g %23g %23g\n", 34, 6, 0.78663, aspa_cdf_K(34,6./34), 0.78663-aspa_cdf_K(34,6./34));
  printf("%5d %5d %23g %23g %23g\n", 45, 13, 0.99919, aspa_cdf_K(45,13./45), 0.99919-aspa_cdf_K(45,13./45));
  printf("%5d %5d %23g %23g %23g\n", 56, 8, 0.81552, aspa_cdf_K(56,8./56), 0.81552-aspa_cdf_K(56,8./56));
  printf("%5d %5d %23g %23g %23g\n", 67, 3, 0.00154, aspa_cdf_K(67,3./67), 0.00154-aspa_cdf_K(67,3./67));
  printf("%5d %5d %23g %23g %23g\n", 77, 14, 0.98936, aspa_cdf_K(77,14./77), 0.98936-aspa_cdf_K(77,14./77));
  printf("%5d %5d %23g %23g %23g\n", 84, 6, 0.24247, aspa_cdf_K(84,6./84), 0.24247-aspa_cdf_K(84,6./84));
  printf("%5d %5d %23g %23g %23g\n", 98, 10, 0.75771, aspa_cdf_K(98,10./98), 0.75771-aspa_cdf_K(98,10./98));
  printf("----------------------------------------------------------------------------------------------------------\n\n\n");
  printf("----------------------------------------------------------------------------------------------------------\n"
	 "Comparison between aspa_cdf_Kplus output and table 1 of Birnbaum and Tingey (1951)\n\n"
	 "We start with table 1\n\n");
  printf("%5s %7g %7g %7g %g7\n", N, 0.1, 0.05, 0.01, 0.001);
  printf("%5d %7g %7g %7g %g7\n", 5, 0.4470, 0.5094, 0.6271, 0.7480);
  printf("%5d %7g %7g %7g %g7\n", 8, 0.3583, 0.4096, 0.5065, 0.6130);
  printf("%5d %7g %7g %7g %g7\n", 10, 0.3226, 0.3687, 0.4566, 0.5550);
  printf("%5d %7g %7g %7g %g7\n", 20, 0.23155, 0.26473, 0.3285, 0.4018);
  printf("%5d %7g %7g %7g %g7\n", 40, 0.16547, 0.18913, 0.2350, 0.2877);
  printf("%5d %7g %7g %7g %g7\n\n", 50, 0.14840, 0.16959, 0.2107, 0.2581);
  printf("Now the output of aspa_cdf_Kplus (doing the inverse)\n\n");
  printf("%5d %7g %7g %7g %g7\n", 5, aspa_cdf_Kplus(5,0.4470),
	 aspa_cdf_Kplus(5,0.5094), aspa_cdf_Kplus(5,0.6271),
	 aspa_cdf_Kplus(5,0.7480));
  printf("%5d %7g %7g %7g %g7\n", 8, aspa_cdf_Kplus(8,0.3583),
	 aspa_cdf_Kplus(8,0.4096), aspa_cdf_Kplus(8,0.5065),
	 aspa_cdf_Kplus(8,0.6130));
  printf("%5d %7g %7g %7g %g7\n", 10, aspa_cdf_Kplus(10,0.3226),
	 aspa_cdf_Kplus(10,0.3687), aspa_cdf_Kplus(10,0.4566),
	 aspa_cdf_Kplus(10,0.5550));
  printf("%5d %7g %7g %7g %g7\n", 20, aspa_cdf_Kplus(20,0.23155),
	 aspa_cdf_Kplus(20,0.26473), aspa_cdf_Kplus(20,0.3285),
	 aspa_cdf_Kplus(20,0.4018));
  printf("%5d %7g %7g %7g %g7\n", 40, aspa_cdf_Kplus(40,0.16547),
	 aspa_cdf_Kplus(40,0.18913), aspa_cdf_Kplus(40,0.2350),
	 aspa_cdf_Kplus(40,0.2877));
  printf("%5d %7g %7g %7g %g7\n\n", 50, aspa_cdf_Kplus(50,0.14840),
	 aspa_cdf_Kplus(50,0.16959), aspa_cdf_Kplus(50,0.2107),
	 aspa_cdf_Kplus(50,0.2581));
  return 0;
}
