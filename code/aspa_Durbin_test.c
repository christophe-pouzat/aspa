/** @file aspa_Durbin_test.c
 *  @brief User program for testing functions aspa_durbin_modification, aspa_Kolmogorov_D, aspa_cdf_Kplus 
 *
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/
#include "aspa.h"

int main()
{
  FILE *fp;
  fp = fopen("../data/traffic.dat","r");
  if (fp == NULL)
  {
    fprintf(stderr,"Can't open file ../data/traffic.dat\n");
    exit(EXIT_FAILURE);
  }
  size_t n=128;
  gsl_vector *traffic=gsl_vector_alloc(n);
  float datum;
  size_t idx=0;
  while (fscanf(fp,"%g",&datum)==1)
  {
    gsl_vector_set(traffic,idx,(double) datum);
    idx++;
  }
  fclose(fp);
  double CV=gsl_stats_sd(traffic->data,1,n)/gsl_stats_mean(traffic->data,1,n);
  for (size_t i=1; i<n; i++)
    gsl_vector_set(traffic,i,gsl_vector_get(traffic,i-1)+gsl_vector_get(traffic,i));
  gsl_vector_scale(traffic,1./gsl_vector_get(traffic,n-1));
  gsl_vector_view traf = gsl_vector_subvector(traffic,0,n-1);
  double K=aspa_Kolmogorov_D(&traf.vector,true,"D");
  double W2=aspa_AndersonDarling_W2(&traf.vector,true);
  printf("The coefficient of variation is: %g; the scaled Kolmogorov's statistic is: %g (%g);"
	 " the Anderson-Darling statistic is: %g (%g).\n",
	 CV,sqrt((double)n-1.)*K,aspa_cdf_K(n-1,K),W2, aspa_cdf_AD_P(n-1,W2));
  gsl_vector_free(traffic);

  fp = fopen("../data/nerve.dat","r");
  n=799;
  gsl_vector *nerve=gsl_vector_alloc(n);
  idx=0;
  while (fscanf(fp,"%g",&datum)==1)
  {
    gsl_vector_set(nerve,idx,(double) datum);
    idx++;
  }
  fclose(fp);
  CV=gsl_stats_sd(nerve->data,1,n)/gsl_stats_mean(nerve->data,1,n);
  for (size_t i=1; i<n; i++)
    gsl_vector_set(nerve,i,gsl_vector_get(nerve,i-1)+gsl_vector_get(nerve,i));
  gsl_vector_scale(nerve,1./gsl_vector_get(nerve,n-1));
  traf = gsl_vector_subvector(nerve,0,n-1);
  K=aspa_Kolmogorov_D(&traf.vector,true,"D");
  W2=aspa_AndersonDarling_W2(&traf.vector,true);
  printf("The coefficient of variation is: %g; the scaled Kolmogorov's statistic is: %g (%g);"
	 " the Anderson-Darling statistic is: %g (%g).\n",
	 CV,sqrt((double)n-1.)*K,aspa_cdf_K(n-1,K),W2, aspa_cdf_AD_P(n-1,W2));
  gsl_vector_free(nerve);

  fp = fopen("../data/mining_disasters.dat","r");
  n=111;
  int interval;
  gsl_vector *coal=gsl_vector_alloc(n);
  for (size_t i=0; i<n; i++)
  {
    if (fscanf(fp,"%d",&interval)==1)
      gsl_vector_set(coal,i,(double) interval);
  }
  fclose(fp);
  CV=gsl_stats_sd(coal->data,1,n)/gsl_stats_mean(coal->data,1,n);
  for (size_t i=1; i<n; i++)
    gsl_vector_set(coal,i,gsl_vector_get(coal,i-1)+gsl_vector_get(coal,i));
  gsl_vector_scale(coal,1./gsl_vector_get(coal,n-1));
  traf = gsl_vector_subvector(coal,0,n-1);
  K=aspa_Kolmogorov_D(&traf.vector,true,"D");
  W2=aspa_AndersonDarling_W2(&traf.vector,true);
  printf("The coefficient of variation is: %g; the scaled Kolmogorov's statistic is: %g (%g);"
	 " the Anderson-Darling statistic is: %g (%g).\n",
	 CV,sqrt((double)n-1.)*K,aspa_cdf_K(n-1,K),W2, aspa_cdf_AD_P(n-1,W2));
  gsl_vector_free(coal);
  return 0;
}
