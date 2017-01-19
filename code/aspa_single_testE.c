#include "aspa.h"

int main()
{
  // Get spike train
  FILE *fp = fopen("locust20010214_Spontaneous_1_tetB_u1.txt","r");
  gsl_vector * st_flat = aspa_raw_fscanf(fp,15000);
  //printf("Just read %d times.\n",(int) st_flat->size);
  fclose(fp);
  // make an aspa_sta structure
  aspa_sta * sta = aspa_sta_from_raw(st_flat, 30, 0, 0, 29);
  gsl_vector_free(st_flat);
  gsl_vector * isi = aspa_sta_isi(sta);
  aspa_fns isi_fns = aspa_fns_get(isi);
  printf("The mean rate is: %4g Hz.\n", aspa_sta_rate(sta));
  printf("The inter spike interval statistics are:\n"),
  aspa_fns_fprintf(stdout,&isi_fns);
  double src = aspa_lagged_spearman(isi, 1);
  printf("A 95%% confidence interval for the lag 1 Spearman rank correlation is: [%g,%g].\n",
	 src-1.96*0.6325/sqrt(isi->size-1),src+1.96*0.6325/sqrt(isi->size-1));
  gsl_vector_free(isi);
  printf("\n\nDoing now the same thing on the aggregated train.\n");
  aspa_sta * asta = aspa_sta_aggregate(sta);
  isi = aspa_sta_isi(asta);
  isi_fns = aspa_fns_get(isi);
  printf("The mean rate is: %4g Hz.\n", aspa_sta_rate(asta));
  printf("The inter spike interval statistics are:\n"),
  aspa_fns_fprintf(stdout,&isi_fns);
  src = aspa_lagged_spearman(isi, 1);
  printf("A 99%% confidence interval for the lag 1 Spearman rank correlation is: [%g,%g].\n",
	 src-2.56*0.6325/sqrt(isi->size-1),src+2.56*0.6325/sqrt(isi->size-1));
  gsl_vector_free(isi);
  aspa_sta_free(sta);
  aspa_sta_free(asta);
  return 0;
}

