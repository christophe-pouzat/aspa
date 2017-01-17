#include "aspa.h"

int main()
{
  // Get spike train
  FILE *fp = fopen("locust20010214_Spontaneous_1_tetB_u1.txt","r");
  gsl_vector * st_flat = aspa_raw_fscanf(fp,15000);
  fclose(fp);
  // make an aspa_sta structure
  aspa_sta * sta = aspa_sta_from_raw(st_flat, 30, 0, 0, 29);
  gsl_vector_free(st_flat);
  // aggregate
  aspa_sta * asta = aspa_sta_aggregate(sta);
  aspa_sta_free(sta);
  // write aggregated structure to file
  aspa_sta_fprintf(stdout,asta,false);
  aspa_sta_free(asta);
  return 0;
}
