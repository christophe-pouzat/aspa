#include "aspa.h"

int main()
{
  // Get spike train
  //FILE *fp = fopen("locust20010214_Spontaneous_1_tetB_u1.txt","r");
  //gsl_vector * st_flat = aspa_raw_fscanf(fp,15000);
  //fclose(fp);
  gsl_vector * st_flat = aspa_raw_fscanf(stdin,15000);
  aspa_sta * sta = aspa_sta_from_raw(st_flat, 30, 0, 0, 29);
  gsl_vector_free(st_flat);
  aspa_sta_fprintf(stdout, sta, true);
  aspa_sta_free(sta);
  return 0;
}
