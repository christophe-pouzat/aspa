#include "aspa.h"

int main()
{
  // Get spike train
  aspa_spike_train_data * st=aspa_get_spike_train_data(15000,30);
  printf("# The spike train contains %d spikes.\n",(int) st->n_spikes);
  printf("# The first spike time is: %g (s),\n",st->spike_train[0]);
  printf("# The last spike time is: %g (s).\n",st->spike_train[st->n_spikes-1]);
  /* // Get isi and number of trials */
  /* size_t n_t; */
  /* gsl_matrix *isir=get_isi_rank(st,30,&n_t); */
  /* printf("# There are %d trials.\n",(int) n_t); */
  /* double isi[isir->size2]; */
  /* for (size_t i=0; i<isir->size2; i++) */
  /*   isi[i] = gsl_matrix_get(isir,1,i); */
  aspa_isi_data * isid=aspa_get_isi_data(st);
  printf("# There are %d trials.\n",(int) isid->n_trials);
  printf("# The shortest isi is: %g (s),\n",gsl_stats_min(isid->isi,1,isid->n_isi));
  printf("# The largest isi is: %g (s),\n",gsl_stats_max(isid->isi,1,isid->n_isi));
  printf("# The mean isi is: %g (s) and its SD is %g (s).\n",
  	 gsl_stats_mean(isid->isi,1,isid->n_isi),
  	 gsl_stats_sd(isid->isi,1,isid->n_isi));
  double rank[isid->n_isi];
  for (size_t i=0; i<isid->n_isi; i++)
    rank[i] = (double) isid->rank[i];
  printf("# The lag 1 rank autocorrelation (+/- se) of the isi is: %g +/- %g.\n",
  	 gsl_stats_lag1_autocorrelation(rank,1,isid->n_isi),
  	 0.6325/sqrt((double) isid->n_isi-1.0));
  aspa_spike_train_data * ast=aspa_aggregate_spike_train(st);
  aspa_spike_train_data_free(st);
  printf("#\n");
  for (size_t i=0; i<ast->n_spikes; i++)
    printf("%g %d\n",ast->spike_train[i],(int) i);
  aspa_spike_train_data_free(ast);
  aspa_isi_data_free(isid);
  return 0;
}
