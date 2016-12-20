#include "aspa.h"

int main()
{
  // Get spike train
  gsl_vector * st=read_spike_train(stdin,15000);
  printf("The spike train contains %d spikes.\n",(int) st->size);
  printf("The first spike time is: %g (s),\n",gsl_vector_get(st,0));
  printf("The last spike time is: %g (s).\n",gsl_vector_get(st,st->size-1));
  // Get isi and number of trials
  size_t n_t;
  gsl_vector *isi=get_isi(st,30,&n_t);
  printf("There are %d trials.\n",(int) n_t);
  printf("The shortest isi is: %g (s),\n",gsl_vector_min(isi));
  printf("The largest isi is: %g (s),\n",gsl_vector_max(isi));
  printf("The mean isi is: %g (s) and its SD is %g (s).\n",
	 gsl_stats_mean(isi->data,1,isi->size),
	 gsl_stats_sd(isi->data,1,isi->size));
  // Get lag 1 autocorrelation of the ranks
  // Check Sec. 12.4 of the GSL manual
  gsl_permutation * perm = gsl_permutation_alloc(isi->size);
  gsl_permutation * rank = gsl_permutation_alloc(isi->size);
  gsl_sort_vector_index(perm,isi);
  gsl_permutation_inverse(rank,perm);
  double rank_d[isi->size]; // a double version of rank
  for (size_t i=0; i<isi->size; i++)
    rank_d[i] = (double) gsl_permutation_get(rank,i);
  printf("The lag 1 rank autocorrelation (+/- se) of the isi is: %g +/- %g.\n",
	 gsl_stats_lag1_autocorrelation(rank_d,1,isi->size),
	 0.6325/sqrt((double) isi->size-1.0));
  gsl_vector_free(st);
  gsl_vector_free(isi);
  gsl_permutation_free(perm);
  gsl_permutation_free(rank);
  return 0;
}
