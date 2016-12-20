#include "aspa.h"

int main()
{
  gsl_vector * st=read_spike_train(stdin,15000);
  printf("The spike train contains %d spikes.\n",(int) st->size);
  printf("The first spike time is: %g (s),\n",gsl_vector_get(st,0));
  printf("The last spike time is: %g (s).\n",gsl_vector_get(st,st->size-1));
  size_t n_t;
  gsl_vector *isi=get_isi(st,30,&n_t);
  printf("There are %d trials.\n",(int) n_t);
  printf("The shortest isi is: %g (s),\n",gsl_vector_min(isi));
  printf("The largest isi is: %g (s),\n",gsl_vector_max(isi));
  printf("The mean isi is: %g (s) and its SD is %g (s).\n",
	 gsl_stats_mean(isi->data,1,isi->size),
	 gsl_stats_sd(isi->data,1,isi->size));
  gsl_vector_free(st);
  gsl_vector_free(isi);
  return 0;
}
