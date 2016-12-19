#include "aspa.h"

int main()
{
  gsl_vector * st=read_spike_train(stdin,15000);
  printf("The spike train contains %d spikes.\n",(int) st->size);
  printf("The first spike time is: %g (s),\n",gsl_vector_get(st,0));
  printf("The last spike time is: %g (s).\n",gsl_vector_get(st,st->size-1));
  gsl_vector_free(st);
  return 0;
}
