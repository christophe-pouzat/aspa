#include "aspa.h"
#define default_length 1000

/** @brief Reads data from STREAM and allocates a gsl_vector
 *
 *  The data pointed to by STREAM are assumed to be organized
 *  in a single column (one spike time per line).
 *
 *  @param[in/out] STREAM pointer to an opened file or stdin
 *  @param[in] sampling_frequency as its name says (in Hz)
 *  @returns a pointer ot an initialized gsl_vector
*/
gsl_vector * read_spike_train(FILE * STREAM, double sampling_frequency)
{
  size_t buffer_length = default_length;
  double *buffer=calloc(buffer_length, sizeof(double));
  size_t counter=0;
  char line[BUFSIZ];
  while (fgets (line, BUFSIZ, STREAM))
  {
    buffer[counter] = atof(line)/sampling_frequency;
    counter++;
    if (counter>=buffer_length)
    {
      buffer_length*=2;
      buffer=realloc(buffer,buffer_length*sizeof(double));
    }
  }
  if (!feof(STREAM))
  {
    fprintf (stderr, "Reading problem\n");
    exit (EXIT_FAILURE);
  }
  gsl_vector * res = gsl_vector_alloc(counter);
  for (size_t i=0; i<counter; i++)
    gsl_vector_set(res,i,buffer[i]);
  free(buffer);
  return res;
}
