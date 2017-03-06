/** @file aspa_dist.c
 *  @brief Function definitions for distribution functions and tests
 *         not included in the GSL
 *
 *  @author Christophe Pouzat <christophe.pouzat@parisdescartes.fr>
*/

#include "aspa.h"

/** @brief Multiplies m x m matrices A and B and return
 *         result in m x m matrix C
 *
 *  Code `mMultiply` of G. Marsaglia, Wai Wan 
 *  Tsang and Jingbo Wong, [J.Stat.Software. 8(18): 1-4](https://www.jstatsoft.org/article/view/v008i18).
 *
 *  @param[in] A pointer to an m x m matrix
 *  @param[in] B pointer to an m x m matrix
 *  @param[out] C pointer to an m x m matrix
 *  @param[in] m int the matrices dimension
 *  @returns nothing C is modified
*/
void mMultiply(const double *A,const double *B,double *C,int m)
{
  int i,j,k;
  double s;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
    {
      s=0.;
      for(k=0;k<m;k++)
	s+=A[i*m+k]*B[k*m+j];
      C[i*m+j]=s;
    }
}

/** @brief Computes nth power of mxm matrix A and stores result
 *         in mxm matrix V
 *
 *  Code `mPower` of G. Marsaglia, Wai Wan 
 *  Tsang and Jingbo Wong, [J.Stat.Software. 8(18): 1-4](https://www.jstatsoft.org/article/view/v008i18).
 *  
 *  @param[in] A pointer to mxm matrix whose nth power is looked for
 *  @param[in] eA an integer 
 *  @param[out] V pointer to mxm matrix containing nth power of A
 *  @param[out] eV pointer to an integer
 *  @param[in] m matrices size
 *  @param[in] n an integer (the sample size)
 *  @returns nothing
*/
void mPower(const double *A,int eA,double *V,int *eV,int m,int n)
{
  double *B;int eB,i;
  if(n==1)
  {
    for(i=0;i<m*m;i++)
      V[i]=A[i];
    *eV=eA;
    return;
  }
  mPower(A,eA,V,eV,m,n/2);
  B=(double*)malloc((m*m)*sizeof(double));
  mMultiply(V,V,B,m);
  eB=2*(*eV);
  if(n%2==0)
  {
    for(i=0;i<m*m;i++)
      V[i]=B[i];
    *eV=eB;
  }
  else
  {
    mMultiply(A,B,V,m);
    *eV=eA+eB;
  }

  if(V[(m/2)*m+(m/2)]>1e140)
  {
    for(i=0;i<m*m;i++)
      V[i]=V[i]*1e-140;
    *eV+=140;
  }
  free(B);
}

/** @brief Returns the Kolmogorov distribution function Prod{D_n <= d}
 *         where D_n is the Kolmogorov statistic and n the sample size
 *
 *  Adapation of code `K` of G. Marsaglia, Wai Wan 
 *  Tsang and Jingbo Wong, [J.Stat.Software. 8(18): 1-4](https://www.jstatsoft.org/article/view/v008i18).
 *  
 *  @param[in] n an integer, the sample size
 *  @param[in] d a double the maximal deviation
 *  @results Prod{D_n <= d}
*/
double aspa_cdf_K(int n,double d)
{
  int i,j,g,eH,eQ;
  double s=d*d*n;
  if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
  int k=(int)ceil(n*d);
  int m=2*k-1;
  double h=k-n*d;
  double *H=malloc((m*m)*sizeof(double));
  double *Q=malloc((m*m)*sizeof(double));
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      if(i-j+1<0) H[i*m+j]=0;
      else     H[i*m+j]=1;
  for(i=0;i<m;i++)
  {
    H[i*m]-=pow(h,i+1);
    H[(m-1)*m+i]-=pow(h,(m-i));
  }
  H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      if(i-j+1>0)
        for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
  eH=0;
  mPower(H,eH,Q,&eQ,m,n);
  s=Q[(k-1)*m+k-1];
  for(i=1;i<=n;i++)
  {
    s=s*i/n;
    if(s<1e-140){s*=1e140; eQ-=140;}
  }
  s*=pow(10.,eQ);
  free(H);
  free(Q);
  return s;
}

/** @brief Returns the Kolmogorov distribution function Prod{D_n_plus <= d}
 *         or Prod{D_n_minus <= d} where D_n_plus/minus are the one sided 
 *         Kolmogorov statistic and n the sample size
 *
 *  The probability is given by equation 3 of Birnbaum and Tingey (1951)
 *  One-sided confidence contours for probability distribution functions
 *  _The Annals of Mathematical Statistics_ __22__: 592-596.
 *  
 *  @param[in] n an integer, the sample size
 *  @param[in] d a double the maximal one sided deviation
 *  @results Prod{D_n_plus/minus <= d}
*/
double aspa_cdf_Kplus(int n,double d)
{
  if (d <= 0.)
    return 0.;
  if (d >= 1.)
    return 1.;
  unsigned int k = (unsigned int) floor(n*(1-d));
  double s=0.;
  double n_inv=1./n;
  for (int j=0; j<=k; j++)
    s += exp(gsl_sf_lnchoose(n,j)+
	     (n-j)*log(1.-d-j*n_inv)+
	     (j-1)*log(d+j*n_inv));
  return 1.-d*s;
}

/** @brief Returns the Kolmogorov statistics
 *
 *  The data are contained in the `gsl_vector` pointed to
 *  by `data`. If the content is not sorted (`sorted==false`)
 *  the data are first copied before being sorted.
 *  Argument `what` can be one of "D", "D+", "D-". If "D" the
 *  two sided statistics is returned, if "D+" the maximal distance
 *  of the dominating part of the empirical cdf to the theoretical 
 *  one is returned, if "D-" the maximal distance of the dominated
 *  part of the empirical cdf to the theoretical one is returned.
 *  If what is given a "wrong" value, -1.0 is returned.
 *
 *  @param[in] data pointer to a `gsl_vector` containing the data
 *  @param[in] sorted a boolean indicated if the `data` content is
 *             already sorted (`true`) or not (`false`)
 *  @param[in] what a character string, "D", "D+" or "D-"
 *  @returns a double with the Kolomogorov statistics of -1.0 if
 *           a wrong value for `what` was given
*/
double aspa_Kolmogorov_D(gsl_vector * data, bool sorted, char * what)
{
  char * choices[] = {"D","D+","D-"};
  gsl_vector * data_s;
  if (sorted == false)
  {
    data_s = gsl_vector_alloc(data->size);
    gsl_vector_memcpy(data_s,data);
    gsl_sort_vector(data_s);
  }
  else
  {
    data_s = data;
  }
  double D_p=0.;
  double D_m=0.;
  double inv_n = 1./data->size;
  for (size_t i=0; i<data->size; i++)
  {
    double x = gsl_vector_get(data_s,i);
    double diff = x-i*inv_n; 
    if (diff > D_p)
      D_p = diff;
    diff = inv_n - diff;
    if (diff > D_m)
      D_m = diff;
  }
  if (sorted == false)
    gsl_vector_free(data_s);
  if (strcmp(what,choices[0])==0)
    return GSL_MAX_DBL(D_p,D_m);
  if (strcmp(what,choices[1])==0)
    return D_p;
  if (strcmp(what,choices[2])==0)
    return D_m;
  return -1.;
}


/** @brief Returns the standard normal distribution function Prod{X <= x}
 *
 *  Code `Phi` of G. Marsaglia [J.Stat.Software. 11(4): 1-11](https://www.jstatsoft.org/article/view/v011i04).
 *  
 *  @param[in] x a double the observec value
 *  @results Prod{X <= x} where X is N(0,1)
*/
double aspa_cdf_norm_P(double x)
{
  long double s=x,t=0,b=x,q=x*x,i=1;
  while(s!=t)
    s=(t=s)+(b*=q/(i+=2));
  return .5+s*exp(-.5*q-.91893853320467274178L);
}


/** @brief Returns the complementary standard normal distribution function Prod{X > x}
 *
 *  Code `cPhi` of G. Marsaglia [J.Stat.Software. 11(4): 1-11](https://www.jstatsoft.org/article/view/v011i04).
 *  
 *  @param[in] x a double the observec value
 *  @results Prod{X > x} where X is N(0,1)
*/
double aspa_cdf_norm_Q(double x)
{
  long double R[9]={1.25331413731550025L,
		    .421369229288054473L,
		    .236652382913560671L,
		    .162377660896867462L,
		    .123131963257932296L,
		    .0990285964717319214L,
		    .0827662865013691773L,
		    .0710695805388521071L,
		    .0622586659950261958L};
       int i,j=.5*(fabs(x)+1);
       long double pwr=1,a=R[j],z=2*j,b=a*z-1,h=fabs(x)-z,s=a+h*b,t=a,q=h*h;
       for(i=2;s!=t;i+=2)
       {
	 a=(a+z*b)/i;
	 b=(b+z*a)/(i+1);
	 pwr*=q;
	 s=(t=s)+pwr*(a+h*b);
       }
       s=s*exp(-.5*x*x-.91893853320467274178L);
       if(x>=0)
	 return (double) s;
       return (double) (1.-s);
}


/** @brief Returns the Anderson-Darling statistics
 *
 *  The data are contained in the `gsl_vector` pointed to
 *  by `data`. If the content is not sorted (`sorted==false`)
 *  the data are first copied before being sorted.
 *
 *  @param[in] data pointer to a `gsl_vector` containing the data
 *  @param[in] sorted a boolean indicated if the `data` content is
 *             already sorted (`true`) or not (`false`)
 *  @returns a double with the Anderson-Darling statistics
*/
double aspa_AndersonDarling_W2(gsl_vector * data, bool sorted)
{
  gsl_vector * data_s;
  if (sorted == false)
  {
    data_s = gsl_vector_alloc(data->size);
    gsl_vector_memcpy(data_s,data);
    gsl_sort_vector(data_s);
  }
  else
  {
    data_s = data;
  }
  size_t n = data->size;
  double A=0.;
  for (size_t i=0; i<n; i++)
  {
    double t = gsl_vector_get(data_s,i)*(1.-gsl_vector_get(data_s,n-1-i));
    A += (i+i+1)*log(t);
  }
  A *= -1./n;
  if (sorted == false)
    gsl_vector_free(data_s);
  return A-(double)n;
}

/** @brief Utility function called by `aspa_cdf_ADinf_P`
 *
 *  Adaptation of function `ADf` of Marsaglia & Marsaglia (2004)
 *  [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 *
 *  @param[in] z a double
 *  @param[in] i an int
 *  @returns a double
*/
double aspa_ADf(double z,int j)
{ 
  double t=(4*j+1)*(4*j+1)*1.23370055013617/z;
  if(t>150.)
    return 0.;
  double a=2.22144146907918*exp(-t)/sqrt(t);
  double b=3.93740248643060*2.*aspa_cdf_norm_Q(sqrt(2*t));
  double r=z*.125;
  double f=a+b*r;
  for(size_t i=1; i<200; i++)
  {
    double c=((i-.5-t)*b+t*a)/i;
    a=b;
    b=c;
    r*=z/(8*i+8);
    if(fabs(r)<1e-40 || fabs(c)<1.e-40)
      return f;
    double fnew=f+c*r;
    if(f==fnew)
      return f;
    f=fnew;
  }
  return f;
}

/** @brief Returns the asymptotic cdf of the Anderson-Darling
 *         statistics
 *
 *  Adaptation of function `ADinf` of Marsaglia & Marsaglia (2004)
 *  [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 *
 *  @param[in] z a double the observed statistics value
 *  @returns a double Prob{W2 <= z}
*/
double aspa_cdf_ADinf_P(double z)
{
  if(z<.01)
    return 0.; /* avoids exponent limits; ADinf(.01)=.528e-52 */
  double r=1./z;
  double ad=r*aspa_ADf(z,0);
  for(size_t j=1; j<100; j++)
  {
    r*=(.5-j)/j;
    double adnew=ad+(4*j+1)*r*aspa_ADf(z,j);
    if(ad==adnew) {
      return ad;
    }
    ad=adnew;
  }
  return ad;
}

/** @brief Short, practical version of full aspa_cdf_ADinf_P(z), z>0.
 *
 *  Adaptation of function `adinf` of Marsaglia & Marsaglia (2004)
 *  [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 *
 *  @param[in] z a double the observed statistics value
 *  @returns a double Prob{W2 <= z}
*/
double aspa_adinf(double z)
{
  if(z<2.)
    return exp(-1.2337141/z)/sqrt(z)*(2.00012+(.247105-(.0649821-(.0347962-(.011672-.00168691*z)*z)*z)*z)*z);
  /* max |error| < .000002 for z<2, (p=.90816...) */
 return exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)*z)*z)*z)*z));
 /* max |error|<.0000008 for 4<z<infinity */
}

/** @brief Corrects the error caused by using the asymptotic 
 *         approximation, x=aspa_adinf(z).
 *
 *  Thus x+errfix(n,x) is uniform in [0,1) for practical purposes;
 *  accuracy may be off at the 5th, rarely at the 4th, digit.
 *  Adaptation of function `errfix` of Marsaglia & Marsaglia (2004)
 *  [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 *  @param[in] n an int
 *  @param[in] x a double
 *  @returns a double
*/
double aspa_errfix(int n, double x)
{
  double t;
  if(x>.8)
    return (-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
  double c=.01265+.1757/n;
  if(x<c)
  {
    t=x/c;
    t=sqrt(t)*(1.-t)*(49*t-102);
    return t*(.0037/(n*n)+.00078/n+.00006)/n;
  }
  t=(x-c)/(.8-c);
  t=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*t)*t)*t)*t)*t;
  return t*(.04213+.01365/n)/n;
}

double aspa_cdf_AD_P(int n,double z)
{
  double v;
  double x=aspa_adinf(z);
  if(x>.8)
  {
    v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
    return x+v;
  }
  double c=.01265+.1757/n;
  if(x<c)
  {
    v=x/c;
    v=sqrt(v)*(1.-v)*(49*v-102);
    return x+v*(.0037/(n*n)+.00078/n+.00006)/n;
  }
  v=(x-c)/(.8-c);
  v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
  return x+v*(.04213+.01365/n)/n;
}

/** @brief Perform "Durbin's modification" on data contained in `seq`
 *
 *  The data are supposed to between 0 and 1. If such is not the case
 *  the function returns -1 and prints an error message to the stderr.
 *
 *  @param[in] seq a pointer to a `gsl_vector` containing the data
 *  @param[out] res a pointer to a `gsl_vector` where the result is
 *              stored
 *  @results 0 if everything goes fine, -1 otherwise
*/
int aspa_durbin_modification(const gsl_vector * seq, gsl_vector * res)
{
  gsl_vector_memcpy(res,seq);
  gsl_sort_vector(res);
  size_t n=res->size;
  if (gsl_vector_get(res,0) < 0)
  {
    fprintf(stderr,"The elements of seq should all be >= 0.\n");
    return -1;
  }
  if (gsl_vector_get(res,n-1) > 1)
  {
    fprintf(stderr,"The elements of seq should all be <= 0.\n");
    return -1;
  }
  gsl_vector *iei = gsl_vector_alloc(n+1);
  gsl_vector_set(iei,0,gsl_vector_get(res,0));
  for (size_t i=1; i < n; i++)
    gsl_vector_set(iei,i,gsl_vector_get(res,i)-gsl_vector_get(res,i-1));
  gsl_vector_set(iei,n,1.-gsl_vector_get(res,n-1));
  gsl_sort_vector(iei);
  for (size_t i=n+1; i>1; i--)
    gsl_vector_set(iei,i,(n+2-i)*(gsl_vector_get(iei,i-1)-gsl_vector_get(iei,i-2)));
  gsl_vector_set(res,0,gsl_vector_get(iei,0));
  for (size_t i=1; i<n; i++)
    gsl_vector_set(res,i,gsl_vector_get(res,i-1)+gsl_vector_get(iei,i));
  gsl_vector_free(iei);
  return 0;
}
