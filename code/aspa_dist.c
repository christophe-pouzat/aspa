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
  int k,m,i,j,g,eH,eQ;
  double h,s,*H,*Q;
  s=d*d*n;
  if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
  k=(int)ceil(n*d);
  m=2*k-1;
  h=k-n*d;
  H=(double*)malloc((m*m)*sizeof(double));
  Q=(double*)malloc((m*m)*sizeof(double));
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

/** @brief Returns the (two-sided) Kolmogorov statistics
 *
 *  The data are contained in the `gsl_vector` pointed to
 *  by `data`. If the content is not sorted (`sorted==false`)
 *  the data are first copied before being sorted.
 *
 *  @param[in] data pointer to a `gsl_vector` containing the data
 *  @param[in] sorted a boolean indicated if the `data` content is
 *             already sorted (`true`) or not (`false`)
 *  @returns a double with the Kolomogorov statistics
*/
double aspa_Kolmogorov_D(gsl_vector * data, bool sorted)
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
  double D=0.;
  double inv_n = 1./data->size;
  for (size_t i=0; i<data->size; i++)
  {
    double x = gsl_vector_get(data_s,i);
    double diff = x-i*inv_n; 
    if (diff > D)
      D = diff;
    diff = inv_n - diff;
    if (diff > D)
      D = diff;
  }
  if (sorted == false)
    gsl_vector_free(data_s);
  return D;
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
