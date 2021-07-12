#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


double x0=0.25;
double sigma=0.05;
double L=1;
FILE *gpipe;

double k0=70*3.14;
double a=0.032;
int n=1000;
double h=1;
#define deltax (0.001*h)
#define deltat (2*deltax*deltax/h)
#define b (-h*h/(4*deltax*deltax))
#define V0 (50*3.14*50*3.14)
#define E_wave (h*h*k0*k0/2)


typedef struct  //Statistic Structure
    {
        double mediaX, VarX, integral, E,mediaP, VarP;
    } Statistics;

typedef struct  //Wave Structure
    {
        complex *fi;
        double k,sigma,x0,c,C;
    } Wave;
void InitializeWave(Wave *o)
{
    int i ;
    o->fi=(complex*) malloc(n*sizeof(complex));
    double sum=0;
    o->k=k0;
    o->c=0;
    o->sigma=sigma;
    o->x0=x0;
    for (i=0; i<n ; i++ )
    {
        sum=sum+cexp(-(i*deltax-o->x0)*(i*deltax-o->x0)/(2*o->sigma*o->sigma));
    }
    o->C=pow(deltax*sum, -0.5);
    for (i=0 ; i<n ; i++)
    {
        o->fi[i]=o->C*cexp(I*o->k*(i*deltax))*cexp(-(i*deltax-o->x0)*(i*deltax-o->x0)/(4*o->sigma*o->sigma));
    }
}
void InitializeWellPotential (double *V)
{
     int i;
     for (i=0; i<n ; i++ )
     {
         if (fabs(i*deltax-L/2)<a)
         {
             V[i]=V0;
         }
         else
            V[i]=0;
     }
}


void InitializeZeroPotential (double *V)
{
     int i;
     for (i=0; i<n ; i++ )
     {
         V[i]=0;
     }
}
 void cranknicholson (Wave *o, double *V)// Evolve wave in time
 {
     complex A[n],x[n],alfa[n],beta[n];
     int i;
     alfa[n-1]=0;
     beta[n-1]=0;
     alfa[0]=0;
     beta[0]=0;
     o->fi[n-1]=0;
     o->fi[0]=0;
     x[0]=0;
     x[n-1]=0;
      for (i=0; i<n ; i++ )
      {
          A[i]=(-h*h/(2*deltax*deltax)-V[i]/2+I*h/deltat);
      }

    for (i=1; i<n-1 ; i++ )
    {
        x[i]=b*o->fi[i+1]+b*o->fi[i-1]+(h*h/(2*deltax*deltax)+V[i]/2+I*h/deltat)*o->fi[i];
    }
    for (i=n-1; 0<i ; i-- )
    {
        alfa[i-1]=b/(A[i-1]-b*alfa[i]);
        beta[i-1]=(x[i-1]+b*beta[i])/(A[i-1]-b*alfa[i]);

    }
    for (i=1; i<n-1 ; i++ )
    {
        o->fi[i]=beta[i]+alfa[i]*o->fi[i-1];

    }
 }
Statistics RunStatistics (Wave *o, double *V)// Calculates some statistics
{
    double sum=0, sumx=0, sumx2=0, sump2=0, sumE=0;
    complex sump=0;
    int i;
    Statistics stats;
    for (i=1; i<n-1 ; i++ )
    {

        sum+=deltax*cabs(o->fi[i])*cabs(o->fi[i]);
        sumx+=deltax*deltax*i*cabs(o->fi[i])*cabs(o->fi[i]);
        sumx2+=deltax*(deltax*i)*(deltax*i)*cabs(o->fi[i])*cabs(o->fi[i]);
        sumE+=conj(o->fi[i])*(2*b*o->fi[i+1]+2*b*o->fi[i-1]+(h*h/(deltax*deltax)+V[i])*o->fi[i])*deltax;
        sump+=conj(o->fi[i])*o->fi[i+1]-conj(o->fi[i])*o->fi[i-1];
        sump2+=deltax*conj(o->fi[i])*(o->fi[i+1]-2*o->fi[i]+o->fi[i-1])/(deltax*deltax);
    }
    i=0;
    sumE=sumE+conj(o->fi[i])*(2*b*o->fi[i+1]+(h*h/(2*deltax*deltax)+V[i])*o->fi[i])*deltax;
    i=n-1;
    sumE+=conj(o->fi[i])*(2*b*o->fi[i-1]+(h*h/(2*deltax*deltax)+V[i])*o->fi[i])*deltax;

    stats.integral=sum;
    stats.mediaX=sumx;
    stats.VarX=sumx2-stats.mediaX*stats.mediaX;
    stats.mediaP=-I*h/2*sump;
    stats.VarP=-h*h*sump2-stats.mediaP*stats.mediaP;
    stats.E=sumE;
    return stats;
}
void Plot (Wave *o, double *V)// Plot the evolution of the wave in time
{
    int i;
    fprintf(gpipe,"set label 'time = %lf' at graph 0.7,0.8  \n",deltat*o->c);
    fprintf(gpipe,"p [] '-' title 'PDF(t)' w l axis x1y1, '-' w l axis x1y2 title 'Potential', '-' w l axis x1y2 title 'Wave Energy'\n");
    for (i=0; i<n; i++)
        fprintf(gpipe, "%lf %lf \n", i*deltax, cabs(o->fi[i])*cabs(o->fi[i]));

    fprintf(gpipe, "e\n");
    
    for (i=0; i<n; i++)
        fprintf(gpipe, "%lf %lf \n", i*deltax, V[i]);    
    
    fprintf(gpipe, "e\n");
    fprintf(gpipe, "0 %lf \n %lf %lf \n e\n",E_wave,L,E_wave);
    fprintf(gpipe,"unset label\n");
    o->c++;
}

int main()
{
  Wave o;
  int i;
  FILE *fid;
  fid=fopen("Statistics.txt", "w");
  Statistics stats;
  gpipe=popen("gnuplot", "w");
  float Tics_Num;
  if (E_wave-V0>=0)
    Tics_Num = E_wave;
  else
    Tics_Num = V0;

  fprintf(gpipe,"set term gif animate \n set output 'Single-wall_potential.gif' \n set xlabel 'x' \n set ylabel 'Probability density function (PDF)' \n set y2label 'Energy [J]'\n set ytics nomirror\n set yrange [0:20]\n set y2range [0:%lf]\n set format y2 '%%.0e' \n set y2tics 0,%lf,%lf \n set title 'Quantum Wave' \n", Tics_Num*1.1,Tics_Num/2,Tics_Num);//arma un gif
  double V[n];
  InitializeWave(&o);
  InitializeZeroPotential(V);

  for (i=0; i<10000; i++)
  {
        stats=RunStatistics(&o,V);
        fprintf(fid,"%lf %lf %lf %lf %lf %lf %lf \n", i*deltat, stats.E, stats.integral, stats.mediaX, stats.VarX, stats.mediaP, stats.VarP);
    if (i%20==0)
        Plot(&o,V);

      cranknicholson(&o,V);
  }

 fclose(fid);
 fprintf(gpipe,"unset output");
 free(o.fi);
 return 0;

}
