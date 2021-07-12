#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


double x0=4;

double L=10;
FILE *gpipe;

double k0=0;
double h=1;
int n = 1000; // L/deltax
#define deltax (0.01*h)
#define deltat (4*deltax*deltax/h)
#define b (-h*h/(4*deltax*deltax))



typedef struct  //Statistic Structure
    {
        double MeanX, VarX, integral, E,MeanP, VarP;
    } Statistics;

typedef struct  //Wave Structure
    {
        complex *fi;
        double k,sigma,x0,c,C;
    } Wave;
    
void InitializeWave(Wave *o)
{
    int i ;
    float sigma = sqrt(h/2);
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

void InitializeDoubleWellPotential (double *V)
{
    int i;
    for (i=0; i<n ; i++ )
    {
        V[i]=-(i*deltax-L/2)*(i*deltax-L/2)/4+(i*deltax-L/2)*(i*deltax-L/2)*(i*deltax-L/2)*(i*deltax-L/2)/8;
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
        if ((L/2)/deltax<i)
        sum+=deltax*cabs(o->fi[i])*cabs(o->fi[i]);
        sumx+=deltax*(i*deltax-L/2)*cabs(o->fi[i])*cabs(o->fi[i]);

    }
    i=0;
    sumE=sumE+conj(o->fi[i])*(2*b*o->fi[i+1]+(h*h/(2*deltax*deltax)+V[i])*o->fi[i])*deltax;
    i=n-1;
    sumE+=conj(o->fi[i])*(2*b*o->fi[i-1]+(h*h/(2*deltax*deltax)+V[i])*o->fi[i])*deltax;

    stats.integral=sum;
    stats.MeanX=sumx;
    return stats;
}

void Plot (Wave *o, double *V)// Plot the evolution of the wave in time
{
    int i;
    fprintf(gpipe,"set label 'Time = %lf' at graph 0.74,0.85  \n",deltat*o->c);
    fprintf(gpipe,"set label 'Wave Energy = 0' at graph 0.74,0.80  \n");
    fprintf(gpipe,"p [] '-' title 'PDF(t)' w l axis x1y1, '-' w l axis x1y2 title 'Potential'\n");
    for (i=0; i<n; i++)
        fprintf(gpipe, "%lf %lf \n", i*deltax-L/2, cabs(o->fi[i])*cabs(o->fi[i]));

    fprintf(gpipe, "e\n");
    
    for (i=0; i<n; i++)
        fprintf(gpipe, "%lf %lf \n", i*deltax-L/2, V[i]);    
    
    fprintf(gpipe, "e\n");

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


  fprintf(gpipe,"set term gif animate \n set output 'Double-well_potential.gif' \n set xlabel 'x' \n set ylabel 'Probability density function (PDF)' \n set ytics nomirror\n set yrange [0:1] \n set xrange [-3:3] \n set y2range [-0.125:0]\n set title 'Quantum Wave / Double-well potential' \n");//arma un gif
  
  double V[n];
  InitializeWave(&o);
  InitializeDoubleWellPotential(V);

  for (i=0; i<10500; i++)
  {
        stats=RunStatistics(&o,V);
        fprintf(fid,"%lf %lf %lf %lf %lf %lf %lf \n", i*deltat, stats.E, stats.integral, stats.MeanX, stats.VarX, stats.MeanP, stats.VarP);
    if (i%20==0)
        Plot(&o,V);

      cranknicholson(&o,V);
  }

 fclose(fid);
 fprintf(gpipe,"unset output");
    free(o.fi);
    return 0;
}
