#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

typedef struct{
    double  x, y, px, py, dt,fx,fy,it,U ;
} particle;

typedef struct{
    double *hist;
    double dx;
} Tuple;

typedef double real_t;

#define Nc 30
#define N Nc*Nc
#define L sqrt(Nc*Nc/0.3)
//    double dt2=0.01*0.5;
//    double dt=0.01;
   double dt2=0.005*0.5;
   double dt=0.005;
   
void initialize (particle* p , int n, int m )// initialize the particle positions 
{
    double a;
    p->py=0;
    double random=(double)rand()/RAND_MAX;
    if (random<0.5)
        p->px=-1.1;
    else
        p->px=1.1;
    a=L/(double)(Nc+1);
    p->x=(double)n*a;
    p->y=(double)m*a;
    p->dt=0.005;



}

void BoundaryConditions(particle* s)// it confines the particles to move within a box of LxL dimension
{
    int i;
    for (i=0; i<= Nc*Nc-1 ; i++)
    {
        if (0>=s[i].x || s[i].x>=L)
        s[i].px=-s[i].px;
        if (0>=s[i].y || s[i].y>=L)
        s[i].py=-s[i].py;
    }

}
real_t V(particle* p1, particle* p2){// potential
real_t r=sqrt(pow(p1->x-p2->x,2)+pow(p1->y-p2->y,2));/////p1->x da el valor de la direccion de memoria p1.x
return ((24./(r*r))*(-1./pow(r,6)+2./pow(r,12)));
}

real_t fx(particle* p1 , particle* p2) {// x Force between two particles.
    return ((p1->x-p2->x))*V(p1,p2);
}
real_t fy(particle* p1 , particle* p2) {// y Force between two particles.
    return ((p1->y-p2->y))*V(p1,p2);
}


void AssessForce(particle* s)//calculates the interaction between all the particles
{
    int k, j;
    real_t Xforce,Yforce;

      for (k=0; k<= Nc*Nc-1 ; k++)
    {
        Xforce=0;
        Yforce=0;
        for (j=0; j<= Nc*Nc-1 ; j++)
        {
            if (k!=j){
            Xforce=Xforce+fx(&s[k],&s[j]);
            Yforce=Yforce+fy(&s[k],&s[j]);
            }
        }
        s[k].fx=Xforce;//es s[k].fx porque s[k] ya es el valor de la direccion de memoria de s
        s[k].fy=Yforce;

    }
}

void Verlet_Iterate(particle *s)// Updates the position and velocity for each particle
{
    int i;

    for (i=0; i<= Nc*Nc-1 ; i++)
    {
            s[i].px=s[i].px+s[i].fx*dt2;
            s[i].py=s[i].py+s[i].fy*dt2;
            s[i].x=s[i].x+s[i].px*dt;
            s[i].y=s[i].y+s[i].py*dt;
    }
        AssessForce(s);

      for (i=0; i<= Nc*Nc-1 ; i++)
        {
            s[i].px+=s[i].fx*dt2;
            s[i].py+=s[i].fy*dt2;
            s[i].it++;
        }

}

real_t KineticEnergy (particle *s)// Calculates the kinetic energy of the system.
{
    int i;
    real_t E=0;
    for (i=0; i<= Nc*Nc-1 ; i++)
    {
        E=E+(pow(s[i].px,2)+pow(s[i].py,2))/2;
    }
return E;
}

real_t Pot(particle* p1, particle* p2) // Calculates the potential between two particles
{
    real_t r=sqrt(pow(p1->x-p2->x,2)+pow(p1->y-p2->y,2));/////p1->x da el valor de la direccion de memoria p1.x
    return (4*(-pow(r,-6)+pow(r,-12)));
}

real_t potential(particle* s)// Calculates the potential of the system.
{
    int k, j;
    real_t U=0;
      for (k=0; k<= Nc*Nc-1 ; k++)
    {
        U=0;
        for (j=k+1; j<= Nc*Nc-1 ; j++)
        {

            U=U+Pot(&s[k],&s[j]);
        }
        s[k].U=U;

    }
    U=0;
   for (k=0; k<= Nc*Nc-1 ; k++)
   {
       U=U+s[k].U;
   }
return U;
}

Tuple histogram (particle *s, int p, bool print)//Calculates the histogram velocity.
{

    double a=-1000000, b=1000000, dx,h ;
    int j,i;

    for (j=0; j<= Nc*Nc-1 ; j++)
        {
            if (sqrt(pow(s[j].px,2)+pow(s[j].py,2))>=a)
            {
                a=sqrt(pow(s[j].px,2)+pow(s[j].py,2));
            }
        }
    for (j=0; j<= Nc*Nc-1 ; j++)
    {
        if (sqrt(pow(s[j].px,2)+pow(s[j].py,2))<=b)
            {
                b=sqrt(pow(s[j].px,2)+pow(s[j].py,2));
            }
    }
    dx=(a-b)/(double)p;
    double *hist=(double*) malloc(p*sizeof (double));
    for (j=0; j<p ; j++)
    {
        hist[j]=0;
    }
    for (j=0; j<= Nc*Nc-1 ; j++)
    {
        i=(int)((sqrt(pow(s[j].px,2)+pow(s[j].py,2))-b)/dx);
        hist[i]++;

    }

    if (print)
    {
        FILE* fid;
        fid= fopen("histogram.txt","w");
        for(j=0; j<p; j++)
        {
            fprintf(fid,"%lf  %lf \n",b+dx*j,hist[j]);
        }
        fclose(fid);
    } 
    // free(hist);
    Tuple tuple = {hist, dx};
    return tuple;
}



real_t Entropy(Tuple tuple, int p)
{
    int j;
    double dx,h;
    double *histogram;
    
    histogram = tuple.hist;
    dx = tuple.dx;
    real_t entropy = 0;
    double c= Nc*Nc*dx;
    for (j=0; j<p; j++)
    {
        h = histogram[j]/c;
        if (h!=0)
            entropy = entropy - h*log(h)*dx; // Boltzman entropy
    }
    return entropy;
}



void ReverseVelocity(particle *s)// reverse velocity
{
    int j;
    for (j=0; j<= Nc*Nc-1 ; j++)
    {
        s[j].px=-s[j].px;
        s[j].py=-s[j].py;
    }
}

void perturbation (particle *s)// Perturbate the positions and velocities of particles.
{
     int j;
     double random;
    for (j=0; j<= Nc*Nc-1 ; j++)
    {
        random=(double)rand()/RAND_MAX;
        if (random<0.5)
        s[j].px=s[j].px+s[j].px*0.01/100;
        else
          s[j].px=s[j].px-s[j].px*0.01/100;
        random=(double)rand()/RAND_MAX;
        if (random<0.5)
        s[j].py=s[j].py+s[j].py*0.01/100;
        else
          s[j].py=s[j].py-s[j].py*0.01/100;
        random=(double)rand()/RAND_MAX;
        if (random<0.5)
        s[j].x=s[j].x+s[j].x*0.01/100;
        else
          s[j].x=s[j].x-s[j].x*0.01/100;
        random=(double)rand()/RAND_MAX;
        if (random<0.5)
        s[j].y=s[j].y+s[j].y*0.01/100;
        else
          s[j].y=s[j].y-s[j].y*0.01/100;
    }
}


void PhysicalProperties (particle *s) // Creates a file with the K, U and total E for a time interval.
{
    FILE* fid;
    fid= fopen("PhysicalProperties.txt","w");
    int j;
    int iteration = 100;

    for (j=0; j<=iteration; j++)
    {
        BoundaryConditions(s);
        Verlet_Iterate(s);
        fprintf(fid, " %lf  %lf  %lf %lf  \n",s[1].dt*s[1].it, KineticEnergy(s), potential(s) ,KineticEnergy(s)+potential(s));
    }
    fclose(fid);
    printf("Hello world!\n");
}

void VelocityDistribution (particle *s) // Creates a histogram with the velocity distribution and a file with the system entropy for a time interval. Check the Maxwell-Boltzman distribution in the histogram.
{
    FILE* fid;
    fid= fopen("entropy.txt","w");
    int j;
    int iteration = 501;
    bool print = {true};
    bool Dontprint = {false};
    real_t entropy;
    for (j=0; j<=iteration; j++)
    {
        BoundaryConditions(s);
        Verlet_Iterate(s);
        entropy = Entropy(histogram(s,30,Dontprint),30);
        fprintf(fid,"%lf %lf \n",dt*s[1].it, entropy);
        if (j==500)
            histogram(s,30,print);
    }
    fclose(fid);
}

void LoschmidtParadox (particle *s) // reverse the velocity of all the particles a time and the entropy decrease.
{
    FILE* fid;
    fid= fopen("LoschmidtParadox.txt","w");
    int j;
    int iteration = 1300;
    bool Dontprint = {false};
    real_t entropy;
    for (j=0; j<=iteration; j++)
    {
        BoundaryConditions(s);
        Verlet_Iterate(s);
        // if (j==350)
        //     perturbation(s);
        if (j==650)
            ReverseVelocity(s);
        entropy = Entropy(histogram(s,30,Dontprint),30);
        fprintf(fid,"%lf %lf \n",dt*s[1].it, entropy);
    }
    fclose(fid);
}


void MolecularDynamics (particle *p)  // Creates a GIF with the movement of the particles
{   
    FILE *gpipe;
    gpipe=popen("gnuplot", "w");
    fprintf(gpipe,"set term gif animate delay 1 \n set output 'MolecularDynamics.gif' \n set style line 1 lc rgb 'black' pt 7 ps 1 \n set xlabel 'x' \n set ylabel 'y' \n set yrange [1:20] \n set xrange [1:20] \n set title 'Molecular Dynamics' \n");//arma un gif
    int j,i;
    int iteration = 2000;
    for (j=0; j<=iteration; j++)
    {   
        BoundaryConditions(p);
        Verlet_Iterate(p);
        fprintf(gpipe,"p [-0.5:%lf] [0:%lf] [] '-' title 't=%lf' w p ls 1\n",L,L,dt*p[1].it);
        for (i=0; i<= Nc*Nc-1 ; i++)
        {
            fprintf(gpipe, "%lf %lf \n", p[i].x, p[i].y);
        }
        fprintf(gpipe, "e\n");
        fflush(gpipe);
    }
    fprintf(gpipe,"unset output");
    fclose(gpipe);
}


int main()
{   int n=0,j,m=0,i;
    srand(time(NULL));
    particle s[Nc*Nc];
    for (i=0; i<= Nc*Nc-1 ; i++)//inicializa la posicion de las particulas
    {
        if(m==Nc)
        {
            m=0;
            n++;
        }
        initialize(&s[i],n,m);
        m++;
    }
    AssessForce(s);

    // VelocityDistribution(s);
    // LoschmidtParadox(s);
    // PhysicalProperties(s);
    MolecularDynamics(s);

    printf("Finished!!! \n");
    return 0;
}
