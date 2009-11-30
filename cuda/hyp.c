/* Linear convection equation with periodic BC
 * solved using MUSCL scheme
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define minmod(a,b)   ( (fabs(a) < fabs(b)) ? (a) : (b) )

//function prototypes
void initCond(int, float, float*);
void fluxFun(float*, float*);

float airk[3] = {0.0, 3.0/4.0, 1.0/3.0};

int main(){
   float *uold, *u, *uper, *fl, res;
   int np = 101;
   float dx = 1.0/(np-1);
   float dt, cfl;
   int niter, maxiter;
   int nirk, rkmax=3;
   int i;
   FILE *fpt;

   uold = (float*)malloc(np*sizeof(float));
   u    = (float*)malloc(np*sizeof(float));
   uper = (float*)malloc((np+2+2)*sizeof(float));
   fl   = (float*)malloc((np+1)*sizeof(float));

   cfl = 0.9;
   dt  = cfl*dx;
   maxiter = 1.0/dt + 1;

   //set initial conditions
   initCond(np, dx, u);

   fpt = fopen("init.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i]);
   fclose(fpt);


   //time-step loop
   for(niter=0; niter<maxiter; niter++){

      for(i=0; i<np; i++) uold[i] = u[i];

      //RK stages
      for(nirk=0; nirk<rkmax; nirk++){

         //periodic array for MUSCL scheme
         uper[0] = u[np-2];
         uper[1] = u[np-1];
         for(i=2; i<=np+1; i++) uper[i] = u[i-2];
         uper[np+2] = u[0];
         uper[np+3] = u[1];

         //flux computation
         for(i=0; i<np+1; i++) fluxFun(&uper[i+1], &fl[i]);

         //update conserved variable
         for(i=0; i<np; i++){
            res = fl[i+1] - fl[i];
            u[i] = airk[nirk]*uold[i] + 
                   (1.0-airk[nirk])*(u[i] - (dt/dx)*res);
         }

      }

   }

   fpt = fopen("final.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i]);
   fclose(fpt);

   free(uold);
   free(u);
   free(uper);
   free(fl);

}

//set initial condition
void initCond(int np, float dx, float *u){
   int i;
   float x;

   for(i=0; i<np; i++){
      x = dx*i;
      u[i] = sin(2.0*M_PI*x);
   }
}

//flux function
void fluxFun(float *u, float *fl){
   float uj, ujp1, ujm1, ujp2;
   float ul, ur;

   uj   = *(u);
   ujp1 = *(u+1);
   ujm1 = *(u-1);
   ujp2 = *(u+2);

   ul = uj   + 0.5*minmod( (uj-ujm1), (ujp1-uj) );
   ur = ujp1 - 0.5*minmod( (ujp1-uj), (ujp2-ujp1) );

   *fl = ul;

}
