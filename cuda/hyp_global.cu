/* Linear convection equation with periodic BC
 * solved using MUSCL scheme
 * CUDA implementation of hyp.c using only global memory
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define minmod(a,b)   ( (fabs(a) < fabs(b)) ? (a) : (b) )


//function prototypes
void initCond(int, float, float*);
__global__ void fluxFun(float*, float*);
__global__ void update(int, float, float, float*, float*, float*);

int main(){
   float *u, *uper;
   float *uold_d, *u_d, *uper_d, *fl_d;
   int np = 101;
   float dx = 1.0/(np-1);
   float dt, cfl;
   int niter, maxiter;
   int nirk, rkmax=3;
   int i;
   FILE *fpt;
   dim3 grid, block;

   u    = (float*)malloc(np*sizeof(float));
   uper = (float*)malloc((np+2+2)*sizeof(float));

   cudaMalloc((void**)&uold_d, (np)*sizeof(float));
   cudaMalloc((void**)&u_d, (np)*sizeof(float));
   cudaMalloc((void**)&uper_d, (np+2+2)*sizeof(float));
   cudaMalloc((void**)&fl_d,   (np+1)*sizeof(float));

   cfl = 0.9;
   dt  = cfl*dx;
   maxiter = 1.0/dt + 1;

   //set initial conditions
   initCond(np, dx, u);

   fpt = fopen("init.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i]);
   fclose(fpt);

   cudaMemcpy(uold_d, u, (np)*sizeof(float), 
              cudaMemcpyHostToDevice);
   cudaMemcpy(u_d, u, (np)*sizeof(float), 
              cudaMemcpyHostToDevice);

   //time-step loop
   for(niter=0; niter<maxiter; niter++){

      //RK stages
      for(nirk=0; nirk<rkmax; nirk++){

         //periodic array for MUSCL scheme
         uper[0] = u[np-2];
         uper[1] = u[np-1];
         for(i=2; i<=np+1; i++) uper[i] = u[i-2];
         uper[np+2] = u[0];
         uper[np+3] = u[1];

         cudaMemcpy(uper_d, uper, (np+2+2)*sizeof(float), 
                    cudaMemcpyHostToDevice);

         //flux computation
         block.x = 3;
         grid.x  = (np+1)/block.x;
         fluxFun<<<grid,block>>>(uper_d, fl_d);

         //update conserved variable
         block.x = 1;
         grid.x  = (np)/block.x;
         update<<<grid,block>>>(nirk, dt, dx, uold_d, fl_d, u_d);
         cudaMemcpy(u, u_d, (np)*sizeof(float), 
                    cudaMemcpyDeviceToHost);
      }
      cudaMemcpy(uold_d, u_d, (np)*sizeof(float), 
                 cudaMemcpyDeviceToDevice);

   }

   fpt = fopen("final.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i]);
   fclose(fpt);

   free(u);
   free(uper);

   cudaFree(uold_d);
   cudaFree(u_d);
   cudaFree(uper_d);
   cudaFree(fl_d);

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
__global__ void fluxFun(float *u, float *fl){
   float uj, ujp1, ujm1, ujp2;
   float ul, ur;

   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   ujm1 = *(u+idx);
   uj   = *(u+idx+1);
   ujp1 = *(u+idx+2);
   ujp2 = *(u+idx+3);

   ul = uj   + 0.5*minmod( (uj-ujm1), (ujp1-uj) );
   ur = ujp1 - 0.5*minmod( (ujp1-uj), (ujp2-ujp1) );

   fl[idx] = ul;

}

__global__ void update(int nirk, float dt, float dx, float *uold, float *fl, 
                       float *u){
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   float res;
   float airk[3] = {0.0, 3.0/4.0, 1.0/3.0};

   res = fl[idx+1] - fl[idx];
   u[idx] = airk[nirk]*uold[idx] + 
          (1.0-airk[nirk])*(u[idx] - (dt/dx)*res);
}
