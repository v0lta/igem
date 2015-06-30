
/***************************************************

Author: Teemu Leppanen (tileppan@lce.hut.fi)

This template is for numerically solving the Brysselator 
model 

du/dt = Du*lap(u)+a-(b+1)*u+u^2*v,
dv/dt = Dv*lap(v)+b*u-u^2*v,

and obtaining nice spatial patterns.

The calculation of the laplacian employs the finite difference method
with periodic boundary conditions.

The data (concentration fields u and v) is saved into a file named
"data.m". The data can be retrieved by MATLAB (you can use the
attached script visual.m).

Enter four command line arguments (e.g. simu 1234 1234 4.5 5.0):
1) First seed for the random number generator (between 0 and 31328)
2) Second seed for the random number generator (between 0 and 30081)
3) Parameter a value
4) Parameter b value

****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#define M 100  // length of the domain in x-direction
#define N 100  // length of the domain in y-direction
#define iterations 50000  // number of iterations 

#define dt .01  // time step
#define dx 1    // lattice constant 

#define D_u 1   // diffusion coefficient of chemical u
#define D_v 8  // diffusion coefficient of chemical v

int main(int argc, char **argv){
  
  int rmarin(int ij, int kl);
  int ranmar(float rvec[], int len);

  double **u, **v;        /* arrays for the concentration fields */
  
  double **ff;          /* array for the reaction kinetics */
  
  double **lap;         /* array for the laplacian */

  float *temp;          /* array for the random number generator */
  
  int seed1, seed2;     /* user given parameters for the random number generator*/
  
  int i, j, p, index=0; /* some indeces for the loops */

  float a, b;             /* the reaction parameters (input from command line) */  

  double f,g;

  double u_0, v_0, var;  /* variables for the initial conditions */

  FILE *tied;           /* file pointer for the output stream*/
  
  if(argc < 5){
    printf("\nEnter four command line arguments (e.g. simu 1234 1234 4.5 5.0):\n1)First seed for the random number generator (between 0 and 31328)\n2)Second seed for the random number generator (between 0 and 30081)\n3)Parameter A value\n 4) Parameter B value\n\n");
    return(1);
  }
  
  /* Allocating space for the arrays */

  u = malloc(M*sizeof(double));
  if(u==NULL)
    {
      printf("u out of memory\n");
      exit(1);
    }

  for(i=0;i<M;i++){
    u[i]=malloc(N*sizeof(double));
    if(u[i]==NULL)
      {
	printf("u[%d] out of memory\n",i);
	exit(1);
      }
  }

  v = malloc(M*sizeof(double));
  if(v==NULL)
    {
      printf("v out of memory\n");
      exit(1);
    }
 
  for(i=0;i<M;i++){
    v[i]=malloc(N*sizeof(double));
    if(v[i]==NULL)
      {
	printf("v[%d] out of memory\n",i);
	exit(1);
      }
  } 
  
  ff = malloc(M*sizeof(double *));
  if(ff==NULL)
    {
      printf("ff out of memory\n");
      exit(1);
    }

  for(i=0;i<M;i++){
    ff[i]=malloc(N*sizeof(double));
    if(ff[i]==NULL)
      {
	printf("ff[%d] out of memory\n",i);
	exit(1);
      }
  }

  lap = malloc(M*sizeof(double *));
  if(lap == NULL)
    {
      printf("lap out of memory\n");
      exit(1);
    }

  for(i=0;i<M;i++)
    {
      lap[i]=malloc(N*sizeof(double));
      if(lap[i]==NULL)
	{
	  printf("lap[%d] out of memory\n",i);
	  exit(1);
	}   
    }  

  temp = calloc(N*M*2,sizeof(float));

  /* Allocation completed */


  /* Read the input line parameters (two seeds for the random number
     generator and the parameters a and b) */

  seed1 = atoi(*++argv);
  seed2 = atoi(*++argv);
  
  a = atof(*++argv);
  b = atof(*++argv);

  
  /* Initialize the random number generator */

  if(rmarin(seed1,seed2)==1)
    return 0;

  
  /* Fill the array temp with random numbers between 0 and 1 */

  ranmar(temp,M*N*2);


  /* The initial values of the concentration fields are random
     deviations around the stationary state */

  /************************

  ADD HERE THE VALUES OF THE STATIONARY STATES YOU HAVE CALCULATED
  (U_0 and V_0) IN THE PROBLEM 1 AND AN AMPLITUDE FOR THE INITIAL
  RANDOM DEVIATIONS FROM THESE STATES (var).

  ************************/

  u_0 = a;
  
  v_0 = b/a;
  
  var = 1.5;
  
  for(i=0;i<M;i++)     
    for(j=0;j<N;j++){
      u[i][j] = u_0 + var*(temp[index]-.5);
      index++;
      v[i][j] = v_0 + var*(temp[index]-.5);
      index++;
    }

  
  /***************************** 

  Here begins the iteration loop

  *****************************/
  
  for(p=0;p<iterations;p++){


    /* calculation of the laplacian with respect to u with finite
       difference method and periodic boundary conditions */
    
    for(i=0;i<M;i++)
      for(j=0;j<N;j++){
	lap[i][j] = 0;
	if(j!=N-1)
	  lap[i][j] += u[i][j+1]-u[i][j];
	if(j!=0)
	  lap[i][j] += u[i][j-1]-u[i][j];
	if(i!=0)
	  lap[i][j] += u[i-1][j]-u[i][j];
	if(i!=M-1)
	  lap[i][j] += u[i+1][j]-u[i][j];
	if(j==N-1)
	  lap[i][j] += u[i][0]-u[i][j];
	if(j==0)
	  lap[i][j] += u[i][N-1]-u[i][j];
	if(i==M-1)
	  lap[i][j] += u[0][j]-u[i][j];
	if(i==0)
	  lap[i][j] += u[M-1][j]-u[i][j];
	lap[i][j] = lap[i][j]/(dx*dx);
      }
    
        
    for(i=0;i<M;i++)
      for(j=0;j<N;j++){
	
        /************************

        ADD HERE THE ITERATION STEP OF THE EULER'S METHOD FOR CONCENTRATION DATA u

        HINT: u THAT IS WRITTEN HERE IS THE u^{t+dt} OF THE ITERATION FORMULA

	************************/
	
	f = a - (b+1)*u[i][j] + u[i][j]*u[i][j] * v[i][j];
	u[i][j] = u[i][j] + dt * (D_u*lap[i][j] + f) ;

      }


    /* calculation of the laplacian with respect to v with finite
       difference method and periodic boundary conditions */  
    
    for(i=0;i<M;i++)
      for(j=0;j<N;j++){	
	lap[i][j] = 0;
	if(j!=N-1)
	  lap[i][j] += v[i][j+1]-v[i][j];
	if(j!=0)
	  lap[i][j] += v[i][j-1]-v[i][j];
	if(i!=0)
	  lap[i][j] += v[i-1][j]-v[i][j];
	if(i!=M-1)
	  lap[i][j] += v[i+1][j]-v[i][j];
	if(j==N-1)
	  lap[i][j] += v[i][0]-v[i][j];
	if(j==0)
	  lap[i][j] += v[i][N-1]-v[i][j];
	if(i==M-1)
	  lap[i][j] += v[0][j]-v[i][j];
	if(i==0)
	  lap[i][j] += v[M-1][j]-v[i][j];
	lap[i][j] = lap[i][j]/(dx*dx);
      }
    
    
    for(i=0;i<M;i++)
      for(j=0;j<N;j++){

        /************************

        ADD HERE THE ITERATION STEP OF THE EULER'S METHOD FOR CONCENTRATION DATA v

        HINT: v THAT IS WRITTEN HERE IS THE v^{t+dt} OF THE ITERATION FORMULA

	************************/

  g = b*u[i][j] - u[i][j]*u[i][j]*v[i][j];
	v[i][j] = v[i][j] + dt *(D_v * lap[i][j] + g);

      }
  
    
    if(p%100 == 0)
      fprintf(stderr,"%d iterations\n",p);
    
  }
  

  /* WRITE THE OUTPUT DATA */
  
  tied = fopen("data.m","w");
  
  fprintf(tied,"u = [");
  
  for(i=0;i<M;i++){
    for(j=0;j<N;j++)
      fprintf(tied,"%f ",u[i][j]);
    fprintf(tied,";");
  }
  fprintf(tied,"];\n");
  
  fprintf(tied,"v = [");
  
  for(i=0;i<M;i++){
    for(j=0;j<N;j++)
      fprintf(tied,"%f ",v[i][j]);
    fprintf(tied,";");
  }
  fprintf(tied,"];\n");
  
  fclose(tied);
  
  return (1);
  
}








