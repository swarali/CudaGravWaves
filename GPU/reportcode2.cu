#include<stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cufft.h>
#include <math.h>
#include<stdlib.h>

/**Correlation coefficient GPU code by taking the inner product of fourier tranform
*
	Calculates the correlation coefficient of a sine wave of given frequency with a given number of sine waves with 
	frequency fr +/- somevalue... -0.5<somevalue<0.5 

	Uses cuda fft library for fft on the device
	Input : sampling frequency, time period, number of samples and frequency(samp * t should be power of 2)
*/

#define PI 3.14159265358979323846

float signal(double fr, double t);
void cudaErrorCheck(char * errormsg);
__global__ void cordev(cufftComplex * x, cufftComplex * y, int size, float *ans);

#define PI 3.14159265358979323846
#define Rse 149600000000	//distance between earth and sun
#define Re 6371000		//earth's radius
#define alpha PI/180.0 * 90.0		//colatitude of the detector
#define beta0 PI/180.0 * 0.0		//local sidereal time
#define theta PI/180.0 * 90.0		//celestial colatitude
#define phi   PI/180.0 * 0.0		//celestial colongitude
#define c_light 299792458.0  //velocity of light
#define omega_rot (7.29211/100000)	//rotational angular velpcity of earth
#define epsilon PI/180.0 * 23.45		//obliquity of the elliptic
#define H0 1.0


float signal(double fr, double t)
{
	float ans, Z, N, xi_rot, P, Q, delta, R;

	Z=2*PI*fr*Rse*sin(theta)/c_light;
	xi_rot = omega_rot*t;
	
	P= 2*PI*fr*Re*sin(alpha)*(cos(beta0)* (sin(theta)*cos(epsilon)*sin(phi) + cos(theta)*sin(epsilon))  - (sin(beta0)*sin(theta)*cos(phi)))/c_light;
	Q= 2*PI*fr*Re*sin(alpha)*(sin(beta0)* (sin(theta)*cos(epsilon)*sin(phi) + cos(theta)*sin(epsilon))  + (cos(beta0)*sin(theta)*cos(phi)))/c_light;
	N= sqrt(P*P + Q*Q);
	delta = atan(P/Q);
	
	R=Z *cos(phi);

//	ans = 2*PI*fr*t - (2*PI*fr*Re/c_light)*(1 - cos(omega_rot*t));
//	printf("\n%lf",(double) P);
	ans = 2*PI*fr*t +Z*cos((xi_rot/365.26) - phi) + N*cos(xi_rot - delta) - R -Q; // omega_orb/omega_rot = 1/365.26
//	ans = 2*PI*fr*t + N*cos(xi_rot - delta) - R -Q; // omega_orb/omega_rot = 1/365.26
//	ans = 2*PI*fr*t;
//	printf("%f..", H0+cos(ans));
	return (H0 + cos(ans));
}



int main()
{
	
		int i, j;
		int samp, t, batch;
		int nx;
		float fr, tempfr;
		float  *d_device_cor, *h_device_cor;
		cufftHandle plan;
		cufftComplex *d_temp1, *h_temp1;	
		cufftComplex *h_data1, *d_data1;
		float temp, sum;
	
	printf("Sampling frequency : ");
	scanf("%d", &samp);

	printf("Time period : ");
	scanf("%d", &t);

	printf("Number of Samples : ");
	scanf("%d", &batch);

	printf("Frequency : ");
	scanf("%f", &fr);
	
	nx = samp *t;
	
	h_device_cor= (float *)malloc(sizeof(float) * batch);
	
	h_temp1 = (cufftComplex *)malloc(sizeof(cufftComplex) * nx * batch);
	h_data1 = (cufftComplex *)malloc(sizeof(cufftComplex) * nx );

	for(i=0;i<nx;i++)
	{
//		h_data1[i].x = sin(2*PI*i*fr/samp)/sqrt((float)nx/2);
		h_data1[i].x = signal((double)fr, (double)i/samp);
		h_data1[i].y = 0;
	}
	
	for(j=0;j<batch;j++)
	{	
		sum=0;
		tempfr = (fr + (float)(j +1)/batch -0.5);
		for(i=0;i<nx;i++)
		{
//			temp = sin(2*PI*i*(fr -0.5 + (float)j/batch)/samp)/sqrt((float)nx/2);			
			temp = signal((double)tempfr, (double)i/samp);
			h_temp1[j*nx + i].x = temp;
			h_temp1[j*nx + i].y = 0;
			sum += temp*temp;			
		}

	}
	
	cudaMalloc((void**)&d_temp1, nx * batch * sizeof(cufftComplex));
	cudaMalloc((void**)&d_data1,  nx * sizeof(cufftComplex));
	
	cudaMemcpy(d_temp1, h_temp1, nx * batch * sizeof(cufftComplex), cudaMemcpyHostToDevice);
	cudaMemcpy(d_data1, h_data1, nx * sizeof(cufftComplex), cudaMemcpyHostToDevice);

	cudaErrorCheck("Memory Copy Error");
	
	//FFT of temp
	cufftPlan1d(&plan, nx, CUFFT_C2C, batch);
	cufftExecC2C(plan, d_temp1, d_temp1,CUFFT_FORWARD);
	cufftDestroy(plan);
	
	//FFT of data
	cufftPlan1d(&plan, nx, CUFFT_C2C, 1);
	cufftExecC2C(plan, d_data1, d_data1,CUFFT_FORWARD);
	cufftDestroy(plan);
	
	dim3 dimGrid((int)sqrt((float)batch) +1, (int)sqrt((float)batch) +1);
	dim3 dimBlk(nx);
	
	cudaMalloc((void**)&d_device_cor,  batch * sizeof(float));
	cudaErrorCheck("Memory Allocation Error");
	
	cordev<<<dimGrid, dimBlk, 3*nx*sizeof(float)>>>(d_data1, d_temp1, nx, d_device_cor);
	
	cudaMemcpy(h_device_cor, d_device_cor, (batch * sizeof(float)), cudaMemcpyDeviceToHost);
	
	for(i=0;i<batch;i++)
		printf("fr : %f : %f \n", (fr + (float)(i +1)/batch -0.5), h_device_cor[i]);

	cudaFree(d_temp1);
	free(h_temp1);
	cudaFree(d_data1);
	free(h_data1);
	free(h_device_cor);
	cudaFree(d_device_cor);
	
	return 0;
	
}
	
void cudaErrorCheck(char * errormsg)
{
	cudaError_t e = cudaGetLastError();
	if(e != cudaSuccess)
	{
		printf("Error %s : %s", errormsg, cudaGetErrorString(e));
		exit(EXIT_FAILURE);
		
	}
}

__global__ void cordev(cufftComplex * x, cufftComplex * y, int nx, float *ans)
{
	extern float __shared__ buffer[];
	float *innr =&buffer[0];
	float *nrm1 =&buffer[nx] ;
	float *nrm2 =&buffer[2*nx];
	int tempsize= nx/2, index_x , index_y;

	
	index_y = (gridDim.x * blockIdx.y + blockIdx.x ) * (blockDim.x * blockDim.y) + threadIdx.y * blockDim.x + threadIdx.x;

	index_x = threadIdx.y * blockDim.x + threadIdx.x;
	
//	innr[index_x] = (x[index_x].x * y[index_y].x  + x[index_x].y * y[index_y].y)/(nx);
	innr[index_x] = (x[index_x].x * y[index_y].x  + x[index_x].y * y[index_y].y);
	nrm1[index_x] = (x[index_x].x * x[index_x].x) +(x[index_x].y *x[index_x].y);
	nrm2[index_x] = (y[index_y].x * y[index_y].x) +(y[index_y].y *y[index_y].y);

	__syncthreads();
	
		while(tempsize > 0)
		{
			if(index_x < tempsize)
			{
				innr[index_x] += innr[index_x + tempsize];
				nrm1[index_x] += nrm1[index_x + tempsize];
				nrm2[index_x] += nrm2[index_x + tempsize];

			}
			
			tempsize=tempsize/2;
			__syncthreads();		
			
		}
	__syncthreads();
		
	if(index_x == 0)
		ans[ blockIdx.y * gridDim.x + blockIdx.x] = innr[0]/sqrt(nrm1[0]*nrm2[0]);
//	  	 ans[ blockIdx.y * gridDim.x + blockIdx.x] = innr[0];

}



