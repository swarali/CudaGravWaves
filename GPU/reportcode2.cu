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

void cudaErrorCheck(char * errormsg);
__global__ void cordev(cufftComplex * x, cufftComplex * y, int size, float *ans);


int main()
{
	
		int i, j;
		int samp, t, batch;
		int nx;
		float fr;
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
		h_data1[i].x = sin(2*PI*i*fr/samp)/sqrt((float)nx/2);
		h_data1[i].y = 0;
	}
	
	for(j=0;j<batch;j++)
	{	
		sum=0;
		for(i=0;i<nx;i++)
		{
			temp = sin(2*PI*i*(fr -0.5 + (float)j/batch)/samp)/sqrt((float)nx/2);			
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
	
	cordev<<<dimGrid, dimBlk, nx*sizeof(float)>>>(d_data1, d_temp1, nx, d_device_cor);
	
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
	int tempsize= nx/2, index_x , index_y;
	
	index_y = (gridDim.x * blockIdx.y + blockIdx.x ) * (blockDim.x * blockDim.y) + threadIdx.y * blockDim.x + threadIdx.x;

	index_x = threadIdx.y * blockDim.x + threadIdx.x;
	
	innr[index_x] = (x[index_x].x * y[index_y].x  + x[index_x].y * y[index_y].y)/(nx);
	
	__syncthreads();
	
		while(tempsize > 0)
		{
			if(index_x < tempsize)
			{
				innr[index_x] += innr[index_x + tempsize];
			}
			
			tempsize=tempsize/2;
			__syncthreads();		
			
		}
	__syncthreads();
		
	if(index_x == 0)
		ans[ blockIdx.y * gridDim.x + blockIdx.x] = innr[0];

}



