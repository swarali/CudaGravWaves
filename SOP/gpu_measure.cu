#include<stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cufft.h>
#include <math.h>
#include<stdlib.h>
#include "cublas_v2.h"

/**Correlation coefficient GPU code by using CUBLAS library function - cublasSgemv : which matrix and vector product
*
	Calculates the correlation coefficient of a sine wave of given frequency with a given number of sine waves with 
	frequency fr +/- somevalue... -0.5<somevalue<0.5 

	Uses cuda fft library for fft on the device
	Input : sampling frequency, time period, number of samples and frequency(samp * t should be power of 2)
*/

#define PI 3.14159265358979323846

float signal(double fr, double t);
void cudaErrorCheck(char * errormsg);


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
	int deviceCount;
	int i, j, k=0, temp1;
	int t, batch;
	int nx, NUM;
	float samp, fr, tempfr;
	float  *d_device_cor, *h_device_cor;
	cufftHandle plan;
//	cufftComplex *d_temp1, *h_temp1;	
//	cufftComplex *h_data1, *d_data1;
	float temp, sum;
//	float memsettime;
//	cudaEvent_t start,stop;
	FILE *file;

	cudaGetDeviceCount(&deviceCount);
	printf("Number of devices on this node = %d\n", deviceCount);
	file = fopen("cor_time1gpu.txt","w"); 

	printf("Frequency : ");
	scanf("%f", &fr);
//	fr =10;
		
//	printf("Sampling frequency : ");
//	scanf("%f", &samp);
//	samp =64;

	printf("Time period : ");
	scanf("%d", &t);
//	t=8;

	printf("Number of tests : ");
	scanf("%d", &NUM);
//	NUM=10;

	printf("Starting number of Samples : ");
	scanf("%d", &batch);
//	batch =100;

//	batch-=10;
	k= NUM;
	while(k--)
	{
	
		cufftComplex **d_temp1, **h_temp1;	
		cufftComplex *h_data1, **d_data1;
	
		d_temp1 = (cufftComplex * *)malloc(sizeof(cufftComplex *)*deviceCount);
		h_temp1 = (cufftComplex * *)malloc(sizeof(cufftComplex *)*deviceCount);
	
		d_data1 = (cufftComplex * *)malloc(sizeof(cufftComplex *)*deviceCount);
	

		temp1 = 2*(int)fr + 1;
		temp1 = t*temp1;
	
		int count =0;
		count++;
		while(temp1=temp1>>1)
			count++;
	
		samp = (float)(1<<count)/(float)t;
		nx = (1<<count);
	
	//	nx = samp *t;
		h_device_cor= (float *)malloc(sizeof(float) * batch);	
	
		h_data1 = (cufftComplex *)malloc(sizeof(cufftComplex) * nx );
	
		for(i=0;i<deviceCount;i++)
		{
			h_temp1[i] = (cufftComplex *)malloc(sizeof(cufftComplex) * nx *( (i+1)*batch/deviceCount - i*batch/deviceCount));			
		}
			
		sum =0.0f;
		for(i=0;i<nx;i++)
		{
	//		h_data1[i].x = sin(2*PI*i*fr/samp)/sqrt((float)nx/2);
			temp = signal((double)fr, (double)i/samp);
			h_data1[i].x = temp;
			h_data1[i].y = 0;
			sum+=temp*temp;	
		}
		
		//normalising the data
		for(i=0;i<nx;i++)
			h_data1[i].x /=(sqrt(nx*sum));
	
		int dev;
			
		for(dev=0;dev<deviceCount;dev++)
		//for(j=0;j<batch;j++)
		{	
			j=0;
			while(j<( (dev+1)*batch/deviceCount - dev*batch/deviceCount))
			{
				tempfr = (fr + (float)(j +1)/batch -0.5);
				sum =0.0f;
			
				for(i=0;i<nx;i++)
				{
	//				temp = sin(2*PI*i*(fr -0.5 + (float)j/batch)/samp)/sqrt((float)nx/2);				
					temp = signal((double)tempfr, (double)i/samp);	
					h_temp1[dev][j*nx + i].x = temp;
					h_temp1[dev][j*nx + i].y = 0;
					sum += temp*temp;
				}
	//				normalising the template values	
				for(i=0;i<nx;i++)
					h_temp1[dev][j*nx+i].x/=(sqrt(nx*sum));
	
				j++;
			}
	
		}

//		cudaEventCreate(&start);
//		cudaEventCreate(&stop);

		for(dev=0;dev<deviceCount;dev++)
		{
		
		cudaSetDevice(dev);
		cudaThreadSynchronize();
//		cudaEventRecord(start,0);
	
		cudaMalloc((void**)&d_temp1[dev], ( (dev+1)*batch/deviceCount - dev*batch/deviceCount) * nx * sizeof(cufftComplex));
		cudaMalloc((void**)&d_data1[dev],  nx * sizeof(cufftComplex));
	
		cudaMemcpy(d_temp1[dev], h_temp1[dev],( (dev+1)*batch/deviceCount - dev*batch/deviceCount) * nx * sizeof(cufftComplex), cudaMemcpyHostToDevice);
		cudaMemcpy(d_data1[dev], h_data1, nx * sizeof(cufftComplex), cudaMemcpyHostToDevice);

		cudaErrorCheck("Memory Copy Error");
	
		//FFT of temp
		cufftPlan1d(&plan, nx, CUFFT_C2C, batch);
		cufftExecC2C(plan, d_temp1[dev], d_temp1[dev],CUFFT_FORWARD);
		cufftDestroy(plan);
		
		//FFT of data
		cufftPlan1d(&plan, nx, CUFFT_C2C, 1);
		cufftExecC2C(plan, d_data1[dev], d_data1[dev],CUFFT_FORWARD);
		cufftDestroy(plan);
	
		//Correlation by taking the vector-matrix product
		float al =1.0f, bet =0.0f;
		cublasHandle_t handle ; // CUBLAS context
		cublasCreate (& handle );
		cudaMalloc((void**)&d_device_cor,  batch * sizeof(float));	
		
		cublasSgemv(handle,CUBLAS_OP_T,2*nx,( (dev+1)*batch/deviceCount - dev*batch/deviceCount),&al,(float *)d_temp1[dev],2*nx,(float *)d_data1[dev],1,&bet,&d_device_cor[dev*batch/deviceCount],1);

		cudaMemcpy(&h_device_cor[dev*batch/deviceCount], &d_device_cor[dev*batch/deviceCount], (( (dev+1)*batch/deviceCount - dev*batch/deviceCount) * sizeof(float)), cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();	
		}
//		cudaEventRecord(stop,0);
//		cudaThreadSynchronize();	
		
//		cudaEventElapsedTime(&memsettime, start, stop);

		//printing correlation coefficients
		printf("\n Batch : %d\n", batch);
		for(i=0;i<batch;i++)
			printf("%d : fr : %f : %f \n",i, (fr + (float)(i +1)/batch -0.5), h_device_cor[i]);

//		printf("%d  %f\n",batch, memsettime/1000);
//		fprintf(file,"%d\t%lf\n",batch,memsettime/1000);

//		cufftComplex **d_temp1, **h_temp1;	
//		cufftComplex *h_data1, **d_data1;
	
		cudaErrorCheck("On line 235\n");
	
		for(i=0;i<deviceCount;i++)
		{
			cudaFree(d_temp1[i]);
			free(h_temp1[i]);
			cudaFree(d_data1[i]);
					
		}
		
		cudaErrorCheck("On line 245\n");
		
		free(d_temp1);
		free(h_temp1);
		free(d_data1);
		free(h_data1);
		free(h_device_cor);
		cudaFree(d_device_cor);
		
		batch+=10;	
	}

	fclose(file);
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
