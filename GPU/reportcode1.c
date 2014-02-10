
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define PI 3.14159265358979323846

/** Correlation coefficient serial code by taking the inner product of fourier tranform
	
	Calculates the correlation coefficient of a sine wave of given frequency with a given number of sine waves with 
	frequency fr +/- somevalue... -0.5<somevalue<0.5

	Uses four1 external function for fft
Input : sampling frequency, time period, number of samples and frequency(samp * t should be power of 2)
*/
	
float signal(double fr, double t); 
void four1(float data[], unsigned long nn, int isign);
void corhost(float* temp, float * data, float *ans, int nx, int batch)
{
	int i, j;
	float sum;
	
	four1(data, nx, 1);	
	
	for(i=0;i<batch;i++)
	{	
		four1(&temp[i*(2*nx+1)], nx, 1);
		ans[i] =0;
		for(j=0;j<nx;j++)
		{
			ans[i]+= ((temp[i*(2*nx +1)+ 2*j + 1]* data[2*j +1])+ (temp[i*(2*nx +1)+ 2*j + 2] * data[2*j +2])); 
		}
	}

}

int main()
{	
	int i, j, k=0;
	int samp = 32, t = 8, batch = 100;
	float step=1/batch;
	int nx, temp1;
	float fr, tempfr, sum;
	float *serial_temp, *serial_data, *serial_cor ;

	printf("Frequency : ");
	scanf("%f", &fr);

//	printf("Sampling frequency : ");
//	scanf("%d", &samp);

	printf("Time period : ");
	scanf("%d", &t);
	
	printf("Number of Samples : ");
	scanf("%d", &batch);

	temp1 = 2*fr + 1;
	temp1 = t*temp1;
	
	int count =0;
	count++;
	while(temp1=temp1>>1)
		count++;
		
	samp = (double)(1<<count)/t; 
	nx = (1<<count);
//	nx = samp*t;

	serial_temp = (float *)malloc(sizeof(float) * (2*nx +1) * batch);
	serial_data = (float *)malloc(sizeof(float) * (2*nx +1));
	serial_cor	= (float *)malloc(sizeof(float) * batch);
	
	//initialising data
	
	sum =0;
	for(i=0;i<=2*nx;i++)
	{
		if(i<=nx)
//			serial_data[i+1] = sin(2*PI*i*fr/samp)/sqrt((float)nx/2); 
			serial_data[i+1] = signal((double)fr, (double)i/samp);
		else
			serial_data[i+1] = 0;	

			sum+= serial_data[i+1]*serial_data[i+1];
	}
	
	for(i=0;i<=2*nx;i++)	 	
			serial_data[i+1]/=(sqrt(nx*sum));

			
	//initialising temp
	for(j=0;j<batch;j++)
	{
		tempfr = (fr + (float)(j +1)/batch -0.5);	 	
		sum =0;
		for(i=0;i<2*nx;i++)
		{

			if(i<=nx)
//		serial_temp[j*(2*nx +1) +i+1] = sin(2*PI*i*tempfr/samp)/sqrt((float)nx/2); 
		serial_temp[j*(2*nx +1) + i +1] = signal((double)tempfr, (double)i/samp);
			else	
			serial_temp[j*(2*nx +1) +i+1] =0;
			
			sum+=serial_temp[j*(2*nx +1) + i+1]*serial_temp[j*(2*nx +1) + i+1];
		}

		for(i=0;i<2*nx;i++)	 	
			serial_temp[j*(2*nx +1) +i+1]/=(sqrt(nx*sum));
			
	}
	
	corhost(serial_temp, serial_data, serial_cor, nx, batch);
	
	for(i=0;i<batch;i++)
		printf("fr : %f : %f \n", (fr + (float)(i +1)/batch -0.5), serial_cor[i]);

	free(serial_temp );	 	
	free(serial_data);	 	
	free(serial_cor);

	return 0;
}


