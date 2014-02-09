#include <stdio.h>
#include <math.h>

/**
	Gnu plotting of sine wave of frequency 20Hz and its fourier transform and inverse transform of the fourier transform 
	
*/


#define T__SIZE 2
#define SAMP_F 256
#define SAMP_SIZE 513		//size should be a power of 2 + 1

extern void plotdata(char fname[]);
extern void four1(float data[], unsigned long nn, int isign);

void main()
{
	int i;
	float dat1[SAMP_SIZE], fr1, j;
	char fname1[] = "fourtest1.txt";
	char fname2[] = "fourtest2.txt";
	char fname3[] = "fourtest3.txt";

	FILE *file= NULL;
	printf("Frequency : ");
	scanf("%f", &fr1);
	file = fopen("fourtest1.txt","w"); 

	if(file == NULL)
	{
		printf("Could not open the file");
//		exit(0);
	}
	
	for(i=0;i<SAMP_SIZE;i++)
	{
		j=i;
		dat1[i] = sin(2*M_PI *i*fr1/SAMP_F);
		//printf("\n%d\t%f",i,dat1[i]);
		fprintf(file,"%f\t%f\n", (j/SAMP_F),dat1[i]);  
	}

	
	plotdata(fname1);	

	fclose(file); 

	
	file = fopen("fourtest2.txt","w");
	if(file == NULL)
	{
		printf("Could not open the file");
//		exit(0);
	}

// Fourier Transform of dat1	
	four1(dat1, SAMP_SIZE /2, 1);

	for(i=2;i<SAMP_SIZE/2;i++)
	{	j=i;
		//printf("\n%d\t%f",i,dat1[i]);
		fprintf(file,"%f\t%f\n",(j-2)/4,dat1[i]); 
	}
	plotdata(fname2);	
	
	fclose(file); 	


	file = fopen("fourtest3.txt","w");
	if(file == NULL)
	{
		printf("Could not open the file");
//		exit(0);
	}
	

// Inverse fourier transform of the dat1
	four1(dat1, SAMP_SIZE /2, -1);

	for(i=0;i<SAMP_SIZE;i++)
	{	j=i;
		//printf("\n%d\t%f",i,dat1[i]);
		fprintf(file,"%f\t%f\n",j/SAMP_F,dat1[i]/(SAMP_F)); 
	}
	plotdata(fname3);	
	
	fclose(file); 	
}

