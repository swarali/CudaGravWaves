#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/**
* generates continuos gravitational wave data in time and frequency domain from the input frequency and time
* Uses four1 function for calculating the fft
*Input : frequency and time
*Sampling frequency is >= 2*fr +1, adjusted to get number of datapoints in form 2^n
*/

/*
#define T 32
#define SAMP 64
#define FR 10.00
*/

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
#define H0 0

void plotdata(char fname[], char imgname[]);

extern void four1(double data[], unsigned long nn, int isign);

int Datapoints, FR, Datapoints;
double SAMP, T;

long double signal(double fr, double t)
{
	long double ans, Z, N, xi_rot, P, Q, delta, R;

	Z=2*PI*fr*Rse*sin(theta)/c_light;
	xi_rot = omega_rot*t;
	
	P= 2*PI*fr*Re*sin(alpha)*(cos(beta0)* (sin(theta)*cos(epsilon)*sin(phi) + cos(theta)*sin(epsilon))  - (sin(beta0)*sin(theta)*cos(phi)))/c_light;
	Q= 2*PI*fr*Re*sin(alpha)*(sin(beta0)* (sin(theta)*cos(epsilon)*sin(phi) + cos(theta)*sin(epsilon))  + (cos(beta0)*sin(theta)*cos(phi)))/c_light;
	N= sqrt(P*P + Q*Q);
	delta = atan(P/Q);
	
	R=Z *cos(phi);

	ans = 2*PI*fr*t - (2*PI*fr*Re/c_light)*(1 - cos(omega_rot*t));
//	printf("\n%lf",(double) P);
//	ans = 2*PI*fr*t +Z*cos((xi_rot/365.26) - phi) + N*cos(xi_rot - delta) - R -Q; // omega_orb/omega_rot = 1/365.26
//	ans = 2*PI*fr*t + N*cos(xi_rot - delta) - R -Q; // omega_orb/omega_rot = 1/365.26
//	ans = 2*PI*fr*t;
	return (H0 + cos(ans));
}

int main()
{

	double *dat;
	long double temp;
	int i, count=0, temp1;

	printf("Frequency : \n");
	scanf("%d", &FR);
	
	printf("Time-period : \n");
	scanf("%lf", &T);

	temp1 = 2*FR + 1;
	temp1 = T*temp1;

	count++;
	while(temp1=temp1>>1)
		count++;

	SAMP = (double)(1<<count)/T;
	Datapoints = (1<<count);

	dat = (double *)malloc(sizeof(double) * Datapoints);

	FILE *file= NULL;

	printf("Starting to compute values for fr = %d and t = %lf with sampling fr = %lf Number of datapoints = %d...\n", FR,T, SAMP, Datapoints);
	
	file = fopen("wavedata1.txt","w"); 

	for(i=0;i<Datapoints;i++)
	{
		dat[i] = signal(FR,(long double)i/SAMP);
		fprintf(file,"\n%lf\t%le",(double)i/SAMP,(double)dat[i]);


	}
	fclose(file);


	printf("Starting fourier transform ....\n");
	four1(dat, Datapoints/2, 1);
		
	file = fopen("wavedata2.txt","w"); 

	for(i=0;i<Datapoints/2 -1;i++)
		{
		fprintf(file,"\n%lf\t%le\t%le",(double)(i)/(T),dat[2*i+1], (double)dat[2*i+2]);
		}
	fclose(file);


	free(dat);
	return 0;
}

