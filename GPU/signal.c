#include<stdio.h>
#include<math.h>
#include<stdlib.h>

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


float signal_cont(double fr, double t)
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

