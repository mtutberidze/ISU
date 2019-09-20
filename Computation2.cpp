// Computation 2
// Author: Mikheil Tutberidze, Associate Professor at Ilia State University
// Application for numerical computation of the approximate solution of one nonlinear parabolic equation.

#include "stdafx.h"

#include <math.h>

const long double pi=4*atan(1.0L);

long double q(long double x, long double t)
{
	return 2*t*x*(1-x);
}

long double g(long double du)
{
	return -0.5 / (du*du + 1) + 0.5;
}

long double U(long double x, long double t)
{
 return exp(t)*cos((pi)*x);//exp(t)*sin((pi)*x);
}

long double f(long double x, long double t)
{
 long double dU=-(pi)*exp(t)*sin((pi)*x);//(pi)*exp(t)*cos((pi)*x);
 long double d2U=(-(pi)*(pi)*exp(t)*cos((pi)*x));//(-(pi)*(pi)*exp(t)*sin((pi)*x));
 long double UU=U(x,t);
		
 return UU-d2U-q(x,t)*dU*dU;
}

long double phi (long double x) 
{
 return cos((pi)*x);
}

long double psi0 (long double t) 
{
 return exp(t);
}

long double psi1 (long double t) 
{
 return -exp(t);
}


const int M=200;
const int N=200;
const long double T=1.0;
const long double eps=0.00001;
const long double L=2;
const long double lambda = 0.5;

FILE *fout;

long double u[M+1];

int _tmain(int argc, _TCHAR* argv[])
{
	const long double h=1.0/M;
	const long double tau=T/N;
	const long double theta=h*h/(h*h+2*tau*L); 
	long double x;
	long double t;
	long double e;
	long double emax;
	long double y[M+1];
	long double y0[M+1];
	int i,j;

	fopen_s(&fout, "C:\\misha.txt", "w");

	x=0;
	for (i=0; i<=M; i++)
	{
	 u[i]=phi(x);
	 fprintf (fout, "%0.17lf\t", u[i]);
	 x=x+h;
	}

	fprintf (fout,"\n"); 

	t=0;
	for (j=1; j<=N; j++)
	{
		t=t+tau;
		for (i=0; i<=M; i++)
			y0[i]=u[i];

		y[0]=psi0(t); 
		y[M]=psi1(t); 
		
		
		do
		{
			emax=0;
			x=0;
			for (i=1; i<M; i++)
			{
				x=x+h;
				y[i]=(1-theta)*y0[i]+theta*(tau*( ( y0[i+1]-2*y0[i]+y0[i-1]+q(x,t)*2*lambda*lambda*g((y0[i+1]-y0[i-1])/(2*lambda)) )/(h*h) + f(x,t)) + u[i]);
				e=abs(y[i]-y0[i]);
				if (e>emax) emax=e;
			}
			for (i=0; i<=M; i++)
				y0[i]=y[i];
		}
		while (emax>=eps);
		
		for (i=0; i<=M; i++)
		{
			u[i]=y[i];
			fprintf (fout, "%0.17lf\t", u[i]-U(i*h,t));
		}
		fprintf (fout,"\n");
		printf("%d\n", j);
	}
	fclose(fout);
	return 0;
}