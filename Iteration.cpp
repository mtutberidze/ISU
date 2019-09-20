// Iteration.cpp : Defines the entry point for the console application.
// Author: Mikheil Tutberidze, Associate Professor at Ilia State University
// Application for numerical computation of the approximate solution of one nonlinear parabolic equation.

#include "stdafx.h"

#include <math.h>

long double k(long double du)
{
	if (du==0) return 1;
	return atan(du)/du;
}

long double U(long double x, long double t)
{
 return exp(t)*cos((4*atan(1.0l))*x);//exp(t)*sin((4*atan(1.0l))*x);
}

long double f(long double x, long double t, long double u)
{
 long double dU=-(4*atan(1.0l))*exp(t)*sin((4*atan(1.0l))*x);//(4*atan(1.0l))*exp(t)*cos((4*atan(1.0l))*x);
 long double d2U=(-(4*atan(1.0l))*(4*atan(1.0l))*exp(t)*cos((4*atan(1.0l))*x));//(-(4*atan(1.0l))*(4*atan(1.0l))*exp(t)*sin((4*atan(1.0l))*x));
 long double UU=U(x,t);
		
 return exp(t)*cos((4*atan(1.0l))*x)-log(u+sqrt(u*u+1))+log(UU+sqrt(UU*UU+1))-
	 d2U/(dU*dU+1)    + x*t;
}



long double phi (long double x) 
{
 return cos((4*atan(1.0l))*x)+x*(1-x);//sin((4*atan(1.0l))*x);;
}

long double psi0 (long double t) 
{
 return exp(t)+t; //0;
}

long double psi1 (long double t) 
{
 return -exp(t)+t; //0;
}


const int M=50;
const int N=100;
const long double T=2.0;
const long double eps=0.0000001;
const long double L=2.0;

FILE *fout;

long double u[M+1];






int _tmain(int argc, _TCHAR* argv[])
{
	const long double h=1.0/M;
	const long double tau=T/N;
	const long double theta=h*h/(h*h+tau*L*(2+h*h)); 
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
				y[i]=(1-theta)*y0[i]+theta*(tau*((k((y0[i+1]-y0[i])/h)*(y0[i+1]-y0[i])-k((y0[i]-y0[i-1])/h)*(y0[i]-y0[i-1]))/(h*h) + f(x,t,y0[i])) + u[i]);
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
			fprintf (fout, "%0.17lf\t", u[i]);
		}
		fprintf (fout,"\n");
	}
	fclose(fout);
	return 0;
}