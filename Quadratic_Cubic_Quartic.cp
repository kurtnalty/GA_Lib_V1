// Routines for Solving Quadratic, Cubic and Quartic polynomials 
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 31 July 2018
// License: Freeware
// Alternative License: BSD
//
//////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

//////////////////////////////////////////////////////

// ************************************************************************

//	https://en.wikipedia.org/wiki/Cubic_function

//	https://en.wikipedia.org/wiki/Quartic_function

// ************************************************************************

void Program_Failed (const char* desc)	//	Report error and quit
{
	cout <<  desc << "\n";
	exit(0);
}

// ************************************************************************

int NearZero(ex a)
{
	double epsilon = 1.0e-13; // KN see how tight this should be.
	return((abs(a)>-epsilon) && (abs(a)<epsilon));

}

// ************************************************************************

ex ccbrt(ex a)	// return ex cube root
{
	return exp(log(a)/3.0);
}

// ************************************************************************

void Build_Quadratic(ex* a, ex* b, ex* c, ex r1, ex r2)
{
	*a = 1.0;
	*b = -(r1 + r2);
	*c = r1*r2;

}

// ************************************************************************

void Solve_Quadratic(ex a, ex b, ex c, ex* r1, ex* r2)
{
	ex B, C, tmp;

//	if((real_part(a) == 0.0) && (imag_part(a) == 0.0)) Program_Failed ("\n\nSolve_Quadratic failed due to 'a == 0'\n");
	if(NearZero(abs(a))  ) Program_Failed ("\n\nSolve_Quadratic failed due to 'a == 0'\n");
	B = 0.5*b/a;	C = c/a;
	*r1 = -B;	*r2 = -B;	// offset
	*r1 += +sqrt(B*B - C);
	*r2 += -sqrt(B*B - C);

}


// ************************************************************************

void Build_Cubic(ex* a, ex* b, ex*c, ex *d, ex r1, ex r2, ex r3)
{
	*a = 1.0;
	*b = -(r1 + r2 + r3);
	*c = r1*r2 + r1*r3 + r2*r3;
	*d = -r1*r2*r3;

}

// ************************************************************************

void Solve_Cubic(ex a, ex b, ex c, ex d, ex* r1, ex* r2, ex *r3)
{
	ex B, C, D;
	ex p, q, zeta1, zeta2, tmp;

	if(NearZero(real_part(a)) && NearZero(imag_part(a)) ) Program_Failed ("\n\nSolve_Cubic failed due to 'a == 0'\n");
	B = b/a;	C = c/a;	D = d/a;	// normalize
	*r1 = -B/3.0;	*r2 = -B/3.0;	*r3 = -B/3.0;	// offset solution to create reduced cubic: r = (t - B/3)

//	r^3 + B*r^2 + C*r + D = 0 becomes t^3 + p*t + q = 0
	
	p = C - ((B*B)/3.0);
	q = (2.0*B*B*B/27.0) - (B*C/3.0) + D;

//	check for trivial case where p = 0, q = 0 (triple root, t = 0, all r = -B/3)
	if( (real_part(p) == 0.0) && (imag_part(p) == 0.0) && (real_part(q) == 0.0) && (imag_part(q) == 0.0) ) return;
	if( NearZero(real_part(p)) &&  NearZero(imag_part(p)) &&  NearZero(real_part(q)) &&  NearZero(imag_part(q)) ) return;

	zeta1 = -0.5 + I*0.5*sqrt(3.0);	// first ex cube root of 1
	zeta2 = -0.5 - I*0.5*sqrt(3.0); // second ex cube root of 1

//	Vieta's substitution method 

//	t = w - (p/(3*w)) yields quadratic equation   w^6 + q*w^3 - p^3/27 = 0

	ex w, w3_1, w3_2, p_3_27;	// w3 = w^3
	p_3_27 = -p*p*p/27.0;

	Solve_Quadratic(1, q, p_3_27, &w3_1, &w3_2);

	w = ccbrt(w3_1);	// get simple ex cube root
	*r1 += w - (p/(3.0*w));
	*r2 += w*zeta1 - (p/(3.0*w*zeta1));
	*r3 += w*zeta2 - (p/(3.0*w*zeta2));


}


// ************************************************************************

void Build_Quartic(ex* a, ex* b, ex*c, ex *d, ex *e, ex r1, ex r2, ex r3, ex r4)
{
	*a = 1.0;
	*b = -(r1 + r2 + r3 + r4);
	*c = r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4;
	*d = -r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4;
	*e = r1*r2*r3*r4;

}

// ************************************************************************

void Solve_Quartic(ex a, ex b, ex c, ex d, ex e, ex* r1, ex* r2, ex *r3, ex *r4)
{
	ex B, C, D, E;
	ex p, q, r, m, m1, m2, m3, tmp;
	ex t1,t2;
	ex y1,y2,y3,y4,sqrt2m;

	if(NearZero(real_part(a)) && NearZero(imag_part(a))) Program_Failed ("\n\nSolve_Quartic failed due to 'a == 0'\n");
	B = b/a;	C = c/a;	D = d/a;	E = e/a;	// normalize
	*r1 = -B/4.0;	*r2 = -B/4.0;	*r3 = -B/4.0;	*r4 = -B/4.0;	// offset solution to create reduced quartic: r = (t - B/4)

//	r^4 + B*r^4 + C*r^2 + D*r + E = 0 becomes t^4 + p*t^2 + q*t + r = 0
	
	p = C - ((0.375*B*B));	// C - 3B^2/8
	q = (0.125*B*B*B) - (0.5*B*C) + D;
	r = (-3.0*B*B*B*B/256.0) + E + (-0.25*B*D) + (0.0625*B*B*C);

//		printf("debug . . . p = (%lg + I %lg) \n",   real_part(p), imag_part(p));
//		printf("debug . . . q = (%lg + I %lg) \n",   real_part(q), imag_part(q));
//		printf("debug . . . r = (%lg + I %lg) \n\n", real_part(r), imag_part(r));


//	check for trivial case where p = 0, q = 0, r = 0 (quad root, t = 0, all r = -B/4)
	if( NearZero(real_part(p)) &&  NearZero(imag_part(p)) 
		&&  NearZero(real_part(q)) &&  NearZero(imag_part(q))
		&&  NearZero(real_part(r)) &&  NearZero(imag_part(r)) ) return;

//	check for the case where q = 0 (biquadratic)  t^4 + p*t^2 + r = 0
	if( NearZero(real_part(q)) &&  NearZero(imag_part(q)) ) {
		Solve_Quadratic(1, p, r, &t1, &t2);	// solve for t^2
		y1 =  sqrt(t1);  y2 = -sqrt(t1);  y3 =  sqrt(t2);  y4 = -sqrt(t2);
		*r1 += y1;	*r2 += y2;	*r3 += y3;	*r4 += y4;
	//	sort roots ascending in real component
	
		if(real_part(*r2) < real_part(*r1)) {tmp = *r1; *r1 = *r2; *r2 = tmp; };
		if(real_part(*r3) < real_part(*r1)) {tmp = *r1; *r1 = *r3; *r3 = tmp; };
		if(real_part(*r4) < real_part(*r1)) {tmp = *r1; *r1 = *r4; *r4 = tmp; };	// should have smallest in *r1
	
		if(real_part(*r3) < real_part(*r2)) {tmp = *r2; *r2 = *r3; *r3 = tmp; };
		if(real_part(*r4) < real_part(*r2)) {tmp = *r2; *r2 = *r4; *r4 = tmp; };
	
		if(real_part(*r4) < real_part(*r3)) {tmp = *r3; *r3 = *r4; *r4 = tmp; };
		return;
	} 

//	Lodovico Ferrari's method. Solve for term m
// 	m^3 + p*m^2 + (-r + 0.25*p^2)*m + (-0.125*q^2) = 0 // resolvent cubic equation
	Solve_Cubic(1, p, (-r + 0.25*p*p), (-0.125*q*q), &m1, &m2, &m3);

//		printf("debug . . . m1 = (%lg + I %lg) \n",   real_part(m1), imag_part(m1));
//		printf("debug . . . m2 = (%lg + I %lg) \n",   real_part(m2), imag_part(m2));
//		printf("debug . . . m3 = (%lg + I %lg) \n\n", real_part(m3), imag_part(m3));

	if(!NearZero(real_part(m1)) || !NearZero(imag_part(m1)) ) m = m1; 
		else if(!NearZero(real_part(m2)) || !NearZero(imag_part(m2)) ) m = m2; 
		else m = m3; // use a NON-ZERO root m
//		printf("debug . . . m  = (%lf + I %lf) \n\n", real_part(m), imag_part(m));
	
	sqrt2m = sqrt(2.0*m);	// convenience variable
//
//	The reduced quartic is now a product of two quadratic equations
//	y^2 + y*( sqrt2m) + (0.5*p + m - (0.5*q/sqrt2m)) = 0  and
//	y^2 + y*(-sqrt2m) + (0.5*p + m + (0.5*q/sqrt2m)) = 0 

	Solve_Quadratic(1, ( sqrt2m), (0.5*p + m - (0.5*q/sqrt2m)), &y1, &y2 );
	Solve_Quadratic(1, (-sqrt2m), (0.5*p + m + (0.5*q/sqrt2m)), &y3, &y4 );

	*r1 += y1;	*r2 += y2;	*r3 += y3;	*r4 += y4;

//	sort roots ascending in real component

	if(real_part(*r2) < real_part(*r1)) {tmp = *r1; *r1 = *r2; *r2 = tmp; };
	if(real_part(*r3) < real_part(*r1)) {tmp = *r1; *r1 = *r3; *r3 = tmp; };
	if(real_part(*r4) < real_part(*r1)) {tmp = *r1; *r1 = *r4; *r4 = tmp; };	// should have smallest in *r1

	if(real_part(*r3) < real_part(*r2)) {tmp = *r2; *r2 = *r3; *r3 = tmp; };
	if(real_part(*r4) < real_part(*r2)) {tmp = *r2; *r2 = *r4; *r4 = tmp; };

	if(real_part(*r4) < real_part(*r3)) {tmp = *r3; *r3 = *r4; *r4 = tmp; };

}

