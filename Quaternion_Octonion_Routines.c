// Routines for Vector, Quaternion and Octonion Algebra
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 9 August 2018
// License: Freeware
// Alternative License: BSD
//
//////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//////////////// Global  //////////////////////////////////////

	int Baez = 0;	// Selects Cayley Table, rather than John Baez' 

//////////////// Structures //////////////////////////////////////


typedef struct
{
	double x;
	double y;
	double z;
} Vector;

typedef struct
{
	double q;
	double i;
	double j;
	double k;
} Quaternion;

typedef struct
{	
	double q;
	double i;
	double j;
	double k;
	double E;
	double I;
	double J;
	double K;
} Octonion;

//////////////////////////////////////////////////////

void Program_Failed (char desc[255])	// Report error and quit
{
	printf("%s \n", desc);
	exit(0);
}

///////////// Vector Routines /////////////////////////////////////////

Vector Zero_Vector(void)		// initializer
{
	Vector a;

	a.x = 0.0; a.y = 0.0; a.z = 0.0;

	return a;

}

//////////////////////////////////////////////////////

Vector Set_Vector(	// initializer
	double x, double y, double z)		
{
	Vector a;

	a.x = x; a.y = y; a.z = z;

	return a;
}

//////////////////////////////////////////////////////

void Println_Vector(Vector v) 
{
	printf("(%10.3e,%10.3e,%10.3e) \n", v.x,v.y,v.z);
}

//////////////////////////////////////////////////////

void Print_Vector(Vector v) 
{
	printf("(%10.3e,%10.3e,%10.3e)", v.x,v.y,v.z);
}

//////////////////////////////////////////////////////

Vector Add_Vector(Vector u, Vector v)
{	
	Vector w;
	w.x = u.x + v.x;
	w.y = u.y + v.y;
	w.z = u.z + v.z;

	return w; 
}

/////////// Vector Routines ///////////////////////////

Vector Sub_Vector(Vector u, Vector v)
{	
	Vector w;
	w.x = u.x - v.x;
	w.y = u.y - v.y;
	w.z = u.z - v.z;

	return w; 
}

//////////////////////////////////////////////////////

Vector Const_Times_Vector(double c, Vector v)
{	
	Vector w;
	w.x = c*v.x;
	w.y = c*v.y;
	w.z = c*v.z;

	return w; 
}

//////////////////////////////////////////////////////

Vector Vector_Times_Const(Vector v, double c)
{	
	Vector w;
	w.x = c*v.x;
	w.y = c*v.y;
	w.z = c*v.z;

	return w; 
}

//////////////////////////////////////////////////////

double Mag_Vector(Vector u) 
{	return sqrt(u.x*u.x + u.y*u.y + u.z*u.z); }

//////////////////////////////////////////////////////

Vector	Cross(Vector a, Vector b)
{
	Vector c;
	
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;
	return (c);
}

//////////////////////////////////////////////////////

double Dot_Vector(Vector a, Vector b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

//////////////////////////////////////////////////////

int EQ_Vector(Vector u, Vector v)
{
	int result;
	result = 	(u.x == v.x )&&
			(u.y == v.y )&&
			(u.z == v.z );
	return result;
}

//////////////////////////////////////////////////////

int NE_Vector(Vector u, Vector v)
{
	int result;
	result = 	(u.x != v.x )||
			(u.y != v.y )||
			(u.z != v.z );
	return result;
}

///////////// Quaternion Routines ////////////////////

Quaternion Zero_Quaternion(void)   // initializer
{
	Quaternion a;

	a.q = 0.0;
	a.i = 0.0; a.j = 0.0; a.k = 0.0;

	return a;
}

//////////////////////////////////////////////////////

Quaternion Set_Quaternion(	// initializer
	double q,
	double i, double j, double k)
{
	Quaternion a;

	a.q = q;
	a.i = i; a.j = j; a.k = k;

	return a;
}

//////////////////////////////////////////////////////

void Println_Quaternion(Quaternion v) 
{
	printf("(%10.3e,  %10.3e,%10.3e,%10.3e) \n ",v.q, v.i,v.j,v.k);
}

//////////////////////////////////////////////////////

void Print_Quaternion(Quaternion v) 
{
	printf("(%10.3e,  %10.3e,%10.3e,%10.3e) ",v.q, v.i,v.j,v.k);
}

//////////////////////////////////////////////////////

Quaternion Add_Quaternion(Quaternion u, Quaternion v)
{
	Quaternion w;
	w.q = u.q + v.q;
	w.i = u.i + v.i;
	w.j = u.j + v.j;
	w.k = u.k + v.k;
	return w;
}

//////////////////////////////////////////////////////

Quaternion Sub_Quaternion(Quaternion u, Quaternion v)
{
	Quaternion w;
	w.q = u.q - v.q;
	w.i = u.i - v.i;
	w.j = u.j - v.j;
	w.k = u.k - v.k;
	return w;
}

//////////////////////////////////////////////////////

double Mag_Quaternion(Quaternion u){
	return sqrt(u.q*u.q + u.i*u.i + u.j*u.j + u.k*u.k);
}

//////////////////////////////////////////////////////

Quaternion Conj_Quaternion(Quaternion u) {
	Quaternion w;
	w.q =  u.q;
	w.i = -u.i;
	w.j = -u.j;
	w.k = -u.k;
	return w;
}

//////////////////////////////////////////////////////

Quaternion Product_Quaternion(Quaternion u, Quaternion v) {
	Quaternion w;
	w.q = (+u.q*v.q-u.i*v.i-u.j*v.j-u.k*v.k);
	w.i = (+u.q*v.i+u.i*v.q+u.j*v.k-u.k*v.j);
	w.j = (+u.q*v.j-u.i*v.k+u.j*v.q+u.k*v.i);
	w.k = (+u.q*v.k+u.i*v.j-u.j*v.i+u.k*v.q);
	return w;
}

//////////////////////////////////////////////////////

Quaternion Reciprocol_Quaternion(Quaternion v) {  // (1/v) is well defined
	double mag2;
	Quaternion w;
	mag2 = v.q*v.q + v.i*v.i + v.j*v.j + v.k*v.k;
	if (mag2 > 1.0e-28) {
		w.q = v.q/mag2;
		w.i = -v.i/mag2;
		w.j = -v.j/mag2;
		w.k = -v.k/mag2;
	} else Program_Failed("Reciprocol_Quaternion failed due to almost zero magnitude\n");
	return w;
}

//////////////////////////////////////////////////////

Quaternion Const_Times_Quaternion(double c, Quaternion  v) {
	Quaternion w;
	w.q = v.q*c;
	w.i = v.i*c;
	w.j = v.j*c;
	w.k = v.k*c;
	return w;
}

//////////////////////////////////////////////////////

Quaternion Quaternion_Times_Const(Quaternion  v, double c) {
	Quaternion w;
	w.q = v.q*c;
	w.i = v.i*c;
	w.j = v.j*c;
	w.k = v.k*c;
	return w;
}

//////////////////////////////////////////////////////

Quaternion Quaternion_Divided_Const(const Quaternion v, const double c) {
	Quaternion w;
	w.q = v.q/c;
	w.i = v.i/c;
	w.j = v.j/c;
	w.k = v.k/c;
	return w;
}

//////////////////////////////////////////////////////

Quaternion Exp_Quaternion(const Quaternion u) {	

	double Mag, theta,c,s;
	Vector v;
	Mag = exp(u.q);
	theta = sqrt(u.i*u.i + u.j*u.j + u.k*u.k);  //really v = |u|sin(phi)
	if(theta > 1.0e-20) {
		v.x = u.i/theta;
		v.y = u.j/theta;
		v.z = u.k/theta;
	} else
		v.x = v.y = v.z = 0;
	c = cos(theta);
	s = sin(theta);
	return Set_Quaternion(Mag*c,Mag*s*v.x,Mag*s*v.y,Mag*s*v.z);
}

//////////////////////////////////////////////////////

Quaternion Log_Quaternion(Quaternion u) {
	double A, theta;
	A = sqrt(u.i*u.i + u.j*u.j + u.k*u.k);
	theta = atan2(A,u.q)/A;	// put normalization of u by A into theta
	return Set_Quaternion(log(Mag_Quaternion(u)), u.i*theta, u.j*theta, u.k*theta);
}

//////////////////////////////////////////////////////

void Q2P(Quaternion u, Vector* v, double* theta, double* Magnitude) {
	double A;
	*Magnitude = Mag_Quaternion(u);
	A = sqrt(u.i*u.i + u.j*u.j + u.k*u.k);
	*theta = atan2(A,u.q);
	(*v).x = u.i/A;
	(*v).y = u.j/A;
	(*v).z = u.k/A;
}

//////////////////////////////////////////////////////

double Dot_Quaternion(Quaternion a, Quaternion b)
{
	return (a.q*b.q + a.i*b.i + a.j*b.j + a.k*b.k);
}

//////////////////////////////////////////////////////

int EQ_Quaternion(Quaternion u, Quaternion v)
{
	int result;
	result = 	(u.q == v.q )&&

			(u.i == v.i )&&
			(u.j == v.j )&&
			(u.k == v.k );
	return result;
}

//////////////////////////////////////////////////////

int NE_Quaternion(Quaternion u, Quaternion v)
{
	int result;
	result = 	(u.q != v.q )||

			(u.i != v.i )||
			(u.j != v.j )||
			(u.k != v.k );
	return result;
}

/////////// Octonion Routines ///////////////////////////////////

Octonion Zero_Octonion(void)		// initializer
{
	Octonion a;

	a.q = 0.0;
	a.i = 0.0; a.j = 0.0; a.k = 0.0;
	a.E = 0.0;
	a.I = 0.0; a.J = 0.0; a.K = 0.0;

	return a;
}
//////////////////////////////////////////////////////

Octonion Set_Octonion(	// initializer
	double q,
	double i, double j, double k,
	double E,
	double I, double J, double K)
{
	Octonion a;

	a.q = q;
	a.i = i; a.j = j; a.k = k;
	a.E = E;
	a.I = I; a.J = J; a.K = K;

	return a;
}

//////////////////////////////////////////////////////

void Println_Octonion(Octonion v) 
{
	printf("(%10.3e,\n  %10.3e,%10.3e,%10.3e,\n  %10.3e,\n %10.3e,%10.3e,%10.3e) \n",v.q, v.i,v.j,v.k, v.E,  v.I,v.J,v.K);
}

//////////////////////////////////////////////////////

void Print_Octonion(Octonion v) 
{
	printf("(%10.3e,  %10.3e,%10.3e,%10.3e,  %10.3e,  %10.3e,%10.3e,%10.3e) ",v.q, v.i,v.j,v.k, v.E,  v.I,v.J,v.K);
}

//////////////////////////////////////////////////////

Octonion Add_Octonion(Octonion u, Octonion v)
{
	Octonion w;
	w.q = (u.q + v.q);
	w.i = (u.i + v.i);
	w.j = (u.j + v.j);
	w.k = (u.k + v.k);
	w.E = (u.E + v.E);
	w.I = (u.I + v.I);
	w.J = (u.J + v.J);
	w.K = (u.K + v.K);
	return w;
}

//////////////////////////////////////////////////////

Octonion Sub_Octonion(Octonion u, Octonion v)
{
	Octonion w;
	w.q = (u.q - v.q);
	w.i = (u.i - v.i);
	w.j = (u.j - v.j);
	w.k = (u.k - v.k);
	w.E = (u.E - v.E);
	w.I = (u.I - v.I);
	w.J = (u.J - v.J);
	w.K = (u.K - v.K);
	return w;
}

//////////////////////////////////////////////////////

double Mag_Octonion(Octonion u){
	return sqrt(u.q*u.q + u.i*u.i + u.j*u.j + u.k*u.k + 
		    u.I*u.I + u.J*u.J + u.K*u.K + u.E*u.E);
}

//////////////////////////////////////////////////////

Octonion Conj_Octonion(const Octonion u) {
	return Set_Octonion(u.q, -u.i, -u.j, -u.k, -u.E, -u.I, -u.J, -u.K);
}

//////////////////////////////////////////////////////

Octonion Product_Octonion(const Octonion u, const Octonion v) { // Cayley, not Baez
	Octonion w;
	double a,b,c,d,e,f,g,h,A,B,C,D,E,F,G,H;

	a = u.q;	b = u.i;	c = u.j;	d = u.k;
	e = u.E;	f = u.I;	g = u.J;	h = u.K;

	A = v.q;	B = v.i;	C = v.j;	D = v.k;
	E = v.E;	F = v.I;	G = v.J;	H = v.K;	

	if(Baez) {
		w.q = (a*A - b*B - c*C - d*D - e*E - f*F - g*G - h*H) ; // Baez
		w.i = (a*B + b*A + c*E + d*H - e*C + f*G - g*F - h*D) ;
		w.j = (a*C - b*E + c*A + d*F + e*B - f*D + g*H - h*G) ;
		w.k = (a*D - b*H - c*F + d*A + e*G + f*C - g*E + h*B) ;
		w.E = (a*E + b*C - c*B - d*G + e*A + f*H + g*D - h*F) ;
		w.I = (a*F - b*G + c*D - d*C - e*H + f*A + g*B + h*E) ;
		w.J = (a*G + b*F - c*H + d*E - e*D - f*B + g*A + h*C) ;
		w.K = (a*H + b*D + c*G - d*B + e*F - f*E - g*C + h*A) ;
	} else {		
		w.q = (a*A - b*B - c*C - d*D - e*E - f*F - g*G - h*H) ; // Cayley
		w.i = (a*B + b*A + c*D - d*C + e*F - f*E - g*H + h*G) ;
		w.j = (a*C - b*D + c*A + d*B + e*G + f*H - g*E - h*F) ;
		w.k = (a*D + b*C - c*B + d*A + e*H - f*G + g*F - h*E) ;
		w.E = (a*E - b*F - c*G - d*H + e*A + f*B + g*C + h*D) ;
		w.I = (a*F + b*E - c*H + d*G - e*B + f*A - g*D + h*C) ;
		w.J = (a*G + b*H + c*E - d*F - e*C + f*D + g*A - h*B) ;
		w.K = (a*H - b*G + c*F + d*E - e*D - f*C + g*B + h*A) ;
	}

	return w;
}

//////////////////////////////////////////////////////

Octonion Reciprocol_Octonion(Octonion v) { // (1/v) is well defined
	double mag2;
	Octonion w;
	mag2 = v.q*v.q + v.i*v.i + v.j*v.j + v.k*v.k + v.I*v.I + v.J*v.J + v.K*v.K + v.E*v.E ;
	if (mag2 > 1.0e-28) {
		w.q =  (v.q/mag2);
		w.i = -(v.i/mag2);
		w.j = -(v.j/mag2);
		w.k = -(v.k/mag2);
		w.E = -(v.E/mag2);
		w.I = -(v.I/mag2);
		w.J = -(v.J/mag2);
		w.K = -(v.K/mag2);
	} else Program_Failed("Reciprocol_Octonion failed due to almost zero magnitude\n");
	return w;
}

//////////////////////////////////////////////////////

Octonion Octonion_Times_Const(Octonion v, const double c) {

	Octonion w;

	w.q = (v.q*c);
	w.i = (v.i*c);
	w.j = (v.j*c);
	w.k = (v.k*c);
	w.E = (v.E*c);
	w.I = (v.I*c);
	w.J = (v.J*c);
	w.K = (v.K*c);

	return w;
}

//////////////////////////////////////////////////////

Octonion Const_Times_Octonion(double c, const Octonion v){

	Octonion w;

	w.q = (v.q*c);
	w.i = (v.i*c);
	w.j = (v.j*c);
	w.k = (v.k*c);
	w.E = (v.E*c);
	w.I = (v.I*c);
	w.J = (v.J*c);
	w.K = (v.K*c);

	return w;
}

//////////////////////////////////////////////////////

Octonion Octonion_Divided_Const(Octonion v, double c) {

	Octonion w;

	w.q = (v.q/c);
	w.i = (v.i/c);
	w.j = (v.j/c);
	w.k = (v.k/c);
	w.E = (v.E/c);
	w.I = (v.I/c);
	w.J = (v.J/c);
	w.K = (v.K/c);

	return w;
}

//////////////////////////////////////////////////////

double Dot_Octonion(Octonion a, Octonion b)
{
	return a.q*b.q + a.i*b.i + a.j*b.j + a.k*b.k + a.E*b.E + a.I*b.I + a.J*b.J + a.K*b.K;
}

//////////////////////////////////////////////////////

int EQ_Octonion(Octonion u, Octonion v)
{
	int result;
	result = 	(u.q == v.q )&&

			(u.i == v.i )&&
			(u.j == v.j )&&
			(u.k == v.k )&&

			(u.E == v.E)&&

			(u.I == v.I)&&
			(u.J == v.J)&&
			(u.K == v.K);
	return result;
}

//////////////////////////////////////////////////////

int NE_Octonion(Octonion u, Octonion v)
{
	int result;
	result = 	(u.q != v.q )||

			(u.i != v.i )||
			(u.j != v.j )||
			(u.k != v.k )||

			(u.E != v.E)||

			(u.I != v.I)||
			(u.J != v.J)||
			(u.K != v.K);
	return result;
}

//////////////////////////////////////////////////////

	Octonion Commutator_Octonion(Octonion a, Octonion b)
{
	Octonion w,x,y;
	w = Product_Octonion(a,b);
	x = Product_Octonion(b,a);
	y = Sub_Octonion(w,x);
	return y;
}

//////////////////////////////////////////////////////

	Octonion Jacobi_Octonion(Octonion a, Octonion b, Octonion c)
{
	Octonion d,e,f,g;

// [A,[B,C]]+[B,[C,A]]+[C,[A,B]]=0

	d = Commutator_Octonion(a, Commutator_Octonion(b,c));
	e = Commutator_Octonion(b, Commutator_Octonion(c,a));
	f = Commutator_Octonion(c, Commutator_Octonion(a,b));

	g = Add_Octonion(Add_Octonion(d ,e), f);

	return g;
}


//////////////////////////////////////////////////////

	Octonion Associator_Octonion(Octonion a, Octonion b, Octonion c)
{
	Octonion d,e,f;
	d = Product_Octonion(Product_Octonion(a,b),c);
	e = Product_Octonion(a,Product_Octonion(b,c));
	f = Sub_Octonion(d,e);
	return f;
}

//////////////////////////////////////////////////////
/*
int	main(void)
{


	return 0;
}

*/


