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
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

//////////////// Global  //////////////////////////////////////

	int Baez = 0;	// Selects Cayley Table, rather than John Baez' 

//////////////// Structures //////////////////////////////////////

struct Vector
{
	ex x,y,z;
	Vector(ex xx, ex yy, ex zz) {x = xx; y = yy; z = zz;}
	Vector() {x=0; y=0; z=0;}
};

Vector origin;	// (0,0,0)

struct Quaternion
{
	ex q,i,j,k;
	Quaternion(ex qq, ex ii, ex jj, ex kk) {q = qq; i = ii; j = jj; k = kk;}
	Quaternion(ex qq) {q = qq; i = 0; j = 0; k = 0;}
	Quaternion(Vector &u) {q = 0; i = u.x; j = u.y; k = u.z;}
	Quaternion() {q = 0; i = 0; j = 0; k = 0;}
};

struct Octonion{	// 0 . 1 2 4 . 3 5 6 . 7 
	ex q,  i,j,k, E, I,J,K;
	Octonion() {q = 0; i = 0; j = 0; k = 0; E = 0;	I = 0; J = 0; K = 0;}
	Octonion(ex qq, ex ii, ex jj, ex kk,  ex EE, ex II, ex JJ, ex KK) 
	 {q = qq; i = ii; j = jj; k = kk; E = EE; I = II; J = JJ; K = KK;}
};

///////////// Vector Routines /////////////////////////////////////////

istream &operator>>(istream &ff, Vector &P)
{	return ff >> P.x >> P.y >> P.z; }

ostream &operator<<(ostream &ff, Vector &P) {
	return ff << "(" << P.x << " , " << P.y << " , " << P.z << ")";
}

Vector operator+(Vector u, Vector v)
{	return Vector(u.x + v.x, u.y + v.y, u.z + v.z); }

Vector operator-(Vector u, Vector v)
{	return Vector(u.x - v.x, u.y - v.y, u.z - v.z);  }

Vector operator*(ex c, Vector v)
{	return Vector(c*v.x, c*v.y, c*v.z); 	}

Vector &operator+=(Vector &u, Vector v)
{	u.x += v.x; u.y += v.y; u.z += v.z; return u;	}

Vector &operator-=(Vector &u, Vector v)
{	u.x -= v.x; u.y -= v.y; u.z -= v.z; return u;	}

Vector &operator*=(Vector &v, ex c)
{	v.x *= c; v.y *= c; v.z *= c; return v; }

//Vector &operator/=(Vector &v, ex c)
//{	v.x /= c; v.y /= c; v.z /= c; return v; }

ex mag(Vector u) 
{	return sqrt(u.x*u.x + u.y*u.y + u.z*u.z); }

Vector	Cross(Vector a, Vector b)
{
	Vector c;
	
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;
	return (c);
}

ex Dot(Vector a, Vector b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

///////////// Quaternion Routines /////////////////////////////////////////

istream &operator>>(istream &ff, Quaternion &v) {
	return ff >> v.q >> v.i >> v.j >> v.k;
}

ostream &operator<<(ostream &ff, Quaternion &v) {
	return ff << "(" << v.q << " , " << v.i << " , " << v.j << " , " << v.k << ")";
}

Quaternion operator+(const Quaternion &u, const Quaternion &v)
{
	Quaternion w;
	w.q = expand(u.q+v.q);
	w.i = expand(u.i+v.i);
	w.j = expand(u.j+v.j);
	w.k = expand(u.k+v.k);
	return w;
}

Quaternion operator+(const Quaternion &u, const ex v)
{
	return Quaternion(expand(u.q+v),u.i,u.j,u.k);
}

Quaternion operator+(const ex u, const Quaternion &v)
{
	return Quaternion(expand(u+v.q),v.i,v.j,v.k);
}

Quaternion operator-(const Quaternion &u, const Quaternion &v)
{
	Quaternion w;
	w.q = expand(u.q-v.q);
	w.i = expand(u.i-v.i);
	w.j = expand(u.j-v.j);
	w.k = expand(u.k-v.k);
	return w;
}

ex mag(const Quaternion &u){
	return sqrt(u.q*u.q + u.i*u.i + u.j*u.j + u.k*u.k);
}

Quaternion conj(const Quaternion &u) {
	return Quaternion(u.q, -u.i, -u.j, -u.k);
}

Quaternion operator*(const Quaternion &u, const Quaternion &v) {
	Quaternion w;
	w.q = expand(+u.q*v.q-u.i*v.i-u.j*v.j-u.k*v.k);
	w.i = expand(+u.q*v.i+u.i*v.q+u.j*v.k-u.k*v.j);
	w.j = expand(+u.q*v.j-u.i*v.k+u.j*v.q+u.k*v.i);
	w.k = expand(+u.q*v.k+u.i*v.j-u.j*v.i+u.k*v.q);
	return w;
}

Quaternion operator/(const Quaternion &u, const Quaternion &v) {
	ex mag2;
	Quaternion w,z;
	mag2 = v.q*v.q + v.i*v.i + v.j*v.j + v.k*v.k;
	if (mag2 > 1.0e-28) {
		z.q = v.q/mag2;
		z.i = -v.i/mag2;
		z.j = -v.j/mag2;
		z.k = -v.k/mag2;
		w = u*z;
	}
	return w;
}

Quaternion operator*(const Quaternion &v, const ex c) {
	return Quaternion(v.q*c, v.i*c, v.j*c, v.k*c);
}

Quaternion operator*(const ex c, const Quaternion &v){
	return Quaternion(v.q*c, v.i*c, v.j*c, v.k*c);
}

Quaternion operator/(const Quaternion &v, const ex c) {
	return Quaternion(v.q/c, v.i/c, v.j/c, v.k/c);
}

Quaternion operator/(const ex c, const Quaternion &v){
	return Quaternion(c)/v;
}

Quaternion &operator+=(Quaternion &u, const Quaternion &v) {
	u.q += v.q; u.i += v.i; u.j += v.j; u.k += v.k;
	return u;
}

Quaternion &operator-=(Quaternion &u, const Quaternion &v) {
	u.q -= v.q; u.i -= v.i; u.j -= v.j; u.k -= v.k;
	return u;
}

Quaternion &operator*=(Quaternion &u, const Quaternion &v) {
	u = u*v;   // lame...
	return u;
}

Quaternion exp(const Quaternion &u) {	
// kn - this looks wrong! - actually okay... 	
	ex Mag, theta,c,s;
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
	return Quaternion(Mag*c,Mag*s*v.x,Mag*s*v.y,Mag*s*v.z);
}

Quaternion log(Quaternion &u) {
	ex A, theta;
	A = sqrt(u.i*u.i + u.j*u.j + u.k*u.k);
	theta = atan2(A,u.q)/A;	// put normalization of u by A into theta
	return Quaternion(log(mag(u)), u.i*theta, u.j*theta, u.k*theta);
}

void Q2P(Quaternion &u, Vector &v, ex &theta, ex &Magnitude) {
	ex A;
	Magnitude = mag(u);
	A = sqrt(u.i*u.i + u.j*u.j + u.k*u.k);
	theta = atan2(A,u.q);
	v.x = u.i/A;
	v.y = u.j/A;
	v.z = u.k/A;
}
	
ex Dot(Quaternion a, Quaternion b)
{
	return expand(a.q*b.q + a.i*b.i + a.j*b.j + a.k*b.k);
}

/////////// Octonion Routines ///////////////////////////////////


istream &operator>>(istream &ff, Octonion &v) {
	return ff >> v.q >> v.i >> v.j >> v.k >> v.E >> v.I >> v.J >> v.K;
}

ostream &operator<<(ostream &ff, Octonion &v) {
	return ff << "(" << v.q << " , " << v.i << " , " << v.j << " , " << v.k  << " , " 
		  << v.E << " , " << v.I << " , " << v.J << " , " << v.K << ")";
}

Octonion operator+(const Octonion &u, const Octonion &v)
{
	Octonion w;
	w.q = expand(u.q + v.q);
	w.i = expand(u.i + v.i);
	w.j = expand(u.j + v.j);
	w.k = expand(u.k + v.k);
	w.E = expand(u.E + v.E);
	w.I = expand(u.I + v.I);
	w.J = expand(u.J + v.J);
	w.K = expand(u.K + v.K);
	return w;
}

Octonion operator+(const Octonion &u, const ex v)
{
	Octonion w;
	w = u;
	w.q = u.q + v;
	return w;
}

Octonion operator+(const ex u, const Octonion &v)
{
	Octonion w;
	w = v;
	w.q = v.q + u;
	return w;
}

Octonion operator-(const Octonion &u, const Octonion &v)
{
	Octonion w;
	w.q = expand(u.q - v.q);
	w.i = expand(u.i - v.i);
	w.j = expand(u.j - v.j);
	w.k = expand(u.k - v.k);
	w.E = expand(u.E - v.E);
	w.I = expand(u.I - v.I);
	w.J = expand(u.J - v.J);
	w.K = expand(u.K - v.K);
	return w;
}

ex mag(const Octonion &u){
	return sqrt(u.q*u.q + u.i*u.i + u.j*u.j + u.k*u.k + 
		    u.I*u.I + u.J*u.J + u.K*u.K + u.E*u.E);
}

Octonion conj(const Octonion &u) {
	return Octonion(u.q, -u.i, -u.j, -u.k, -u.E, -u.I, -u.J, -u.K);
}

Octonion operator*(const Octonion &u, const Octonion &v) { // Cayley, not Baez
	Octonion w;
	ex a,b,c,d,e,f,g,h,A,B,C,D,E,F,G,H;

	a = u.q;	b = u.i;	c = u.j;	d = u.k;
	e = u.E;	f = u.I;	g = u.J;	h = u.K;

	A = v.q;	B = v.i;	C = v.j;	D = v.k;
	E = v.E;	F = v.I;	G = v.J;	H = v.K;	

	if(Baez) {
		w.q = expand(a*A - b*B - c*C - d*D - e*E - f*F - g*G - h*H) ; // Baez
		w.i = expand(a*B + b*A + c*E + d*H - e*C + f*G - g*F - h*D) ;
		w.j = expand(a*C - b*E + c*A + d*F + e*B - f*D + g*H - h*G) ;
		w.k = expand(a*D - b*H - c*F + d*A + e*G + f*C - g*E + h*B) ;
		w.E = expand(a*E + b*C - c*B - d*G + e*A + f*H + g*D - h*F) ;
		w.I = expand(a*F - b*G + c*D - d*C - e*H + f*A + g*B + h*E) ;
		w.J = expand(a*G + b*F - c*H + d*E - e*D - f*B + g*A + h*C) ;
		w.K = expand(a*H + b*D + c*G - d*B + e*F - f*E - g*C + h*A) ;
	} else {		
		w.q = expand(a*A - b*B - c*C - d*D - e*E - f*F - g*G - h*H) ; // Cayley
		w.i = expand(a*B + b*A + c*D - d*C + e*F - f*E - g*H + h*G) ;
		w.j = expand(a*C - b*D + c*A + d*B + e*G + f*H - g*E - h*F) ;
		w.k = expand(a*D + b*C - c*B + d*A + e*H - f*G + g*F - h*E) ;
		w.E = expand(a*E - b*F - c*G - d*H + e*A + f*B + g*C + h*D) ;
		w.I = expand(a*F + b*E - c*H + d*G - e*B + f*A - g*D + h*C) ;
		w.J = expand(a*G + b*H + c*E - d*F - e*C + f*D + g*A - h*B) ;
		w.K = expand(a*H - b*G + c*F + d*E - e*D - f*C + g*B + h*A) ;
	}

	return w;
}

Octonion operator/(const Octonion &u, const Octonion &v) {
	ex mag2;
	Octonion w,z;
	mag2 = v.q*v.q + v.i*v.i + v.j*v.j + v.k*v.k + v.I*v.I + v.J*v.J + v.K*v.K + v.E*v.E ;
	if (mag2 > 1.0e-28) {
		z.q =  expand(v.q/mag2);
		z.i = -expand(v.i/mag2);
		z.j = -expand(v.j/mag2);
		z.k = -expand(v.k/mag2);
		z.E = -expand(v.E/mag2);
		z.I = -expand(v.I/mag2);
		z.J = -expand(v.J/mag2);
		z.K = -expand(v.K/mag2);
		w = u*z;
	}
	return w;
}

Octonion operator*(const Octonion &v, const ex c) {

	Octonion w;

	w.q = expand(v.q*c);
	w.i = expand(v.i*c);
	w.j = expand(v.j*c);
	w.k = expand(v.k*c);
	w.E = expand(v.E*c);
	w.I = expand(v.I*c);
	w.J = expand(v.J*c);
	w.K = expand(v.K*c);

	return w;
}

Octonion operator*(const ex c, const Octonion &v){

	Octonion w;

	w.q = expand(v.q*c);
	w.i = expand(v.i*c);
	w.j = expand(v.j*c);
	w.k = expand(v.k*c);
	w.E = expand(v.E*c);
	w.I = expand(v.I*c);
	w.J = expand(v.J*c);
	w.K = expand(v.K*c);

	return w;
}

Octonion operator/(const Octonion &v, const ex c) {

	Octonion w;

	w.q = expand(v.q/c);
	w.i = expand(v.i/c);
	w.j = expand(v.j/c);
	w.k = expand(v.k/c);
	w.E = expand(v.E/c);
	w.I = expand(v.I/c);
	w.J = expand(v.J/c);
	w.K = expand(v.K/c);

	return w;
}

Octonion operator/(const ex c, const Octonion &v){
	Octonion w,z;
	w.q = c;
	z = w/v;
	return z;
}

Octonion &operator+=(Octonion &u, const Octonion &v) {
	u.q += v.q; u.i += v.i; u.j += v.j; u.k += v.k;
	u.E += v.E; u.I += v.I; u.J += v.J; u.K += v.K;
	return u;
}

Octonion &operator-=(Octonion &u, const Octonion &v) {
	u.q -= v.q; u.i -= v.i; u.j -= v.j; u.k -= v.k;
	u.E -= v.E; u.I -= v.I; u.J -= v.J; u.K -= v.K;
	return u;
}

Octonion &operator*=(Octonion &u, const Octonion &v) {
	u = u*v;  
	return u;
}
	
ex Dot(Octonion a, Octonion b)
{
	return a.q*b.q + a.i*b.i + a.j*b.j + a.k*b.k + a.E*b.E + a.I*b.I + a.J*b.J + a.K*b.K;
}


//////////////////////////////////////////////////////

	Octonion Commutator(Octonion a, Octonion b)
{
	return (a*b - b*a);
}

//////////////////////////////////////////////////////

	Octonion Jacobi(Octonion a, Octonion b, Octonion c)
{
	Octonion d,e,f,g;

//	g = [[a,b],c] + [[b,c],a] + [[c,a],b] ;
	
//	d = Commutator(Commutator(a,b),c);
//	e = Commutator(Commutator(b,c),a);
//	f = Commutator(Commutator(c,a),b);

// [A,[B,C]]+[B,[C,A]]+[C,[A,B]]=0

	d = Commutator(a, Commutator(b,c));
	e = Commutator(b, Commutator(c,a));
	f = Commutator(c, Commutator(a,b));


	g = d + e + f;

	return g;
}


//////////////////////////////////////////////////////

	Octonion Associator(Octonion a, Octonion b, Octonion c)
{
	Octonion d;
	d = (a*b)*c - a*(b*c);
	return d;
}

//////////////////////////////////////////////////////

int operator==(const Octonion &u, const Octonion &v)
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

int operator!=(const Octonion &u, const Octonion &v)
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

int operator==(const Quaternion &u, const Quaternion &v)
{
	int result;
	result = 	(u.q == v.q )&&

			(u.i == v.i )&&
			(u.j == v.j )&&
			(u.k == v.k );
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const Quaternion &u, const Quaternion &v)
{
	int result;
	result = 	(u.q != v.q )||

			(u.i != v.i )||
			(u.j != v.j )||
			(u.k != v.k );
	return result;
}


//////////////////////////////////////////////////////

int operator==(const Vector &u, const Vector &v)
{
	int result;
	result = 	(u.x == v.x )&&
			(u.y == v.y )&&
			(u.z == v.z );
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const Vector &u, const Vector &v)
{
	int result;
	result = 	(u.x != v.x )||
			(u.y != v.y )||
			(u.z != v.z );
	return result;
}



