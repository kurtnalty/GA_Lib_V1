// Routines for Geometric Algebra in Four Dimensional Euclidean spacetime
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 12 March 2018
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

struct GA4E{	// 0 . 1 2 4 8 . 3 5 6 9 10 12 . 7 11 13 14 . 15
	ex q,  x,y,z,w, xy,xz,yz,xw,yw,zw,  xyz,xyw,xzw,yzw,  xyzw;
	GA4E() {q = 0; x = 0; y = 0; z = 0; w = 0; 
		xy = 0; xz = 0; yz = 0;  xw = 0; yw = 0; zw = 0;
		xyz = 0; xyw = 0; xzw = 0; yzw = 0; xyzw=0;}
	GA4E(ex qq, ex xx, ex yy, ex zz, ex ww, 
		ex xxyy, ex xxzz, ex yyzz, ex xxww, ex yyww, ex zzww,
		ex xxyyzz, ex xxyyww, ex xxzzww, ex yyzzww,   ex xxyyzzww) 

		{q = qq; x = xx; y = yy; z = zz; w = ww;
		 xy = xxyy; xz = xxzz; yz = yyzz; xw = xxww; yw = yyww; zw = zzww; 
		 xyz = xxyyzz;  xyw = xxyyww; xzw = xxzzww; yzw = yyzzww; xyzw=xxyyzzww;}
};


//////////////////////////////////////////////////////

// Necessary forward declarations

GA4E LeftContraction (const GA4E &a, const GA4E &b) ;
GA4E Product(const GA4E &a, const GA4E &b) ;

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, GA4E &v) 
{
	return ff << "\n(" 
		<< v.q << ", \n" 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.w << ", \n" 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << "," 
		<< v.xw << "," 
		<< v.yw << "," 
		<< v.zw << ", \n" 
		<< v.xyz << ","
		<< v.xyw << ","
		<< v.xzw << ","
		<< v.yzw << ", \n"
		<< v.xyzw << ")\n";
}


//////////////////////////////////////////////////////

void PrintMV(GA4E &v) 
{
	cout 	<< "( "
		<< v.q << ", " 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.w << ", " 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << "," 
		<< v.xw << "," 
		<< v.yw << "," 
		<< v.zw << ", " 
		<< v.xyz << ","
		<< v.xyw << ","
		<< v.xzw << ","
		<< v.yzw << ", "
		<< v.xyzw << ")";
}

//////////////////////////////////////////////////////

GA4E OverBar(GA4E a)
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA4E b;		// a blade wedge b blade = pseudovector blade
// OverBar(Blade[i]) = Reciprocal(Blade[i])*wxyz ;

	b.q     =  a.xyzw;	// xyzw

	b.x     = -a.yzw;
	b.y     =  a.xzw;
	b.z     = -a.xyw;
	b.w     =  a.xyz;

	b.xy   =   a.zw;
	b.xz   =  -a.yw;
	b.yz   =   a.xw;
	b.xw   =   a.yz;
	b.yw   =  -a.xz;
	b.zw   =   a.xy;

	b.xyz  = -a.w;
	b.xyw  =  a.z;
	b.xzw  = -a.y;
	b.yzw  =  a.x;

	b.xyzw = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA4E UnderBar(GA4E a)   
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA4E b;		// b blade wedge a blade = pseudovector blade
// UnderBar(Blade[i]) = xyzw*Reciprocal(Blade[i]) ;

	b.q     =  a.xyzw;	// xyzw

	b.x     =  a.yzw;
	b.y     = -a.xzw;
	b.z     =  a.xyw;
	b.w     = -a.xyz;

	b.xy   =   a.zw;
	b.xz   =  -a.yw;
	b.yz   =   a.xw;
	b.xw   =   a.yz;
	b.yw   =  -a.xz;
	b.zw   =   a.xy;

	b.xyz  =  a.w;
	b.xyw  = -a.z;
	b.xzw  =  a.y;
	b.yzw  = -a.x;

	b.xyzw = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA4E Reverse(GA4E a)	
// (A,  B, C, D, E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R,  S) score = 10 N( 2046) Reverse
{
	GA4E b;
	b.q =  a.q;

	b.x = a.x;
	b.y = a.y;
	b.z = a.z;
	b.w = a.w;

	b.xy = -a.xy;
	b.xz = -a.xz;
	b.yz = -a.yz;
	b.xw = -a.xw;
	b.yw = -a.yw;
	b.zw = -a.zw;

	b.xyz = -a.xyz;
	b.xyw = -a.xyw;
	b.xzw = -a.xzw;
	b.yzw = -a.yzw;

	b.xyzw =  a.xyzw;

	return b;
}

//////////////////////////////////////////////////////

GA4E Involution(GA4E a)	
// Corresponds to PT parity transform: x -> -x, y -> -y, z ->-z, t -> -t
// (A, -B,-C,-D,-E,  F, G, H, J, K, L, -M,-N,-P,-R,  S) score =  0 N(30750) Involution
{
	GA4E b;
	b.q =  a.q;

	b.x = -a.x;
	b.y = -a.y;
	b.z = -a.z;
	b.w = -a.w;

	b.xy = a.xy;
	b.xz = a.xz;
	b.yz = a.yz;
	b.xw = a.xw;
	b.yw = a.yw;
	b.zw = a.zw;

	b.xyz = -a.xyz;
	b.xyw = -a.xyw;
	b.xzw = -a.xzw;
	b.yzw = -a.yzw;

	b.xyzw =  a.xyzw;

	return b;
}

//////////////////////////////////////////////////////

GA4E Transpose(GA4E a)	
//(A,  B, C, D,-E, -F,-G,-H, J, K, L, -M, N, P, R, -S) score =  6 N( 3857) Transpose  
{
	GA4E b;
/*	b.q =  a.q;   Minkowski

	b.x =  a.x;
	b.y =  a.y;
	b.z =  a.z;
	b.w = -a.w;

	b.xy = -a.xy;
	b.xz = -a.xz;
	b.yz = -a.yz;
	b.xw =  a.xw;
	b.yw =  a.yw;
	b.zw =  a.zw;

	b.xyz = -a.xyz;
	b.xyw =  a.xyw;
	b.xzw =  a.xzw;
	b.yzw =  a.yzw;

	b.xyzw = -a.xyzw;
*/
	b.q =  a.q;

	b.x =  a.x;
	b.y =  a.y;
	b.z =  a.z;
	b.w =  a.w;

	b.xy = -a.xy;
	b.xz = -a.xz;
	b.yz = -a.yz;
	b.xw = -a.xw;
	b.yw = -a.yw;
	b.zw = -a.zw;

	b.xyz = -a.xyz;
	b.xyw = -a.xyw;
	b.xzw = -a.xzw;
	b.yzw = -a.yzw;

	b.xyzw =  a.xyzw;

	return b;
}

//////////////////////////////////////////////////////

GA4E Conjugation(GA4E a)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R, -S)
{
	GA4E b;
	b.q =  a.q;

	b.x = -a.x;
	b.y = -a.y;
	b.z = -a.z;
	b.w = -a.w;

	b.xy = -a.xy;
	b.xz = -a.xz;
	b.yz = -a.yz;
	b.xw = -a.xw;
	b.yw = -a.yw;
	b.zw = -a.zw;

	b.xyz = -a.xyz;
	b.xyw = -a.xyw;
	b.xzw = -a.xzw;
	b.yzw = -a.yzw;

	b.xyzw =  -a.xyzw;

	return b;
}

//////////////////////////////////////////////////////

GA4E CliffordConjugation(GA4E a)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L,  M, N, P, R,  S) score = 10 N(32736) Clifford Conjugation
{
	GA4E b;
	b.q =  a.q;

	b.x = -a.x;
	b.y = -a.y;
	b.z = -a.z;
	b.w = -a.w;

	b.xy = -a.xy;
	b.xz = -a.xz;
	b.yz = -a.yz;
	b.xw = -a.xw;
	b.yw = -a.yw;
	b.zw = -a.zw;

	b.xyz = a.xyz;
	b.xyw = a.xyw;
	b.xzw = a.xzw;
	b.yzw = a.yzw;

	b.xyzw =  a.xyzw;

	return b;
}

//////////////////////////////////////////////////////

GA4E Dual(GA4E a)   // return w*I_inv 
// Dual(r) = u = (p, o,-n,m,-l, -k,j,-i,-h,g,-f, -e,d,-c,b, a)

{
	GA4E b;
//	GA4E I_inv;
//	I_inv.xyzw = 1;

	b.q =  a.xyzw;

	b.x =  a.yzw;
	b.y = -a.xzw;
	b.z =  a.xyw;
	b.w = -a.xyz;

	b.xy = -a.zw;
	b.xz =  a.yw;
	b.yz = -a.xw;
	b.xw = -a.yz;
	b.yw =  a.xz;
	b.zw = -a.xy;

	b.xyz = -a.w;
	b.xyw =  a.z;
	b.xzw = -a.y;
	b.yzw =  a.x;

	b.xyzw =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA4E DorstDual(GA4E a)
{
	GA4E b, I_inv;

	I_inv.xyzw = 1;

	b = LeftContraction(a,I_inv);
	return b;
}

//////////////////////////////////////////////////////

GA4E DorstUnDual(GA4E a)
{
	GA4E b, I;

	I.xyzw = 1;

	b = LeftContraction(a,I);
	return b;
}

//////////////////////////////////////////////////////

GA4E operator+(const GA4E &u, const GA4E &v)
{
	GA4E w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.w   = u.w   + v.w  ;
	w.xy  = u.xy  + v.xy ;
	w.xz  = u.xz  + v.xz ;
	w.yz  = u.yz  + v.yz ;
	w.xw  = u.xw  + v.xw ;
	w.yw  = u.yw  + v.yw ;
	w.zw  = u.zw  + v.zw ;
	w.xyz = u.xyz + v.xyz;
	w.xyw = u.xyw + v.xyw;
	w.xzw = u.xzw + v.xzw;
	w.yzw = u.yzw + v.yzw;
	w.xyzw = u.xyzw + v.xyzw;

	return w;
}

//////////////////////////////////////////////////////

GA4E operator-(const GA4E &u, const GA4E &v)
{
	GA4E w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.w   = u.w   - v.w  ;
	w.xy  = u.xy  - v.xy ;
	w.xz  = u.xz  - v.xz ;
	w.yz  = u.yz  - v.yz ;
	w.xw  = u.xw  - v.xw ;
	w.yw  = u.yw  - v.yw ;
	w.zw  = u.zw  - v.zw ;
	w.xyz = u.xyz - v.xyz;
	w.xyw = u.xyw - v.xyw;
	w.xzw = u.xzw - v.xzw;
	w.yzw = u.yzw - v.yzw;
	w.xyzw = u.xyzw - v.xyzw;

	return w;
}

//////////////////////////////////////////////////////

int operator==(const GA4E &u, const GA4E &v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&
			(u.w ==v.w )&&

			(u.xy==v.xy)&&
			(u.xz==v.xz)&&
			(u.yz==v.yz)&&
			(u.xw==v.xw)&&
			(u.yw==v.yw)&&
			(u.zw==v.zw)&&

			(u.xyz==v.xyz)&&
			(u.xyw==v.xyw)&&
			(u.xzw==v.xzw)&&
			(u.yzw==v.yzw)&&

			(u.xyzw==v.xyzw);
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const GA4E &u, const GA4E &v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||
			(u.w !=v.w )||

			(u.xy!=v.xy)||
			(u.xz!=v.xz)||
			(u.yz!=v.yz)||
			(u.xw!=v.xw)||
			(u.yw!=v.yw)||
			(u.zw!=v.zw)||

			(u.xyz!=v.xyz)||
			(u.xyw!=v.xyw)||
			(u.xzw!=v.xzw)||
			(u.yzw!=v.yzw)||

			(u.xyzw!=v.xyzw);
	return result;
}

//////////////////////////////////////////////////////

GA4E operator*(const GA4E &a, const GA4E &b) {
	GA4E c;

c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.w   *b.w    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xw  *b.xw   - a.yw  *b.yw   - a.zw  *b.zw   - a.xyz *b.xyz  - a.xyw *b.xyw  - a.xzw *b.xzw  - a.yzw *b.yzw  + a.xyzw*b.xyzw;
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.y   *b.xy   - a.z   *b.xz   - a.w   *b.xw   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  + a.xw  *b.w    - a.yw  *b.xyw  - a.zw  *b.xzw  - a.xyz *b.yz   - a.xyw *b.yw   - a.xzw *b.zw   + a.yzw *b.xyzw - a.xyzw*b.yzw ;
c.y    =  + a.q   *b.y    + a.x   *b.xy   + a.y   *b.q    - a.z   *b.yz   - a.w   *b.yw   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    + a.xw  *b.xyw  + a.yw  *b.w    - a.zw  *b.yzw  + a.xyz *b.xz   + a.xyw *b.xw   - a.xzw *b.xyzw - a.yzw *b.zw   + a.xyzw*b.xzw ;
c.z    =  + a.q   *b.z    + a.x   *b.xz   + a.y   *b.yz   + a.z   *b.q    - a.w   *b.zw   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    + a.xw  *b.xzw  + a.yw  *b.yzw  + a.zw  *b.w    - a.xyz *b.xy   + a.xyw *b.xyzw + a.xzw *b.xw   + a.yzw *b.yw   - a.xyzw*b.xyw ;
c.w    =  + a.q   *b.w    + a.x   *b.xw   + a.y   *b.yw   + a.z   *b.zw   + a.w   *b.q    - a.xy  *b.xyw  - a.xz  *b.xzw  - a.yz  *b.yzw  - a.xw  *b.x    - a.yw  *b.y    - a.zw  *b.z    - a.xyz *b.xyzw - a.xyw *b.xy   - a.xzw *b.xz   - a.yzw *b.yz   + a.xyzw*b.xyz ;
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.z   *b.xyz  + a.w   *b.xyw  + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   - a.xw  *b.yw   + a.yw  *b.xw   - a.zw  *b.xyzw + a.xyz *b.z    + a.xyw *b.w    - a.xzw *b.yzw  + a.yzw *b.xzw  - a.xyzw*b.zw  ;
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.y   *b.xyz  - a.z   *b.x    + a.w   *b.xzw  + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   - a.xw  *b.zw   + a.yw  *b.xyzw + a.zw  *b.xw   - a.xyz *b.y    + a.xyw *b.yzw  + a.xzw *b.w    - a.yzw *b.xyw  + a.xyzw*b.yw  ;
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.y   *b.z    - a.z   *b.y    + a.w   *b.yzw  - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    - a.xw  *b.xyzw - a.yw  *b.zw   + a.zw  *b.yw   + a.xyz *b.x    - a.xyw *b.xzw  + a.xzw *b.xyw  + a.yzw *b.w    - a.xyzw*b.xw  ;
c.xw   =  + a.q   *b.xw   + a.x   *b.w    - a.y   *b.xyw  - a.z   *b.xzw  - a.w   *b.x    + a.xy  *b.yw   + a.xz  *b.zw   - a.yz  *b.xyzw + a.xw  *b.q    - a.yw  *b.xy   - a.zw  *b.xz   - a.xyz *b.yzw  - a.xyw *b.y    - a.xzw *b.z    + a.yzw *b.xyz  - a.xyzw*b.yz  ;
c.yw   =  + a.q   *b.yw   + a.x   *b.xyw  + a.y   *b.w    - a.z   *b.yzw  - a.w   *b.y    - a.xy  *b.xw   + a.xz  *b.xyzw + a.yz  *b.zw   + a.xw  *b.xy   + a.yw  *b.q    - a.zw  *b.yz   + a.xyz *b.xzw  + a.xyw *b.x    - a.xzw *b.xyz  - a.yzw *b.z    + a.xyzw*b.xz  ;
c.zw   =  + a.q   *b.zw   + a.x   *b.xzw  + a.y   *b.yzw  + a.z   *b.w    - a.w   *b.z    - a.xy  *b.xyzw - a.xz  *b.xw   - a.yz  *b.yw   + a.xw  *b.xz   + a.yw  *b.yz   + a.zw  *b.q    - a.xyz *b.xyw  + a.xyw *b.xyz  + a.xzw *b.x    + a.yzw *b.y    - a.xyzw*b.xy  ;
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   - a.w   *b.xyzw + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  + a.xyz *b.q    - a.xyw *b.zw   + a.xzw *b.yw   - a.yzw *b.xw   + a.xyzw*b.w   ;
c.xyw  =  + a.q   *b.xyw  + a.x   *b.yw   - a.y   *b.xw   + a.z   *b.xyzw + a.w   *b.xy   + a.xy  *b.w    - a.xz  *b.yzw  + a.yz  *b.xzw  - a.xw  *b.y    + a.yw  *b.x    - a.zw  *b.xyz  + a.xyz *b.zw   + a.xyw *b.q    - a.xzw *b.yz   + a.yzw *b.xz   - a.xyzw*b.z   ;
c.xzw  =  + a.q   *b.xzw  + a.x   *b.zw   - a.y   *b.xyzw - a.z   *b.xw   + a.w   *b.xz   + a.xy  *b.yzw  + a.xz  *b.w    - a.yz  *b.xyw  - a.xw  *b.z    + a.yw  *b.xyz  + a.zw  *b.x    - a.xyz *b.yw   + a.xyw *b.yz   + a.xzw *b.q    - a.yzw *b.xy   + a.xyzw*b.y   ;
c.yzw  =  + a.q   *b.yzw  + a.x   *b.xyzw + a.y   *b.zw   - a.z   *b.yw   + a.w   *b.yz   - a.xy  *b.xzw  + a.xz  *b.xyw  + a.yz  *b.w    - a.xw  *b.xyz  - a.yw  *b.z    + a.zw  *b.y    + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   + a.yzw *b.q    - a.xyzw*b.x   ;
c.xyzw =  + a.q   *b.xyzw + a.x   *b.yzw  - a.y   *b.xzw  + a.z   *b.xyw  - a.w   *b.xyz  + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   + a.xyz *b.w    - a.xyw *b.z    + a.xzw *b.y    - a.yzw *b.x    + a.xyzw*b.q   ;


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

GA4E Product(const GA4E &a, const GA4E &b)
{

	GA4E c;
/*
 *    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
 x    |  x     q     xy    xz    xw    y     z     xyz   w     xyw   xzw   yz    yw    zw    xyzw  yzw  
 y    |  y    -xy    q     yz    yw   -x    -xyz   z    -xyw   w     yzw  -xz   -xw   -xyzw  zw   -xzw  
 z    |  z    -xz   -yz    q     zw    xyz  -x    -y    -xzw  -yzw   w     xy    xyzw -xw   -yw    xyw  
 w    |  w    -xw   -yw   -zw    q     xyw   xzw   yzw  -x    -y    -z    -xyzw  xy    xz    yz   -xyz  
 xy   |  xy   -y     x     xyz   xyw  -q    -yz    xz   -yw    xw    xyzw -z    -w    -yzw   xzw  -zw   
 xz   |  xz   -z    -xyz   x     xzw   yz   -q    -xy   -zw   -xyzw  xw    y     yzw  -w    -xyw   yw   
 yz   |  yz    xyz  -z     y     yzw  -xz    xy   -q     xyzw -zw    yw   -x    -xzw   xyw  -w    -xw   
 xw   |  xw   -w    -xyw  -xzw   x     yw    zw    xyzw -q    -xy   -xz   -yzw   y     z     xyz  -yz   
 yw   |  yw    xyw  -w    -yzw   y    -xw   -xyzw  zw    xy   -q    -yz    xzw  -x    -xyz   z     xz   
 zw   |  zw    xzw   yzw  -w     z     xyzw -xw   -yw    xz    yz   -q    -xyw   xyz  -x    -y    -xy   
 xyz  |  xyz   yz   -xz    xy    xyzw -z     y    -x     yzw  -xzw   xyw  -q    -zw    yw   -xw   -w    
 xyw  |  xyw   yw   -xw   -xyzw  xy   -w    -yzw   xzw   y    -x    -xyz   zw   -q    -yz    xz    z    
 xzw  |  xzw   zw    xyzw -xw    xz    yzw  -w    -xyw   z     xyz  -x    -yw    yz   -q    -xy   -y    
 yzw  |  yzw  -xyzw  zw   -yw    yz   -xzw   xyw  -w    -xyz   z    -y     xw   -xz    xy   -q     x    
 xyzw |  xyzw -yzw   xzw  -xyw   xyz  -zw    yw   -xw   -yz    xz   -xy    w    -z     y    -x     q    

*/

c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.w   *b.w    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xw  *b.xw   - a.yw  *b.yw   - a.zw  *b.zw   - a.xyz *b.xyz  - a.xyw *b.xyw  - a.xzw *b.xzw  - a.yzw *b.yzw  + a.xyzw*b.xyzw;
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.y   *b.xy   - a.z   *b.xz   - a.w   *b.xw   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  + a.xw  *b.w    - a.yw  *b.xyw  - a.zw  *b.xzw  - a.xyz *b.yz   - a.xyw *b.yw   - a.xzw *b.zw   + a.yzw *b.xyzw - a.xyzw*b.yzw ;
c.y    =  + a.q   *b.y    + a.x   *b.xy   + a.y   *b.q    - a.z   *b.yz   - a.w   *b.yw   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    + a.xw  *b.xyw  + a.yw  *b.w    - a.zw  *b.yzw  + a.xyz *b.xz   + a.xyw *b.xw   - a.xzw *b.xyzw - a.yzw *b.zw   + a.xyzw*b.xzw ;
c.z    =  + a.q   *b.z    + a.x   *b.xz   + a.y   *b.yz   + a.z   *b.q    - a.w   *b.zw   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    + a.xw  *b.xzw  + a.yw  *b.yzw  + a.zw  *b.w    - a.xyz *b.xy   + a.xyw *b.xyzw + a.xzw *b.xw   + a.yzw *b.yw   - a.xyzw*b.xyw ;
c.w    =  + a.q   *b.w    + a.x   *b.xw   + a.y   *b.yw   + a.z   *b.zw   + a.w   *b.q    - a.xy  *b.xyw  - a.xz  *b.xzw  - a.yz  *b.yzw  - a.xw  *b.x    - a.yw  *b.y    - a.zw  *b.z    - a.xyz *b.xyzw - a.xyw *b.xy   - a.xzw *b.xz   - a.yzw *b.yz   + a.xyzw*b.xyz ;
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.z   *b.xyz  + a.w   *b.xyw  + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   - a.xw  *b.yw   + a.yw  *b.xw   - a.zw  *b.xyzw + a.xyz *b.z    + a.xyw *b.w    - a.xzw *b.yzw  + a.yzw *b.xzw  - a.xyzw*b.zw  ;
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.y   *b.xyz  - a.z   *b.x    + a.w   *b.xzw  + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   - a.xw  *b.zw   + a.yw  *b.xyzw + a.zw  *b.xw   - a.xyz *b.y    + a.xyw *b.yzw  + a.xzw *b.w    - a.yzw *b.xyw  + a.xyzw*b.yw  ;
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.y   *b.z    - a.z   *b.y    + a.w   *b.yzw  - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    - a.xw  *b.xyzw - a.yw  *b.zw   + a.zw  *b.yw   + a.xyz *b.x    - a.xyw *b.xzw  + a.xzw *b.xyw  + a.yzw *b.w    - a.xyzw*b.xw  ;
c.xw   =  + a.q   *b.xw   + a.x   *b.w    - a.y   *b.xyw  - a.z   *b.xzw  - a.w   *b.x    + a.xy  *b.yw   + a.xz  *b.zw   - a.yz  *b.xyzw + a.xw  *b.q    - a.yw  *b.xy   - a.zw  *b.xz   - a.xyz *b.yzw  - a.xyw *b.y    - a.xzw *b.z    + a.yzw *b.xyz  - a.xyzw*b.yz  ;
c.yw   =  + a.q   *b.yw   + a.x   *b.xyw  + a.y   *b.w    - a.z   *b.yzw  - a.w   *b.y    - a.xy  *b.xw   + a.xz  *b.xyzw + a.yz  *b.zw   + a.xw  *b.xy   + a.yw  *b.q    - a.zw  *b.yz   + a.xyz *b.xzw  + a.xyw *b.x    - a.xzw *b.xyz  - a.yzw *b.z    + a.xyzw*b.xz  ;
c.zw   =  + a.q   *b.zw   + a.x   *b.xzw  + a.y   *b.yzw  + a.z   *b.w    - a.w   *b.z    - a.xy  *b.xyzw - a.xz  *b.xw   - a.yz  *b.yw   + a.xw  *b.xz   + a.yw  *b.yz   + a.zw  *b.q    - a.xyz *b.xyw  + a.xyw *b.xyz  + a.xzw *b.x    + a.yzw *b.y    - a.xyzw*b.xy  ;
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   - a.w   *b.xyzw + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  + a.xyz *b.q    - a.xyw *b.zw   + a.xzw *b.yw   - a.yzw *b.xw   + a.xyzw*b.w   ;
c.xyw  =  + a.q   *b.xyw  + a.x   *b.yw   - a.y   *b.xw   + a.z   *b.xyzw + a.w   *b.xy   + a.xy  *b.w    - a.xz  *b.yzw  + a.yz  *b.xzw  - a.xw  *b.y    + a.yw  *b.x    - a.zw  *b.xyz  + a.xyz *b.zw   + a.xyw *b.q    - a.xzw *b.yz   + a.yzw *b.xz   - a.xyzw*b.z   ;
c.xzw  =  + a.q   *b.xzw  + a.x   *b.zw   - a.y   *b.xyzw - a.z   *b.xw   + a.w   *b.xz   + a.xy  *b.yzw  + a.xz  *b.w    - a.yz  *b.xyw  - a.xw  *b.z    + a.yw  *b.xyz  + a.zw  *b.x    - a.xyz *b.yw   + a.xyw *b.yz   + a.xzw *b.q    - a.yzw *b.xy   + a.xyzw*b.y   ;
c.yzw  =  + a.q   *b.yzw  + a.x   *b.xyzw + a.y   *b.zw   - a.z   *b.yw   + a.w   *b.yz   - a.xy  *b.xzw  + a.xz  *b.xyw  + a.yz  *b.w    - a.xw  *b.xyz  - a.yw  *b.z    + a.zw  *b.y    + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   + a.yzw *b.q    - a.xyzw*b.x   ;
c.xyzw =  + a.q   *b.xyzw + a.x   *b.yzw  - a.y   *b.xzw  + a.z   *b.xyw  - a.w   *b.xyz  + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   + a.xyz *b.w    - a.xyw *b.z    + a.xzw *b.y    - a.yzw *b.x    + a.xyzw*b.q   ;


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

GA4E operator*(const GA4E &v, const ex &u) {
	GA4E w;

	w.q    = expand(u*v.q   ) ;
	w.w    = expand(u*v.w   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;

	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xw   = expand(u*v.xw  ) ;
	w.yw   = expand(u*v.yw  ) ;
	w.zw   = expand(u*v.zw  ) ;

	w.xyz  = expand(u*v.xyz ) ;
	w.xyw  = expand(u*v.xyw ) ;
	w.xzw  = expand(u*v.xzw ) ;
	w.yzw  = expand(u*v.yzw ) ;
	w.xyzw = expand(u*v.xyzw) ;

	return w;
}

//////////////////////////////////////////////////////

GA4E operator*(const ex &u, const GA4E &v) {
	GA4E w;

	w.q    = expand(u*v.q   ) ;
	w.w    = expand(u*v.w   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;

	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xw   = expand(u*v.xw  ) ;
	w.yw   = expand(u*v.yw  ) ;
	w.zw   = expand(u*v.zw  ) ;

	w.xyz  = expand(u*v.xyz ) ;
	w.xyw  = expand(u*v.xyw ) ;
	w.xzw  = expand(u*v.xzw ) ;
	w.yzw  = expand(u*v.yzw ) ;
	w.xyzw = expand(u*v.xyzw) ;

	return w;
}


//////////////////////////////////////////////////////

GA4E operator/(const GA4E &u, const int i)
{
	GA4E w;
	w.q = u.q/i;

	w.x = u.x/i;
	w.y = u.y/i;
	w.z = u.z/i;
	w.w = u.w/i;

	w.xy = u.xy/i;
	w.xz = u.xz/i;
	w.yz = u.yz/i;
	w.xw = u.xw/i;
	w.yw = u.yw/i;
	w.zw = u.zw/i;

	w.xyz = u.xyz/i;
	w.xyw = u.xyw/i;
	w.xzw = u.xzw/i;
	w.yzw = u.yzw/i;

	w.xyzw = u.xyzw/i;

	return w;
}

//////////////////////////////////////////////////////

GA4E operator^(const GA4E &u, const GA4E &v) {
	GA4E w;

/*

Wedge product, four dimensional spacetime

 ^    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
 x    |  x     0     xy    xz    xw    0     0     xyz   0     xyw   xzw   0     0     0     xyzw  0    
 y    |  y    -xy    0     yz    yw    0    -xyz   0    -xyw   0     yzw   0     0    -xyzw  0     0    
 z    |  z    -xz   -yz    0     zw    xyz   0     0    -xzw  -yzw   0     0     xyzw  0     0     0    
 w    |  w    -xw   -yw   -zw    0     xyw   xzw   yzw   0     0     0    -xyzw  0     0     0     0    
 xy   |  xy    0     0     xyz   xyw   0     0     0     0     0     xyzw  0     0     0     0     0    
 xz   |  xz    0    -xyz   0     xzw   0     0     0     0    -xyzw  0     0     0     0     0     0    
 yz   |  yz    xyz   0     0     yzw   0     0     0     xyzw  0     0     0     0     0     0     0    
 xw   |  xw    0    -xyw  -xzw   0     0     0     xyzw  0     0     0     0     0     0     0     0    
 yw   |  yw    xyw   0    -yzw   0     0    -xyzw  0     0     0     0     0     0     0     0     0    
 zw   |  zw    xzw   yzw   0     0     xyzw  0     0     0     0     0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     xyzw  0     0     0     0     0     0     0     0     0     0     0    
 xyw  |  xyw   0     0    -xyzw  0     0     0     0     0     0     0     0     0     0     0     0    
 xzw  |  xzw   0     xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0    
 yzw  |  yzw  -xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xyzw |  xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

w.q    =  + u.q   *v.q   ;
w.x    =  + u.q   *v.x    + u.x   *v.q   ;
w.y    =  + u.q   *v.y    + u.y   *v.q   ;
w.z    =  + u.q   *v.z    + u.z   *v.q   ;
w.w    =  + u.q   *v.w    + u.w   *v.q   ;
w.xy   =  + u.q   *v.xy   + u.x   *v.y    - u.y   *v.x    + u.xy  *v.q   ;
w.xz   =  + u.q   *v.xz   + u.x   *v.z    - u.z   *v.x    + u.xz  *v.q   ;
w.yz   =  + u.q   *v.yz   + u.y   *v.z    - u.z   *v.y    + u.yz  *v.q   ;
w.xw   =  + u.q   *v.xw   + u.x   *v.w    - u.w   *v.x    + u.xw  *v.q   ;
w.yw   =  + u.q   *v.yw   + u.y   *v.w    - u.w   *v.y    + u.yw  *v.q   ;
w.zw   =  + u.q   *v.zw   + u.z   *v.w    - u.w   *v.z    + u.zw  *v.q   ;
w.xyz  =  + u.q   *v.xyz  + u.x   *v.yz   - u.y   *v.xz   + u.z   *v.xy   + u.xy  *v.z    - u.xz  *v.y    + u.yz  *v.x    + u.xyz *v.q   ;
w.xyw  =  + u.q   *v.xyw  + u.x   *v.yw   - u.y   *v.xw   + u.w   *v.xy   + u.xy  *v.w    - u.xw  *v.y    + u.yw  *v.x    + u.xyw *v.q   ;
w.xzw  =  + u.q   *v.xzw  + u.x   *v.zw   - u.z   *v.xw   + u.w   *v.xz   + u.xz  *v.w    - u.xw  *v.z    + u.zw  *v.x    + u.xzw *v.q   ;
w.yzw  =  + u.q   *v.yzw  + u.y   *v.zw   - u.z   *v.yw   + u.w   *v.yz   + u.yz  *v.w    - u.yw  *v.z    + u.zw  *v.y    + u.yzw *v.q   ;
w.xyzw =  + u.q   *v.xyzw + u.x   *v.yzw  - u.y   *v.xzw  + u.z   *v.xyw  - u.w   *v.xyz  + u.xy  *v.zw   - u.xz  *v.yw   + u.yz  *v.xw   + u.xw  *v.yz   - u.yw  *v.xz   + u.zw  *v.xy   + u.xyz *v.w    - u.xyw *v.z    + u.xzw *v.y    - u.yzw *v.x    + u.xyzw*v.q   ;

	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.w = expand(w.w);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xw = expand(w.xw);
	w.yw = expand(w.yw);
	w.zw = expand(w.zw);
	w.xyz = expand(w.xyz);
	w.xyw = expand(w.xyw);
	w.xzw = expand(w.xzw);
	w.yzw = expand(w.yzw);
	w.xyzw = expand(w.xyzw);

	return w;
}

//////////////////////////////////////////////////////

GA4E Wedge(const GA4E &u, const GA4E &v) {
	GA4E w;

/*

Wedge product, four dimensional spacetime

 ^    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
 x    |  x     0     xy    xz    xw    0     0     xyz   0     xyw   xzw   0     0     0     xyzw  0    
 y    |  y    -xy    0     yz    yw    0    -xyz   0    -xyw   0     yzw   0     0    -xyzw  0     0    
 z    |  z    -xz   -yz    0     zw    xyz   0     0    -xzw  -yzw   0     0     xyzw  0     0     0    
 w    |  w    -xw   -yw   -zw    0     xyw   xzw   yzw   0     0     0    -xyzw  0     0     0     0    
 xy   |  xy    0     0     xyz   xyw   0     0     0     0     0     xyzw  0     0     0     0     0    
 xz   |  xz    0    -xyz   0     xzw   0     0     0     0    -xyzw  0     0     0     0     0     0    
 yz   |  yz    xyz   0     0     yzw   0     0     0     xyzw  0     0     0     0     0     0     0    
 xw   |  xw    0    -xyw  -xzw   0     0     0     xyzw  0     0     0     0     0     0     0     0    
 yw   |  yw    xyw   0    -yzw   0     0    -xyzw  0     0     0     0     0     0     0     0     0    
 zw   |  zw    xzw   yzw   0     0     xyzw  0     0     0     0     0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     xyzw  0     0     0     0     0     0     0     0     0     0     0    
 xyw  |  xyw   0     0    -xyzw  0     0     0     0     0     0     0     0     0     0     0     0    
 xzw  |  xzw   0     xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0    
 yzw  |  yzw  -xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xyzw |  xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

w.q    =  + u.q   *v.q   ;
w.x    =  + u.q   *v.x    + u.x   *v.q   ;
w.y    =  + u.q   *v.y    + u.y   *v.q   ;
w.z    =  + u.q   *v.z    + u.z   *v.q   ;
w.w    =  + u.q   *v.w    + u.w   *v.q   ;
w.xy   =  + u.q   *v.xy   + u.x   *v.y    - u.y   *v.x    + u.xy  *v.q   ;
w.xz   =  + u.q   *v.xz   + u.x   *v.z    - u.z   *v.x    + u.xz  *v.q   ;
w.yz   =  + u.q   *v.yz   + u.y   *v.z    - u.z   *v.y    + u.yz  *v.q   ;
w.xw   =  + u.q   *v.xw   + u.x   *v.w    - u.w   *v.x    + u.xw  *v.q   ;
w.yw   =  + u.q   *v.yw   + u.y   *v.w    - u.w   *v.y    + u.yw  *v.q   ;
w.zw   =  + u.q   *v.zw   + u.z   *v.w    - u.w   *v.z    + u.zw  *v.q   ;
w.xyz  =  + u.q   *v.xyz  + u.x   *v.yz   - u.y   *v.xz   + u.z   *v.xy   + u.xy  *v.z    - u.xz  *v.y    + u.yz  *v.x    + u.xyz *v.q   ;
w.xyw  =  + u.q   *v.xyw  + u.x   *v.yw   - u.y   *v.xw   + u.w   *v.xy   + u.xy  *v.w    - u.xw  *v.y    + u.yw  *v.x    + u.xyw *v.q   ;
w.xzw  =  + u.q   *v.xzw  + u.x   *v.zw   - u.z   *v.xw   + u.w   *v.xz   + u.xz  *v.w    - u.xw  *v.z    + u.zw  *v.x    + u.xzw *v.q   ;
w.yzw  =  + u.q   *v.yzw  + u.y   *v.zw   - u.z   *v.yw   + u.w   *v.yz   + u.yz  *v.w    - u.yw  *v.z    + u.zw  *v.y    + u.yzw *v.q   ;
w.xyzw =  + u.q   *v.xyzw + u.x   *v.yzw  - u.y   *v.xzw  + u.z   *v.xyw  - u.w   *v.xyz  + u.xy  *v.zw   - u.xz  *v.yw   + u.yz  *v.xw   + u.xw  *v.yz   - u.yw  *v.xz   + u.zw  *v.xy   + u.xyz *v.w    - u.xyw *v.z    + u.xzw *v.y    - u.yzw *v.x    + u.xyzw*v.q   ;

	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.w = expand(w.w);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xw = expand(w.xw);
	w.yw = expand(w.yw);
	w.zw = expand(w.zw);
	w.xyz = expand(w.xyz);
	w.xyw = expand(w.xyw);
	w.xzw = expand(w.xzw);
	w.yzw = expand(w.yzw);
	w.xyzw = expand(w.xyzw);

	return w;
}

//////////////////////////////////////////////////////

GA4E AntiWedge(const GA4E a, const GA4E b)
{
	GA4E c;

/*
Lengyel's AntiWedge product

 V    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0    -q     0     y    
 z    |  0     0     0     0     0     0     0     0     0     0     0     0     q     0     0     z    
 w    |  0     0     0     0     0     0     0     0     0     0     0    -q     0     0     0     w    
 xy   |  0     0     0     0     0     0     0     0     0     0     q     0     0     x     y     xy   
 xz   |  0     0     0     0     0     0     0     0     0    -q     0     0    -x     0     z     xz   
 yz   |  0     0     0     0     0     0     0     0     q     0     0     0    -y    -z     0     yz   
 xw   |  0     0     0     0     0     0     0     q     0     0     0     x     0     0     w     xw   
 yw   |  0     0     0     0     0     0    -q     0     0     0     0     y     0    -w     0     yw   
 zw   |  0     0     0     0     0     q     0     0     0     0     0     z     w     0     0     zw   
 xyz  |  0     0     0     0     q     0     0     0     x     y     z     0     xy    xz    yz    xyz  
 xyw  |  0     0     0    -q     0     0    -x    -y     0     0     w    -xy    0     xw    yw    xyw  
 xzw  |  0     0     q     0     0     x     0    -z     0    -w     0    -xz   -xw    0     zw    xzw  
 yzw  |  0    -q     0     0     0     y     z     0     w     0     0    -yz   -yw   -zw    0     yzw  
 xyzw |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 

*/

c.q    =  + a.q   *b.xyzw + a.x   *b.yzw  - a.y   *b.xzw  + a.z   *b.xyw  - a.w   *b.xyz  + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   + a.xyz *b.w    - a.xyw *b.z    + a.xzw *b.y    - a.yzw *b.x    + a.xyzw*b.q    ; 
c.x    =  + a.x   *b.xyzw + a.xy  *b.xzw  - a.xz  *b.xyw  + a.xw  *b.xyz  + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   + a.xyzw*b.x    ; 
c.y    =  + a.y   *b.xyzw + a.xy  *b.yzw  - a.yz  *b.xyw  + a.yw  *b.xyz  + a.xyz *b.yw   - a.xyw *b.yz   + a.yzw *b.xy   + a.xyzw*b.y    ; 
c.z    =  + a.z   *b.xyzw + a.xz  *b.yzw  - a.yz  *b.xzw  + a.zw  *b.xyz  + a.xyz *b.zw   - a.xzw *b.yz   + a.yzw *b.xz   + a.xyzw*b.z    ; 
c.w    =  + a.w   *b.xyzw + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  + a.xyw *b.zw   - a.xzw *b.yw   + a.yzw *b.xw   + a.xyzw*b.w    ; 
c.xy   =  + a.xy  *b.xyzw + a.xyz *b.xyw  - a.xyw *b.xyz  + a.xyzw*b.xy   ; 
c.xz   =  + a.xz  *b.xyzw + a.xyz *b.xzw  - a.xzw *b.xyz  + a.xyzw*b.xz   ; 
c.yz   =  + a.yz  *b.xyzw + a.xyz *b.yzw  - a.yzw *b.xyz  + a.xyzw*b.yz   ; 
c.xw   =  + a.xw  *b.xyzw + a.xyw *b.xzw  - a.xzw *b.xyw  + a.xyzw*b.xw   ; 
c.yw   =  + a.yw  *b.xyzw + a.xyw *b.yzw  - a.yzw *b.xyw  + a.xyzw*b.yw   ; 
c.zw   =  + a.zw  *b.xyzw + a.xzw *b.yzw  - a.yzw *b.xzw  + a.xyzw*b.zw   ; 
c.xyz  =  + a.xyz *b.xyzw + a.xyzw*b.xyz  ; 
c.xyw  =  + a.xyw *b.xyzw + a.xyzw*b.xyw  ; 
c.xzw  =  + a.xzw *b.xyzw + a.xyzw*b.xzw  ; 
c.yzw  =  + a.yzw *b.xyzw + a.xyzw*b.yzw  ; 
c.xyzw =  + a.xyzw*b.xyzw ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////


GA4E Regressive(GA4E a, GA4E b) 
{
	GA4E c;

//	GA4E I, I_inv;
//	I.xyzw = 1;
//	I_inv.xyzw = 1;
//	c = ((a*I_inv)^(b*I_inv))*I;

/*

Hestenes' Regressive product differs from Lengyel's AntiWedge product in sign details

 V    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0    -q     x    
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     0     0     0     0     0     0     0     0    -q     0     0     z    
 w    |  0     0     0     0     0     0     0     0     0     0     0     q     0     0     0     w    
 xy   |  0     0     0     0     0     0     0     0     0     0     q     0     0     x     y     xy   
 xz   |  0     0     0     0     0     0     0     0     0    -q     0     0    -x     0     z     xz   
 yz   |  0     0     0     0     0     0     0     0     q     0     0     0    -y    -z     0     yz   
 xw   |  0     0     0     0     0     0     0     q     0     0     0     x     0     0     w     xw   
 yw   |  0     0     0     0     0     0    -q     0     0     0     0     y     0    -w     0     yw   
 zw   |  0     0     0     0     0     q     0     0     0     0     0     z     w     0     0     zw   
 xyz  |  0     0     0     0    -q     0     0     0     x     y     z     0    -xy   -xz   -yz    xyz  
 xyw  |  0     0     0     q     0     0    -x    -y     0     0     w     xy    0    -xw   -yw    xyw  
 xzw  |  0     0    -q     0     0     x     0    -z     0    -w     0     xz    xw    0    -zw    xzw  
 yzw  |  0     q     0     0     0     y     z     0     w     0     0     yz    yw    zw    0     yzw  
 xyzw |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 

*/

c.q    =  + a.q   *b.xyzw - a.x   *b.yzw  + a.y   *b.xzw  - a.z   *b.xyw  + a.w   *b.xyz  + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   - a.xyz *b.w    + a.xyw *b.z    - a.xzw *b.y    + a.yzw *b.x    + a.xyzw*b.q   ;
c.x    =  + a.x   *b.xyzw + a.xy  *b.xzw  - a.xz  *b.xyw  + a.xw  *b.xyz  + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   + a.xyzw*b.x   ;
c.y    =  + a.y   *b.xyzw + a.xy  *b.yzw  - a.yz  *b.xyw  + a.yw  *b.xyz  + a.xyz *b.yw   - a.xyw *b.yz   + a.yzw *b.xy   + a.xyzw*b.y   ;
c.z    =  + a.z   *b.xyzw + a.xz  *b.yzw  - a.yz  *b.xzw  + a.zw  *b.xyz  + a.xyz *b.zw   - a.xzw *b.yz   + a.yzw *b.xz   + a.xyzw*b.z   ;
c.w    =  + a.w   *b.xyzw + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  + a.xyw *b.zw   - a.xzw *b.yw   + a.yzw *b.xw   + a.xyzw*b.w   ;
c.xy   =  + a.xy  *b.xyzw - a.xyz *b.xyw  + a.xyw *b.xyz  + a.xyzw*b.xy  ;
c.xz   =  + a.xz  *b.xyzw - a.xyz *b.xzw  + a.xzw *b.xyz  + a.xyzw*b.xz  ;
c.yz   =  + a.yz  *b.xyzw - a.xyz *b.yzw  + a.yzw *b.xyz  + a.xyzw*b.yz  ;
c.xw   =  + a.xw  *b.xyzw - a.xyw *b.xzw  + a.xzw *b.xyw  + a.xyzw*b.xw  ;
c.yw   =  + a.yw  *b.xyzw - a.xyw *b.yzw  + a.yzw *b.xyw  + a.xyzw*b.yw  ;
c.zw   =  + a.zw  *b.xyzw - a.xzw *b.yzw  + a.yzw *b.xzw  + a.xyzw*b.zw  ;
c.xyz  =  + a.xyz *b.xyzw + a.xyzw*b.xyz ;
c.xyw  =  + a.xyw *b.xyzw + a.xyzw*b.xyw ;
c.xzw  =  + a.xzw *b.xyzw + a.xyzw*b.xzw ;
c.yzw  =  + a.yzw *b.xyzw + a.xyzw*b.yzw ;
c.xyzw =  + a.xyzw*b.xyzw ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

GA4E RegressiveViaFormula(GA4E a, GA4E b)
{
	GA4E c;

	GA4E I,I_inv,d,e,f;
	
	I.xyzw = 1;	I_inv.xyzw = 1;

	d = a*I_inv^b*I_inv;
	c = d*I;

	return c;
}

//////////////////////////////////////////////////////

GA4E LowerRightViaFormula(GA4E a, GA4E b)  // Like a mirror of Wedge along rising diagonal
{

/*

LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     xyzw 
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     xyzw  yzw  
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0    -xyzw  0    -xzw  
 z    |  0     0     0     0     0     0     0     0     0     0     0     0     xyzw  0     0     xyw  
 w    |  0     0     0     0     0     0     0     0     0     0     0    -xyzw  0     0     0     xyz  
 xy   |  0     0     0     0     0     0     0     0     0     0     xyzw  0     0    -yzw   xzw  -zw   
 xz   |  0     0     0     0     0     0     0     0     0    -xyzw  0     0     yzw   0    -xyw   yw   
 yz   |  0     0     0     0     0     0     0     0     xyzw  0     0     0    -xzw   xyw   0    -xw   
 xw   |  0     0     0     0     0     0     0     xyzw  0     0     0    -yzw   0     0    -xyz   yz   
 yw   |  0     0     0     0     0     0    -xyzw  0     0     0     0     xzw   0     xyz   0    -xz   
 zw   |  0     0     0     0     0     xyzw  0     0     0     0     0    -xyw  -xyz   0     0     xy   
 xyz  |  0     0     0     0     xyzw  0     0     0     yzw  -xzw   xyw   0    -zw    yw   -xw   -w    
 xyw  |  0     0     0    -xyzw  0     0    -yzw   xzw   0     0     xyz   zw    0     yz   -xz   -z    
 xzw  |  0     0     xyzw  0     0     yzw   0    -xyw   0    -xyz   0    -yw   -yz    0     xy    y    
 yzw  |  0    -xyzw  0     0     0    -xzw   xyw   0     xyz   0     0     xw    xz   -xy    0    -x    
 xyzw |  xyzw -yzw   xzw  -xyw  -xyz  -zw    yw   -xw    yz   -xz    xy    w     z    -y     x    -q    

*/

// -Wedge(Blade[i]*xyzt,xyzt*Blade[j])

	GA4E c;

	GA4E I,I_inv,d,e,f;
	
	I.xyzw = 1;	I_inv.xyzw = 1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

GA4E Expander(const GA4E &a, const GA4E &b)
{
/*

Terms with increased rank
 >    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    xw    0     0     xyz   0     xyw   xzw   0     0     0     xyzw  0    
 y    |  0    -xy    0     yz    yw    0    -xyz   0    -xyw   0     yzw   0     0    -xyzw  0     0    
 z    |  0    -xz   -yz    0     zw    xyz   0     0    -xzw  -yzw   0     0     xyzw  0     0     0    
 w    |  0    -xw   -yw   -zw    0     xyw   xzw   yzw   0     0     0    -xyzw  0     0     0     0    
 xy   |  0     0     0     xyz   xyw   0     0     0     0     0     xyzw  0     0     0     0     0    
 xz   |  0     0    -xyz   0     xzw   0     0     0     0    -xyzw  0     0     0     0     0     0    
 yz   |  0     xyz   0     0     yzw   0     0     0     xyzw  0     0     0     0     0     0     0    
 xw   |  0     0    -xyw  -xzw   0     0     0     xyzw  0     0     0     0     0     0     0     0    
 yw   |  0     xyw   0    -yzw   0     0    -xyzw  0     0     0     0     0     0     0     0     0    
 zw   |  0     xzw   yzw   0     0     xyzw  0     0     0     0     0     0     0     0     0     0    
 xyz  |  0     0     0     0     xyzw  0     0     0     0     0     0     0     0     0     0     0    
 xyw  |  0     0     0    -xyzw  0     0     0     0     0     0     0     0     0     0     0     0    
 xzw  |  0     0     xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0    
 yzw  |  0    -xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xyzw |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/
	GA4E c;

c.q    =  0 ; 
c.x    =  0 ; 
c.y    =  0 ; 
c.z    =  0 ; 
c.w    =  0 ; 
c.xy   =  + a.x   *b.y    - a.y   *b.x    ; 
c.xz   =  + a.x   *b.z    - a.z   *b.x    ; 
c.yz   =  + a.y   *b.z    - a.z   *b.y    ; 
c.xw   =  + a.x   *b.w    - a.w   *b.x    ; 
c.yw   =  + a.y   *b.w    - a.w   *b.y    ; 
c.zw   =  + a.z   *b.w    - a.w   *b.z    ; 
c.xyz  =  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    ; 
c.xyw  =  + a.x   *b.yw   - a.y   *b.xw   + a.w   *b.xy   + a.xy  *b.w    - a.xw  *b.y    + a.yw  *b.x    ; 
c.xzw  =  + a.x   *b.zw   - a.z   *b.xw   + a.w   *b.xz   + a.xz  *b.w    - a.xw  *b.z    + a.zw  *b.x    ; 
c.yzw  =  + a.y   *b.zw   - a.z   *b.yw   + a.w   *b.yz   + a.yz  *b.w    - a.yw  *b.z    + a.zw  *b.y    ; 
c.xyzw =  + a.x   *b.yzw  - a.y   *b.xzw  + a.z   *b.xyw  - a.w   *b.xyz  + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   + a.xyz *b.w    - a.xyw *b.z    + a.xzw *b.y    - a.yzw *b.x    ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);


	return c;
}

//////////////////////////////////////////////////////

GA4E Conserver(const GA4E &a, const GA4E &b)
{
/*

Terms with preserved rank
 =    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
 x    |  x     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 w    |  w     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0     0    -yz    xz   -yw    xw    0     0     0    -yzw   xzw   0    
 xz   |  xz    0     0     0     0     yz    0    -xy   -zw    0     xw    0     yzw   0    -xyw   0    
 yz   |  yz    0     0     0     0    -xz    xy    0     0    -zw    yw    0    -xzw   xyw   0     0    
 xw   |  xw    0     0     0     0     yw    zw    0     0    -xy   -xz   -yzw   0     0     xyz   0    
 yw   |  yw    0     0     0     0    -xw    0     zw    xy    0    -yz    xzw   0    -xyz   0     0    
 zw   |  zw    0     0     0     0     0    -xw   -yw    xz    yz    0    -xyw   xyz   0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0     yzw  -xzw   xyw   0     0     0     0     0    
 xyw  |  xyw   0     0     0     0     0    -yzw   xzw   0     0    -xyz   0     0     0     0     0    
 xzw  |  xzw   0     0     0     0     yzw   0    -xyw   0     xyz   0     0     0     0     0     0    
 yzw  |  yzw   0     0     0     0    -xzw   xyw   0    -xyz   0     0     0     0     0     0     0    
 xyzw |  xyzw  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

	GA4E c;

c.q    =  + a.q   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   - a.xw  *b.yw   + a.yw  *b.xw   ; 
c.xz   =  + a.q   *b.xz   + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   - a.xw  *b.zw   + a.zw  *b.xw   ; 
c.yz   =  + a.q   *b.yz   - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    - a.yw  *b.zw   + a.zw  *b.yw   ; 
c.xw   =  + a.q   *b.xw   + a.xy  *b.yw   + a.xz  *b.zw   + a.xw  *b.q    - a.yw  *b.xy   - a.zw  *b.xz   ; 
c.yw   =  + a.q   *b.yw   - a.xy  *b.xw   + a.yz  *b.zw   + a.xw  *b.xy   + a.yw  *b.q    - a.zw  *b.yz   ; 
c.zw   =  + a.q   *b.zw   - a.xz  *b.xw   - a.yz  *b.yw   + a.xw  *b.xz   + a.yw  *b.yz   + a.zw  *b.q    ; 
c.xyz  =  + a.q   *b.xyz  + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  + a.xyz *b.q    - a.xyw *b.zw   + a.xzw *b.yw   - a.yzw *b.xw   ; 
c.xyw  =  + a.q   *b.xyw  - a.xz  *b.yzw  + a.yz  *b.xzw  - a.zw  *b.xyz  + a.xyz *b.zw   + a.xyw *b.q    - a.xzw *b.yz   + a.yzw *b.xz   ; 
c.xzw  =  + a.q   *b.xzw  + a.xy  *b.yzw  - a.yz  *b.xyw  + a.yw  *b.xyz  - a.xyz *b.yw   + a.xyw *b.yz   + a.xzw *b.q    - a.yzw *b.xy   ; 
c.yzw  =  + a.q   *b.yzw  - a.xy  *b.xzw  + a.xz  *b.xyw  - a.xw  *b.xyz  + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   + a.yzw *b.q    ; 
c.xyzw =  + a.q   *b.xyzw + a.xyzw*b.q    ; 


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);



	return c;
}

//////////////////////////////////////////////////////

GA4E Shrinker(const GA4E &a, const GA4E &b)
{

	GA4E c;

/*

Terms with reduced rank
 <    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     0     y     z     0     w     0     0     yz    yw    zw    0     yzw  
 y    |  0     0     q     0     0    -x     0     z     0     w     0    -xz   -xw    0     zw   -xzw  
 z    |  0     0     0     q     0     0    -x    -y     0     0     w     xy    0    -xw   -yw    xyw  
 w    |  0     0     0     0     q     0     0     0    -x    -y    -z     0     xy    xz    yz   -xyz  
 xy   |  0    -y     x     0     0    -q     0     0     0     0     0    -z    -w     0     0    -zw   
 xz   |  0    -z     0     x     0     0    -q     0     0     0     0     y     0    -w     0     yw   
 yz   |  0     0    -z     y     0     0     0    -q     0     0     0    -x     0     0    -w    -xw   
 xw   |  0    -w     0     0     x     0     0     0    -q     0     0     0     y     z     0    -yz   
 yw   |  0     0    -w     0     y     0     0     0     0    -q     0     0    -x     0     z     xz   
 zw   |  0     0     0    -w     z     0     0     0     0     0    -q     0     0    -x    -y    -xy   
 xyz  |  0     yz   -xz    xy    0    -z     y    -x     0     0     0    -q    -zw    yw   -xw   -w    
 xyw  |  0     yw   -xw    0     xy   -w     0     0     y    -x     0     zw   -q    -yz    xz    z    
 xzw  |  0     zw    0    -xw    xz    0    -w     0     z     0    -x    -yw    yz   -q    -xy   -y    
 yzw  |  0     0     zw   -yw    yz    0     0    -w     0     z    -y     xw   -xz    xy   -q     x    
 xyzw |  0    -yzw   xzw  -xyw   xyz  -zw    yw   -xw   -yz    xz   -xy    w    -z     y    -x     q    

*/


// Shrinker equation set for test purposes
c.q    =  + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.w   *b.w    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xw  *b.xw   - a.yw  *b.yw   - a.zw  *b.zw   - a.xyz *b.xyz  - a.xyw *b.xyw  - a.xzw *b.xzw  - a.yzw *b.yzw  + a.xyzw*b.xyzw ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   - a.w   *b.xw   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  + a.xw  *b.w    - a.yw  *b.xyw  - a.zw  *b.xzw  - a.xyz *b.yz   - a.xyw *b.yw   - a.xzw *b.zw   + a.yzw *b.xyzw - a.xyzw*b.yzw  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   - a.w   *b.yw   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    + a.xw  *b.xyw  + a.yw  *b.w    - a.zw  *b.yzw  + a.xyz *b.xz   + a.xyw *b.xw   - a.xzw *b.xyzw - a.yzw *b.zw   + a.xyzw*b.xzw  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   - a.w   *b.zw   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    + a.xw  *b.xzw  + a.yw  *b.yzw  + a.zw  *b.w    - a.xyz *b.xy   + a.xyw *b.xyzw + a.xzw *b.xw   + a.yzw *b.yw   - a.xyzw*b.xyw  ; 
c.w    =  + a.x   *b.xw   + a.y   *b.yw   + a.z   *b.zw   - a.xy  *b.xyw  - a.xz  *b.xzw  - a.yz  *b.yzw  - a.xw  *b.x    - a.yw  *b.y    - a.zw  *b.z    - a.xyz *b.xyzw - a.xyw *b.xy   - a.xzw *b.xz   - a.yzw *b.yz   + a.xyzw*b.xyz  ; 
c.xy   =  + a.z   *b.xyz  + a.w   *b.xyw  - a.zw  *b.xyzw + a.xyz *b.z    + a.xyw *b.w    - a.xzw *b.yzw  + a.yzw *b.xzw  - a.xyzw*b.zw   ; 
c.xz   =  - a.y   *b.xyz  + a.w   *b.xzw  + a.yw  *b.xyzw - a.xyz *b.y    + a.xyw *b.yzw  + a.xzw *b.w    - a.yzw *b.xyw  + a.xyzw*b.yw   ; 
c.yz   =  + a.x   *b.xyz  + a.w   *b.yzw  - a.xw  *b.xyzw + a.xyz *b.x    - a.xyw *b.xzw  + a.xzw *b.xyw  + a.yzw *b.w    - a.xyzw*b.xw   ; 
c.xw   =  - a.y   *b.xyw  - a.z   *b.xzw  - a.yz  *b.xyzw - a.xyz *b.yzw  - a.xyw *b.y    - a.xzw *b.z    + a.yzw *b.xyz  - a.xyzw*b.yz   ; 
c.yw   =  + a.x   *b.xyw  - a.z   *b.yzw  + a.xz  *b.xyzw + a.xyz *b.xzw  + a.xyw *b.x    - a.xzw *b.xyz  - a.yzw *b.z    + a.xyzw*b.xz   ; 
c.zw   =  + a.x   *b.xzw  + a.y   *b.yzw  - a.xy  *b.xyzw - a.xyz *b.xyw  + a.xyw *b.xyz  + a.xzw *b.x    + a.yzw *b.y    - a.xyzw*b.xy   ; 
c.xyz  =  - a.w   *b.xyzw + a.xyzw*b.w    ; 
c.xyw  =  + a.z   *b.xyzw - a.xyzw*b.z    ; 
c.xzw  =  - a.y   *b.xyzw + a.xyzw*b.y    ; 
c.yzw  =  + a.x   *b.xyzw - a.xyzw*b.x    ; 
c.xyzw = 0 ; 



	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);


	return c;
}

///////////////////////////////////////////////////////

GA4E Symmetric(const GA4E &a, const GA4E &b)
{
/*
Symmetric Product
 ?    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
 x    |  x     q     0     0     0     0     0     xyz   0     xyw   xzw   yz    yw    zw    0     0    
 y    |  y     0     q     0     0     0    -xyz   0    -xyw   0     yzw  -xz   -xw    0     zw    0    
 z    |  z     0     0     q     0     xyz   0     0    -xzw  -yzw   0     xy    0    -xw   -yw    0    
 w    |  w     0     0     0     q     xyw   xzw   yzw   0     0     0     0     xy    xz    yz    0    
 xy   |  xy    0     0     xyz   xyw  -q     0     0     0     0     xyzw -z    -w     0     0    -zw   
 xz   |  xz    0    -xyz   0     xzw   0    -q     0     0    -xyzw  0     y     0    -w     0     yw   
 yz   |  yz    xyz   0     0     yzw   0     0    -q     xyzw  0     0    -x     0     0    -w    -xw   
 xw   |  xw    0    -xyw  -xzw   0     0     0     xyzw -q     0     0     0     y     z     0    -yz   
 yw   |  yw    xyw   0    -yzw   0     0    -xyzw  0     0    -q     0     0    -x     0     z     xz   
 zw   |  zw    xzw   yzw   0     0     xyzw  0     0     0     0    -q     0     0    -x    -y    -xy   
 xyz  |  xyz   yz   -xz    xy    0    -z     y    -x     0     0     0    -q     0     0     0     0    
 xyw  |  xyw   yw   -xw    0     xy   -w     0     0     y    -x     0     0    -q     0     0     0    
 xzw  |  xzw   zw    0    -xw    xz    0    -w     0     z     0    -x     0     0    -q     0     0    
 yzw  |  yzw   0     zw   -yw    yz    0     0    -w     0     z    -y     0     0     0    -q     0    
 xyzw |  xyzw  0     0     0     0    -zw    yw   -xw   -yz    xz   -xy    0     0     0     0     q    

*/
	GA4E c;

//	c = (a*b + b*a)/2;

// Symmetric Product Equations
c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.w   *b.w    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xw  *b.xw   - a.yw  *b.yw   - a.zw  *b.zw   - a.xyz *b.xyz  - a.xyw *b.xyw  - a.xzw *b.xzw  - a.yzw *b.yzw  + a.xyzw*b.xyzw ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.yz  *b.xyz  - a.yw  *b.xyw  - a.zw  *b.xzw  - a.xyz *b.yz   - a.xyw *b.yw   - a.xzw *b.zw   ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    + a.xz  *b.xyz  + a.xw  *b.xyw  - a.zw  *b.yzw  + a.xyz *b.xz   + a.xyw *b.xw   - a.yzw *b.zw   ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    - a.xy  *b.xyz  + a.xw  *b.xzw  + a.yw  *b.yzw  - a.xyz *b.xy   + a.xzw *b.xw   + a.yzw *b.yw   ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    - a.xy  *b.xyw  - a.xz  *b.xzw  - a.yz  *b.yzw  - a.xyw *b.xy   - a.xzw *b.xz   - a.yzw *b.yz   ; 
c.xy   =  + a.q   *b.xy   + a.z   *b.xyz  + a.w   *b.xyw  + a.xy  *b.q    - a.zw  *b.xyzw + a.xyz *b.z    + a.xyw *b.w    - a.xyzw*b.zw   ; 
c.xz   =  + a.q   *b.xz   - a.y   *b.xyz  + a.w   *b.xzw  + a.xz  *b.q    + a.yw  *b.xyzw - a.xyz *b.y    + a.xzw *b.w    + a.xyzw*b.yw   ; 
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.w   *b.yzw  + a.yz  *b.q    - a.xw  *b.xyzw + a.xyz *b.x    + a.yzw *b.w    - a.xyzw*b.xw   ; 
c.xw   =  + a.q   *b.xw   - a.y   *b.xyw  - a.z   *b.xzw  - a.yz  *b.xyzw + a.xw  *b.q    - a.xyw *b.y    - a.xzw *b.z    - a.xyzw*b.yz   ; 
c.yw   =  + a.q   *b.yw   + a.x   *b.xyw  - a.z   *b.yzw  + a.xz  *b.xyzw + a.yw  *b.q    + a.xyw *b.x    - a.yzw *b.z    + a.xyzw*b.xz   ; 
c.zw   =  + a.q   *b.zw   + a.x   *b.xzw  + a.y   *b.yzw  - a.xy  *b.xyzw + a.zw  *b.q    + a.xzw *b.x    + a.yzw *b.y    - a.xyzw*b.xy   ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyw  =  + a.q   *b.xyw  + a.x   *b.yw   - a.y   *b.xw   + a.w   *b.xy   + a.xy  *b.w    - a.xw  *b.y    + a.yw  *b.x    + a.xyw *b.q    ; 
c.xzw  =  + a.q   *b.xzw  + a.x   *b.zw   - a.z   *b.xw   + a.w   *b.xz   + a.xz  *b.w    - a.xw  *b.z    + a.zw  *b.x    + a.xzw *b.q    ; 
c.yzw  =  + a.q   *b.yzw  + a.y   *b.zw   - a.z   *b.yw   + a.w   *b.yz   + a.yz  *b.w    - a.yw  *b.z    + a.zw  *b.y    + a.yzw *b.q    ; 
c.xyzw =  + a.q   *b.xyzw + a.xy  *b.zw   - a.xz  *b.yw   + a.yz  *b.xw   + a.xw  *b.yz   - a.yw  *b.xz   + a.zw  *b.xy   + a.xyzw*b.q    ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);


	return c;
}

//////////////////////////////////////////////////////

GA4E AntiSymmetric(const GA4E &a, const GA4E &b)
{
/*

AntiSymmetric Product
 ?    |  q     x     y     z     w     xy    xz    yz    xw    yw    zw    xyz   xyw   xzw   yzw   xyzw 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    xw    y     z     0     w     0     0     0     0     0     xyzw  yzw  
 y    |  0    -xy    0     yz    yw   -x     0     z     0     w     0     0     0    -xyzw  0    -xzw  
 z    |  0    -xz   -yz    0     zw    0    -x    -y     0     0     w     0     xyzw  0     0     xyw  
 w    |  0    -xw   -yw   -zw    0     0     0     0    -x    -y    -z    -xyzw  0     0     0    -xyz  
 xy   |  0    -y     x     0     0     0    -yz    xz   -yw    xw    0     0     0    -yzw   xzw   0    
 xz   |  0    -z     0     x     0     yz    0    -xy   -zw    0     xw    0     yzw   0    -xyw   0    
 yz   |  0     0    -z     y     0    -xz    xy    0     0    -zw    yw    0    -xzw   xyw   0     0    
 xw   |  0    -w     0     0     x     yw    zw    0     0    -xy   -xz   -yzw   0     0     xyz   0    
 yw   |  0     0    -w     0     y    -xw    0     zw    xy    0    -yz    xzw   0    -xyz   0     0    
 zw   |  0     0     0    -w     z     0    -xw   -yw    xz    yz    0    -xyw   xyz   0     0     0    
 xyz  |  0     0     0     0     xyzw  0     0     0     yzw  -xzw   xyw   0    -zw    yw   -xw   -w    
 xyw  |  0     0     0    -xyzw  0     0    -yzw   xzw   0     0    -xyz   zw    0    -yz    xz    z    
 xzw  |  0     0     xyzw  0     0     yzw   0    -xyw   0     xyz   0    -yw    yz    0    -xy   -y    
 yzw  |  0    -xyzw  0     0     0    -xzw   xyw   0    -xyz   0     0     xw   -xz    xy    0     x    
 xyzw |  0    -yzw   xzw  -xyw   xyz   0     0     0     0     0     0     w    -z     y    -x     0    

*/

	GA4E c;

//	c = (a*b - b*a)/2;

// AntiSymmetric Product Equations
c.q    =  0 ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   - a.w   *b.xw   + a.xy  *b.y    + a.xz  *b.z    + a.xw  *b.w    + a.yzw *b.xyzw - a.xyzw*b.yzw  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   - a.w   *b.yw   - a.xy  *b.x    + a.yz  *b.z    + a.yw  *b.w    - a.xzw *b.xyzw + a.xyzw*b.xzw  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   - a.w   *b.zw   - a.xz  *b.x    - a.yz  *b.y    + a.zw  *b.w    + a.xyw *b.xyzw - a.xyzw*b.xyw  ; 
c.w    =  + a.x   *b.xw   + a.y   *b.yw   + a.z   *b.zw   - a.xw  *b.x    - a.yw  *b.y    - a.zw  *b.z    - a.xyz *b.xyzw + a.xyzw*b.xyz  ; 
c.xy   =  + a.x   *b.y    - a.y   *b.x    - a.xz  *b.yz   + a.yz  *b.xz   - a.xw  *b.yw   + a.yw  *b.xw   - a.xzw *b.yzw  + a.yzw *b.xzw  ; 
c.xz   =  + a.x   *b.z    - a.z   *b.x    + a.xy  *b.yz   - a.yz  *b.xy   - a.xw  *b.zw   + a.zw  *b.xw   + a.xyw *b.yzw  - a.yzw *b.xyw  ; 
c.yz   =  + a.y   *b.z    - a.z   *b.y    - a.xy  *b.xz   + a.xz  *b.xy   - a.yw  *b.zw   + a.zw  *b.yw   - a.xyw *b.xzw  + a.xzw *b.xyw  ; 
c.xw   =  + a.x   *b.w    - a.w   *b.x    + a.xy  *b.yw   + a.xz  *b.zw   - a.yw  *b.xy   - a.zw  *b.xz   - a.xyz *b.yzw  + a.yzw *b.xyz  ; 
c.yw   =  + a.y   *b.w    - a.w   *b.y    - a.xy  *b.xw   + a.yz  *b.zw   + a.xw  *b.xy   - a.zw  *b.yz   + a.xyz *b.xzw  - a.xzw *b.xyz  ; 
c.zw   =  + a.z   *b.w    - a.w   *b.z    - a.xz  *b.xw   - a.yz  *b.yw   + a.xw  *b.xz   + a.yw  *b.yz   - a.xyz *b.xyw  + a.xyw *b.xyz  ; 
c.xyz  =  - a.w   *b.xyzw + a.xw  *b.yzw  - a.yw  *b.xzw  + a.zw  *b.xyw  - a.xyw *b.zw   + a.xzw *b.yw   - a.yzw *b.xw   + a.xyzw*b.w    ; 
c.xyw  =  + a.z   *b.xyzw - a.xz  *b.yzw  + a.yz  *b.xzw  - a.zw  *b.xyz  + a.xyz *b.zw   - a.xzw *b.yz   + a.yzw *b.xz   - a.xyzw*b.z    ; 
c.xzw  =  - a.y   *b.xyzw + a.xy  *b.yzw  - a.yz  *b.xyw  + a.yw  *b.xyz  - a.xyz *b.yw   + a.xyw *b.yz   - a.yzw *b.xy   + a.xyzw*b.y    ; 
c.yzw  =  + a.x   *b.xyzw - a.xy  *b.xzw  + a.xz  *b.xyw  - a.xw  *b.xyz  + a.xyz *b.xw   - a.xyw *b.xz   + a.xzw *b.xy   - a.xyzw*b.x    ; 
c.xyzw =  + a.x   *b.yzw  - a.y   *b.xzw  + a.z   *b.xyw  - a.w   *b.xyz  + a.xyz *b.w    - a.xyw *b.z    + a.xzw *b.y    - a.yzw *b.x    ; 


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);


	return c;
}

//////////////////////////////////////////////////////

GA4E Inner(const GA4E &a, const GA4E &b)
{

	GA4E c;

	c.q    = a.w*b.w + a.x*b.x - a.xw*b.xw - a.xy*b.xy - a.xyw*b.xyw - a.xyz*b.xyz + a.xyzw*b.xyzw
		 - a.xz*b.xz - a.xzw*b.xzw + a.y*b.y - a.yw*b.yw - a.yz*b.yz - a.yzw*b.yzw + a.z*b.z - a.zw*b.zw;

	c.x    =  + (-a.w*b.xw + a.xw*b.w + a.xy*b.y - a.xyw*b.yw - a.xyz*b.yz - a.xyzw*b.yzw + a.xz*b.z
		 - a.xzw*b.zw - a.y*b.xy - a.yw*b.xyw - a.yz*b.xyz + a.yzw*b.xyzw - a.z*b.xz - a.zw*b.xzw);
	c.y    =  + (-a.w*b.yw + a.x*b.xy + a.xw*b.xyw - a.xy*b.x + a.xyw*b.xw + a.xyz*b.xz + a.xyzw*b.xzw
		 + a.xz*b.xyz - a.xzw*b.xyzw + a.yw*b.w + a.yz*b.z - a.yzw*b.zw - a.z*b.yz - a.zw*b.yzw);
	c.z    =  + (-a.w*b.zw + a.x*b.xz + a.xw*b.xzw - a.xy*b.xyz + a.xyw*b.xyzw - a.xyz*b.xy - a.xyzw*b.xyw
		 - a.xz*b.x + a.xzw*b.xw + a.y*b.yz + a.yw*b.yzw - a.yz*b.y + a.yzw*b.yw + a.zw*b.w);
	c.w    =  + (a.x*b.xw - a.xw*b.x - a.xy*b.xyw - a.xyw*b.xy - a.xyz*b.xyzw + a.xyzw*b.xyz
		 - a.xz*b.xzw - a.xzw*b.xz + a.y*b.yw - a.yw*b.y - a.yz*b.yzw - a.yzw*b.yz + a.z*b.zw - a.zw*b.z);


	c.xy = + ( a.w*b.xyw + a.xyw*b.w + a.xyz*b.z - a.xyzw*b.zw + a.z*b.xyz - a.zw*b.xyzw);
	c.xz = + ( a.w*b.xzw - a.xyz*b.y + a.xyzw*b.yw + a.xzw*b.w - a.y*b.xyz + a.yw*b.xyzw);
	c.yz = + ( a.w*b.yzw + a.x*b.xyz - a.xw*b.xyzw + a.xyz*b.x - a.xyzw*b.xw + a.yzw*b.w);
	c.xw = + (-a.xyw*b.y - a.xyzw*b.yz - a.xzw*b.z - a.y*b.xyw - a.yz*b.xyzw - a.z*b.xzw);
	c.yw = + ( a.x*b.xyw + a.xyw*b.x + a.xyzw*b.xz + a.xz*b.xyzw - a.yzw*b.z - a.z*b.yzw);
	c.zw = + ( a.x*b.xzw - a.xy*b.xyzw - a.xyzw*b.xy + a.xzw*b.x + a.y*b.yzw + a.yzw*b.y);

	c.xyz = + (-a.w*b.xyzw + a.xyzw*b.w);
	c.xyw = + (-a.xyzw*b.z + a.z*b.xyzw);
	c.xzw = + ( a.xyzw*b.y - a.y*b.xyzw);
	c.yzw = + ( a.x*b.xyzw - a.xyzw*b.x);

	c.xyzw = 0;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

GA4E LeftContraction (const GA4E &a, const GA4E &b)
{

	GA4E c;

	c.q = a.q*b.q + a.w*b.w + a.x*b.x - a.xw*b.xw - a.xy*b.xy - a.xyw*b.xyw - a.xyz*b.xyz + a.xyzw*b.xyzw
		 - a.xz*b.xz - a.xzw*b.xzw + a.y*b.y - a.yw*b.yw - a.yz*b.yz - a.yzw*b.yzw + a.z*b.z - a.zw*b.zw ;

	c.x = + (a.q*b.x - a.w*b.xw - a.y*b.xy - a.yw*b.xyw - a.yz*b.xyz + a.yzw*b.xyzw - a.z*b.xz - a.zw*b.xzw) ;
	c.y = + (a.q*b.y - a.w*b.yw + a.x*b.xy + a.xw*b.xyw + a.xz*b.xyz - a.xzw*b.xyzw - a.z*b.yz - a.zw*b.yzw) ;
	c.z = + (a.q*b.z - a.w*b.zw + a.x*b.xz + a.xw*b.xzw - a.xy*b.xyz + a.xyw*b.xyzw + a.y*b.yz + a.yw*b.yzw) ;
	c.w = + (a.q*b.w + a.x*b.xw - a.xy*b.xyw - a.xyz*b.xyzw - a.xz*b.xzw + a.y*b.yw - a.yz*b.yzw + a.z*b.zw) ;

	c.xy = + (a.q*b.xy + a.w*b.xyw + a.z*b.xyz - a.zw*b.xyzw) ;
	c.xz = + (a.q*b.xz + a.w*b.xzw - a.y*b.xyz + a.yw*b.xyzw) ;
	c.yz = + (a.q*b.yz + a.w*b.yzw + a.x*b.xyz - a.xw*b.xyzw) ;
	c.xw = + (a.q*b.xw - a.y*b.xyw - a.yz*b.xyzw - a.z*b.xzw) ;
	c.yw = + (a.q*b.yw + a.x*b.xyw + a.xz*b.xyzw - a.z*b.yzw) ;
	c.zw = + (a.q*b.zw + a.x*b.xzw - a.xy*b.xyzw + a.y*b.yzw) ;

	c.xyz = + (a.q*b.xyz - a.w*b.xyzw) ;
	c.xyw = + (a.q*b.xyw + a.z*b.xyzw) ;
	c.xzw = + (a.q*b.xzw - a.y*b.xyzw) ;
	c.yzw = + (a.q*b.yzw + a.x*b.xyzw) ;

	c.xyzw = + a.q*b.xyzw ;


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

GA4E RightContraction (const GA4E &a, const GA4E &b)
{

	GA4E c;

	c.q = a.q*b.q + a.w*b.w + a.x*b.x - a.xw*b.xw - a.xy*b.xy - a.xyw*b.xyw - a.xyz*b.xyz + a.xyzw*b.xyzw
		 - a.xz*b.xz - a.xzw*b.xzw + a.y*b.y - a.yw*b.yw - a.yz*b.yz - a.yzw*b.yzw + a.z*b.z - a.zw*b.zw ;

	c.x = + ( a.x*b.q + a.xw*b.w + a.xy*b.y - a.xyw*b.yw - a.xyz*b.yz - a.xyzw*b.yzw + a.xz*b.z - a.xzw*b.zw) ;
	c.y = + (-a.xy*b.x + a.xyw*b.xw + a.xyz*b.xz + a.xyzw*b.xzw + a.y*b.q + a.yw*b.w + a.yz*b.z - a.yzw*b.zw) ;
	c.z = + (-a.xyz*b.xy - a.xyzw*b.xyw - a.xz*b.x + a.xzw*b.xw - a.yz*b.y + a.yzw*b.yw + a.z*b.q + a.zw*b.w) ;
	c.w = + ( a.w*b.q - a.xw*b.x - a.xyw*b.xy + a.xyzw*b.xyz - a.xzw*b.xz - a.yw*b.y - a.yzw*b.yz - a.zw*b.z) ;

	c.xy = + ( a.xy*b.q + a.xyw*b.w + a.xyz*b.z - a.xyzw*b.zw) ;
	c.xz = + (-a.xyz*b.y + a.xyzw*b.yw + a.xz*b.q + a.xzw*b.w) ;
	c.yz = + ( a.xyz*b.x - a.xyzw*b.xw + a.yz*b.q + a.yzw*b.w) ;
	c.xw = + ( a.xw*b.q - a.xyw*b.y - a.xyzw*b.yz - a.xzw*b.z) ;
	c.yw = + ( a.xyw*b.x + a.xyzw*b.xz + a.yw*b.q - a.yzw*b.z) ;
	c.zw = + (-a.xyzw*b.xy + a.xzw*b.x + a.yzw*b.y + a.zw*b.q) ;

	c.xyz = + ( a.xyz*b.q + a.xyzw*b.w) ;
	c.xyw = + ( a.xyw*b.q - a.xyzw*b.z) ;
	c.xzw = + ( a.xyzw*b.y + a.xzw*b.q) ;
	c.yzw = + (-a.xyzw*b.x + a.yzw*b.q) ;

	c.xyzw = + a.xyzw*b.q ;


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.w = expand(c.w);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xw = expand(c.xw);
	c.yw = expand(c.yw);
	c.zw = expand(c.zw);
	c.xyz = expand(c.xyz);
	c.xyw = expand(c.xyw);
	c.xzw = expand(c.xzw);
	c.yzw = expand(c.yzw);
	c.xyzw = expand(c.xyzw);

	return c;
}

//////////////////////////////////////////////////////

ex Determinant(GA4E A) {

	GA4E B, C;
	ex a, b, c, d, e, s, det;

	B = Reverse(A);
	C = Product(B,A);

	a = C.q;
	b = C.x;
	c = C.y;
	d = C.z;
	e = C.w;

	s = C.xyzw;
	det = expand(a*a - b*b - c*c - d*d - e*e - s*s);

	return(det);
}

//////////////////////////////////////////////////////

GA4E Adjugate(GA4E a)
{

	GA4E u;

//Product(Reverse(r),Conjugation(Product(r,Reverse(r)))) = Adjugate(r);

	u = Product(Reverse(a),Conjugation(Product(a,Reverse(a))));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.w = expand(u.w);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xw = expand(u.xw);
	u.yw = expand(u.yw);
	u.zw = expand(u.zw);
	u.xyz = expand(u.xyz);
	u.xyw = expand(u.xyw);
	u.xzw = expand(u.xzw);
	u.yzw = expand(u.yzw);
	u.xyzw = expand(u.xyzw);

	return u;

}

//////////////////////////////////////////////////////

GA4E Reciprocal(GA4E a)
{
	ex b;
	GA4E c,d;

	b = Determinant(a);

	c = Adjugate(a);

	d.q = expand(c.q/b);
	d.x = expand(c.x/b);
	d.y = expand(c.y/b);
	d.z = expand(c.z/b);
	d.w = expand(c.w/b);
	d.xy = expand(c.xy/b);
	d.xz = expand(c.xz/b);
	d.yz = expand(c.yz/b);
	d.xw = expand(c.xw/b);
	d.yw = expand(c.yw/b);
	d.zw = expand(c.zw/b);
	d.xyz = expand(c.xyz/b);
	d.xyw = expand(c.xyw/b);
	d.xzw = expand(c.xzw/b);
	d.yzw = expand(c.yzw/b);
	d.xyzw = expand(c.xyzw/b);

	return d;
}

//////////////////////////////////////////////////////

GA4E	Comp(GA4E u, int i) {

	GA4E MV;
	ex a, b,c,d,e, f,g,h,j,k,l, m,n,p,r, s;
	int Mask, Sign[16];

	Mask = 32;
	if((Mask & i) == 0) Sign[5] = +1; else Sign[5] = -1;  Mask >>= 1;  // q
	if((Mask & i) == 0) Sign[4] = +1; else Sign[4] = -1;  Mask >>= 1;  // x
	if((Mask & i) == 0) Sign[3] = +1; else Sign[3] = -1;  Mask >>= 1;  // y
	if((Mask & i) == 0) Sign[2] = +1; else Sign[2] = -1;  Mask >>= 1;  // z
	if((Mask & i) == 0) Sign[1] = +1; else Sign[1] = -1;  Mask >>= 1;  // w
	if((Mask & i) == 0) Sign[0] = +1; else Sign[0] = -1;  Mask >>= 1;  // xy

	Sign[6]  = Sign[0]*Sign[2]*Sign[3];	 
	Sign[7]  = Sign[0]*Sign[1]*Sign[3];	
	Sign[8]  = Sign[0]*Sign[2]*Sign[4];	
	Sign[9]  = Sign[0]*Sign[1]*Sign[4];	
	Sign[10] = Sign[0]*Sign[1]*Sign[2]*Sign[3]*Sign[4];	
	
	Sign[11] = Sign[0]*Sign[2]*Sign[5];	
	Sign[12] = Sign[0]*Sign[1]*Sign[5];	
	Sign[13] = Sign[0]*Sign[1]*Sign[2]*Sign[3]*Sign[5]; 
	Sign[14] = Sign[0]*Sign[1]*Sign[2]*Sign[4]*Sign[5];
	
	Sign[15] = Sign[1]*Sign[2]*Sign[3]*Sign[4]*Sign[5];

	a = Sign[5]*u.q;   b = Sign[4]*u.x;   c = Sign[3]*u.y;   d = Sign[2]*u.z;   e = Sign[1]*u.w;
	f = Sign[0]*u.xy;  g = Sign[6]*u.xz;  j = Sign[7]*u.xw;  h = Sign[8]*u.yz;  k = Sign[9]*u.yw;  l = Sign[10]*u.zw;
	m = Sign[11]*u.xyz;  n = Sign[12]*u.xyw;  p = Sign[13]*u.xzw;   r = Sign[14]*u.yzw;  s = Sign[15]*u.xyzw;

	MV = GA4E(a, b,c,d,e, f,g,h,j,k,l, m,n,p,r, s) ; 
	
	return MV;
}

//////////////////////////////////////////////////////

GA4E Magic(const GA4E &u, int i) {  

	ex  a,  b, c, d, e,  f, g, h, j, k, l,  m, n, p, r,  s;

	a = u.q;	b = u.x;	c = u.y;	d = u.z;	e = u.w;
	f = u.xy;	g = u.xz;	h = u.yz;	j = u.xw;	k = u.yw;	l = u.zw;
	m = u.xyz;	n = u.xyw;	p = u.xzw;	r = u.yzw;	s = u.xyzw;

	GA4E MV;
	MV = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);

	if(i ==   0) MV = GA4E( a, +b,+c,+d,+e, +f,+g,+h,+j,+k,+l, +m,+n,+p,+r, +s) ; 
	if(i ==   1) MV = GA4E( a, +b,+c,+d,-s, +f,+g,+h,-r,+p,-n, +m,+l,-k,+j, +e) ; 
	if(i ==   2) MV = GA4E( a, +b,+c,+e,+d, +f,+j,+k,+g,+h,-l, +n,+m,-p,-r, -s) ; 
	if(i ==   3) MV = GA4E( a, +b,+c,+s,+d, +f,+r,-p,+g,+h,-n, -l,+m,-k,+j, +e) ; 
	if(i ==   4) MV = GA4E( a, +b,+c,+e,+s, +f,+j,+k,+r,-p,-m, +n,-l,-h,+g, +d) ; 
	if(i ==   5) MV = GA4E( a, +b,+c,-s,+e, +f,-r,+p,+j,+k,-m, +l,+n,-h,+g, +d) ; 
	if(i ==   6) MV = GA4E( a, +b,+d,+c,+e, +g,+f,-h,+j,+l,+k, -m,+p,+n,-r, -s) ; 
	if(i ==   7) MV = GA4E( a, +b,+d,+c,+s, +g,+f,-h,+r,+n,-p, -m,+k,-l,+j, +e) ; 
	if(i ==   8) MV = GA4E( a, +b,+e,+c,+d, +j,+f,-k,+g,-l,+h, -n,-p,+m,+r, +s) ; 
	if(i ==   9) MV = GA4E( a, +b,-s,+c,+d, -r,+f,-p,+g,+n,+h, -l,+k,+m,+j, +e) ; 
	if(i ==  10) MV = GA4E( a, +b,+e,+c,-s, +j,+f,-k,-r,+m,+p, -n,+h,+l,+g, +d) ; 
	if(i ==  11) MV = GA4E( a, +b,+s,+c,+e, +r,+f,+p,+j,+m,+k, +l,+h,+n,+g, +d) ; 
	if(i ==  12) MV = GA4E( a, +b,+d,+e,+c, +g,+j,+l,+f,-h,-k, +p,-m,-n,+r, +s) ; 
	if(i ==  13) MV = GA4E( a, +b,+d,-s,+c, +g,-r,-n,+f,-h,-p, -k,-m,-l,+j, +e) ; 
	if(i ==  14) MV = GA4E( a, +b,+e,+d,+c, +j,+g,-l,+f,-k,-h, -p,-n,-m,-r, -s) ; 
	if(i ==  15) MV = GA4E( a, +b,+s,+d,+c, +r,+g,-n,+f,+p,-h, -k,+l,-m,+j, +e) ; 
	if(i ==  16) MV = GA4E( a, +b,+e,+s,+c, +j,+r,-m,+f,-k,+p, -h,-n,+l,+g, +d) ; 
	if(i ==  17) MV = GA4E( a, +b,-s,+e,+c, -r,+j,-m,+f,-p,-k, -h,-l,-n,+g, +d) ; 
	if(i ==  18) MV = GA4E( a, +b,+d,+e,-s, +g,+j,+l,-r,-n,+m, +p,-k,+h,+f, +c) ; 
	if(i ==  19) MV = GA4E( a, +b,+d,+s,+e, +g,+r,+n,+j,+l,+m, +k,+p,+h,+f, +c) ; 
	if(i ==  20) MV = GA4E( a, +b,+e,+d,+s, +j,+g,-l,+r,-m,+n, -p,-h,+k,+f, +c) ; 
	if(i ==  21) MV = GA4E( a, +b,-s,+d,+e, -r,+g,+n,+j,-m,+l, +k,-h,+p,+f, +c) ; 
	if(i ==  22) MV = GA4E( a, +b,+e,-s,+d, +j,-r,+m,+g,-l,+n, +h,-p,+k,+f, +c) ; 
	if(i ==  23) MV = GA4E( a, +b,+s,+e,+d, +r,+j,+m,+g,-n,-l, +h,-k,-p,+f, +c) ; 
	if(i ==  24) MV = GA4E( a, +c,+b,+d,+e, -f,+h,+g,+k,+j,+l, -m,-n,+r,+p, -s) ; 
	if(i ==  25) MV = GA4E( a, +c,+b,+d,+s, -f,+h,+g,-p,+r,+n, -m,+l,-j,+k, +e) ; 
	if(i ==  26) MV = GA4E( a, +c,+b,+e,+d, -f,+k,+j,+h,+g,-l, -n,-m,-r,-p, +s) ; 
	if(i ==  27) MV = GA4E( a, +c,+b,-s,+d, -f,+p,-r,+h,+g,+n, -l,-m,-j,+k, +e) ; 
	if(i ==  28) MV = GA4E( a, +c,+b,+e,-s, -f,+k,+j,+p,-r,+m, -n,-l,-g,+h, +d) ; 
	if(i ==  29) MV = GA4E( a, +c,+b,+s,+e, -f,-p,+r,+k,+j,+m, +l,-n,-g,+h, +d) ; 
	if(i ==  30) MV = GA4E( a, +d,+b,+c,+e, -g,-h,+f,+l,+j,+k, +m,-p,-r,+n, +s) ; 
	if(i ==  31) MV = GA4E( a, +d,+b,+c,-s, -g,-h,+f,-n,-r,+p, +m,+k,-j,+l, +e) ; 
	if(i ==  32) MV = GA4E( a, +e,+b,+c,+d, -j,-k,+f,-l,+g,+h, +n,+p,+r,+m, -s) ; 
	if(i ==  33) MV = GA4E( a, +s,+b,+c,+d, -r,+p,+f,-n,+g,+h, -l,+k,-j,+m, +e) ; 
	if(i ==  34) MV = GA4E( a, +e,+b,+c,+s, -j,-k,+f,-m,+r,-p, +n,+h,-g,-l, +d) ; 
	if(i ==  35) MV = GA4E( a, -s,+b,+c,+e, +r,-p,+f,-m,+j,+k, +l,+h,-g,+n, +d) ; 
	if(i ==  36) MV = GA4E( a, +d,+b,+e,+c, -g,+l,+j,-h,+f,-k, -p,+m,+r,-n, -s) ; 
	if(i ==  37) MV = GA4E( a, +d,+b,+s,+c, -g,+n,+r,-h,+f,+p, -k,+m,-j,+l, +e) ; 
	if(i ==  38) MV = GA4E( a, +e,+b,+d,+c, -j,-l,+g,-k,+f,-h, +p,+n,-r,-m, +s) ; 
	if(i ==  39) MV = GA4E( a, -s,+b,+d,+c, +r,+n,+g,-p,+f,-h, -k,+l,-j,-m, +e) ; 
	if(i ==  40) MV = GA4E( a, +e,+b,-s,+c, -j,+m,-r,-k,+f,-p, -h,+n,-g,-l, +d) ; 
	if(i ==  41) MV = GA4E( a, +s,+b,+e,+c, -r,+m,+j,+p,+f,-k, -h,-l,-g,-n, +d) ; 
	if(i ==  42) MV = GA4E( a, +d,+b,+e,+s, -g,+l,+j,+n,+r,-m, -p,-k,-f,-h, +c) ; 
	if(i ==  43) MV = GA4E( a, +d,+b,-s,+e, -g,-n,-r,+l,+j,-m, +k,-p,-f,-h, +c) ; 
	if(i ==  44) MV = GA4E( a, +e,+b,+d,-s, -j,-l,+g,+m,-r,-n, +p,-h,-f,-k, +c) ; 
	if(i ==  45) MV = GA4E( a, +s,+b,+d,+e, -r,-n,+g,+m,+j,+l, +k,-h,-f,+p, +c) ; 
	if(i ==  46) MV = GA4E( a, +e,+b,+s,+d, -j,-m,+r,-l,+g,-n, +h,+p,-f,-k, +c) ; 
	if(i ==  47) MV = GA4E( a, -s,+b,+e,+d, +r,-m,+j,+n,+g,-l, +h,-k,-f,-p, +c) ; 
	if(i ==  48) MV = GA4E( a, +c,+d,+b,+e, +h,-f,-g,+k,+l,+j, +m,+r,-n,-p, +s) ; 
	if(i ==  49) MV = GA4E( a, +c,+d,+b,-s, +h,-f,-g,+p,-n,-r, +m,+j,-l,+k, +e) ; 
	if(i ==  50) MV = GA4E( a, +c,+e,+b,+d, +k,-f,-j,+h,-l,+g, +n,-r,-m,+p, -s) ; 
	if(i ==  51) MV = GA4E( a, +c,+s,+b,+d, -p,-f,-r,+h,-n,+g, -l,+j,-m,+k, +e) ; 
	if(i ==  52) MV = GA4E( a, +c,+e,+b,+s, +k,-f,-j,-p,-m,+r, +n,+g,+l,+h, +d) ; 
	if(i ==  53) MV = GA4E( a, +c,-s,+b,+e, +p,-f,+r,+k,-m,+j, +l,+g,-n,+h, +d) ; 
	if(i ==  54) MV = GA4E( a, +d,+c,+b,+e, -h,-g,-f,+l,+k,+j, -m,-r,-p,-n, -s) ; 
	if(i ==  55) MV = GA4E( a, +d,+c,+b,+s, -h,-g,-f,+n,-p,+r, -m,+j,-k,+l, +e) ; 
	if(i ==  56) MV = GA4E( a, +e,+c,+b,+d, -k,-j,-f,-l,+h,+g, -n,+r,+p,-m, +s) ; 
	if(i ==  57) MV = GA4E( a, -s,+c,+b,+d, -p,+r,-f,+n,+h,+g, -l,+j,-k,-m, +e) ; 
	if(i ==  58) MV = GA4E( a, +e,+c,+b,-s, -k,-j,-f,+m,+p,-r, -n,+g,-h,-l, +d) ; 
	if(i ==  59) MV = GA4E( a, +s,+c,+b,+e, +p,-r,-f,+m,+k,+j, +l,+g,-h,-n, +d) ; 
	if(i ==  60) MV = GA4E( a, +d,+e,+b,+c, +l,-g,-j,-h,-k,+f, +p,+r,+m,+n, +s) ; 
	if(i ==  61) MV = GA4E( a, +d,-s,+b,+c, -n,-g,+r,-h,-p,+f, -k,+j,+m,+l, +e) ; 
	if(i ==  62) MV = GA4E( a, +e,+d,+b,+c, -l,-j,-g,-k,-h,+f, -p,-r,+n,+m, -s) ; 
	if(i ==  63) MV = GA4E( a, +s,+d,+b,+c, -n,-r,-g,+p,-h,+f, -k,+j,-l,+m, +e) ; 
	if(i ==  64) MV = GA4E( a, +e,+s,+b,+c, -m,-j,-r,-k,+p,+f, -h,+g,+n,-l, +d) ; 
	if(i ==  65) MV = GA4E( a, -s,+e,+b,+c, -m,+r,-j,-p,-k,+f, -h,+g,+l,+n, +d) ; 
	if(i ==  66) MV = GA4E( a, +d,+e,+b,-s, +l,-g,-j,-n,+m,-r, +p,+f,+k,-h, +c) ; 
	if(i ==  67) MV = GA4E( a, +d,+s,+b,+e, +n,-g,-r,+l,+m,+j, +k,+f,-p,-h, +c) ; 
	if(i ==  68) MV = GA4E( a, +e,+d,+b,+s, -l,-j,-g,-m,+n,+r, -p,+f,+h,-k, +c) ; 
	if(i ==  69) MV = GA4E( a, -s,+d,+b,+e, +n,+r,-g,-m,+l,+j, +k,+f,+h,-p, +c) ; 
	if(i ==  70) MV = GA4E( a, +e,-s,+b,+d, +m,-j,+r,-l,+n,+g, +h,+f,+p,-k, +c) ; 
	if(i ==  71) MV = GA4E( a, +s,+e,+b,+d, +m,-r,-j,-n,-l,+g, +h,+f,+k,+p, +c) ; 
	if(i ==  72) MV = GA4E( a, +c,+d,+e,+b, +h,+k,+l,-f,-g,-j, +r,+m,+n,+p, -s) ; 
	if(i ==  73) MV = GA4E( a, +c,+d,+s,+b, +h,-p,+n,-f,-g,-r, -j,+m,-l,+k, +e) ; 
	if(i ==  74) MV = GA4E( a, +c,+e,+d,+b, +k,+h,-l,-f,-j,-g, -r,+n,+m,-p, +s) ; 
	if(i ==  75) MV = GA4E( a, +c,-s,+d,+b, +p,+h,+n,-f,+r,-g, -j,+l,+m,+k, +e) ; 
	if(i ==  76) MV = GA4E( a, +c,+e,-s,+b, +k,+p,+m,-f,-j,+r, -g,+n,+l,+h, +d) ; 
	if(i ==  77) MV = GA4E( a, +c,+s,+e,+b, -p,+k,+m,-f,-r,-j, -g,-l,+n,+h, +d) ; 
	if(i ==  78) MV = GA4E( a, +d,+c,+e,+b, -h,+l,+k,-g,-f,-j, -r,-m,+p,+n, +s) ; 
	if(i ==  79) MV = GA4E( a, +d,+c,-s,+b, -h,-n,+p,-g,-f,+r, -j,-m,-k,+l, +e) ; 
	if(i ==  80) MV = GA4E( a, +e,+c,+d,+b, -k,-l,+h,-j,-f,-g, +r,-n,-p,+m, -s) ; 
	if(i ==  81) MV = GA4E( a, +s,+c,+d,+b, +p,-n,+h,-r,-f,-g, -j,+l,-k,+m, +e) ; 
	if(i ==  82) MV = GA4E( a, +e,+c,+s,+b, -k,-m,-p,-j,-f,-r, -g,-n,-h,-l, +d) ; 
	if(i ==  83) MV = GA4E( a, -s,+c,+e,+b, -p,-m,+k,+r,-f,-j, -g,-l,-h,+n, +d) ; 
	if(i ==  84) MV = GA4E( a, +d,+e,+c,+b, +l,-h,-k,-g,-j,-f, +r,+p,-m,-n, -s) ; 
	if(i ==  85) MV = GA4E( a, +d,+s,+c,+b, +n,-h,+p,-g,-r,-f, -j,+k,-m,+l, +e) ; 
	if(i ==  86) MV = GA4E( a, +e,+d,+c,+b, -l,-k,-h,-j,-g,-f, -r,-p,-n,-m, +s) ; 
	if(i ==  87) MV = GA4E( a, -s,+d,+c,+b, +n,-p,-h,+r,-g,-f, -j,+k,-l,-m, +e) ; 
	if(i ==  88) MV = GA4E( a, +e,-s,+c,+b, +m,-k,-p,-j,+r,-f, -g,+h,-n,-l, +d) ; 
	if(i ==  89) MV = GA4E( a, +s,+e,+c,+b, +m,+p,-k,-r,-j,-f, -g,+h,+l,-n, +d) ; 
	if(i ==  90) MV = GA4E( a, +d,+e,+s,+b, +l,+n,-m,-g,-j,-r, -f,+p,+k,-h, +c) ; 
	if(i ==  91) MV = GA4E( a, +d,-s,+e,+b, -n,+l,-m,-g,+r,-j, -f,-k,+p,-h, +c) ; 
	if(i ==  92) MV = GA4E( a, +e,+d,-s,+b, -l,+m,-n,-j,-g,+r, -f,-p,+h,-k, +c) ; 
	if(i ==  93) MV = GA4E( a, +s,+d,+e,+b, -n,+m,+l,-r,-g,-j, -f,-k,+h,+p, +c) ; 
	if(i ==  94) MV = GA4E( a, +e,+s,+d,+b, -m,-l,-n,-j,-r,-g, -f,-h,-p,-k, +c) ; 
	if(i ==  95) MV = GA4E( a, -s,+e,+d,+b, -m,+n,-l,+r,-j,-g, -f,-h,+k,-p, +c) ; 
	if(i ==  96) MV = GA4E( a, +c,+d,+e,+s, +h,+k,+l,-p,+n,-m, +r,-j,+g,-f, +b) ; 
	if(i ==  97) MV = GA4E( a, +c,+d,-s,+e, +h,+p,-n,+k,+l,-m, +j,+r,+g,-f, +b) ; 
	if(i ==  98) MV = GA4E( a, +c,+e,+d,-s, +k,+h,-l,+p,+m,-n, -r,-g,+j,-f, +b) ; 
	if(i ==  99) MV = GA4E( a, +c,+s,+d,+e, -p,+h,-n,+k,+m,+l, +j,-g,+r,-f, +b) ; 
	if(i == 100) MV = GA4E( a, +c,+e,+s,+d, +k,-p,-m,+h,-l,-n, +g,-r,+j,-f, +b) ; 
	if(i == 101) MV = GA4E( a, +c,-s,+e,+d, +p,+k,-m,+h,+n,-l, +g,-j,-r,-f, +b) ; 
	if(i == 102) MV = GA4E( a, +d,+c,+e,-s, -h,+l,+k,-n,+p,+m, -r,-j,+f,-g, +b) ; 
	if(i == 103) MV = GA4E( a, +d,+c,+s,+e, -h,+n,-p,+l,+k,+m, +j,-r,+f,-g, +b) ; 
	if(i == 104) MV = GA4E( a, +e,+c,+d,+s, -k,-l,+h,-m,-p,+n, +r,-g,+f,-j, +b) ; 
	if(i == 105) MV = GA4E( a, -s,+c,+d,+e, -p,+n,+h,-m,+k,+l, +j,-g,+f,+r, +b) ; 
	if(i == 106) MV = GA4E( a, +e,+c,-s,+d, -k,+m,+p,-l,+h,+n, +g,+r,+f,-j, +b) ; 
	if(i == 107) MV = GA4E( a, +s,+c,+e,+d, +p,+m,+k,-n,+h,-l, +g,-j,+f,-r, +b) ; 
	if(i == 108) MV = GA4E( a, +d,+e,+c,+s, +l,-h,-k,+n,-m,-p, +r,-f,+j,-g, +b) ; 
	if(i == 109) MV = GA4E( a, +d,-s,+c,+e, -n,-h,-p,+l,-m,+k, +j,-f,-r,-g, +b) ; 
	if(i == 110) MV = GA4E( a, +e,+d,+c,-s, -l,-k,-h,+m,-n,+p, -r,-f,+g,-j, +b) ; 
	if(i == 111) MV = GA4E( a, +s,+d,+c,+e, -n,+p,-h,+m,+l,+k, +j,-f,+g,-r, +b) ; 
	if(i == 112) MV = GA4E( a, +e,+s,+c,+d, -m,-k,+p,-l,-n,+h, +g,-f,+r,-j, +b) ; 
	if(i == 113) MV = GA4E( a, -s,+e,+c,+d, -m,-p,-k,+n,-l,+h, +g,-f,+j,+r, +b) ; 
	if(i == 114) MV = GA4E( a, +d,+e,-s,+c, +l,-n,+m,-h,-k,-p, +f,+r,+j,-g, +b) ; 
	if(i == 115) MV = GA4E( a, +d,+s,+e,+c, +n,+l,+m,-h,+p,-k, +f,-j,+r,-g, +b) ; 
	if(i == 116) MV = GA4E( a, +e,+d,+s,+c, -l,-m,+n,-k,-h,+p, +f,-r,+g,-j, +b) ; 
	if(i == 117) MV = GA4E( a, -s,+d,+e,+c, +n,-m,+l,-p,-h,-k, +f,-j,+g,+r, +b) ; 
	if(i == 118) MV = GA4E( a, +e,-s,+d,+c, +m,-l,+n,-k,-p,-h, +f,-g,-r,-j, +b) ; 
	if(i == 119) MV = GA4E( a, +s,+e,+d,+c, +m,-n,-l,+p,-k,-h, +f,-g,+j,-r, +b) ; 
	
	return MV;
}

////////////////////////////////////


