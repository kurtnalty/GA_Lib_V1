// Routines for Geometric Algebra in Three Dimensional Euclidean spacetime
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

struct GA3E{	// 0 . 1 2 4 . 3 5 6 . 7 
	ex q,  x,y,z, xy,xz,yz,  xyz;
	GA3E() {q = 0; x = 0; y = 0; z = 0;
		xy = 0; xz = 0; yz = 0;
		xyz = 0;}
	GA3E(ex qq, ex xx, ex yy, ex zz, 
		ex xxyy, ex xxzz, ex yyzz,
		ex xxyyzz) 

		{q = qq; x = xx; y = yy; z = zz;
		 xy = xxyy; xz = xxzz; yz = yyzz;
		 xyz = xxyyzz;}
};


//////////////////////////////////////////////////////

// Necessary forward declarations

GA3E LeftContraction (const GA3E &a, const GA3E &b) ;
GA3E Product(const GA3E &a, const GA3E &b) ;

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, GA3E &v) 
{
	return ff << "\n(" 
		<< v.q << ", \n" 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << ", \n" 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << ", \n" 
		<< v.xyz << ")\n";
}


//////////////////////////////////////////////////////

void PrintMV(GA3E &v) 
{
	cout 	<< "( "
		<< v.q << ", " 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << ", " 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << ", " 
		<< v.xyz << ")";
}

//////////////////////////////////////////////////////

GA3E OverBar(GA3E a)	 // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3E b;		// a blade wedge b blade = pseudovector blade
// OverBar(Blade[i]) = Reciprocal(Blade[i])*xyz ;

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     = -a.xz;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.xz   =  -a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3E UnderBar(GA3E a)   // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3E b;		// b blade wedge a blade = pseudovector blade
// UnderBar(Blade[i]) = xyz*Reciprocal(Blade[i]) ;

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     = -a.xz;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.xz   =  -a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3E Reverse(GA3E w)	
// (A,  B, C, D,  -E,-F,-G, -H) 
{
	GA3E v;
	v.q =  w.q;

	v.x = w.x;
	v.y = w.y;
	v.z = w.z;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Involution(GA3E w)	
// Corresponds to P parity transform: x -> -x, y -> -y, z ->-z
// (A, -B,-C,-D,  E, F, G,  -H) 
{
	GA3E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = w.xy;
	v.xz = w.xz;
	v.yz = w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Hermitian(GA3E w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(A,  B, C, D, -E,-F,-G, -H)  
{
	GA3E v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Hermitian_Component(GA3E w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(A,  B, C, D, 0, 0, 0,  0)  
{
	GA3E v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = 0;
	v.xz = 0;
	v.yz = 0;

	v.xyz = 0;

	return v;
}

//////////////////////////////////////////////////////

GA3E Anti_Hermitian_Component(GA3E w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(0, 0, 0, 0,  E, F, G,  H)  
{
	GA3E v;

	v.q =  0;

	v.x =  0;
	v.y =  0;
	v.z =  0;

	v.xy = w.xy;
	v.xz = w.xz;
	v.yz = w.yz;

	v.xyz = w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Transpose(GA3E w)	// only correct for 4x4 real matrix implementation
//(A,  B, C, D, -E,-F,-G, -H)  
{
	GA3E v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Conjugation(GA3E w)
// (A, -B,-C,-D, -E,-F,-G, -H)
{
	GA3E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E CliffordConjugation(GA3E w)
// (A, -B,-C,-D, -E,-F,-G,  H) 
{
	GA3E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;

	v.xyz = w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3E Dual(GA3E w)   // return w*I_inv 
// Dual(r) = u = (h, g,-f, e, -d, c,-b, -a)

{
	GA3E v;
//	GA3E I_inv;
//	I_inv.xyz = -1;
//	v = Product(w,I_inv);

	v.q =  w.xyz;

	v.x =  w.yz;
	v.y = -w.xz;
	v.z =  w.xy;

	v.xy = -w.z;
	v.xz =  w.y;
	v.yz = -w.x;

	v.xyz = -w.q;

	return v;
}



//////////////////////////////////////////////////////

GA3E DorstDual(GA3E a)  //DorstDual(r) = u = (h, g,-f,e, -d,c,-b, -a)
{
	GA3E b;
//	GA3E I_inv;	I_inv.xyz = -1;

//	b = LeftContraction(a,I_inv);

	b.q =  a.xyz;

	b.x =  a.yz;
	b.y = -a.xz;
	b.z =  a.xy;

	b.xy = -a.z;
	b.xz =  a.y;
	b.yz = -a.x;

	b.xyz = -a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3E DorstUnDual(GA3E a)  //DorstUnDual(r) = u = (-h,-g, f,-e,  d,-c, b,  a)
{
	GA3E b;
//	GA3E I;	I.xyz = 1;

//	b = LeftContraction(a,I);

	b.q = -a.xyz;

	b.x = -a.yz;
	b.y =  a.xz;
	b.z = -a.xy;

	b.xy =  a.z;
	b.xz = -a.y;
	b.yz =  a.x;

	b.xyz =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3E operator+(const GA3E &u, const GA3E &v)
{
	GA3E w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.xy  = u.xy  + v.xy ;
	w.xz  = u.xz  + v.xz ;
	w.yz  = u.yz  + v.yz ;
	w.xyz = u.xyz + v.xyz;

	return w;
}

//////////////////////////////////////////////////////

GA3E operator-(const GA3E &u, const GA3E &v)
{
	GA3E w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.xy  = u.xy  - v.xy ;
	w.xz  = u.xz  - v.xz ;
	w.yz  = u.yz  - v.yz ;
	w.xyz = u.xyz - v.xyz;

	return w;
}

//////////////////////////////////////////////////////

int operator==(const GA3E &u, const GA3E &v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&

			(u.xy==v.xy)&&
			(u.xz==v.xz)&&
			(u.yz==v.yz)&&

			(u.xyz==v.xyz);
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const GA3E &u, const GA3E &v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||

			(u.xy!=v.xy)||
			(u.xz!=v.xz)||
			(u.yz!=v.yz)||

			(u.xyz!=v.xyz);
	return result;
}

//////////////////////////////////////////////////////

GA3E operator*(const GA3E &a, const GA3E &b) {
	GA3E c;

/*

 *    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     q     xy    xz    y     z     xyz   yz   
 y    |  y    -xy    q     yz   -x    -xyz   z    -xz   
 z    |  z    -xz   -yz    q     xyz  -x    -y     xy   
 xy   |  xy   -y     x     xyz  -q    -yz    xz   -z    
 xz   |  xz   -z    -xyz   x     yz   -q    -xy    y    
 yz   |  yz    xyz  -z     y    -xz    xy   -q    -x    
 xyz  |  xyz   yz   -xz    xy   -z     y    -x    -q    

*/
	c.q   =  + a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.xz*b.xz  - a.yz*b.yz  - a.xyz*b.xyz ;
	c.x   =  + a.q*b.x   + a.x*b.q   - a.y*b.xy  - a.z*b.xz  + a.xy*b.y   + a.xz*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   =  + a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   + a.xz*b.xyz + a.yz*b.z   + a.xyz*b.xz ;
	c.z   =  + a.q*b.z   + a.x*b.xz  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz - a.xz*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy  =  + a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   - a.xz*b.yz  + a.yz*b.xz  + a.xyz*b.z  ;
	c.xz  =  + a.q*b.xz  + a.x*b.z   - a.y*b.xyz - a.z*b.x   + a.xy*b.yz  + a.xz*b.q   - a.yz*b.xy  - a.xyz*b.y  ;
	c.yz  =  + a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   - a.xy*b.xz  + a.xz*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz =  + a.q*b.xyz + a.x*b.yz  - a.y*b.xz  + a.z*b.xy  + a.xy*b.z   - a.xz*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E Product(const GA3E &a, const GA3E &b)
{

/*

 *    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     q     xy    xz    y     z     xyz   yz   
 y    |  y    -xy    q     yz   -x    -xyz   z    -xz   
 z    |  z    -xz   -yz    q     xyz  -x    -y     xy   
 xy   |  xy   -y     x     xyz  -q    -yz    xz   -z    
 xz   |  xz   -z    -xyz   x     yz   -q    -xy    y    
 yz   |  yz    xyz  -z     y    -xz    xy   -q    -x    
 xyz  |  xyz   yz   -xz    xy   -z     y    -x    -q    

*/
	GA3E c;

	c.q   =  + a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.xz*b.xz  - a.yz*b.yz  - a.xyz*b.xyz ;
	c.x   =  + a.q*b.x   + a.x*b.q   - a.y*b.xy  - a.z*b.xz  + a.xy*b.y   + a.xz*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   =  + a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   + a.xz*b.xyz + a.yz*b.z   + a.xyz*b.xz ;
	c.z   =  + a.q*b.z   + a.x*b.xz  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz - a.xz*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy  =  + a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   - a.xz*b.yz  + a.yz*b.xz  + a.xyz*b.z  ;
	c.xz  =  + a.q*b.xz  + a.x*b.z   - a.y*b.xyz - a.z*b.x   + a.xy*b.yz  + a.xz*b.q   - a.yz*b.xy  - a.xyz*b.y  ;
	c.yz  =  + a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   - a.xy*b.xz  + a.xz*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz =  + a.q*b.xyz + a.x*b.yz  - a.y*b.xz  + a.z*b.xy  + a.xy*b.z   - a.xz*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E operator*(const GA3E &v, const ex &u) {
	GA3E w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xyz  = expand(u*v.xyz ) ;

	return w;
}

//////////////////////////////////////////////////////

GA3E operator*(const ex &u, const GA3E &v) {
	GA3E w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xyz  = expand(u*v.xyz ) ;

	return w;
}


//////////////////////////////////////////////////////

GA3E operator/(const GA3E &u, const int i)
{
	GA3E w;

	w.q =   u.q/i;

	w.x   = u.x/i;
	w.y   = u.y/i;
	w.z   = u.z/i;

	w.xy  = u.xy/i;
	w.xz  = u.xz/i;
	w.yz  = u.yz/i;

	w.xyz = u.xyz/i;

	return w;
}

//////////////////////////////////////////////////////

GA3E operator^(const GA3E &a, const GA3E &b) {

	GA3E c;

/*

 ^    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     0     xy    xz    0     0     xyz   0    
 y    |  y    -xy    0     yz    0    -xyz   0     0    
 z    |  z    -xz   -yz    0     xyz   0     0     0    
 xy   |  xy    0     0     xyz   0     0     0     0    
 xz   |  xz    0    -xyz   0     0     0     0     0    
 yz   |  yz    xyz   0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/

	c.q   =  + a.q*b.q  ;
	c.x   =  + a.q*b.x   + a.x*b.q  ;
	c.y   =  + a.q*b.y   + a.y*b.q  ;
	c.z   =  + a.q*b.z   + a.z*b.q  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q  ;
	c.xz  =  + a.q*b.xz  + a.x*b.z  - a.z*b.x  + a.xz*b.q  ;
	c.yz  =  + a.q*b.yz  + a.y*b.z  - a.z*b.y  + a.yz*b.q  ;
	c.xyz =  + a.q*b.xyz + a.x*b.yz - a.y*b.xz + a.z*b.xy + a.xy*b.z - a.xz*b.y + a.yz*b.x + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E Wedge(const GA3E &a, const GA3E &b) {

	GA3E c;

/*

 ^    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     0     xy    xz    0     0     xyz   0    
 y    |  y    -xy    0     yz    0    -xyz   0     0    
 z    |  z    -xz   -yz    0     xyz   0     0     0    
 xy   |  xy    0     0     xyz   0     0     0     0    
 xz   |  xz    0    -xyz   0     0     0     0     0    
 yz   |  yz    xyz   0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/

	c.q   =  + a.q*b.q  ;
	c.x   =  + a.q*b.x   + a.x*b.q  ;
	c.y   =  + a.q*b.y   + a.y*b.q  ;
	c.z   =  + a.q*b.z   + a.z*b.q  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q  ;
	c.xz  =  + a.q*b.xz  + a.x*b.z  - a.z*b.x  + a.xz*b.q  ;
	c.yz  =  + a.q*b.yz  + a.y*b.z  - a.z*b.y  + a.yz*b.q  ;
	c.xyz =  + a.q*b.xyz + a.x*b.yz - a.y*b.xz + a.z*b.xy + a.xy*b.z - a.xz*b.y + a.yz*b.x + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;

}

//////////////////////////////////////////////////////

GA3E AntiWedge(const GA3E a, const GA3E b)
{
	GA3E c;
/*

Lengyel's AntiWedge

 V    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0    -q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0     x     y     xy   
 xz   |  0     0    -q     0    -x     0     z     xz   
 yz   |  0     q     0     0    -y    -z     0     yz   
 xyz  |  q     x     y     z     xy    xz    yz    xyz  

*/

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

// AntiWedge
	c.xyz =  + a.xyz*b.xyz ; 
	c.yz  =  + a.yz *b.xyz + a.xyz*b.yz  ; 
	c.xz  =  + a.xz *b.xyz + a.xyz*b.xz  ; 
	c.xy  =  + a.xy *b.xyz + a.xyz*b.xy  ; 
	c.z   =  + a.z  *b.xyz + a.xz *b.yz  - a.yz *b.xz  + a.xyz*b.z   ; 
	c.y   =  + a.y  *b.xyz + a.xy *b.yz  - a.yz *b.xy  + a.xyz*b.y   ; 
	c.x  =  + a.x  *b.xyz + a.xy *b.xz  - a.xz *b.xy  + a.xyz*b.x   ; 
	c.q   =  + a.q  *b.xyz + a.x  *b.yz  - a.y  *b.xz  + a.z  *b.xy
                 + a.xy *b.z   - a.xz *b.y   + a.yz *b.x   + a.xyz*b.q   ; 


	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////


GA3E Regressive(GA3E a, GA3E b) 
{
	GA3E c;	

/*

Hestenes' Regressive (differs from AntiWedge in signs)

 V    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0    -q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0    -x    -y     xy   
 xz   |  0     0    -q     0     x     0    -z     xz   
 yz   |  0     q     0     0     y     z     0     yz   
 xyz  |  q     x     y     z     xy    xz    yz    xyz  

*/

//	GA3E I, I_inv;
//	I.xyz = 1;
//	I_inv.xyz = -1;

//	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = a.q*b.xyz + a.x*b.yz + a.xy*b.z + a.xyz*b.q - a.xz*b.y - a.y*b.xz + a.yz*b.x + a.z*b.xy ;

	c.x = + ( a.x*b.xyz - a.xy*b.xz + a.xyz*b.x + a.xz*b.xy) ;
	c.y = + (-a.xy*b.yz + a.xyz*b.y + a.y*b.xyz + a.yz*b.xy) ;
	c.z = + ( a.xyz*b.z - a.xz*b.yz + a.yz*b.xz + a.z*b.xyz) ;

	c.xy = + (a.xy*b.xyz + a.xyz*b.xy) ;
	c.xz = + (a.xyz*b.xz + a.xz*b.xyz) ;
	c.yz = + (a.xyz*b.yz + a.yz*b.xyz) ;

	c.xyz = + a.xyz*b.xyz ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E RegressiveViaFormula(GA3E a, GA3E b)
{
	GA3E c;

	GA3E I,I_inv;
	
	I.xyz = 1;	I_inv.xyz = -1;

	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E LowerRightViaFormula(GA3E a, GA3E b)  // Like a mirror of Wedge along rising diagonal
{

/*


LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     xyz  
 x    |  0     0     0     0     0     0     xyz   yz   
 y    |  0     0     0     0     0    -xyz   0    -xz   
 z    |  0     0     0     0     xyz   0     0     xy   
 xy   |  0     0     0     xyz   0    -yz    xz   -z    
 xz   |  0     0    -xyz   0     yz    0    -xy    y    
 yz   |  0     xyz   0     0    -xz    xy    0    -x    
 xyz  |  xyz   yz   -xz    xy   -z     y    -x    -q    


*/
// Wedge(a*I_inv,I*b);

	GA3E c;

	GA3E I,I_inv,d,e,f;
	
	I.xyz = 1;	I_inv.xyz = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

GA3E Expander(const GA3E &a, const GA3E &b)
{ 
/*

Terms with increased rank
 >    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    0     0     xyz   0    
 y    |  0    -xy    0     yz    0    -xyz   0     0    
 z    |  0    -xz   -yz    0     xyz   0     0     0    
 xy   |  0     0     0     xyz   0     0     0     0    
 xz   |  0     0    -xyz   0     0     0     0     0    
 yz   |  0     xyz   0     0     0     0     0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3E c;

	c.q    =  0 ; 

	c.x    =  0 ; 
	c.y    =  0 ; 
	c.z    =  0 ; 

	c.xy   =  + a.x*b.y - a.y*b.x    ; 
	c.xz   =  + a.x*b.z - a.z*b.x    ; 
	c.yz   =  + a.y*b.z - a.z*b.y    ; 

	c.xyz  =  + a.x*b.yz - a.y*b.xz + a.z*b.xy + a.xy*b.z - a.xz*b.y + a.yz*b.x    ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E Conserver(const GA3E &a, const GA3E &b)
{  

/*

Terms with preserved rank
 =    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0    -yz    xz    0    
 xz   |  xz    0     0     0     yz    0    -xy    0    
 yz   |  yz    0     0     0    -xz    xy    0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/
	GA3E c;

	c.q    =  + a.q*b.q    ; 

	c.x    =  + a.q*b.x    + a.x*b.q    ; 
	c.y    =  + a.q*b.y    + a.y*b.q    ; 	
	c.z    =  + a.q*b.z    + a.z*b.q    ; 

	c.xy   =  + a.q*b.xy   + a.xy*b.q    - a.xz*b.yz   + a.yz*b.xz ; 
	c.xz   =  + a.q*b.xz   + a.xy*b.yz   + a.xz*b.q    - a.yz*b.xy ; 
	c.yz   =  + a.q*b.yz   - a.xy*b.xz   + a.xz*b.xy   + a.yz*b.q  ; 

	c.xyz  =  + a.q*b.xyz  + a.xyz *b.q  ; 


	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);



	return c;
}

//////////////////////////////////////////////////////

GA3E Shrinker(const GA3E &a, const GA3E &b)
{

/*

Terms with reduced rank
 <    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y     z     0     yz   
 y    |  0     0     q     0    -x     0     z    -xz   
 z    |  0     0     0     q     0    -x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 xz   |  0    -z     0     x     0    -q     0     y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz   -xz    xy   -z     y    -x    -q    

*/
	GA3E c;

// Shrinker equation set
	c.q    =  + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.xz*b.xz - a.yz*b.yz - a.xyz *b.xyz; 

	c.x    =  - a.y*b.xy - a.z*b.xz + a.xy*b.y + a.xz*b.z - a.yz*b.xyz - a.xyz*b.yz ; 
	c.y    =  + a.x*b.xy - a.z*b.yz - a.xy*b.x + a.xz*b.xyz + a.yz*b.z + a.xyz*b.xz ; 
	c.z    =  + a.x*b.xz + a.y*b.yz - a.xy*b.xyz - a.xz*b.x - a.yz*b.y - a.xyz*b.xy ; 

	c.xy   =  + a.z*b.xyz + a.xyz*b.z  ; 
	c.xz   =  - a.y*b.xyz - a.xyz*b.y ; 
	c.yz   =  + a.x*b.xyz + a.xyz*b.x ; 

	c.xyz  =  0  ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

///////////////////////////////////////////////////////

GA3E Symmetric(const GA3E &a, const GA3E &b)
{ 
/*

Symmetric Product
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  x     q     0     0     0     0     xyz   yz   
 y    |  y     0     q     0     0    -xyz   0    -xz   
 z    |  z     0     0     q     xyz   0     0     xy   
 xy   |  xy    0     0     xyz  -q     0     0    -z    
 xz   |  xz    0    -xyz   0     0    -q     0     y    
 yz   |  yz    xyz   0     0     0     0    -q    -x    
 xyz  |  xyz   yz   -xz    xy   -z     y    -x    -q    

*/
	GA3E c;

//	c = (a*b + b*a)/2;

	c.q = a.q*b.q + a.x*b.x - a.xy*b.xy - a.xyz*b.xyz - a.xz*b.xz + a.y*b.y - a.yz*b.yz + a.z*b.z ;

	c.x = + (a.q*b.x + a.x*b.q - a.xyz*b.yz - a.yz*b.xyz) ;
	c.y = + (a.q*b.y + a.xyz*b.xz + a.xz*b.xyz + a.y*b.q) ;
	c.z = + (a.q*b.z - a.xy*b.xyz - a.xyz*b.xy + a.z*b.q) ;

	c.xy = + (a.q*b.xy + a.xy*b.q + a.xyz*b.z + a.z*b.xyz) ;
	c.xz = + (a.q*b.xz - a.xyz*b.y + a.xz*b.q - a.y*b.xyz) ;
	c.yz = + (a.q*b.yz + a.x*b.xyz + a.xyz*b.x + a.yz*b.q) ;

	c.xyz = + (a.q*b.xyz + a.x*b.yz + a.xy*b.z + a.xyz*b.q - a.xz*b.y - a.y*b.xz + a.yz*b.x + a.z*b.xy) ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E AntiSymmetric(const GA3E &a, const GA3E &b)
{ 
/*

AntiSymmetric Product
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    y     z     0     0    
 y    |  0    -xy    0     yz   -x     0     z     0    
 z    |  0    -xz   -yz    0     0    -x    -y     0    
 xy   |  0    -y     x     0     0    -yz    xz    0    
 xz   |  0    -z     0     x     yz    0    -xy    0    
 yz   |  0     0    -z     y    -xz    xy    0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3E c;

//	c = (a*b - b*a)/2;

	c.q = 0;

	c.x = + (a.xy*b.y + a.xz*b.z - a.y*b.xy - a.z*b.xz) ;
	c.y = + (a.x*b.xy - a.xy*b.x + a.yz*b.z - a.z*b.yz) ;
	c.z = + (a.x*b.xz - a.xz*b.x + a.y*b.yz - a.yz*b.y) ;

	c.xy = + ( a.x*b.y - a.xz*b.yz - a.y*b.x + a.yz*b.xz) ;
	c.xz = + ( a.x*b.z + a.xy*b.yz - a.yz*b.xy - a.z*b.x) ;
	c.yz = + (-a.xy*b.xz + a.xz*b.xy + a.y*b.z - a.z*b.y) ;

	c.xyz = 0;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E Inner(const GA3E &a, const GA3E &b)
{ 
/*

Inner Product
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y     z     0     yz   
 y    |  0     0     q     0    -x     0     z    -xz   
 z    |  0     0     0     q     0    -x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 xz   |  0    -z     0     x     0    -q     0     y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz   -xz    xy   -z     y    -x    -q    

*/
	GA3E c;

	c.q = a.x*b.x - a.xy*b.xy - a.xyz*b.xyz - a.xz*b.xz + a.y*b.y - a.yz*b.yz + a.z*b.z ;

	c.x = + (a.xy*b.y - a.xyz*b.yz + a.xz*b.z - a.y*b.xy - a.yz*b.xyz - a.z*b.xz) ;
	c.y = + (a.x*b.xy - a.xy*b.x + a.xyz*b.xz + a.xz*b.xyz + a.yz*b.z - a.z*b.yz) ;
	c.z = + (a.x*b.xz - a.xy*b.xyz - a.xyz*b.xy - a.xz*b.x + a.y*b.yz - a.yz*b.y) ;

	c.xy = + ( a.xyz*b.z + a.z*b.xyz) ;
	c.xz = + (-a.xyz*b.y - a.y*b.xyz) ;
	c.yz = + ( a.x*b.xyz + a.xyz*b.x) ;

	c.xyz = 0;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E LeftContraction (const GA3E &a, const GA3E &b)
{ 
/*

LeftContraction Product
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    xz    yz    xyz  
 x    |  0     q     0     0     y     z     0     yz   
 y    |  0     0     q     0    -x     0     z    -xz   
 z    |  0     0     0     q     0    -x    -y     xy   
 xy   |  0     0     0     0    -q     0     0    -z    
 xz   |  0     0     0     0     0    -q     0     y    
 yz   |  0     0     0     0     0     0    -q    -x    
 xyz  |  0     0     0     0     0     0     0    -q    


*/
	GA3E c;

	c.q = a.q*b.q + a.x*b.x - a.xy*b.xy - a.xyz*b.xyz - a.xz*b.xz + a.y*b.y - a.yz*b.yz + a.z*b.z ;

	c.x = + (a.q*b.x - a.y*b.xy - a.yz*b.xyz - a.z*b.xz) ;
	c.y = + (a.q*b.y + a.x*b.xy + a.xz*b.xyz - a.z*b.yz) ;
	c.z = + (a.q*b.z + a.x*b.xz - a.xy*b.xyz + a.y*b.yz) ;

	c.xy = + (a.q*b.xy + a.z*b.xyz) ;
	c.xz = + (a.q*b.xz - a.y*b.xyz) ;
	c.yz = + (a.q*b.yz + a.x*b.xyz) ;

	c.xyz = + a.q*b.xyz ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3E RightContraction (const GA3E &a, const GA3E &b)
{ 
/*

RightContraction Product
 ?    |  q     x     y     z     xy    xz    yz    xyz  
-------------------------------------------------------
 q    |  q     0     0     0     0     0     0     0    
 x    |  x     q     0     0     0     0     0     0    
 y    |  y     0     q     0     0     0     0     0    
 z    |  z     0     0     q     0     0     0     0    
 xy   |  xy   -y     x     0    -q     0     0     0    
 xz   |  xz   -z     0     x     0    -q     0     0    
 yz   |  yz    0    -z     y     0     0    -q     0    
 xyz  |  xyz   yz   -xz    xy   -z     y    -x    -q    


*/
	GA3E c;

	c.q = a.q*b.q + a.x*b.x - a.xy*b.xy - a.xyz*b.xyz - a.xz*b.xz + a.y*b.y - a.yz*b.yz + a.z*b.z ;

	c.x = + ( a.x*b.q + a.xy*b.y - a.xyz*b.yz + a.xz*b.z) ;
	c.y = + (-a.xy*b.x + a.xyz*b.xz + a.y*b.q + a.yz*b.z) ;
	c.z = + (-a.xyz*b.xy - a.xz*b.x - a.yz*b.y + a.z*b.q) ;

	c.xy = + ( a.xy*b.q + a.xyz*b.z) ;
	c.xz = + (-a.xyz*b.y + a.xz*b.q) ;
	c.yz = + ( a.xyz*b.x + a.yz*b.q) ;

	c.xyz = + a.xyz*b.q ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

ex Determinant(GA3E A) {

	GA3E B, C;
	ex a, b, c, d, e, s, det;

	B = CliffordConjugation(A);
	C = Product(B,A);

	a = C.q;
	b = C.xyz;

	det = expand(a*a + b*b);

	return(det);
}

//////////////////////////////////////////////////////

GA3E Adjugate(GA3E a)
{

	GA3E u;

//Product(CliffordConjugation(a),Conjugation(Product(a,CliffordConjugation(a)))) = Adjugate(a);

	u = Product(CliffordConjugation(a),Conjugation(Product(a,CliffordConjugation(a))));

	u.q = expand(u.q);

	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);

	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);

	u.xyz = expand(u.xyz);

	return u;

}

//////////////////////////////////////////////////////

GA3E Reciprocal(GA3E a)
{
	ex b;
	GA3E c,d;

	b = Determinant(a);

	c = Adjugate(a);

	d.q = expand(c.q/b);

	d.x = expand(c.x/b);
	d.y = expand(c.y/b);
	d.z = expand(c.z/b);

	d.xy = expand(c.xy/b);
	d.xz = expand(c.xz/b);
	d.yz = expand(c.yz/b);

	d.xyz = expand(c.xyz/b);

	return d;
}

/////////////////////////////////// 
 
GA3E	Magic(GA3E u, int i) {

	GA3E MV;
	ex a, b,c,d, e,f,g, h;

	a = u.q;	b = u.x;	c = u.y;	d = u.z;
	e = u.xy;	f = u.xz;	g = u.yz;	h = u.xyz;


	if(i == 0) MV = GA3E( a, +b,+c,+d, +e,+f,+g, +h) ; 
	if(i == 1) MV = GA3E( a, +b,+d,+c, +f,+e,-g, -h) ; 
	if(i == 2) MV = GA3E( a, +c,+b,+d, -e,+g,+f, -h) ; 
	if(i == 3) MV = GA3E( a, +c,+d,+b, +g,-e,-f, +h) ; 
	if(i == 4) MV = GA3E( a, +d,+b,+c, -f,-g,+e, +h) ; 
	if(i == 5) MV = GA3E( a, +d,+c,+b, -g,-f,-e, -h) ; 
	
	return MV;
}

//////////////////////////////////////////////////////

GA3E	Comp(GA3E u, int i) {

	GA3E MV;
	ex a, b,c,d, e,f,g, h;
	int Mask, Sign[8];

	Mask = 16;
	if((Mask & i) == 0) Sign[4] = +1; else Sign[4] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[3] = +1; else Sign[3] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[2] = +1; else Sign[2] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[1] = +1; else Sign[1] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[0] = +1; else Sign[0] = -1;  Mask >>= 1;
	Sign[5] = Sign[2]*Sign[1]*Sign[0];
	Sign[6] = Sign[3]*Sign[1]*Sign[0];
	Sign[7] = Sign[4]*Sign[1]*Sign[0];

	a = Sign[4]*u.q;   b = Sign[3]*u.x;   c = Sign[2]*u.y;   d = Sign[1]*u.z;
	e = Sign[0]*u.xy;  f = Sign[5]*u.xz;  g = Sign[6]*u.yz;  h = Sign[7]*u.xyz;

	MV = GA3E(a, b,c,d, e,f,g, h) ; 
	
	return MV;
}


//////////////////////////////////////////////////////
/*
matrix GA3E_To_Matrix(GA3E MV) // newer ginac calling sequence
{
	matrix q, x,y,z, xy,yz,xz, xyz;
	matrix M;


// q = [ 1, 0]    x = [ 0, 1]     y = [ 1, 0]     xy = [ 0,-1]
//     [ 0, 1]        [ 1, 0]         [ 0,-1]          [ 1, 0]

// z = [ 0, I]   xz = [-I, 0]    yz = [ 0, I]    xyz = [ I, 0]
//     [-I, 0]        [ 0, I]         [ I, 0]          [ 0, I]

	Zero = matrix({{0, 0},{0, 0}}); 


	q =  matrix({{1, 0},{ 0, 1}}) ;
	x =  matrix({{0, 1},{ 1, 0}}) ;
	y =  matrix({{1, 0},{ 0,-1}}) ;
	z =  matrix({{0, I},{-I, 0}}) ;
	
	xy = matrix({{ 0,-1},{ 1, 0}}) ;
	xz = matrix({{-I, 0},{ 0, I}}) ;
	yz = matrix({{ 0, I},{ I, 0}}) ;

	xyz = matrix({{ I, 0},{ 0, I}}) ;

	M = Zero;
	M = M.add(q.mul_scalar(MV.q));

	M = M.add(x.mul_scalar(MV.x));
	M = M.add(y.mul_scalar(MV.y));
	M = M.add(z.mul_scalar(MV.z));

	M = M.add(xy.mul_scalar(MV.xy));
	M = M.add(xz.mul_scalar(MV.xz));
	M = M.add(yz.mul_scalar(MV.yz));

	M = M.add(xyz.mul_scalar(MV.xyz));

	return M;

}
*/
//////////////////////////////////////////////////////

matrix GA3E_To_Matrix(GA3E MV)
{
	matrix q(2,2), x(2,2),y(2,2),z(2,2), xy(2,2),yz(2,2),xz(2,2), xyz(2,2);
	matrix M(2,2);


// q = [ 1, 0]    x = [ 0, 1]     y = [ 1, 0]     xy = [ 0,-1]
//     [ 0, 1]        [ 1, 0]         [ 0,-1]          [ 1, 0]

// z = [ 0, I]   xz = [-I, 0]    yz = [ 0, I]    xyz = [ I, 0]
//     [-I, 0]        [ 0, I]         [ I, 0]          [ 0, I]

	q =  1, 0, 0, 1 ;
	x =  0, 1, 1, 0 ;
	y =  1, 0, 0,-1 ;
	z =  0, I,-I, 0 ;
	
	xy =  0,-1, 1, 0 ;
	xz = -I, 0, 0, I ;
	yz =  0, I, I, 0 ;

	xyz =  I, 0, 0, I ;

	M = 0, 0, 0, 0 ;

	M = M.add(q.mul_scalar(MV.q));

	M = M.add(x.mul_scalar(MV.x));
	M = M.add(y.mul_scalar(MV.y));
	M = M.add(z.mul_scalar(MV.z));

	M = M.add(xy.mul_scalar(MV.xy));
	M = M.add(xz.mul_scalar(MV.xz));
	M = M.add(yz.mul_scalar(MV.yz));

	M = M.add(xyz.mul_scalar(MV.xyz));

	return M;

}

//////////////////////////////////////////////////////

GA3E  Matrix_To_GA3E(matrix M)
{
	GA3E W;

	ex A, B,C,D, E,F,G, H;

	A = real_part( + M(0,0) + M(1,1) )/2;
	B = real_part( + M(0,1) + M(1,0) )/2;
	C = real_part( + M(0,0) - M(1,1) )/2;
	D = imag_part( + M(0,1) - M(1,0) )/2;
	E = real_part( - M(0,1) + M(1,0) )/2;
	F = imag_part( - M(0,0) + M(1,1) )/2;
	G = imag_part( + M(0,1) + M(1,0) )/2;
	H = imag_part( + M(0,0) + M(1,1) )/2;

	W = GA3E(A, B,C,D, E,F,G, H);

	return W;

}


