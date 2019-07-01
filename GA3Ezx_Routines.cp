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

struct GA3Ezx{	// 0 . 1 2 4 . 3 5 6 . 7 
	ex q,  x,y,z, xy,zx,yz,  xyz;
	GA3Ezx() {q = 0; x = 0; y = 0; z = 0;
		xy = 0; zx = 0; yz = 0;
		xyz = 0;}
	GA3Ezx(ex qq, ex xx, ex yy, ex zz, 
		ex xxyy, ex zzxx, ex yyzz,
		ex xxyyzz) 

		{q = qq; x = xx; y = yy; z = zz;
		 xy = xxyy; zx = zzxx; yz = yyzz;
		 xyz = xxyyzz;}
};


//////////////////////////////////////////////////////

// Necessary forward declarations

GA3Ezx LeftContraction (const GA3Ezx &a, const GA3Ezx &b) ;
GA3Ezx Product(const GA3Ezx &a, const GA3Ezx &b) ;

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, GA3Ezx &v) 
{
	return ff << "\n(" 
		<< v.q << ", \n" 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << ", \n" 
		<< v.xy << "," 
		<< v.zx << "," 
		<< v.yz << ", \n" 
		<< v.xyz << ")\n";
}


//////////////////////////////////////////////////////

void PrintMV(GA3Ezx &v) 
{
	cout 	<< "( "
		<< v.q << ", " 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << ", " 
		<< v.xy << "," 
		<< v.zx << "," 
		<< v.yz << ", " 
		<< v.xyz << ")";
}

//////////////////////////////////////////////////////

GA3Ezx OverBar(GA3Ezx a)	 // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3Ezx b;	
// OverBar(Blade[i]) = Blade[i]*xyz ; 

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     =  a.zx;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.zx   =   a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx UnderBar(GA3Ezx a)   // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3Ezx b;	
// UnderBar(Blade[i]) = xyz*Blade[i] ;

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     =  a.zx;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.zx   =   a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx Reverse(GA3Ezx w)	
// (A,  B, C, D,  -E,-F,-G, -H) 
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = w.x;
	v.y = w.y;
	v.z = w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Involution(GA3Ezx w)	
// Corresponds to P parity transform: x -> -x, y -> -y, z ->-z
// (A, -B,-C,-D,  E, F, G,  -H) 
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = w.xy;
	v.zx = w.zx;
	v.yz = w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Hermitian(GA3Ezx w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(A,  B, C, D, -E,-F,-G, -H)  // KN - Please verify
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Hermitian_Component(GA3Ezx w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(A,  B, C, D, 0, 0, 0,  0)   // KN - Please verify
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = 0;
	v.zx = 0;
	v.yz = 0;

	v.xyz = 0;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Anti_Hermitian_Component(GA3Ezx w)	// correct for both 2x2 complex and 4x4 real matrix implementations
//(0, 0, 0, 0,  E, F, G,  H)   // KN - Please verify
{
	GA3Ezx v;

	v.q =  0;

	v.x =  0;
	v.y =  0;
	v.z =  0;

	v.xy = w.xy;
	v.zx = w.zx;
	v.yz = w.yz;

	v.xyz = w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Transpose(GA3Ezx w)	// only correct for 4x4 real matrix implementation
//(A,  B, C, D, -E,-F,-G, -H)   // KN - Please verify
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx Conjugation(GA3Ezx w)
// (A, -B,-C,-D, -E,-F,-G, -H)
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx CliffordConjugation(GA3Ezx w)
// (A, -B,-C,-D, -E,-F,-G,  H) 
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = w.xyz;

	return v;
}

////////////////////////////////////////////////////// 

GA3Ezx Dual(GA3Ezx w)   // return w*I_inv 
// Dual(r) = u = (h, g, f, e, -d,-c,-b, -a)

{
	GA3Ezx v;
//	GA3Ezx I_inv;
//	I_inv.xyz = -1;
//	v = Product(w,I_inv);

	v.q =  w.xyz;

	v.x =  w.yz;
	v.y =  w.zx;
	v.z =  w.xy;

	v.xy = -w.z;
	v.zx = -w.y;
	v.yz = -w.x;

	v.xyz = -w.q;

	return v;
}



////////////////////////////////////////////////////// 

GA3Ezx DorstDual(GA3Ezx a)  //DorstDual(r) = u = (h, g, f,e, -d,-c,-b, -a)
{
	GA3Ezx b;
//	GA3Ezx I_inv;	I_inv.xyz = -1;

//	b = LeftContraction(a,I_inv); 

	b.q =  a.xyz;

	b.x =  a.yz;
	b.y =  a.zx;
	b.z =  a.xy;

	b.xy = -a.z;
	b.zx = -a.y;
	b.yz = -a.x;

	b.xyz = -a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx DorstUnDual(GA3Ezx a)  //DorstUnDual(r) = u = (-h,-g,-f,-e,  d, c, b,  a)
{
	GA3Ezx b;
//	GA3Ezx I;	I.xyz = 1;

//	b = LeftContraction(a,I);

	b.q = -a.xyz;

	b.x = -a.yz;
	b.y = -a.zx;
	b.z = -a.xy;

	b.xy =  a.z;
	b.zx =  a.y;
	b.yz =  a.x;

	b.xyz =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx operator+(const GA3Ezx &u, const GA3Ezx &v)
{
	GA3Ezx w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.xy  = u.xy  + v.xy ;
	w.zx  = u.zx  + v.zx ;
	w.yz  = u.yz  + v.yz ;
	w.xyz = u.xyz + v.xyz;

	return w;
}

//////////////////////////////////////////////////////

GA3Ezx operator-(const GA3Ezx &u, const GA3Ezx &v)
{
	GA3Ezx w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.xy  = u.xy  - v.xy ;
	w.zx  = u.zx  - v.zx ;
	w.yz  = u.yz  - v.yz ;
	w.xyz = u.xyz - v.xyz;

	return w;
}

//////////////////////////////////////////////////////

int operator==(const GA3Ezx &u, const GA3Ezx &v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&

			(u.xy==v.xy)&&
			(u.zx==v.zx)&&
			(u.yz==v.yz)&&

			(u.xyz==v.xyz);
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const GA3Ezx &u, const GA3Ezx &v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||

			(u.xy!=v.xy)||
			(u.zx!=v.zx)||
			(u.yz!=v.yz)||

			(u.xyz!=v.xyz);
	return result;
}

//////////////////////////////////////////////////////

GA3Ezx operator*(const GA3Ezx &a, const GA3Ezx &b) {
	GA3Ezx c;


/*


 *    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     q     xy   -zx    y    -z     xyz   yz   
 y    |  y    -xy    q     yz   -x     xyz   z     zx   
 z    |  z     zx   -yz    q     xyz   x    -y     xy   
 xy   |  xy   -y     x     xyz  -q     yz   -zx   -z    
 zx   |  zx    z     xyz  -x    -yz   -q     xy   -y    
 yz   |  yz    xyz  -z     y     zx   -xy   -q    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
	c.q   = a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.zx*b.zx  - a.yz*b.yz  - a.xyz*b.xyz;
	c.x   = a.q*b.x   + a.x*b.q   - a.y*b.xy  + a.z*b.zx  + a.xy*b.y   - a.zx*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   = a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   - a.zx*b.xyz + a.yz*b.z   - a.xyz*b.zx ;
	c.z   = a.q*b.z   - a.x*b.zx  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz + a.zx*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy  = a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   + a.zx*b.yz  - a.yz*b.zx  + a.xyz*b.z  ;
	c.zx  = a.q*b.zx  - a.x*b.z   + a.y*b.xyz + a.z*b.x   - a.xy*b.yz  + a.zx*b.q   + a.yz*b.xy  + a.xyz*b.y  ;
	c.yz  = a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   + a.xy*b.zx  - a.zx*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz = a.q*b.xyz + a.x*b.yz  + a.y*b.zx  + a.z*b.xy  + a.xy*b.z   + a.zx*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx Product(const GA3Ezx &a, const GA3Ezx &b)
{

/*  

 *    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     q     xy   -zx    y    -z     xyz   yz   
 y    |  y    -xy    q     yz   -x     xyz   z     zx   
 z    |  z     zx   -yz    q     xyz   x    -y     xy   
 xy   |  xy   -y     x     xyz  -q     yz   -zx   -z    
 zx   |  zx    z     xyz  -x    -yz   -q     xy   -y    
 yz   |  yz    xyz  -z     y     zx   -xy   -q    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

	c.q   = a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.zx*b.zx  - a.yz*b.yz  - a.xyz*b.xyz;
	c.x   = a.q*b.x   + a.x*b.q   - a.y*b.xy  + a.z*b.zx  + a.xy*b.y   - a.zx*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   = a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   - a.zx*b.xyz + a.yz*b.z   - a.xyz*b.zx ;
	c.z   = a.q*b.z   - a.x*b.zx  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz + a.zx*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy  = a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   + a.zx*b.yz  - a.yz*b.zx  + a.xyz*b.z  ;
	c.zx  = a.q*b.zx  - a.x*b.z   + a.y*b.xyz + a.z*b.x   - a.xy*b.yz  + a.zx*b.q   + a.yz*b.xy  + a.xyz*b.y  ;
	c.yz  = a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   + a.xy*b.zx  - a.zx*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz = a.q*b.xyz + a.x*b.yz  + a.y*b.zx  + a.z*b.xy  + a.xy*b.z   + a.zx*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx operator*(const GA3Ezx &v, const ex &u) {
	GA3Ezx w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.xy   = expand(u*v.xy  ) ;
	w.zx   = expand(u*v.zx  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xyz  = expand(u*v.xyz ) ;

	return w;
}

//////////////////////////////////////////////////////

GA3Ezx operator*(const ex &u, const GA3Ezx &v) {
	GA3Ezx w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.xy   = expand(u*v.xy  ) ;
	w.zx   = expand(u*v.zx  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xyz  = expand(u*v.xyz ) ;

	return w;
}


//////////////////////////////////////////////////////

GA3Ezx operator/(const GA3Ezx &u, const int i)
{
	GA3Ezx w;

	w.q =   u.q/i;

	w.x   = u.x/i;
	w.y   = u.y/i;
	w.z   = u.z/i;

	w.xy  = u.xy/i;
	w.zx  = u.zx/i;
	w.yz  = u.yz/i;

	w.xyz = u.xyz/i;

	return w;
}

//////////////////////////////////////////////////////

GA3Ezx operator^(const GA3Ezx &a, const GA3Ezx &b) {

	GA3Ezx c;

/*

 ^    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     0     xy   -zx    0     0     xyz   0    
 y    |  y    -xy    0     yz    0     xyz   0     0    
 z    |  z     zx   -yz    0     xyz   0     0     0    
 xy   |  xy    0     0     xyz   0     0     0     0    
 zx   |  zx    0     xyz   0     0     0     0     0    
 yz   |  yz    xyz   0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/

	c.q   =  a.q*b.q;  
	c.x   =  a.q*b.x   + a.x*b.q;  
	c.y   =  a.q*b.y   + a.y*b.q;  
	c.z   =  a.q*b.z   + a.z*b.q; 
	c.xy  =  a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q;  
	c.zx  =  a.q*b.zx  - a.x*b.z  + a.z*b.x  + a.zx*b.q;  
	c.yz  =  a.q*b.yz  + a.y*b.z  - a.z*b.y  + a.yz*b.q;  
	c.xyz =  a.q*b.xyz + a.x*b.yz + a.y*b.zx + a.xy*b.z + a.z*b.xy + a.zx*b.y + a.yz*b.x + a.xyz*b.q;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx Wedge(const GA3Ezx &a, const GA3Ezx &b) {

	GA3Ezx c;

/*

 ^    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     0     xy   -zx    0     0     xyz   0    
 y    |  y    -xy    0     yz    0     xyz   0     0    
 z    |  z     zx   -yz    0     xyz   0     0     0    
 xy   |  xy    0     0     xyz   0     0     0     0    
 zx   |  zx    0     xyz   0     0     0     0     0    
 yz   |  yz    xyz   0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/

	c.q   =  a.q*b.q;  
	c.x   =  a.q*b.x   + a.x*b.q;  
	c.y   =  a.q*b.y   + a.y*b.q;  
	c.z   =  a.q*b.z   + a.z*b.q; 
	c.xy  =  a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q;  
	c.zx  =  a.q*b.zx  - a.x*b.z  + a.z*b.x  + a.zx*b.q;  
	c.yz  =  a.q*b.yz  + a.y*b.z  - a.z*b.y  + a.yz*b.q;  
	c.xyz =  a.q*b.xyz + a.x*b.yz + a.y*b.zx + a.xy*b.z + a.z*b.xy + a.zx*b.y + a.yz*b.x + a.xyz*b.q;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;

}

//////////////////////////////////////////////////////

GA3Ezx AntiWedge(const GA3Ezx a, const GA3Ezx b) 
{
	GA3Ezx c;
/*

Lengyel's AntiWedge
 V    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0    -x     y     xy   
 zx   |  0     0     q     0     x     0    -z     zx   
 yz   |  0     q     0     0    -y     z     0     yz   
 xyz  |  q     x     y     z     xy    zx    yz    xyz  

AntiWedge Formula is Consistent with OverBar((UnderBar(r)^UnderBar(s)))

AntiWedge Formula is Consistent with UnderBar((OverBar(r)^OverBar(s)))


*/

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

// AntiWedge
	c.q   =  + a.q  *b.xyz + a.x  *b.yz  + a.y  *b.zx  + a.z  *b.xy
                 + a.xy *b.z   + a.zx *b.y   + a.yz *b.x   + a.xyz*b.q   ; 
	c.x   =  + a.x  *b.xyz - a.xy *b.zx  + a.zx *b.xy  + a.xyz*b.x   ; 
	c.y   =  + a.y  *b.xyz + a.xy *b.yz  - a.yz *b.xy  + a.xyz*b.y   ; 
	c.z   =  + a.z  *b.xyz - a.zx *b.yz  + a.yz *b.zx  + a.xyz*b.z   ; 
	c.xy  =  + a.xy *b.xyz + a.xyz*b.xy  ; 
	c.zx  =  + a.zx *b.xyz + a.xyz*b.zx  ; 
	c.yz  =  + a.yz *b.xyz + a.xyz*b.yz  ; 
	c.xyz =  + a.xyz*b.xyz ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////   

GA3Ezx Regressive(GA3Ezx a, GA3Ezx b) 
{
	GA3Ezx c;	

/*

Hestenes' Regressive (differs from AntiWedge in signs)

 V    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0     x    -y     xy   
 zx   |  0     0     q     0    -x     0     z     zx   
 yz   |  0     q     0     0     y    -z     0     yz   
 xyz  |  q     x     y     z     xy    zx    yz    xyz  

*/

//	GA3Ezx I, I_inv;
//	I.xyz = 1;
//	I_inv.xyz = -1;

//	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = a.q*b.xyz + a.x*b.yz + a.y*b.zx + a.z*b.xy + a.xy*b.z + a.zx*b.y + a.yz*b.x + a.xyz*b.q ;

	c.x = + ( a.x*b.xyz + a.xy*b.zx - a.zx*b.xy + a.xyz*b.x) ;
	c.y = + ( a.y*b.xyz - a.xy*b.yz + a.yz*b.xy + a.xyz*b.y) ;
	c.z = + ( a.z*b.xyz + a.zx*b.yz - a.yz*b.zx + a.xyz*b.z) ;

	c.xy = + ( a.xy*b.xyz + a.xyz*b.xy) ;
	c.zx = + ( a.zx*b.xyz + a.xyz*b.zx) ;
	c.yz = + ( a.yz*b.xyz + a.xyz*b.yz) ;

	c.xyz = + a.xyz*b.xyz ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx RegressiveViaFormula(GA3Ezx a, GA3Ezx b)
{
	GA3Ezx c;

	GA3Ezx I,I_inv;
	
	I.xyz = 1;	I_inv.xyz = -1;

	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx LowerRightViaFormula(GA3Ezx a, GA3Ezx b)  // Like a mirror of Wedge along rising diagonal
{

/*


LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     xyz  
 x    |  0     0     0     0     0     0     xyz   yz   
 y    |  0     0     0     0     0     xyz   0     zx   
 z    |  0     0     0     0     xyz   0     0     xy   
 xy   |  0     0     0     xyz   0     yz   -zx   -z    
 zx   |  0     0     xyz   0    -yz    0     xy   -y    
 yz   |  0     xyz   0     0     zx   -xy    0    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
// Wedge(a*I_inv,I*b);

	GA3Ezx c;

	GA3Ezx I,I_inv,d,e,f;
	
	I.xyz = 1;	I_inv.xyz = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx Expander(const GA3Ezx &a, const GA3Ezx &b)
{ 
/*

Terms with increased rank
 >    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy   -zx    0     0     xyz   0    
 y    |  0    -xy    0     yz    0     xyz   0     0    
 z    |  0     zx   -yz    0     xyz   0     0     0    
 xy   |  0     0     0     xyz   0     0     0     0    
 zx   |  0     0     xyz   0     0     0     0     0    
 yz   |  0     xyz   0     0     0     0     0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3Ezx c;

	c.q    =  0 ; 

	c.x    =  0 ; 
	c.y    =  0 ; 
	c.z    =  0 ; 

	c.xy   =  + a.x*b.y - a.y*b.x    ; 
	c.zx   =  - a.x*b.z + a.z*b.x    ; 
	c.yz   =  + a.y*b.z - a.z*b.y    ; 

	c.xyz  =  + a.x*b.yz + a.y*b.zx + a.z*b.xy + a.xy*b.z + a.zx*b.y + a.yz*b.x    ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx Conserver(const GA3Ezx &a, const GA3Ezx &b)
{  

/*

Terms with preserved rank
 =    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0     yz   -zx    0    
 zx   |  zx    0     0     0    -yz    0     xy    0    
 yz   |  yz    0     0     0     zx   -xy    0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/
	GA3Ezx c;

	c.q    =  + a.q*b.q    ; 

	c.x    =  + a.q*b.x    + a.x*b.q    ; 
	c.y    =  + a.q*b.y    + a.y*b.q    ; 	
	c.z    =  + a.q*b.z    + a.z*b.q    ; 

	c.xy   =  + a.q*b.xy   + a.xy*b.q    + a.zx*b.yz   - a.yz*b.zx ; 
	c.zx   =  + a.q*b.zx   - a.xy*b.yz   + a.zx*b.q    + a.yz*b.xy ; 
	c.yz   =  + a.q*b.yz   + a.xy*b.zx   - a.zx*b.xy   + a.yz*b.q  ; 

	c.xyz  =  + a.q*b.xyz  + a.xyz *b.q  ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);



	return c;
}

//////////////////////////////////////////////////////

GA3Ezx Shrinker(const GA3Ezx &a, const GA3Ezx &b)
{

/*

Terms with reduced rank
 <    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 zx   |  0     z     0    -x     0    -q     0    -y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

// Shrinker equation set
	c.q    =  + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz *b.xyz; 

	c.x    =  - a.y*b.xy + a.z*b.zx + a.xy*b.y - a.zx*b.z - a.yz*b.xyz - a.xyz*b.yz ; 
	c.y    =  + a.x*b.xy - a.z*b.yz - a.xy*b.x - a.zx*b.xyz + a.yz*b.z - a.xyz*b.zx ; 
	c.z    =  - a.x*b.zx + a.y*b.yz - a.xy*b.xyz + a.zx*b.x - a.yz*b.y - a.xyz*b.xy ; 

	c.xy   =  + a.z*b.xyz + a.xyz*b.z  ; 
	c.zx   =  + a.y*b.xyz + a.xyz*b.y ; 
	c.yz   =  + a.x*b.xyz + a.xyz*b.x ; 

	c.xyz  =  0  ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

/////////////////////////////////////////////////////// 

GA3Ezx Symmetric(const GA3Ezx &a, const GA3Ezx &b)
{ 
/*

Symmetric Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     q     0     0     0     0     xyz   yz   
 y    |  y     0     q     0     0     xyz   0     zx   
 z    |  z     0     0     q     xyz   0     0     xy   
 xy   |  xy    0     0     xyz  -q     0     0    -z    
 zx   |  zx    0     xyz   0     0    -q     0    -y    
 yz   |  yz    xyz   0     0     0     0    -q    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

//	c = (a*b + b*a)/2;

	c.q   =  + a.q*b.q   + a.x*b.x   + a.y *b.y   + a.z  *b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ; 
	c.x   =  + a.q*b.x   + a.x*b.q   - a.yz*b.xyz - a.xyz*b.yz  ; 
	c.y   =  + a.q*b.y   + a.y*b.q   - a.zx*b.xyz - a.xyz*b.zx  ; 
	c.z   =  + a.q*b.z   + a.z*b.q   - a.xy*b.xyz - a.xyz*b.xy  ; 
	c.xy  =  + a.q*b.xy  + a.z*b.xyz + a.xy*b.q   + a.xyz*b.z   ; 
	c.zx  =  + a.q*b.zx  + a.y*b.xyz + a.zx*b.q   + a.xyz*b.y   ; 
	c.yz  =  + a.q*b.yz  + a.x*b.xyz + a.yz*b.q   + a.xyz*b.x   ; 
	c.xyz =  + a.q*b.xyz + a.x*b.yz  + a.y *b.zx  + a.z*b.xy  + a.xy*b.z + a.zx*b.y + a.yz*b.x + a.xyz*b.q   ; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx AntiSymmetric(const GA3Ezx &a, const GA3Ezx &b)
{ 
/*

AntiSymmetric Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy   -zx    y    -z     0     0    
 y    |  0    -xy    0     yz   -x     0     z     0    
 z    |  0     zx   -yz    0     0     x    -y     0    
 xy   |  0    -y     x     0     0     yz   -zx    0    
 zx   |  0     z     0    -x    -yz    0     xy    0    
 yz   |  0     0    -z     y     zx   -xy    0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3Ezx c;

//	c = (a*b - b*a)/2;

	c.q   =  0; 
	c.x   =  - a.y*b.xy  + a.z*b.zx  + a.xy*b.y   - a.zx*b.z   ; 
	c.y   =  + a.x*b.xy  - a.z*b.yz  - a.xy*b.x   + a.yz*b.z   ; 
	c.z   =  - a.x*b.zx  + a.y*b.yz  + a.zx*b.x   - a.yz*b.y   ; 
	c.xy  =  + a.x*b.y   - a.y*b.x   + a.zx*b.yz  - a.yz*b.zx  ; 
	c.zx  =  - a.x*b.z   + a.z*b.x   - a.xy*b.yz  + a.yz*b.xy  ; 
	c.yz  =  + a.y*b.z   - a.z*b.y   + a.xy*b.zx  - a.zx*b.xy  ; 
	c.xyz =  0; 

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////// 

GA3Ezx Inner(const GA3Ezx &a, const GA3Ezx &b)  // Based Upon Bromborsky (A|B)
{ 
/*

Inner Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 zx   |  0     z     0    -x     0    -q     0    -y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz   -zx    xy   -z    -y    -x    -q    


*/
	GA3Ezx c;


	c.q =  + a.x*b.x  + a.y*b.y  + a.z*b.z  - a.xy*b.xy  - a.zx*b.zx  - a.yz*b.yz  - a.xyz*b.xyz ;
	c.x =             - a.y*b.xy + a.z*b.zx + a.xy*b.y   - a.zx*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y =  + a.x*b.xy            - a.z*b.yz - a.xy*b.x   - a.zx*b.xyz + a.yz*b.z   - a.xyz*b.zx ;
	c.z =  - a.x*b.zx + a.y*b.yz            - a.xy*b.xyz + a.zx*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy = + a.z*b.xyz + a.xyz*b.z  ;
	c.zx = + a.y*b.xyz + a.xyz*b.y  ;
	c.yz = + a.x*b.xyz + a.xyz*b.x ;
	c.xyz = 0;


	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx LeftContraction (const GA3Ezx &a, const GA3Ezx &b)	// Based Upon Bromborsky (A<B)
{ 
/*

LeftContraction Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0     0     0     0    -q     0     0    -z    
 zx   |  0     0     0     0     0    -q     0    -y    
 yz   |  0     0     0     0     0     0    -q    -x    
 xyz  |  0     0     0     0     0     0     0    -q    

*/
	GA3Ezx c;

	c.q = a.q*b.q + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ;

	c.x = + (a.q*b.x - a.y*b.xy + a.z*b.zx - a.yz*b.xyz) ;
	c.y = + (a.q*b.y + a.x*b.xy - a.z*b.yz - a.zx*b.xyz) ;
	c.z = + (a.q*b.z - a.x*b.zx + a.y*b.yz - a.xy*b.xyz) ;

	c.xy = + (a.q*b.xy + a.z*b.xyz) ;
	c.zx = + (a.q*b.zx + a.y*b.xyz) ;
	c.yz = + (a.q*b.yz + a.x*b.xyz) ;

	c.xyz = + a.q*b.xyz ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx RightContraction (const GA3Ezx &a, const GA3Ezx &b) // Based Upon Bromborsky (A>B)
{ 
/*

RightContraction Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     0     0     0     0     0     0     0    
 x    |  x     q     0     0     0     0     0     0    
 y    |  y     0     q     0     0     0     0     0    
 z    |  z     0     0     q     0     0     0     0    
 xy   |  xy   -y     x     0    -q     0     0     0    
 zx   |  zx    z     0    -x     0    -q     0     0    
 yz   |  yz    0    -z     y     0     0    -q     0    
 xyz  |  xyz   yz   -zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

	c.q = a.q*b.q + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ;

	c.x = + ( a.x*b.q + a.xy*b.y - a.zx*b.z - a.xyz*b.yz) ;
	c.y = + ( a.y*b.q - a.xy*b.x + a.yz*b.z - a.xyz*b.zx) ;
	c.z = + ( a.z*b.q + a.zx*b.x - a.yz*b.y - a.xyz*b.xy) ;

	c.xy =  ( a.xy*b.q + a.xyz*b.z) ;
	c.zx =  ( a.zx*b.q - a.xyz*b.y) ; 
	c.yz =  ( a.yz*b.q + a.xyz*b.x) ;

	c.xyz = + a.xyz*b.q ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.zx = expand(c.zx);
	c.yz = expand(c.yz);

	c.xyz = expand(c.xyz);

	return c;
}

//////////////////////////////////////////////////////

ex Determinant(GA3Ezx A) {

	GA3Ezx B, C;
	ex a, b, c, d, e, s, det;

	B = CliffordConjugation(A);
	C = Product(B,A);

	a = C.q;
	b = C.xyz;

	det = expand(a*a + b*b);

	return(det);
}

//////////////////////////////////////////////////////

GA3Ezx Adjugate(GA3Ezx a)
{

	GA3Ezx u;

//Product(CliffordConjugation(a),Conjugation(Product(a,CliffordConjugation(a)))) = Adjugate(a);

	u = Product(CliffordConjugation(a),Conjugation(Product(a,CliffordConjugation(a))));

	u.q = expand(u.q);

	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);

	u.xy = expand(u.xy);
	u.zx = expand(u.zx);
	u.yz = expand(u.yz);

	u.xyz = expand(u.xyz);

	return u;

}

//////////////////////////////////////////////////////

GA3Ezx Reciprocal(GA3Ezx a)
{
	ex b;
	GA3Ezx c,d;

	b = Determinant(a);

	c = Adjugate(a);

	d.q = expand(c.q/b);

	d.x = expand(c.x/b);
	d.y = expand(c.y/b);
	d.z = expand(c.z/b);

	d.xy = expand(c.xy/b);
	d.zx = expand(c.zx/b);
	d.yz = expand(c.yz/b);

	d.xyz = expand(c.xyz/b);

	return d;
}

/////////////////////////////////// 

GA3Ezx	Magic(GA3Ezx u, int i) { 

	GA3Ezx MV;
	ex a, b,c,d, e,f,g, h;

	a = u.q;	b = u.x;	c = u.y;	d = u.z;
	e = u.xy;	f = u.zx;	g = u.yz;	h = u.xyz;

	if(i == 0) MV = GA3Ezx( a, +b,+c,+d, +e,+f,+g, +h) ;
	if(i == 1) MV = GA3Ezx( a, +b,+d,+c, -f,-e,-g, -h) ;
	if(i == 2) MV = GA3Ezx( a, +c,+b,+d, -e,-g,-f, -h) ; 
	if(i == 3) MV = GA3Ezx( a, +c,+d,+b, +g,+e,+f, +h) ; 
	if(i == 4) MV = GA3Ezx( a, +d,+b,+c, +f,+g,+e, +h) ; 
	if(i == 5) MV = GA3Ezx( a, +d,+c,+b, -g,-f,-e, -h) ;
	
	return MV;
}

//////////////////////////////////////////////////////

GA3Ezx	Comp(GA3Ezx u, int i) {

	GA3Ezx MV;
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
	e = Sign[0]*u.xy;  f = Sign[5]*u.zx;  g = Sign[6]*u.yz;  h = Sign[7]*u.xyz;

	MV = GA3Ezx(a, b,c,d, e,f,g, h) ; 
	
	return MV;
}


////////////////////////////////////////////////////// 
/*
matrix GA3Ezx_To_Matrix(GA3Ezx MV) // newer ginac calling sequence
{
	matrix q, x,y,z, xy,zx,yz, xyz;
	matrix M;


// q = [ 1, 0]    x = [ 0, 1]     y = [ 1, 0]     xy = [ 0,-1]
//     [ 0, 1]        [ 1, 0]         [ 0,-1]          [ 1, 0]

// z = [ 0, I]   zx = [ I, 0]    yz = [ 0, I]    xyz = [ I, 0]
//     [-I, 0]        [ 0,-I]         [ I, 0]          [ 0, I]

	Zero = matrix({{0, 0},{0, 0}}); 


	q =  matrix({{1, 0},{ 0, 1}}) ;
	x =  matrix({{0, 1},{ 1, 0}}) ;
	y =  matrix({{1, 0},{ 0,-1}}) ;
	z =  matrix({{0, I},{-I, 0}}) ;
	
	xy = matrix({{ 0,-1},{ 1, 0}}) ;
	zx = matrix({{ I, 0},{ 0,-I}}) ;
	yz = matrix({{ 0, I},{ I, 0}}) ;

	xyz = matrix({{ I, 0},{ 0, I}}) ;

	M = Zero;
	M = M.add(q.mul_scalar(MV.q));

	M = M.add(x.mul_scalar(MV.x));
	M = M.add(y.mul_scalar(MV.y));
	M = M.add(z.mul_scalar(MV.z));

	M = M.add(xy.mul_scalar(MV.xy));
	M = M.add(xz.mul_scalar(MV.zx));
	M = M.add(yz.mul_scalar(MV.yz));

	M = M.add(xyz.mul_scalar(MV.xyz));

	return M;

}
*/
//////////////////////////////////////////////////////

matrix GA3Ezx_To_Matrix(GA3Ezx MV)
{
	matrix q(2,2), x(2,2),y(2,2),z(2,2), xy(2,2),zx(2,2),yz(2,2), xyz(2,2);
	matrix M(2,2);


// q = [ 1, 0]    x = [ 0, 1]     y = [ 1, 0]     xy = [ 0,-1]
//     [ 0, 1]        [ 1, 0]         [ 0,-1]          [ 1, 0]

// z = [ 0, I]   zx = [ I, 0]    yz = [ 0, I]    xyz = [ I, 0]
//     [-I, 0]        [ 0,-I]         [ I, 0]          [ 0, I]

	q =  1, 0, 0, 1 ;
	x =  0, 1, 1, 0 ;
	y =  1, 0, 0,-1 ;
	z =  0, I,-I, 0 ;
	
	xy =  0,-1, 1, 0 ;
	zx =  I, 0, 0,-I ;
	yz =  0, I, I, 0 ;

	xyz =  I, 0, 0, I ;

	M = 0, 0, 0, 0 ;

	M = M.add(q.mul_scalar(MV.q));

	M = M.add(x.mul_scalar(MV.x));
	M = M.add(y.mul_scalar(MV.y));
	M = M.add(z.mul_scalar(MV.z));

	M = M.add(xy.mul_scalar(MV.xy));
	M = M.add(zx.mul_scalar(MV.zx));
	M = M.add(yz.mul_scalar(MV.yz));

	M = M.add(xyz.mul_scalar(MV.xyz));

	return M;

}

//////////////////////////////////////////////////////

GA3Ezx  Matrix_To_GA3Ezx(matrix M)
{
	GA3Ezx W;

	ex A, B,C,D, E,F,G, H;

	A = real_part( + M(0,0) + M(1,1) )/2;
	B = real_part( + M(0,1) + M(1,0) )/2;
	C = real_part( + M(0,0) - M(1,1) )/2;
	D = imag_part( + M(0,1) - M(1,0) )/2;
	E = real_part( - M(0,1) + M(1,0) )/2;
	F = imag_part( + M(0,0) - M(1,1) )/2;
	G = imag_part( + M(0,1) + M(1,0) )/2;
	H = imag_part( + M(0,0) + M(1,1) )/2;

	W = GA3Ezx(A, B,C,D, E,F,G, H);

	return W;

}

/////////////////////////////////////////////////////
/*
int main(void) {

	return 0;
}
*/


