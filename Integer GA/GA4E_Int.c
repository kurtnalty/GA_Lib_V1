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

//////////////////////////////////////////////////////

typedef struct
{
	long q;
	long x,  y,  z, t;
	long xy, xz, yz, xt, yt, zt;
	long xyz, xyt, xzt, yzt;
	long xyzt;
} GA4E;


//////////////////////////////////////////////////////

GA4E Zero(void)		// initializer
{
	GA4E a;

	a.q = 0;
	a.x = 0; a.y = 0; a.z = 0; a.t = 0;
	a.xy = 0; a.xz = 0; a.yz = 0; a.xt = 0; a.yt = 0; a.zt = 0;
	a.xyz = 0; a.xyt = 0; a.xzt = 0; a.yzt = 0;
	a.xyzt = 0;

	return a;

}

//////////////////////////////////////////////////////

GA4E Set(	long q,
	long x,  long y,  long z, long t,
	long xy, long xz, long yz, long xt, long yt, long zt,
	long xyz, long xyt, long xzt, long yzt,
	long xyzt)		// initializer
{
	GA4E a;

	a.q = q;
	a.x = x; a.y = y; a.z = z; a.t = t;
	a.xy = xy; a.xz = xz; a.yz = yz; a.xt = xt; a.yt = yt; a.zt = zt;
	a.xyz = xyz; a.xyt = xyt; a.xzt = xzt;  a.yzt = yzt;
	a.xyzt = xyzt;

	return a;

}

//////////////////////////////////////////////////////

// Necessary forward declarations

GA4E LeftContraction (GA4E a, GA4E b) ;
GA4E Product(GA4E a, GA4E b) ;

//////////////////////////////////////////////////////

void PrintlnMV(GA4E v) 
{
	printf("(%10ld,\n  %10ld,%10ld,%10ld,%10ld,\n  %10ld,%10ld,%10ld, %10ld,%10ld,%10ld,\n  %10ld,%10ld,%10ld,%10ld,\n %10ld) \n",
		v.q, v.x,v.y,v.z,v.t, v.xy,v.xz,v.yz,v.xt,v.yt,v.zt,  v.xyz,v.xyt,v.xzt,v.yzt,  v.xyzt);
}

//////////////////////////////////////////////////////

void PrintMV(GA4E v) 
{
	printf("(%10ld,  %10ld,%10ld,%10ld,%10ld,  %10ld,%10ld,%10ld, %10ld,%10ld,%10ld,  %10ld,%10ld,%10ld,%10ld, %10ld) ",
		v.q, v.x,v.y,v.z,v.t, v.xy,v.xz,v.yz,v.xt,v.yt,v.zt,  v.xyz,v.xyt,v.xzt,v.yzt,  v.xyzt);
}

//////////////////////////////////////////////////////

GA4E OverBar(GA4E a)
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA4E b;		// a blade wedge b blade = pseudovector blade

	b.q     =  a.xyzt;	// xyzt

	b.x     = -a.yzt;
	b.y     =  a.xzt;
	b.z     = -a.xyt;
	b.t     =  a.xyz;

	b.xy   =   a.zt;
	b.xz   =  -a.yt;
	b.yz   =   a.xt;
	b.xt   =   a.yz;
	b.yt   =  -a.xz;
	b.zt   =   a.xy;

	b.xyz  = -a.t;
	b.xyt  =  a.z;
	b.xzt  = -a.y;
	b.yzt  =  a.x;

	b.xyzt = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA4E UnderBar(GA4E a)   
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA4E b;		// b blade wedge a blade = pseudovector blade

	b.q     =  a.xyzt;	// xyzt

	b.x     =  a.yzt;
	b.y     = -a.xzt;
	b.z     =  a.xyt;
	b.t     = -a.xyz;

	b.xy   =   a.zt;
	b.xz   =  -a.yt;
	b.yz   =   a.xt;
	b.xt   =   a.yz;
	b.yt   =  -a.xz;
	b.zt   =   a.xy;

	b.xyz  =  a.t;
	b.xyt  = -a.z;
	b.xzt  =  a.y;
	b.yzt  = -a.x;

	b.xyzt = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA4E Reverse(GA4E w)	
// (A,  B, C, D, E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R,  S) score = 10 N( 2046) Reverse
{
	GA4E v;
	v.q =  w.q;

	v.x = w.x;
	v.y = w.y;
	v.z = w.z;
	v.t = w.t;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;
	v.xt = -w.xt;
	v.yt = -w.yt;
	v.zt = -w.zt;

	v.xyz = -w.xyz;
	v.xyt = -w.xyt;
	v.xzt = -w.xzt;
	v.yzt = -w.yzt;

	v.xyzt =  w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

GA4E Involution(GA4E w)	
// Corresponds to PT parity transform: x -> -x, y -> -y, z ->-z, t -> -t
// (A, -B,-C,-D,-E,  F, G, H, J, K, L, -M,-N,-P,-R,  S) score =  0 N(30750) Involution
{
	GA4E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;
	v.t = -w.t;

	v.xy = w.xy;
	v.xz = w.xz;
	v.yz = w.yz;
	v.xt = w.xt;
	v.yt = w.yt;
	v.zt = w.zt;

	v.xyz = -w.xyz;
	v.xyt = -w.xyt;
	v.xzt = -w.xzt;
	v.yzt = -w.yzt;

	v.xyzt =  w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

GA4E Transpose(GA4E w)	// KN - Check this
//(A,  B, C, D,-E, -F,-G,-H, J, K, L, -M, N, P, R, -S) score =  6 N( 3857) Transpose  
{
	GA4E v;
/*	v.q =  w.q;   Minkowski

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;
	v.t = -w.t;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;
	v.xt =  w.xt;
	v.yt =  w.yt;
	v.zt =  w.zt;

	v.xyz = -w.xyz;
	v.xyt =  w.xyt;
	v.xzt =  w.xzt;
	v.yzt =  w.yzt;

	v.xyzt = -w.xyzt;
*/
	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;
	v.t =  w.t;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;
	v.xt = -w.xt;
	v.yt = -w.yt;
	v.zt = -w.zt;

	v.xyz = -w.xyz;
	v.xyt = -w.xyt;
	v.xzt = -w.xzt;
	v.yzt = -w.yzt;

	v.xyzt =  w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

GA4E Conjugation(GA4E w)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R, -S)
{
	GA4E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;
	v.t = -w.t;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;
	v.xt = -w.xt;
	v.yt = -w.yt;
	v.zt = -w.zt;

	v.xyz = -w.xyz;
	v.xyt = -w.xyt;
	v.xzt = -w.xzt;
	v.yzt = -w.yzt;

	v.xyzt =  -w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

GA4E CliffordConjugation(GA4E w)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L,  M, N, P, R,  S) score = 10 N(32736) Clifford Conjugation
{
	GA4E v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;
	v.t = -w.t;

	v.xy = -w.xy;
	v.xz = -w.xz;
	v.yz = -w.yz;
	v.xt = -w.xt;
	v.yt = -w.yt;
	v.zt = -w.zt;

	v.xyz = w.xyz;
	v.xyt = w.xyt;
	v.xzt = w.xzt;
	v.yzt = w.yzt;

	v.xyzt =  w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

GA4E Dual(GA4E w)   // return w*I_inv 
// Dual(r) = u = (p, o,-n,m,-l, -k,j,-i,-h,g,-f, -e,d,-c,b, a)

{
	GA4E v;
//	GA4E I_inv;
//	I_inv.xyzt = 1;

	v.q =  w.xyzt;

	v.x =  w.yzt;
	v.y = -w.xzt;
	v.z =  w.xyt;
	v.t = -w.xyz;

	v.xy = -w.zt;
	v.xz =  w.yt;
	v.yz = -w.xt;
	v.xt = -w.yz;
	v.yt =  w.xz;
	v.zt = -w.xy;

	v.xyz = -w.t;
	v.xyt =  w.z;
	v.xzt = -w.y;
	v.yzt =  w.x;

	v.xyzt =  w.q;

	return v;
}

//////////////////////////////////////////////////////

GA4E DorstDual(GA4E a)
{
	GA4E b, I_inv;

	I_inv = Zero();
	I_inv.xyzt = 1;

	b = LeftContraction(a,I_inv);
	return b;
}

//////////////////////////////////////////////////////

GA4E DorstUnDual(GA4E a)
{
	GA4E b, I;

	I = Zero();
	I.xyzt = 1;

	b = LeftContraction(a,I);
	return b;
}

//////////////////////////////////////////////////////

GA4E Add(GA4E u, GA4E v)
{
	GA4E w;

	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.t   = u.t   + v.t  ;
	w.xy  = u.xy  + v.xy ;
	w.xz  = u.xz  + v.xz ;
	w.yz  = u.yz  + v.yz ;
	w.xt  = u.xt  + v.xt ;
	w.yt  = u.yt  + v.yt ;
	w.zt  = u.zt  + v.zt ;
	w.xyz = u.xyz + v.xyz;
	w.xyt = u.xyt + v.xyt;
	w.xzt = u.xzt + v.xzt;
	w.yzt = u.yzt + v.yzt;
	w.xyzt = u.xyzt + v.xyzt;

	return w;
}

//////////////////////////////////////////////////////

GA4E Subtract(GA4E u, GA4E v)
{
	GA4E w;

	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.t   = u.t   - v.t  ;
	w.xy  = u.xy  - v.xy ;
	w.xz  = u.xz  - v.xz ;
	w.yz  = u.yz  - v.yz ;
	w.xt  = u.xt  - v.xt ;
	w.yt  = u.yt  - v.yt ;
	w.zt  = u.zt  - v.zt ;
	w.xyz = u.xyz - v.xyz;
	w.xyt = u.xyt - v.xyt;
	w.xzt = u.xzt - v.xzt;
	w.yzt = u.yzt - v.yzt;
	w.xyzt = u.xyzt - v.xyzt;

	return w;
}

//////////////////////////////////////////////////////

int Equal(GA4E u, GA4E v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&
			(u.t ==v.t )&&

			(u.xy==v.xy)&&
			(u.xz==v.xz)&&
			(u.yz==v.yz)&&
			(u.xt==v.xt)&&
			(u.yt==v.yt)&&
			(u.zt==v.zt)&&

			(u.xyz==v.xyz)&&
			(u.xyt==v.xyt)&&
			(u.xzt==v.xzt)&&
			(u.yzt==v.yzt)&&

			(u.xyzt==v.xyzt);
	return result;
}

//////////////////////////////////////////////////////

int Not_Equal(GA4E u, GA4E v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||
			(u.t !=v.t )||

			(u.xy!=v.xy)||
			(u.xz!=v.xz)||
			(u.yz!=v.yz)||
			(u.xt!=v.xt)||
			(u.yt!=v.yt)||
			(u.zt!=v.zt)||

			(u.xyz!=v.xyz)||
			(u.xyt!=v.xyt)||
			(u.xzt!=v.xzt)||
			(u.yzt!=v.yzt)||

			(u.xyzt!=v.xyzt);
	return result;
}

//////////////////////////////////////////////////////

GA4E Product(GA4E a, GA4E b)
{

	GA4E c;
/*
 *    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     q     xy    xz    xt    y     z     xyz   t     xyt   xzt   yz    yt    zt    xyzt  yzt  
 y    |  y    -xy    q     yz    yt   -x    -xyz   z    -xyt   t     yzt  -xz   -xt   -xyzt  zt   -xzt  
 z    |  z    -xz   -yz    q     zt    xyz  -x    -y    -xzt  -yzt   t     xy    xyzt -xt   -yt    xyt  
 t    |  t    -xt   -yt   -zt    q     xyt   xzt   yzt  -x    -y    -z    -xyzt  xy    xz    yz   -xyz  
 xy   |  xy   -y     x     xyz   xyt  -q    -yz    xz   -yt    xt    xyzt -z    -t    -yzt   xzt  -zt   
 xz   |  xz   -z    -xyz   x     xzt   yz   -q    -xy   -zt   -xyzt  xt    y     yzt  -t    -xyt   yt   
 yz   |  yz    xyz  -z     y     yzt  -xz    xy   -q     xyzt -zt    yt   -x    -xzt   xyt  -t    -xt   
 xt   |  xt   -t    -xyt  -xzt   x     yt    zt    xyzt -q    -xy   -xz   -yzt   y     z     xyz  -yz   
 yt   |  yt    xyt  -t    -yzt   y    -xt   -xyzt  zt    xy   -q    -yz    xzt  -x    -xyz   z     xz   
 zt   |  zt    xzt   yzt  -t     z     xyzt -xt   -yt    xz    yz   -q    -xyt   xyz  -x    -y    -xy   
 xyz  |  xyz   yz   -xz    xy    xyzt -z     y    -x     yzt  -xzt   xyt  -q    -zt    yt   -xt   -t    
 xyt  |  xyt   yt   -xt   -xyzt  xy   -t    -yzt   xzt   y    -x    -xyz   zt   -q    -yz    xz    z    
 xzt  |  xzt   zt    xyzt -xt    xz    yzt  -t    -xyt   z     xyz  -x    -yt    yz   -q    -xy   -y    
 yzt  |  yzt  -xyzt  zt   -yt    yz   -xzt   xyt  -t    -xyz   z    -y     xt   -xz    xy   -q     x    
 xyzt |  xyzt -yzt   xzt  -xyt   xyz  -zt    yt   -xt   -yz    xz   -xy    t    -z     y    -x     q    

*/

c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xt  *b.xt   - a.yt  *b.yt   - a.zt  *b.zt   - a.xyz *b.xyz  - a.xyt *b.xyt  - a.xzt *b.xzt  - a.yzt *b.yzt  + a.xyzt*b.xyzt;
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.y   *b.xy   - a.z   *b.xz   - a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  + a.xt  *b.t    - a.yt  *b.xyt  - a.zt  *b.xzt  - a.xyz *b.yz   - a.xyt *b.yt   - a.xzt *b.zt   + a.yzt *b.xyzt - a.xyzt*b.yzt ;
c.y    =  + a.q   *b.y    + a.x   *b.xy   + a.y   *b.q    - a.z   *b.yz   - a.t   *b.yt   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    + a.xt  *b.xyt  + a.yt  *b.t    - a.zt  *b.yzt  + a.xyz *b.xz   + a.xyt *b.xt   - a.xzt *b.xyzt - a.yzt *b.zt   + a.xyzt*b.xzt ;
c.z    =  + a.q   *b.z    + a.x   *b.xz   + a.y   *b.yz   + a.z   *b.q    - a.t   *b.zt   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    + a.xt  *b.xzt  + a.yt  *b.yzt  + a.zt  *b.t    - a.xyz *b.xy   + a.xyt *b.xyzt + a.xzt *b.xt   + a.yzt *b.yt   - a.xyzt*b.xyt ;
c.t    =  + a.q   *b.t    + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   + a.t   *b.q    - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   + a.xyzt*b.xyz ;
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.z   *b.xyz  + a.t   *b.xyt  + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   - a.xt  *b.yt   + a.yt  *b.xt   - a.zt  *b.xyzt + a.xyz *b.z    + a.xyt *b.t    - a.xzt *b.yzt  + a.yzt *b.xzt  - a.xyzt*b.zt  ;
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.y   *b.xyz  - a.z   *b.x    + a.t   *b.xzt  + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   - a.xt  *b.zt   + a.yt  *b.xyzt + a.zt  *b.xt   - a.xyz *b.y    + a.xyt *b.yzt  + a.xzt *b.t    - a.yzt *b.xyt  + a.xyzt*b.yt  ;
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.y   *b.z    - a.z   *b.y    + a.t   *b.yzt  - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    - a.xt  *b.xyzt - a.yt  *b.zt   + a.zt  *b.yt   + a.xyz *b.x    - a.xyt *b.xzt  + a.xzt *b.xyt  + a.yzt *b.t    - a.xyzt*b.xt  ;
c.xt   =  + a.q   *b.xt   + a.x   *b.t    - a.y   *b.xyt  - a.z   *b.xzt  - a.t   *b.x    + a.xy  *b.yt   + a.xz  *b.zt   - a.yz  *b.xyzt + a.xt  *b.q    - a.yt  *b.xy   - a.zt  *b.xz   - a.xyz *b.yzt  - a.xyt *b.y    - a.xzt *b.z    + a.yzt *b.xyz  - a.xyzt*b.yz  ;
c.yt   =  + a.q   *b.yt   + a.x   *b.xyt  + a.y   *b.t    - a.z   *b.yzt  - a.t   *b.y    - a.xy  *b.xt   + a.xz  *b.xyzt + a.yz  *b.zt   + a.xt  *b.xy   + a.yt  *b.q    - a.zt  *b.yz   + a.xyz *b.xzt  + a.xyt *b.x    - a.xzt *b.xyz  - a.yzt *b.z    + a.xyzt*b.xz  ;
c.zt   =  + a.q   *b.zt   + a.x   *b.xzt  + a.y   *b.yzt  + a.z   *b.t    - a.t   *b.z    - a.xy  *b.xyzt - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   + a.zt  *b.q    - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy  ;
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   - a.t   *b.xyzt + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xt  *b.yzt  - a.yt  *b.xzt  + a.zt  *b.xyt  + a.xyz *b.q    - a.xyt *b.zt   + a.xzt *b.yt   - a.yzt *b.xt   + a.xyzt*b.t   ;
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.z   *b.xyzt + a.t   *b.xy   + a.xy  *b.t    - a.xz  *b.yzt  + a.yz  *b.xzt  - a.xt  *b.y    + a.yt  *b.x    - a.zt  *b.xyz  + a.xyz *b.zt   + a.xyt *b.q    - a.xzt *b.yz   + a.yzt *b.xz   - a.xyzt*b.z   ;
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.y   *b.xyzt - a.z   *b.xt   + a.t   *b.xz   + a.xy  *b.yzt  + a.xz  *b.t    - a.yz  *b.xyt  - a.xt  *b.z    + a.yt  *b.xyz  + a.zt  *b.x    - a.xyz *b.yt   + a.xyt *b.yz   + a.xzt *b.q    - a.yzt *b.xy   + a.xyzt*b.y   ;
c.yzt  =  + a.q   *b.yzt  + a.x   *b.xyzt + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   - a.xy  *b.xzt  + a.xz  *b.xyt  + a.yz  *b.t    - a.xt  *b.xyz  - a.yt  *b.z    + a.zt  *b.y    + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.yzt *b.q    - a.xyzt*b.x   ;
c.xyzt =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q   ;

	return c;
}

//////////////////////////////////////////////////////

GA4E Divide_By_Constant(GA4E u, long a)
{
	GA4E w;
	w.q = u.q/a;

	w.x = u.x/a;
	w.y = u.y/a;
	w.z = u.z/a;
	w.t = u.t/a;

	w.xy = u.xy/a;
	w.xz = u.xz/a;
	w.yz = u.yz/a;
	w.xt = u.xt/a;
	w.yt = u.yt/a;
	w.zt = u.zt/a;

	w.xyz = u.xyz/a;
	w.xyt = u.xyt/a;
	w.xzt = u.xzt/a;
	w.yzt = u.yzt/a;

	w.xyzt = u.xyzt/a;

	return w;
}

//////////////////////////////////////////////////////

GA4E Wedge(GA4E u, GA4E v) {
	GA4E w;

/*

Wedge product, four dimensional spacetime

 ^    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     0     xy    xz    xt    0     0     xyz   0     xyt   xzt   0     0     0     xyzt  0    
 y    |  y    -xy    0     yz    yt    0    -xyz   0    -xyt   0     yzt   0     0    -xyzt  0     0    
 z    |  z    -xz   -yz    0     zt    xyz   0     0    -xzt  -yzt   0     0     xyzt  0     0     0    
 t    |  t    -xt   -yt   -zt    0     xyt   xzt   yzt   0     0     0    -xyzt  0     0     0     0    
 xy   |  xy    0     0     xyz   xyt   0     0     0     0     0     xyzt  0     0     0     0     0    
 xz   |  xz    0    -xyz   0     xzt   0     0     0     0    -xyzt  0     0     0     0     0     0    
 yz   |  yz    xyz   0     0     yzt   0     0     0     xyzt  0     0     0     0     0     0     0    
 xt   |  xt    0    -xyt  -xzt   0     0     0     xyzt  0     0     0     0     0     0     0     0    
 yt   |  yt    xyt   0    -yzt   0     0    -xyzt  0     0     0     0     0     0     0     0     0    
 zt   |  zt    xzt   yzt   0     0     xyzt  0     0     0     0     0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     xyzt  0     0     0     0     0     0     0     0     0     0     0    
 xyt  |  xyt   0     0    -xyzt  0     0     0     0     0     0     0     0     0     0     0     0    
 xzt  |  xzt   0     xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0    
 yzt  |  yzt  -xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xyzt |  xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

w.q    =  + u.q   *v.q   ;
w.x    =  + u.q   *v.x    + u.x   *v.q   ;
w.y    =  + u.q   *v.y    + u.y   *v.q   ;
w.z    =  + u.q   *v.z    + u.z   *v.q   ;
w.t    =  + u.q   *v.t    + u.t   *v.q   ;
w.xy   =  + u.q   *v.xy   + u.x   *v.y    - u.y   *v.x    + u.xy  *v.q   ;
w.xz   =  + u.q   *v.xz   + u.x   *v.z    - u.z   *v.x    + u.xz  *v.q   ;
w.yz   =  + u.q   *v.yz   + u.y   *v.z    - u.z   *v.y    + u.yz  *v.q   ;
w.xt   =  + u.q   *v.xt   + u.x   *v.t    - u.t   *v.x    + u.xt  *v.q   ;
w.yt   =  + u.q   *v.yt   + u.y   *v.t    - u.t   *v.y    + u.yt  *v.q   ;
w.zt   =  + u.q   *v.zt   + u.z   *v.t    - u.t   *v.z    + u.zt  *v.q   ;
w.xyz  =  + u.q   *v.xyz  + u.x   *v.yz   - u.y   *v.xz   + u.z   *v.xy   + u.xy  *v.z    - u.xz  *v.y    + u.yz  *v.x    + u.xyz *v.q   ;
w.xyt  =  + u.q   *v.xyt  + u.x   *v.yt   - u.y   *v.xt   + u.t   *v.xy   + u.xy  *v.t    - u.xt  *v.y    + u.yt  *v.x    + u.xyt *v.q   ;
w.xzt  =  + u.q   *v.xzt  + u.x   *v.zt   - u.z   *v.xt   + u.t   *v.xz   + u.xz  *v.t    - u.xt  *v.z    + u.zt  *v.x    + u.xzt *v.q   ;
w.yzt  =  + u.q   *v.yzt  + u.y   *v.zt   - u.z   *v.yt   + u.t   *v.yz   + u.yz  *v.t    - u.yt  *v.z    + u.zt  *v.y    + u.yzt *v.q   ;
w.xyzt =  + u.q   *v.xyzt + u.x   *v.yzt  - u.y   *v.xzt  + u.z   *v.xyt  - u.t   *v.xyz  + u.xy  *v.zt   - u.xz  *v.yt   + u.yz  *v.xt   + u.xt  *v.yz   - u.yt  *v.xz   + u.zt  *v.xy   + u.xyz *v.t    - u.xyt *v.z    + u.xzt *v.y    - u.yzt *v.x    + u.xyzt*v.q   ;

	return w;
}

//////////////////////////////////////////////////////

GA4E AntiWedge(GA4E a, GA4E b)
{
	GA4E c;

/*
Lengyel's AntiWedge product

 V    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0    -q     0     y    
 z    |  0     0     0     0     0     0     0     0     0     0     0     0     q     0     0     z    
 t    |  0     0     0     0     0     0     0     0     0     0     0    -q     0     0     0     t    
 xy   |  0     0     0     0     0     0     0     0     0     0     q     0     0     x     y     xy   
 xz   |  0     0     0     0     0     0     0     0     0    -q     0     0    -x     0     z     xz   
 yz   |  0     0     0     0     0     0     0     0     q     0     0     0    -y    -z     0     yz   
 xt   |  0     0     0     0     0     0     0     q     0     0     0     x     0     0     t     xt   
 yt   |  0     0     0     0     0     0    -q     0     0     0     0     y     0    -t     0     yt   
 zt   |  0     0     0     0     0     q     0     0     0     0     0     z     t     0     0     zt   
 xyz  |  0     0     0     0     q     0     0     0     x     y     z     0     xy    xz    yz    xyz  
 xyt  |  0     0     0    -q     0     0    -x    -y     0     0     t    -xy    0     xt    yt    xyt  
 xzt  |  0     0     q     0     0     x     0    -z     0    -t     0    -xz   -xt    0     zt    xzt  
 yzt  |  0    -q     0     0     0     y     z     0     t     0     0    -yz   -yt   -zt    0     yzt  
 xyzt |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 

*/

c.q    =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q    ; 
c.x    =  + a.x   *b.xyzt + a.xy  *b.xzt  - a.xz  *b.xyt  + a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.xyzt*b.x    ; 
c.y    =  + a.y   *b.xyzt + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  + a.xyz *b.yt   - a.xyt *b.yz   + a.yzt *b.xy   + a.xyzt*b.y    ; 
c.z    =  + a.z   *b.xyzt + a.xz  *b.yzt  - a.yz  *b.xzt  + a.zt  *b.xyz  + a.xyz *b.zt   - a.xzt *b.yz   + a.yzt *b.xz   + a.xyzt*b.z    ; 
c.t    =  + a.t   *b.xyzt + a.xt  *b.yzt  - a.yt  *b.xzt  + a.zt  *b.xyt  + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   + a.xyzt*b.t    ; 
c.xy   =  + a.xy  *b.xyzt + a.xyz *b.xyt  - a.xyt *b.xyz  + a.xyzt*b.xy   ; 
c.xz   =  + a.xz  *b.xyzt + a.xyz *b.xzt  - a.xzt *b.xyz  + a.xyzt*b.xz   ; 
c.yz   =  + a.yz  *b.xyzt + a.xyz *b.yzt  - a.yzt *b.xyz  + a.xyzt*b.yz   ; 
c.xt   =  + a.xt  *b.xyzt + a.xyt *b.xzt  - a.xzt *b.xyt  + a.xyzt*b.xt   ; 
c.yt   =  + a.yt  *b.xyzt + a.xyt *b.yzt  - a.yzt *b.xyt  + a.xyzt*b.yt   ; 
c.zt   =  + a.zt  *b.xyzt + a.xzt *b.yzt  - a.yzt *b.xzt  + a.xyzt*b.zt   ; 
c.xyz  =  + a.xyz *b.xyzt + a.xyzt*b.xyz  ; 
c.xyt  =  + a.xyt *b.xyzt + a.xyzt*b.xyt  ; 
c.xzt  =  + a.xzt *b.xyzt + a.xyzt*b.xzt  ; 
c.yzt  =  + a.yzt *b.xyzt + a.xyzt*b.yzt  ; 
c.xyzt =  + a.xyzt*b.xyzt ; 

	return c;
}

//////////////////////////////////////////////////////


GA4E Regressive(GA4E a, GA4E b) 
{
	GA4E c;

//	GA4E I, I_inv;
//	I.xyzt = 1;
//	I_inv.xyzt = 1;
//	c = ((a*I_inv)^(b*I_inv))*I;

/*

Hestenes' Regressive product differs from Lengyel's AntiWedge product in sign details

 V    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0    -q     x    
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     0     0     0     0     0     0     0     0    -q     0     0     z    
 t    |  0     0     0     0     0     0     0     0     0     0     0     q     0     0     0     t    
 xy   |  0     0     0     0     0     0     0     0     0     0     q     0     0     x     y     xy   
 xz   |  0     0     0     0     0     0     0     0     0    -q     0     0    -x     0     z     xz   
 yz   |  0     0     0     0     0     0     0     0     q     0     0     0    -y    -z     0     yz   
 xt   |  0     0     0     0     0     0     0     q     0     0     0     x     0     0     t     xt   
 yt   |  0     0     0     0     0     0    -q     0     0     0     0     y     0    -t     0     yt   
 zt   |  0     0     0     0     0     q     0     0     0     0     0     z     t     0     0     zt   
 xyz  |  0     0     0     0    -q     0     0     0     x     y     z     0    -xy   -xz   -yz    xyz  
 xyt  |  0     0     0     q     0     0    -x    -y     0     0     t     xy    0    -xt   -yt    xyt  
 xzt  |  0     0    -q     0     0     x     0    -z     0    -t     0     xz    xt    0    -zt    xzt  
 yzt  |  0     q     0     0     0     y     z     0     t     0     0     yz    yt    zt    0     yzt  
 xyzt |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 

*/

c.q    =  + a.q   *b.xyzt - a.x   *b.yzt  + a.y   *b.xzt  - a.z   *b.xyt  + a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   - a.xyz *b.t    + a.xyt *b.z    - a.xzt *b.y    + a.yzt *b.x    + a.xyzt*b.q   ;
c.x    =  + a.x   *b.xyzt + a.xy  *b.xzt  - a.xz  *b.xyt  + a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.xyzt*b.x   ;
c.y    =  + a.y   *b.xyzt + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  + a.xyz *b.yt   - a.xyt *b.yz   + a.yzt *b.xy   + a.xyzt*b.y   ;
c.z    =  + a.z   *b.xyzt + a.xz  *b.yzt  - a.yz  *b.xzt  + a.zt  *b.xyz  + a.xyz *b.zt   - a.xzt *b.yz   + a.yzt *b.xz   + a.xyzt*b.z   ;
c.t    =  + a.t   *b.xyzt + a.xt  *b.yzt  - a.yt  *b.xzt  + a.zt  *b.xyt  + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   + a.xyzt*b.t   ;
c.xy   =  + a.xy  *b.xyzt - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xyzt*b.xy  ;
c.xz   =  + a.xz  *b.xyzt - a.xyz *b.xzt  + a.xzt *b.xyz  + a.xyzt*b.xz  ;
c.yz   =  + a.yz  *b.xyzt - a.xyz *b.yzt  + a.yzt *b.xyz  + a.xyzt*b.yz  ;
c.xt   =  + a.xt  *b.xyzt - a.xyt *b.xzt  + a.xzt *b.xyt  + a.xyzt*b.xt  ;
c.yt   =  + a.yt  *b.xyzt - a.xyt *b.yzt  + a.yzt *b.xyt  + a.xyzt*b.yt  ;
c.zt   =  + a.zt  *b.xyzt - a.xzt *b.yzt  + a.yzt *b.xzt  + a.xyzt*b.zt  ;
c.xyz  =  + a.xyz *b.xyzt + a.xyzt*b.xyz ;
c.xyt  =  + a.xyt *b.xyzt + a.xyzt*b.xyt ;
c.xzt  =  + a.xzt *b.xyzt + a.xyzt*b.xzt ;
c.yzt  =  + a.yzt *b.xyzt + a.xyzt*b.yzt ;
c.xyzt =  + a.xyzt*b.xyzt ;

	return c;
}

//////////////////////////////////////////////////////

GA4E RegressiveViaFormula(GA4E a, GA4E b)
{
	GA4E c;

	GA4E I,I_inv;

	I = Zero();	I_inv = Zero();	
	I.xyzt = 1;	I_inv.xyzt = 1;

	c = Product(Wedge(Product(a,I_inv),Product(b,I_inv)),I);

	return c;
}

//////////////////////////////////////////////////////

GA4E LowerRightViaFormula(GA4E a, GA4E b)  // Like a mirror of Wedge along rising diagonal
{

/*

LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     xyzt 
 x    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     xyzt  yzt  
 y    |  0     0     0     0     0     0     0     0     0     0     0     0     0    -xyzt  0    -xzt  
 z    |  0     0     0     0     0     0     0     0     0     0     0     0     xyzt  0     0     xyt  
 t    |  0     0     0     0     0     0     0     0     0     0     0    -xyzt  0     0     0     xyz  
 xy   |  0     0     0     0     0     0     0     0     0     0     xyzt  0     0    -yzt   xzt  -zt   
 xz   |  0     0     0     0     0     0     0     0     0    -xyzt  0     0     yzt   0    -xyt   yt   
 yz   |  0     0     0     0     0     0     0     0     xyzt  0     0     0    -xzt   xyt   0    -xt   
 xt   |  0     0     0     0     0     0     0     xyzt  0     0     0    -yzt   0     0    -xyz   yz   
 yt   |  0     0     0     0     0     0    -xyzt  0     0     0     0     xzt   0     xyz   0    -xz   
 zt   |  0     0     0     0     0     xyzt  0     0     0     0     0    -xyt  -xyz   0     0     xy   
 xyz  |  0     0     0     0     xyzt  0     0     0     yzt  -xzt   xyt   0    -zt    yt   -xt   -t    
 xyt  |  0     0     0    -xyzt  0     0    -yzt   xzt   0     0     xyz   zt    0     yz   -xz   -z    
 xzt  |  0     0     xyzt  0     0     yzt   0    -xyt   0    -xyz   0    -yt   -yz    0     xy    y    
 yzt  |  0    -xyzt  0     0     0    -xzt   xyt   0     xyz   0     0     xt    xz   -xy    0    -x    
 xyzt |  xyzt -yzt   xzt  -xyt  -xyz  -zt    yt   -xt    yz   -xz    xy    t     z    -y     x    -q    

*/

// -Wedge(Blade[i]*xyzt,xyzt*Blade[j])

	GA4E c;

	GA4E I,I_inv;
	
	I = Zero();	I_inv = Zero();
	I.xyzt = 1;	I_inv.xyzt = 1;

	c = Wedge(Product(a,I_inv),Product(I,b));
	

	return c;
}

//////////////////////////////////////////////////////

GA4E Expander(GA4E a, GA4E b)
{
/*

Terms with increased rank
 >    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    xt    0     0     xyz   0     xyt   xzt   0     0     0     xyzt  0    
 y    |  0    -xy    0     yz    yt    0    -xyz   0    -xyt   0     yzt   0     0    -xyzt  0     0    
 z    |  0    -xz   -yz    0     zt    xyz   0     0    -xzt  -yzt   0     0     xyzt  0     0     0    
 t    |  0    -xt   -yt   -zt    0     xyt   xzt   yzt   0     0     0    -xyzt  0     0     0     0    
 xy   |  0     0     0     xyz   xyt   0     0     0     0     0     xyzt  0     0     0     0     0    
 xz   |  0     0    -xyz   0     xzt   0     0     0     0    -xyzt  0     0     0     0     0     0    
 yz   |  0     xyz   0     0     yzt   0     0     0     xyzt  0     0     0     0     0     0     0    
 xt   |  0     0    -xyt  -xzt   0     0     0     xyzt  0     0     0     0     0     0     0     0    
 yt   |  0     xyt   0    -yzt   0     0    -xyzt  0     0     0     0     0     0     0     0     0    
 zt   |  0     xzt   yzt   0     0     xyzt  0     0     0     0     0     0     0     0     0     0    
 xyz  |  0     0     0     0     xyzt  0     0     0     0     0     0     0     0     0     0     0    
 xyt  |  0     0     0    -xyzt  0     0     0     0     0     0     0     0     0     0     0     0    
 xzt  |  0     0     xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0    
 yzt  |  0    -xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xyzt |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/
	GA4E c;

c.q    =  0 ; 
c.x    =  0 ; 
c.y    =  0 ; 
c.z    =  0 ; 
c.t    =  0 ; 
c.xy   =  + a.x   *b.y    - a.y   *b.x    ; 
c.xz   =  + a.x   *b.z    - a.z   *b.x    ; 
c.yz   =  + a.y   *b.z    - a.z   *b.y    ; 
c.xt   =  + a.x   *b.t    - a.t   *b.x    ; 
c.yt   =  + a.y   *b.t    - a.t   *b.y    ; 
c.zt   =  + a.z   *b.t    - a.t   *b.z    ; 
c.xyz  =  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    ; 
c.xyt  =  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    ; 
c.xzt  =  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    ; 
c.yzt  =  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    ; 
c.xyzt =  + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    ; 

	return c;
}

//////////////////////////////////////////////////////

GA4E Conserver(GA4E a, GA4E b)
{
/*

Terms with preserved rank
 =    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 t    |  t     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0     0    -yz    xz   -yt    xt    0     0     0    -yzt   xzt   0    
 xz   |  xz    0     0     0     0     yz    0    -xy   -zt    0     xt    0     yzt   0    -xyt   0    
 yz   |  yz    0     0     0     0    -xz    xy    0     0    -zt    yt    0    -xzt   xyt   0     0    
 xt   |  xt    0     0     0     0     yt    zt    0     0    -xy   -xz   -yzt   0     0     xyz   0    
 yt   |  yt    0     0     0     0    -xt    0     zt    xy    0    -yz    xzt   0    -xyz   0     0    
 zt   |  zt    0     0     0     0     0    -xt   -yt    xz    yz    0    -xyt   xyz   0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0     yzt  -xzt   xyt   0     0     0     0     0    
 xyt  |  xyt   0     0     0     0     0    -yzt   xzt   0     0    -xyz   0     0     0     0     0    
 xzt  |  xzt   0     0     0     0     yzt   0    -xyt   0     xyz   0     0     0     0     0     0    
 yzt  |  yzt   0     0     0     0    -xzt   xyt   0    -xyz   0     0     0     0     0     0     0    
 xyzt |  xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

	GA4E c;

c.q    =  + a.q   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   - a.xt  *b.yt   + a.yt  *b.xt   ; 
c.xz   =  + a.q   *b.xz   + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   - a.xt  *b.zt   + a.zt  *b.xt   ; 
c.yz   =  + a.q   *b.yz   - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    - a.yt  *b.zt   + a.zt  *b.yt   ; 
c.xt   =  + a.q   *b.xt   + a.xy  *b.yt   + a.xz  *b.zt   + a.xt  *b.q    - a.yt  *b.xy   - a.zt  *b.xz   ; 
c.yt   =  + a.q   *b.yt   - a.xy  *b.xt   + a.yz  *b.zt   + a.xt  *b.xy   + a.yt  *b.q    - a.zt  *b.yz   ; 
c.zt   =  + a.q   *b.zt   - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   + a.zt  *b.q    ; 
c.xyz  =  + a.q   *b.xyz  + a.xt  *b.yzt  - a.yt  *b.xzt  + a.zt  *b.xyt  + a.xyz *b.q    - a.xyt *b.zt   + a.xzt *b.yt   - a.yzt *b.xt   ; 
c.xyt  =  + a.q   *b.xyt  - a.xz  *b.yzt  + a.yz  *b.xzt  - a.zt  *b.xyz  + a.xyz *b.zt   + a.xyt *b.q    - a.xzt *b.yz   + a.yzt *b.xz   ; 
c.xzt  =  + a.q   *b.xzt  + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  - a.xyz *b.yt   + a.xyt *b.yz   + a.xzt *b.q    - a.yzt *b.xy   ; 
c.yzt  =  + a.q   *b.yzt  - a.xy  *b.xzt  + a.xz  *b.xyt  - a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.yzt *b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.xyzt*b.q    ; 

	return c;
}

//////////////////////////////////////////////////////

GA4E Shrinker(GA4E a, GA4E b)
{

	GA4E c;

/*

Terms with reduced rank
 <    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     0     y     z     0     t     0     0     yz    yt    zt    0     yzt  
 y    |  0     0     q     0     0    -x     0     z     0     t     0    -xz   -xt    0     zt   -xzt  
 z    |  0     0     0     q     0     0    -x    -y     0     0     t     xy    0    -xt   -yt    xyt  
 t    |  0     0     0     0     q     0     0     0    -x    -y    -z     0     xy    xz    yz   -xyz  
 xy   |  0    -y     x     0     0    -q     0     0     0     0     0    -z    -t     0     0    -zt   
 xz   |  0    -z     0     x     0     0    -q     0     0     0     0     y     0    -t     0     yt   
 yz   |  0     0    -z     y     0     0     0    -q     0     0     0    -x     0     0    -t    -xt   
 xt   |  0    -t     0     0     x     0     0     0    -q     0     0     0     y     z     0    -yz   
 yt   |  0     0    -t     0     y     0     0     0     0    -q     0     0    -x     0     z     xz   
 zt   |  0     0     0    -t     z     0     0     0     0     0    -q     0     0    -x    -y    -xy   
 xyz  |  0     yz   -xz    xy    0    -z     y    -x     0     0     0    -q    -zt    yt   -xt   -t    
 xyt  |  0     yt   -xt    0     xy   -t     0     0     y    -x     0     zt   -q    -yz    xz    z    
 xzt  |  0     zt    0    -xt    xz    0    -t     0     z     0    -x    -yt    yz   -q    -xy   -y    
 yzt  |  0     0     zt   -yt    yz    0     0    -t     0     z    -y     xt   -xz    xy   -q     x    
 xyzt |  0    -yzt   xzt  -xyt   xyz  -zt    yt   -xt   -yz    xz   -xy    t    -z     y    -x     q    

*/


// Shrinker equation set for test purposes
c.q    =  + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xt  *b.xt   - a.yt  *b.yt   - a.zt  *b.zt   - a.xyz *b.xyz  - a.xyt *b.xyt  - a.xzt *b.xzt  - a.yzt *b.yzt  + a.xyzt*b.xyzt ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   - a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  + a.xt  *b.t    - a.yt  *b.xyt  - a.zt  *b.xzt  - a.xyz *b.yz   - a.xyt *b.yt   - a.xzt *b.zt   + a.yzt *b.xyzt - a.xyzt*b.yzt  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   - a.t   *b.yt   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    + a.xt  *b.xyt  + a.yt  *b.t    - a.zt  *b.yzt  + a.xyz *b.xz   + a.xyt *b.xt   - a.xzt *b.xyzt - a.yzt *b.zt   + a.xyzt*b.xzt  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   - a.t   *b.zt   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    + a.xt  *b.xzt  + a.yt  *b.yzt  + a.zt  *b.t    - a.xyz *b.xy   + a.xyt *b.xyzt + a.xzt *b.xt   + a.yzt *b.yt   - a.xyzt*b.xyt  ; 
c.t    =  + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   + a.xyzt*b.xyz  ; 
c.xy   =  + a.z   *b.xyz  + a.t   *b.xyt  - a.zt  *b.xyzt + a.xyz *b.z    + a.xyt *b.t    - a.xzt *b.yzt  + a.yzt *b.xzt  - a.xyzt*b.zt   ; 
c.xz   =  - a.y   *b.xyz  + a.t   *b.xzt  + a.yt  *b.xyzt - a.xyz *b.y    + a.xyt *b.yzt  + a.xzt *b.t    - a.yzt *b.xyt  + a.xyzt*b.yt   ; 
c.yz   =  + a.x   *b.xyz  + a.t   *b.yzt  - a.xt  *b.xyzt + a.xyz *b.x    - a.xyt *b.xzt  + a.xzt *b.xyt  + a.yzt *b.t    - a.xyzt*b.xt   ; 
c.xt   =  - a.y   *b.xyt  - a.z   *b.xzt  - a.yz  *b.xyzt - a.xyz *b.yzt  - a.xyt *b.y    - a.xzt *b.z    + a.yzt *b.xyz  - a.xyzt*b.yz   ; 
c.yt   =  + a.x   *b.xyt  - a.z   *b.yzt  + a.xz  *b.xyzt + a.xyz *b.xzt  + a.xyt *b.x    - a.xzt *b.xyz  - a.yzt *b.z    + a.xyzt*b.xz   ; 
c.zt   =  + a.x   *b.xzt  + a.y   *b.yzt  - a.xy  *b.xyzt - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy   ; 
c.xyz  =  - a.t   *b.xyzt + a.xyzt*b.t    ; 
c.xyt  =  + a.z   *b.xyzt - a.xyzt*b.z    ; 
c.xzt  =  - a.y   *b.xyzt + a.xyzt*b.y    ; 
c.yzt  =  + a.x   *b.xyzt - a.xyzt*b.x    ; 
c.xyzt = 0 ; 

	return c;
}

///////////////////////////////////////////////////////

GA4E Symmetric(GA4E a, GA4E b)
{
/*
Symmetric Product
 ?    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     q     0     0     0     0     0     xyz   0     xyt   xzt   yz    yt    zt    0     0    
 y    |  y     0     q     0     0     0    -xyz   0    -xyt   0     yzt  -xz   -xt    0     zt    0    
 z    |  z     0     0     q     0     xyz   0     0    -xzt  -yzt   0     xy    0    -xt   -yt    0    
 t    |  t     0     0     0     q     xyt   xzt   yzt   0     0     0     0     xy    xz    yz    0    
 xy   |  xy    0     0     xyz   xyt  -q     0     0     0     0     xyzt -z    -t     0     0    -zt   
 xz   |  xz    0    -xyz   0     xzt   0    -q     0     0    -xyzt  0     y     0    -t     0     yt   
 yz   |  yz    xyz   0     0     yzt   0     0    -q     xyzt  0     0    -x     0     0    -t    -xt   
 xt   |  xt    0    -xyt  -xzt   0     0     0     xyzt -q     0     0     0     y     z     0    -yz   
 yt   |  yt    xyt   0    -yzt   0     0    -xyzt  0     0    -q     0     0    -x     0     z     xz   
 zt   |  zt    xzt   yzt   0     0     xyzt  0     0     0     0    -q     0     0    -x    -y    -xy   
 xyz  |  xyz   yz   -xz    xy    0    -z     y    -x     0     0     0    -q     0     0     0     0    
 xyt  |  xyt   yt   -xt    0     xy   -t     0     0     y    -x     0     0    -q     0     0     0    
 xzt  |  xzt   zt    0    -xt    xz    0    -t     0     z     0    -x     0     0    -q     0     0    
 yzt  |  yzt   0     zt   -yt    yz    0     0    -t     0     z    -y     0     0     0    -q     0    
 xyzt |  xyzt  0     0     0     0    -zt    yt   -xt   -yz    xz   -xy    0     0     0     0     q    

*/
	GA4E c;

//	c = (a*b + b*a)/2;

// Symmetric Product Equations
c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    + a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   - a.xt  *b.xt   - a.yt  *b.yt   - a.zt  *b.zt   - a.xyz *b.xyz  - a.xyt *b.xyt  - a.xzt *b.xzt  - a.yzt *b.yzt  + a.xyzt*b.xyzt ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.yz  *b.xyz  - a.yt  *b.xyt  - a.zt  *b.xzt  - a.xyz *b.yz   - a.xyt *b.yt   - a.xzt *b.zt   ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    + a.xz  *b.xyz  + a.xt  *b.xyt  - a.zt  *b.yzt  + a.xyz *b.xz   + a.xyt *b.xt   - a.yzt *b.zt   ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    - a.xy  *b.xyz  + a.xt  *b.xzt  + a.yt  *b.yzt  - a.xyz *b.xy   + a.xzt *b.xt   + a.yzt *b.yt   ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   ; 
c.xy   =  + a.q   *b.xy   + a.z   *b.xyz  + a.t   *b.xyt  + a.xy  *b.q    - a.zt  *b.xyzt + a.xyz *b.z    + a.xyt *b.t    - a.xyzt*b.zt   ; 
c.xz   =  + a.q   *b.xz   - a.y   *b.xyz  + a.t   *b.xzt  + a.xz  *b.q    + a.yt  *b.xyzt - a.xyz *b.y    + a.xzt *b.t    + a.xyzt*b.yt   ; 
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.t   *b.yzt  + a.yz  *b.q    - a.xt  *b.xyzt + a.xyz *b.x    + a.yzt *b.t    - a.xyzt*b.xt   ; 
c.xt   =  + a.q   *b.xt   - a.y   *b.xyt  - a.z   *b.xzt  - a.yz  *b.xyzt + a.xt  *b.q    - a.xyt *b.y    - a.xzt *b.z    - a.xyzt*b.yz   ; 
c.yt   =  + a.q   *b.yt   + a.x   *b.xyt  - a.z   *b.yzt  + a.xz  *b.xyzt + a.yt  *b.q    + a.xyt *b.x    - a.yzt *b.z    + a.xyzt*b.xz   ; 
c.zt   =  + a.q   *b.zt   + a.x   *b.xzt  + a.y   *b.yzt  - a.xy  *b.xyzt + a.zt  *b.q    + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy   ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    + a.xyt *b.q    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    + a.xzt *b.q    ; 
c.yzt  =  + a.q   *b.yzt  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    + a.yzt *b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyzt*b.q    ; 

	return c;
}

//////////////////////////////////////////////////////

GA4E AntiSymmetric(GA4E a, GA4E b)
{
/*

AntiSymmetric Product
 ?    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     0     xy    xz    xt    y     z     0     t     0     0     0     0     0     xyzt  yzt  
 y    |  0    -xy    0     yz    yt   -x     0     z     0     t     0     0     0    -xyzt  0    -xzt  
 z    |  0    -xz   -yz    0     zt    0    -x    -y     0     0     t     0     xyzt  0     0     xyt  
 t    |  0    -xt   -yt   -zt    0     0     0     0    -x    -y    -z    -xyzt  0     0     0    -xyz  
 xy   |  0    -y     x     0     0     0    -yz    xz   -yt    xt    0     0     0    -yzt   xzt   0    
 xz   |  0    -z     0     x     0     yz    0    -xy   -zt    0     xt    0     yzt   0    -xyt   0    
 yz   |  0     0    -z     y     0    -xz    xy    0     0    -zt    yt    0    -xzt   xyt   0     0    
 xt   |  0    -t     0     0     x     yt    zt    0     0    -xy   -xz   -yzt   0     0     xyz   0    
 yt   |  0     0    -t     0     y    -xt    0     zt    xy    0    -yz    xzt   0    -xyz   0     0    
 zt   |  0     0     0    -t     z     0    -xt   -yt    xz    yz    0    -xyt   xyz   0     0     0    
 xyz  |  0     0     0     0     xyzt  0     0     0     yzt  -xzt   xyt   0    -zt    yt   -xt   -t    
 xyt  |  0     0     0    -xyzt  0     0    -yzt   xzt   0     0    -xyz   zt    0    -yz    xz    z    
 xzt  |  0     0     xyzt  0     0     yzt   0    -xyt   0     xyz   0    -yt    yz    0    -xy   -y    
 yzt  |  0    -xyzt  0     0     0    -xzt   xyt   0    -xyz   0     0     xt   -xz    xy    0     x    
 xyzt |  0    -yzt   xzt  -xyt   xyz   0     0     0     0     0     0     t    -z     y    -x     0    

*/

	GA4E c;

//	c = (a*b - b*a)/2;

// AntiSymmetric Product Equations
c.q    =  0 ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   - a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    + a.xt  *b.t    + a.yzt *b.xyzt - a.xyzt*b.yzt  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   - a.t   *b.yt   - a.xy  *b.x    + a.yz  *b.z    + a.yt  *b.t    - a.xzt *b.xyzt + a.xyzt*b.xzt  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   - a.t   *b.zt   - a.xz  *b.x    - a.yz  *b.y    + a.zt  *b.t    + a.xyt *b.xyzt - a.xyzt*b.xyt  ; 
c.t    =  + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt + a.xyzt*b.xyz  ; 
c.xy   =  + a.x   *b.y    - a.y   *b.x    - a.xz  *b.yz   + a.yz  *b.xz   - a.xt  *b.yt   + a.yt  *b.xt   - a.xzt *b.yzt  + a.yzt *b.xzt  ; 
c.xz   =  + a.x   *b.z    - a.z   *b.x    + a.xy  *b.yz   - a.yz  *b.xy   - a.xt  *b.zt   + a.zt  *b.xt   + a.xyt *b.yzt  - a.yzt *b.xyt  ; 
c.yz   =  + a.y   *b.z    - a.z   *b.y    - a.xy  *b.xz   + a.xz  *b.xy   - a.yt  *b.zt   + a.zt  *b.yt   - a.xyt *b.xzt  + a.xzt *b.xyt  ; 
c.xt   =  + a.x   *b.t    - a.t   *b.x    + a.xy  *b.yt   + a.xz  *b.zt   - a.yt  *b.xy   - a.zt  *b.xz   - a.xyz *b.yzt  + a.yzt *b.xyz  ; 
c.yt   =  + a.y   *b.t    - a.t   *b.y    - a.xy  *b.xt   + a.yz  *b.zt   + a.xt  *b.xy   - a.zt  *b.yz   + a.xyz *b.xzt  - a.xzt *b.xyz  ; 
c.zt   =  + a.z   *b.t    - a.t   *b.z    - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   - a.xyz *b.xyt  + a.xyt *b.xyz  ; 
c.xyz  =  - a.t   *b.xyzt + a.xt  *b.yzt  - a.yt  *b.xzt  + a.zt  *b.xyt  - a.xyt *b.zt   + a.xzt *b.yt   - a.yzt *b.xt   + a.xyzt*b.t    ; 
c.xyt  =  + a.z   *b.xyzt - a.xz  *b.yzt  + a.yz  *b.xzt  - a.zt  *b.xyz  + a.xyz *b.zt   - a.xzt *b.yz   + a.yzt *b.xz   - a.xyzt*b.z    ; 
c.xzt  =  - a.y   *b.xyzt + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  - a.xyz *b.yt   + a.xyt *b.yz   - a.yzt *b.xy   + a.xyzt*b.y    ; 
c.yzt  =  + a.x   *b.xyzt - a.xy  *b.xzt  + a.xz  *b.xyt  - a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   - a.xyzt*b.x    ; 
c.xyzt =  + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    ; 

	return c;
}

//////////////////////////////////////////////////////

GA4E Inner(GA4E a, GA4E b)
{

	GA4E c;

	c.q    = a.t*b.t + a.x*b.x - a.xt*b.xt - a.xy*b.xy - a.xyt*b.xyt - a.xyz*b.xyz + a.xyzt*b.xyzt
		 - a.xz*b.xz - a.xzt*b.xzt + a.y*b.y - a.yt*b.yt - a.yz*b.yz - a.yzt*b.yzt + a.z*b.z - a.zt*b.zt;

	c.x    =  + (-a.t*b.xt + a.xt*b.t + a.xy*b.y - a.xyt*b.yt - a.xyz*b.yz - a.xyzt*b.yzt + a.xz*b.z
		 - a.xzt*b.zt - a.y*b.xy - a.yt*b.xyt - a.yz*b.xyz + a.yzt*b.xyzt - a.z*b.xz - a.zt*b.xzt);
	c.y    =  + (-a.t*b.yt + a.x*b.xy + a.xt*b.xyt - a.xy*b.x + a.xyt*b.xt + a.xyz*b.xz + a.xyzt*b.xzt
		 + a.xz*b.xyz - a.xzt*b.xyzt + a.yt*b.t + a.yz*b.z - a.yzt*b.zt - a.z*b.yz - a.zt*b.yzt);
	c.z    =  + (-a.t*b.zt + a.x*b.xz + a.xt*b.xzt - a.xy*b.xyz + a.xyt*b.xyzt - a.xyz*b.xy - a.xyzt*b.xyt
		 - a.xz*b.x + a.xzt*b.xt + a.y*b.yz + a.yt*b.yzt - a.yz*b.y + a.yzt*b.yt + a.zt*b.t);
	c.t    =  + (a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz
		 - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z);


	c.xy = + ( a.t*b.xyt + a.xyt*b.t + a.xyz*b.z - a.xyzt*b.zt + a.z*b.xyz - a.zt*b.xyzt);
	c.xz = + ( a.t*b.xzt - a.xyz*b.y + a.xyzt*b.yt + a.xzt*b.t - a.y*b.xyz + a.yt*b.xyzt);
	c.yz = + ( a.t*b.yzt + a.x*b.xyz - a.xt*b.xyzt + a.xyz*b.x - a.xyzt*b.xt + a.yzt*b.t);
	c.xt = + (-a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt);
	c.yt = + ( a.x*b.xyt + a.xyt*b.x + a.xyzt*b.xz + a.xz*b.xyzt - a.yzt*b.z - a.z*b.yzt);
	c.zt = + ( a.x*b.xzt - a.xy*b.xyzt - a.xyzt*b.xy + a.xzt*b.x + a.y*b.yzt + a.yzt*b.y);

	c.xyz = + (-a.t*b.xyzt + a.xyzt*b.t);
	c.xyt = + (-a.xyzt*b.z + a.z*b.xyzt);
	c.xzt = + ( a.xyzt*b.y - a.y*b.xyzt);
	c.yzt = + ( a.x*b.xyzt - a.xyzt*b.x);

	c.xyzt = 0;

	return c;
}

//////////////////////////////////////////////////////

GA4E LeftContraction (GA4E a, GA4E b)
{

	GA4E c;

	c.q = a.q*b.q + a.t*b.t + a.x*b.x - a.xt*b.xt - a.xy*b.xy - a.xyt*b.xyt - a.xyz*b.xyz + a.xyzt*b.xyzt
		 - a.xz*b.xz - a.xzt*b.xzt + a.y*b.y - a.yt*b.yt - a.yz*b.yz - a.yzt*b.yzt + a.z*b.z - a.zt*b.zt ;

	c.x = + (a.q*b.x - a.t*b.xt - a.y*b.xy - a.yt*b.xyt - a.yz*b.xyz + a.yzt*b.xyzt - a.z*b.xz - a.zt*b.xzt) ;
	c.y = + (a.q*b.y - a.t*b.yt + a.x*b.xy + a.xt*b.xyt + a.xz*b.xyz - a.xzt*b.xyzt - a.z*b.yz - a.zt*b.yzt) ;
	c.z = + (a.q*b.z - a.t*b.zt + a.x*b.xz + a.xt*b.xzt - a.xy*b.xyz + a.xyt*b.xyzt + a.y*b.yz + a.yt*b.yzt) ;
	c.t = + (a.q*b.t + a.x*b.xt - a.xy*b.xyt - a.xyz*b.xyzt - a.xz*b.xzt + a.y*b.yt - a.yz*b.yzt + a.z*b.zt) ;

	c.xy = + (a.q*b.xy + a.t*b.xyt + a.z*b.xyz - a.zt*b.xyzt) ;
	c.xz = + (a.q*b.xz + a.t*b.xzt - a.y*b.xyz + a.yt*b.xyzt) ;
	c.yz = + (a.q*b.yz + a.t*b.yzt + a.x*b.xyz - a.xt*b.xyzt) ;
	c.xt = + (a.q*b.xt - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt) ;
	c.yt = + (a.q*b.yt + a.x*b.xyt + a.xz*b.xyzt - a.z*b.yzt) ;
	c.zt = + (a.q*b.zt + a.x*b.xzt - a.xy*b.xyzt + a.y*b.yzt) ;

	c.xyz = + (a.q*b.xyz - a.t*b.xyzt) ;
	c.xyt = + (a.q*b.xyt + a.z*b.xyzt) ;
	c.xzt = + (a.q*b.xzt - a.y*b.xyzt) ;
	c.yzt = + (a.q*b.yzt + a.x*b.xyzt) ;

	c.xyzt = + a.q*b.xyzt ;

	return c;
}

//////////////////////////////////////////////////////

GA4E RightContraction (GA4E a, GA4E b)
{

	GA4E c;

	c.q = a.q*b.q + a.t*b.t + a.x*b.x - a.xt*b.xt - a.xy*b.xy - a.xyt*b.xyt - a.xyz*b.xyz + a.xyzt*b.xyzt
		 - a.xz*b.xz - a.xzt*b.xzt + a.y*b.y - a.yt*b.yt - a.yz*b.yz - a.yzt*b.yzt + a.z*b.z - a.zt*b.zt ;

	c.x = + ( a.x*b.q + a.xt*b.t + a.xy*b.y - a.xyt*b.yt - a.xyz*b.yz - a.xyzt*b.yzt + a.xz*b.z - a.xzt*b.zt) ;
	c.y = + (-a.xy*b.x + a.xyt*b.xt + a.xyz*b.xz + a.xyzt*b.xzt + a.y*b.q + a.yt*b.t + a.yz*b.z - a.yzt*b.zt) ;
	c.z = + (-a.xyz*b.xy - a.xyzt*b.xyt - a.xz*b.x + a.xzt*b.xt - a.yz*b.y + a.yzt*b.yt + a.z*b.q + a.zt*b.t) ;
	c.t = + ( a.t*b.q - a.xt*b.x - a.xyt*b.xy + a.xyzt*b.xyz - a.xzt*b.xz - a.yt*b.y - a.yzt*b.yz - a.zt*b.z) ;

	c.xy = + ( a.xy*b.q + a.xyt*b.t + a.xyz*b.z - a.xyzt*b.zt) ;
	c.xz = + (-a.xyz*b.y + a.xyzt*b.yt + a.xz*b.q + a.xzt*b.t) ;
	c.yz = + ( a.xyz*b.x - a.xyzt*b.xt + a.yz*b.q + a.yzt*b.t) ;
	c.xt = + ( a.xt*b.q - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z) ;
	c.yt = + ( a.xyt*b.x + a.xyzt*b.xz + a.yt*b.q - a.yzt*b.z) ;
	c.zt = + (-a.xyzt*b.xy + a.xzt*b.x + a.yzt*b.y + a.zt*b.q) ;

	c.xyz = + ( a.xyz*b.q + a.xyzt*b.t) ;
	c.xyt = + ( a.xyt*b.q - a.xyzt*b.z) ;
	c.xzt = + ( a.xyzt*b.y + a.xzt*b.q) ;
	c.yzt = + (-a.xyzt*b.x + a.yzt*b.q) ;

	c.xyzt = + a.xyzt*b.q ;

	return c;
}

//////////////////////////////////////////////////////

long Determinant(GA4E A) {

	GA4E B, C;
	long a, b, c, d, e, s, det;

	B = Reverse(A);
	C = Product(B,A);

	a = C.q;
	b = C.x;
	c = C.y;
	d = C.z;
	e = C.t;

	s = C.xyzt;
	det = (a*a - b*b - c*c - d*d - e*e - s*s);

	return(det);
}

//////////////////////////////////////////////////////

GA4E Adjugate(GA4E a)
{

	GA4E u;

//Product(Reverse(r),Conjugation(Product(r,Reverse(r)))) = Adjugate(r);

	u = Product(Reverse(a),Conjugation(Product(a,Reverse(a))));

	return u;

}

//////////////////////////////////////////////////////

// int	main(void) {  return 0;}

