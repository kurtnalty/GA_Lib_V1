// Routines for Geometric Algebra in Minkowski spacetime
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

struct Mink{	// 0 . 1 2 4 8 . 3 5 6 9 10 12 . 7 11 13 14 . 15
	ex q,  x,y,z,t, xy,xz,yz,xt,yt,zt,  xyz,xyt,xzt,yzt,  xyzt;
	Mink() {q = 0; x = 0; y = 0; z = 0; t = 0; 
		xy = 0; xz = 0; yz = 0;  xt = 0; yt = 0; zt = 0;
		xyz = 0; xyt = 0; xzt = 0; yzt = 0; xyzt=0;}
	Mink(ex qq, ex xx, ex yy, ex zz, ex tt, 
		ex xxyy, ex xxzz, ex yyzz, ex xxtt, ex yytt, ex zztt,
		ex xxyyzz, ex xxyytt, ex xxzztt, ex yyzztt,   ex xxyyzztt) 

		{q = qq; x = xx; y = yy; z = zz; t = tt;
		 xy = xxyy; xz = xxzz; yz = yyzz; xt = xxtt; yt = yytt; zt = zztt; 
		 xyz = xxyyzz;  xyt = xxyytt; xzt = xxzztt; yzt = yyzztt; xyzt=xxyyzztt;}
};


//////////////////////////////////////////////////////

// Necessary forward declarations

Mink LeftContraction (const Mink &a, const Mink &b);

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, Mink &v) 
{
	return ff << "\n(" 
		<< v.q << ", \n" 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.t << ", \n" 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << "," 
		<< v.xt << "," 
		<< v.yt << "," 
		<< v.zt << ", \n" 
		<< v.xyz << ","
		<< v.xyt << ","
		<< v.xzt << ","
		<< v.yzt << ", \n"
		<< v.xyzt << ")\n";
}


//////////////////////////////////////////////////////

void PrintMV(Mink &v) 
{
	cout 	<< "( "
		<< v.q << ", " 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.t << ", " 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.yz << "," 
		<< v.xt << "," 
		<< v.yt << "," 
		<< v.zt << ", " 
		<< v.xyz << ","
		<< v.xyt << ","
		<< v.xzt << ","
		<< v.yzt << ", "
		<< v.xyzt << ")";
}

//////////////////////////////////////////////////////

Mink OverBar(Mink a)
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	Mink b;		// a blade wedge b blade = pseudovector blade
// OverBar(Blade[i]) = Reciprocal(Blade[i])*xyzt ;

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

Mink UnderBar(Mink a)   
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	Mink b;		// b blade wedge a blade = pseudovector blade
// UnderBar(Blade[i]) = xyzt*Reciprocal(Blade[i]) ;

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

Mink Reverse(Mink w)	
// (A,  B, C, D, E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R,  S) score = 10 N( 2046) Reverse
{
	Mink v;
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

Mink Involution(Mink w)	
// Corresponds to PT parity transform: x -> -x, y -> -y, z ->-z, t -> -t
// (A, -B,-C,-D,-E,  F, G, H, J, K, L, -M,-N,-P,-R,  S) score =  0 N(30750) Involution
{
	Mink v;
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

Mink Hermitian(Mink w)	
//(A,  B, C, D,-E, -F,-G,-H, J, K, L, -M, N, P, R, -S) Reverse and time sign change
{
	Mink v;
	v.q =  w.q;

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

	return v;
}

//////////////////////////////////////////////////////

Mink Hermitian_Component(Mink w)	
//(A,  B, C, D, 0,  0, 0, 0, J, K, L,  0, N, P, R,  0)
{
	Mink v;
	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;
	v.t = 0;

	v.xy = 0;
	v.xz = 0;
	v.yz = 0;
	v.xt =  w.xt;
	v.yt =  w.yt;
	v.zt =  w.zt;

	v.xyz = 0;
	v.xyt =  w.xyt;
	v.xzt =  w.xzt;
	v.yzt =  w.yzt;

	v.xyzt = 0;

	return v;
}

//////////////////////////////////////////////////////

Mink Anti_Hermitian_Component(Mink w)	
//(0,  0, 0, 0, E,  F, G, H, 0, 0, 0, M, 0, 0, 0, S)
{
	Mink v;
	v.q =  0;

	v.x =  0;
	v.y =  0;
	v.z =  0;
	v.t = w.t;

	v.xy = w.xy;
	v.xz = w.xz;
	v.yz = w.yz;
	v.xt =  0;
	v.yt =  0;
	v.zt =  0;

	v.xyz = w.xyz;
	v.xyt =  0;
	v.xzt =  0;
	v.yzt =  0;

	v.xyzt = w.xyzt;

	return v;
}

//////////////////////////////////////////////////////

Mink Transpose(Mink w)	
//(A,  B, C, D,-E, -F,-G,-H, J, K, L, -M, N, P, R, -S) score =  6 N( 3857) Transpose  
{
	Mink v;
	v.q =  w.q;

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

	return v;
}

//////////////////////////////////////////////////////

Mink Conjugation(Mink w)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R, -S)
{
	Mink v;
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

Mink CliffordConjugation(Mink w)
// (A, -B,-C,-D,-E, -F,-G,-H,-J,-K,-L,  M, N, P, R,  S) score = 10 N(32736) Clifford Conjugation
{
	Mink v;
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

Mink Dual(Mink w)   // return w*I_inv 
// u = r*s = (p, o,-n,m,l, -k,j,-i,h,-g,f, -e,-d,c,-b, -a)

{
	Mink v;
	v.q =  w.xyzt;

	v.x =  w.yzt;
	v.y = -w.xzt;
	v.z =  w.xyt;
	v.t =  w.xyz;

	v.xy = -w.zt;
	v.xz =  w.yt;
	v.yz = -w.xt;
	v.xt =  w.yz;
	v.yt = -w.xz;
	v.zt =  w.xy;

	v.xyz = -w.t;
	v.xyt = -w.z;
	v.xzt =  w.y;
	v.yzt = -w.x;

	v.xyzt = -w.q;

	return v;
}

//////////////////////////////////////////////////////

Mink DorstDual(Mink a)
{
	Mink b, I_inv;

	I_inv.xyzt = -1;

	b = LeftContraction(a,I_inv);
	return b;
}

//////////////////////////////////////////////////////

Mink DorstUnDual(Mink a)
{
	Mink b, I;

	I.xyzt = 1;

	b = LeftContraction(a,I);
	return b;
}

//////////////////////////////////////////////////////

Mink operator+(const Mink &u, const Mink &v)
{
	Mink w;
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

Mink operator-(const Mink &u, const Mink &v)
{
	Mink w;
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

int operator==(const Mink &u, const Mink &v)
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

int operator!=(const Mink &u, const Mink &v)
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

Mink operator*(const Mink &u, const Mink &v) {
	Mink w;

/*  

Geometric product, Minkowski spacetime

 *    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     q     xy    xz    xt    y     z     xyz   t     xyt   xzt   yz    yt    zt    xyzt  yzt  
 y    |  y    -xy    q     yz    yt   -x    -xyz   z    -xyt   t     yzt  -xz   -xt   -xyzt  zt   -xzt  
 z    |  z    -xz   -yz    q     zt    xyz  -x    -y    -xzt  -yzt   t     xy    xyzt -xt   -yt    xyt  
 t    |  t    -xt   -yt   -zt   -q     xyt   xzt   yzt   x     y     z    -xyzt -xy   -xz   -yz    xyz  
 xy   |  xy   -y     x     xyz   xyt  -q    -yz    xz   -yt    xt    xyzt -z    -t    -yzt   xzt  -zt   
 xz   |  xz   -z    -xyz   x     xzt   yz   -q    -xy   -zt   -xyzt  xt    y     yzt  -t    -xyt   yt   
 yz   |  yz    xyz  -z     y     yzt  -xz    xy   -q     xyzt -zt    yt   -x    -xzt   xyt  -t    -xt   
 xt   |  xt   -t    -xyt  -xzt  -x     yt    zt    xyzt  q     xy    xz   -yzt  -y    -z    -xyz   yz   
 yt   |  yt    xyt  -t    -yzt  -y    -xt   -xyzt  zt   -xy    q     yz    xzt   x     xyz  -z    -xz   
 zt   |  zt    xzt   yzt  -t    -z     xyzt -xt   -yt   -xz   -yz    q    -xyt  -xyz   x     y     xy   
 xyz  |  xyz   yz   -xz    xy    xyzt -z     y    -x     yzt  -xzt   xyt  -q    -zt    yt   -xt   -t    
 xyt  |  xyt   yt   -xt   -xyzt -xy   -t    -yzt   xzt  -y     x     xyz   zt    q     yz   -xz   -z    
 xzt  |  xzt   zt    xyzt -xt   -xz    yzt  -t    -xyt  -z    -xyz   x    -yt   -yz    q     xy    y    
 yzt  |  yzt  -xyzt  zt   -yt   -yz   -xzt   xyt  -t     xyz  -z     y     xt    xz   -xy    q    -x    
 xyzt |  xyzt -yzt   xzt  -xyt  -xyz  -zt    yt   -xt    yz   -xz    xy    t     z    -y     x    -q    

*/

w.q    = + u.q*v.q + u.x*v.x + u.y*v.y + u.z*v.z - u.t*v.t - u.xy*v.xy - u.xz*v.xz - u.yz*v.yz
         + u.xt*v.xt + u.yt*v.yt + u.zt*v.zt - u.xyz*v.xyz + u.xyt*v.xyt + u.xzt*v.xzt + u.yzt*v.yzt - u.xyzt*v.xyzt;
w.x    = + u.q*v.x + u.x*v.q - u.y*v.xy - u.z*v.xz + u.t*v.xt + u.xy*v.y + u.xz*v.z - u.yz*v.xyz
         - u.xt*v.t + u.yt*v.xyt + u.zt*v.xzt - u.xyz*v.yz + u.xyt*v.yt + u.xzt*v.zt - u.yzt*v.xyzt + u.xyzt*v.yzt ;
w.y    = + u.q*v.y + u.x*v.xy + u.y*v.q - u.z*v.yz + u.t*v.yt - u.xy*v.x + u.xz*v.xyz + u.yz*v.z 
         - u.xt*v.xyt - u.yt*v.t + u.zt*v.yzt + u.xyz*v.xz - u.xyt*v.xt + u.xzt*v.xyzt + u.yzt*v.zt - u.xyzt*v.xzt ;
w.z    = + u.q*v.z + u.x*v.xz + u.y*v.yz + u.z*v.q + u.t*v.zt - u.xy*v.xyz - u.xz*v.x - u.yz*v.y
         - u.xt*v.xzt - u.yt*v.yzt - u.zt*v.t - u.xyz*v.xy - u.xyt*v.xyzt - u.xzt*v.xt - u.yzt*v.yt + u.xyzt*v.xyt ;
w.t    = + u.q*v.t + u.x*v.xt + u.y*v.yt + u.z*v.zt + u.t*v.q - u.xy*v.xyt - u.xz*v.xzt - u.yz*v.yzt 
         - u.xt*v.x - u.yt*v.y - u.zt*v.z - u.xyz*v.xyzt - u.xyt*v.xy - u.xzt*v.xz - u.yzt*v.yz + u.xyzt*v.xyz ;
w.xy   = + u.q*v.xy + u.x*v.y - u.y*v.x + u.z*v.xyz - u.t*v.xyt + u.xy*v.q - u.xz*v.yz + u.yz*v.xz
         + u.xt*v.yt - u.yt*v.xt + u.zt*v.xyzt + u.xyz*v.z - u.xyt*v.t + u.xzt*v.yzt - u.yzt*v.xzt + u.xyzt*v.zt ;
w.xz   = + u.q*v.xz + u.x*v.z - u.y*v.xyz - u.z*v.x - u.t*v.xzt + u.xy*v.yz + u.xz*v.q - u.yz*v.xy
         + u.xt*v.zt - u.yt*v.xyzt - u.zt*v.xt - u.xyz*v.y - u.xyt*v.yzt - u.xzt*v.t + u.yzt*v.xyt - u.xyzt*v.yt ;
w.yz   = + u.q*v.yz + u.x*v.xyz + u.y*v.z - u.z*v.y - u.t*v.yzt - u.xy*v.xz + u.xz*v.xy + u.yz*v.q 
         + u.xt*v.xyzt + u.yt*v.zt - u.zt*v.yt + u.xyz*v.x + u.xyt*v.xzt - u.xzt*v.xyt - u.yzt*v.t + u.xyzt*v.xt ;
w.xt   = + u.q*v.xt + u.x*v.t - u.y*v.xyt - u.z*v.xzt - u.t*v.x + u.xy*v.yt + u.xz*v.zt - u.yz*v.xyzt
         + u.xt*v.q - u.yt*v.xy - u.zt*v.xz - u.xyz*v.yzt - u.xyt*v.y - u.xzt*v.z + u.yzt*v.xyz - u.xyzt*v.yz ;
w.yt   = + u.q*v.yt + u.x*v.xyt + u.y*v.t - u.z*v.yzt - u.t*v.y - u.xy*v.xt + u.xz*v.xyzt + u.yz*v.zt
         + u.xt*v.xy + u.yt*v.q - u.zt*v.yz + u.xyz*v.xzt + u.xyt*v.x - u.xzt*v.xyz - u.yzt*v.z + u.xyzt*v.xz ;
w.zt   = + u.q*v.zt + u.x*v.xzt + u.y*v.yzt + u.z*v.t - u.t*v.z - u.xy*v.xyzt - u.xz*v.xt - u.yz*v.yt
         + u.xt*v.xz + u.yt*v.yz + u.zt*v.q - u.xyz*v.xyt + u.xyt*v.xyz + u.xzt*v.x + u.yzt*v.y - u.xyzt*v.xy ;
w.xyz  = + u.q*v.xyz + u.x*v.yz - u.y*v.xz + u.z*v.xy + u.t*v.xyzt + u.xy*v.z - u.xz*v.y + u.yz*v.x
         - u.xt*v.yzt + u.yt*v.xzt - u.zt*v.xyt + u.xyz*v.q + u.xyt*v.zt - u.xzt*v.yt + u.yzt*v.xt - u.xyzt*v.t ;
w.xyt  = + u.q*v.xyt + u.x*v.yt - u.y*v.xt + u.z*v.xyzt + u.t*v.xy + u.xy*v.t - u.xz*v.yzt + u.yz*v.xzt
         - u.xt*v.y + u.yt*v.x - u.zt*v.xyz + u.xyz*v.zt + u.xyt*v.q - u.xzt*v.yz + u.yzt*v.xz - u.xyzt*v.z ;
w.xzt  = + u.q*v.xzt + u.x*v.zt - u.y*v.xyzt - u.z*v.xt + u.t*v.xz + u.xy*v.yzt + u.xz*v.t - u.yz*v.xyt
         - u.xt*v.z + u.yt*v.xyz + u.zt*v.x - u.xyz*v.yt + u.xyt*v.yz + u.xzt*v.q - u.yzt*v.xy + u.xyzt*v.y ;
w.yzt  = + u.q*v.yzt + u.x*v.xyzt + u.y*v.zt - u.z*v.yt + u.t*v.yz - u.xy*v.xzt + u.xz*v.xyt + u.yz*v.t
         - u.xt*v.xyz - u.yt*v.z + u.zt*v.y + u.xyz*v.xt - u.xyt*v.xz + u.xzt*v.xy + u.yzt*v.q - u.xyzt*v.x ;
w.xyzt = + u.q*v.xyzt + u.x*v.yzt - u.y*v.xzt + u.z*v.xyt - u.t*v.xyz + u.xy*v.zt - u.xz*v.yt + u.yz*v.xt
         + u.xt*v.yz - u.yt*v.xz + u.zt*v.xy + u.xyz*v.t - u.xyt*v.z + u.xzt*v.y - u.yzt*v.x + u.xyzt*v.q ;


	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.t = expand(w.t);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xt = expand(w.xt);
	w.yt = expand(w.yt);
	w.zt = expand(w.zt);
	w.xyz = expand(w.xyz);
	w.xyt = expand(w.xyt);
	w.xzt = expand(w.xzt);
	w.yzt = expand(w.yzt);
	w.xyzt = expand(w.xyzt);

	return w;
}

//////////////////////////////////////////////////////

Mink Product(const Mink &u, const Mink &v)
{

	Mink w;

w.q    = + u.q*v.q + u.x*v.x + u.y*v.y + u.z*v.z - u.t*v.t - u.xy*v.xy - u.xz*v.xz - u.yz*v.yz
         + u.xt*v.xt + u.yt*v.yt + u.zt*v.zt - u.xyz*v.xyz + u.xyt*v.xyt + u.xzt*v.xzt + u.yzt*v.yzt - u.xyzt*v.xyzt;
w.x    = + u.q*v.x + u.x*v.q - u.y*v.xy - u.z*v.xz + u.t*v.xt + u.xy*v.y + u.xz*v.z - u.yz*v.xyz
         - u.xt*v.t + u.yt*v.xyt + u.zt*v.xzt - u.xyz*v.yz + u.xyt*v.yt + u.xzt*v.zt - u.yzt*v.xyzt + u.xyzt*v.yzt ;
w.y    = + u.q*v.y + u.x*v.xy + u.y*v.q - u.z*v.yz + u.t*v.yt - u.xy*v.x + u.xz*v.xyz + u.yz*v.z 
         - u.xt*v.xyt - u.yt*v.t + u.zt*v.yzt + u.xyz*v.xz - u.xyt*v.xt + u.xzt*v.xyzt + u.yzt*v.zt - u.xyzt*v.xzt ;
w.z    = + u.q*v.z + u.x*v.xz + u.y*v.yz + u.z*v.q + u.t*v.zt - u.xy*v.xyz - u.xz*v.x - u.yz*v.y
         - u.xt*v.xzt - u.yt*v.yzt - u.zt*v.t - u.xyz*v.xy - u.xyt*v.xyzt - u.xzt*v.xt - u.yzt*v.yt + u.xyzt*v.xyt ;
w.t    = + u.q*v.t + u.x*v.xt + u.y*v.yt + u.z*v.zt + u.t*v.q - u.xy*v.xyt - u.xz*v.xzt - u.yz*v.yzt 
         - u.xt*v.x - u.yt*v.y - u.zt*v.z - u.xyz*v.xyzt - u.xyt*v.xy - u.xzt*v.xz - u.yzt*v.yz + u.xyzt*v.xyz ;
w.xy   = + u.q*v.xy + u.x*v.y - u.y*v.x + u.z*v.xyz - u.t*v.xyt + u.xy*v.q - u.xz*v.yz + u.yz*v.xz
         + u.xt*v.yt - u.yt*v.xt + u.zt*v.xyzt + u.xyz*v.z - u.xyt*v.t + u.xzt*v.yzt - u.yzt*v.xzt + u.xyzt*v.zt ;
w.xz   = + u.q*v.xz + u.x*v.z - u.y*v.xyz - u.z*v.x - u.t*v.xzt + u.xy*v.yz + u.xz*v.q - u.yz*v.xy
         + u.xt*v.zt - u.yt*v.xyzt - u.zt*v.xt - u.xyz*v.y - u.xyt*v.yzt - u.xzt*v.t + u.yzt*v.xyt - u.xyzt*v.yt ;
w.yz   = + u.q*v.yz + u.x*v.xyz + u.y*v.z - u.z*v.y - u.t*v.yzt - u.xy*v.xz + u.xz*v.xy + u.yz*v.q 
         + u.xt*v.xyzt + u.yt*v.zt - u.zt*v.yt + u.xyz*v.x + u.xyt*v.xzt - u.xzt*v.xyt - u.yzt*v.t + u.xyzt*v.xt ;
w.xt   = + u.q*v.xt + u.x*v.t - u.y*v.xyt - u.z*v.xzt - u.t*v.x + u.xy*v.yt + u.xz*v.zt - u.yz*v.xyzt
         + u.xt*v.q - u.yt*v.xy - u.zt*v.xz - u.xyz*v.yzt - u.xyt*v.y - u.xzt*v.z + u.yzt*v.xyz - u.xyzt*v.yz ;
w.yt   = + u.q*v.yt + u.x*v.xyt + u.y*v.t - u.z*v.yzt - u.t*v.y - u.xy*v.xt + u.xz*v.xyzt + u.yz*v.zt
         + u.xt*v.xy + u.yt*v.q - u.zt*v.yz + u.xyz*v.xzt + u.xyt*v.x - u.xzt*v.xyz - u.yzt*v.z + u.xyzt*v.xz ;
w.zt   = + u.q*v.zt + u.x*v.xzt + u.y*v.yzt + u.z*v.t - u.t*v.z - u.xy*v.xyzt - u.xz*v.xt - u.yz*v.yt
         + u.xt*v.xz + u.yt*v.yz + u.zt*v.q - u.xyz*v.xyt + u.xyt*v.xyz + u.xzt*v.x + u.yzt*v.y - u.xyzt*v.xy ;
w.xyz  = + u.q*v.xyz + u.x*v.yz - u.y*v.xz + u.z*v.xy + u.t*v.xyzt + u.xy*v.z - u.xz*v.y + u.yz*v.x
         - u.xt*v.yzt + u.yt*v.xzt - u.zt*v.xyt + u.xyz*v.q + u.xyt*v.zt - u.xzt*v.yt + u.yzt*v.xt - u.xyzt*v.t ;
w.xyt  = + u.q*v.xyt + u.x*v.yt - u.y*v.xt + u.z*v.xyzt + u.t*v.xy + u.xy*v.t - u.xz*v.yzt + u.yz*v.xzt
         - u.xt*v.y + u.yt*v.x - u.zt*v.xyz + u.xyz*v.zt + u.xyt*v.q - u.xzt*v.yz + u.yzt*v.xz - u.xyzt*v.z ;
w.xzt  = + u.q*v.xzt + u.x*v.zt - u.y*v.xyzt - u.z*v.xt + u.t*v.xz + u.xy*v.yzt + u.xz*v.t - u.yz*v.xyt
         - u.xt*v.z + u.yt*v.xyz + u.zt*v.x - u.xyz*v.yt + u.xyt*v.yz + u.xzt*v.q - u.yzt*v.xy + u.xyzt*v.y ;
w.yzt  = + u.q*v.yzt + u.x*v.xyzt + u.y*v.zt - u.z*v.yt + u.t*v.yz - u.xy*v.xzt + u.xz*v.xyt + u.yz*v.t
         - u.xt*v.xyz - u.yt*v.z + u.zt*v.y + u.xyz*v.xt - u.xyt*v.xz + u.xzt*v.xy + u.yzt*v.q - u.xyzt*v.x ;
w.xyzt = + u.q*v.xyzt + u.x*v.yzt - u.y*v.xzt + u.z*v.xyt - u.t*v.xyz + u.xy*v.zt - u.xz*v.yt + u.yz*v.xt
         + u.xt*v.yz - u.yt*v.xz + u.zt*v.xy + u.xyz*v.t - u.xyt*v.z + u.xzt*v.y - u.yzt*v.x + u.xyzt*v.q ;


	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.t = expand(w.t);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xt = expand(w.xt);
	w.yt = expand(w.yt);
	w.zt = expand(w.zt);
	w.xyz = expand(w.xyz);
	w.xyt = expand(w.xyt);
	w.xzt = expand(w.xzt);
	w.yzt = expand(w.yzt);
	w.xyzt = expand(w.xyzt);

	return w;
}

//////////////////////////////////////////////////////

Mink operator*(const Mink &v, const ex &u) {
	Mink w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.t    = expand(u*v.t   ) ;

	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xt   = expand(u*v.xt  ) ;
	w.yt   = expand(u*v.yt  ) ;
	w.zt   = expand(u*v.zt  ) ;

	w.xyz  = expand(u*v.xyz ) ;
	w.xyt  = expand(u*v.xyt ) ;
	w.xzt  = expand(u*v.xzt ) ;
	w.yzt  = expand(u*v.yzt ) ;
	w.xyzt = expand(u*v.xyzt) ;

	return w;
}

//////////////////////////////////////////////////////

Mink operator*(const ex &u, const Mink &v) {
	Mink w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.z    = expand(u*v.z   ) ;
	w.t    = expand(u*v.t   ) ;

	w.xy   = expand(u*v.xy  ) ;
	w.xz   = expand(u*v.xz  ) ;
	w.yz   = expand(u*v.yz  ) ;
	w.xt   = expand(u*v.xt  ) ;
	w.yt   = expand(u*v.yt  ) ;
	w.zt   = expand(u*v.zt  ) ;

	w.xyz  = expand(u*v.xyz ) ;
	w.xyt  = expand(u*v.xyt ) ;
	w.xzt  = expand(u*v.xzt ) ;
	w.yzt  = expand(u*v.yzt ) ;
	w.xyzt = expand(u*v.xyzt) ;

	return w;
}


//////////////////////////////////////////////////////

Mink operator/(const Mink &u, const int i)
{
	Mink w;
	w.q = u.q/i;

	w.x = u.x/i;
	w.y = u.y/i;
	w.z = u.z/i;
	w.t = u.t/i;

	w.xy = u.xy/i;
	w.xz = u.xz/i;
	w.yz = u.yz/i;
	w.xt = u.xt/i;
	w.yt = u.yt/i;
	w.zt = u.zt/i;

	w.xyz = u.xyz/i;
	w.xyt = u.xyt/i;
	w.xzt = u.xzt/i;
	w.yzt = u.yzt/i;

	w.xyzt = u.xyzt/i;

	return w;
}

//////////////////////////////////////////////////////

Mink operator^(const Mink &u, const Mink &v) {
	Mink w;

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

	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.t = expand(w.t);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xt = expand(w.xt);
	w.yt = expand(w.yt);
	w.zt = expand(w.zt);
	w.xyz = expand(w.xyz);
	w.xyt = expand(w.xyt);
	w.xzt = expand(w.xzt);
	w.yzt = expand(w.yzt);
	w.xyzt = expand(w.xyzt);

	return w;
}

//////////////////////////////////////////////////////

Mink Wedge(const Mink &u, const Mink &v) {
	Mink w;

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

	w.q = expand(w.q);
	w.x = expand(w.x);
	w.y = expand(w.y);
	w.z = expand(w.z);
	w.t = expand(w.t);
	w.xy = expand(w.xy);
	w.xz = expand(w.xz);
	w.yz = expand(w.yz);
	w.xt = expand(w.xt);
	w.yt = expand(w.yt);
	w.zt = expand(w.zt);
	w.xyz = expand(w.xyz);
	w.xyt = expand(w.xyt);
	w.xzt = expand(w.xzt);
	w.yzt = expand(w.yzt);
	w.xyzt = expand(w.xyzt);

	return w;
}

//////////////////////////////////////////////////////

Mink AntiWedge(const Mink a, const Mink b)
{
	Mink c;

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

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);

	return c;
}

//////////////////////////////////////////////////////


Mink Regressive(Mink a, Mink b) 
{
//	Mink I.xyzt = 1;	Mink I_inv.xyzt = -1;

//	c = ((a*I_inv)^(b*I_inv))*I;

	Mink c;
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

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);

	return c;
}

//////////////////////////////////////////////////////

Mink RegressiveViaFormula(Mink a, Mink b)
{
	Mink c;

	Mink I,I_inv,d,e,f;
	
	I.xyzt = 1;	I_inv.xyzt = -1;

	d = a*I_inv^b*I_inv;
	c = d*I;

	return c;
}

//////////////////////////////////////////////////////

Mink LowerRightViaFormula(Mink a, Mink b)  // Like a mirror of Wedge along rising diagonal
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

	Mink c;

	Mink I,I_inv,d,e,f;
	
	I.xyzt = 1;	I_inv.xyzt = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

Mink Expander(const Mink &a, const Mink &b)
{
/*

Terms with increased rank
 <    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
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
	Mink c;

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


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);


	return c;
}

//////////////////////////////////////////////////////

Mink Conserver(const Mink &a, const Mink &b)
{
/*

Terms with preserved rank
 <    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
 x    |  x     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 t    |  t     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0     0    -yz    xz   -yt    xt    0     0     0    -yzt   xzt   0    
 xz   |  xz    0     0     0     0     yz    0    -xy   -zt    0     xt    0     yzt   0    -xyt   0    
 yz   |  yz    0     0     0     0    -xz    xy    0     0    -zt    yt    0    -xzt   xyt   0     0    
 xt   |  xt    0     0     0     0     yt    zt    0     0     xy    xz   -yzt   0     0    -xyz   0    
 yt   |  yt    0     0     0     0    -xt    0     zt   -xy    0     yz    xzt   0     xyz   0     0    
 zt   |  zt    0     0     0     0     0    -xt   -yt   -xz   -yz    0    -xyt  -xyz   0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0     yzt  -xzt   xyt   0     0     0     0     0    
 xyt  |  xyt   0     0     0     0     0    -yzt   xzt   0     0     xyz   0     0     0     0     0    
 xzt  |  xzt   0     0     0     0     yzt   0    -xyt   0    -xyz   0     0     0     0     0     0    
 yzt  |  yzt   0     0     0     0    -xzt   xyt   0     xyz   0     0     0     0     0     0     0    
 xyzt |  xyzt  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    

*/

	Mink c;


// Conserver equation set for test purposes
c.q    =  + a.q   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   + a.xt  *b.yt   - a.yt  *b.xt   ; 
c.xz   =  + a.q   *b.xz   + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   + a.xt  *b.zt   - a.zt  *b.xt   ; 
c.yz   =  + a.q   *b.yz   - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    + a.yt  *b.zt   - a.zt  *b.yt   ; 
c.xt   =  + a.q   *b.xt   + a.xy  *b.yt   + a.xz  *b.zt   + a.xt  *b.q    - a.yt  *b.xy   - a.zt  *b.xz   ; 
c.yt   =  + a.q   *b.yt   - a.xy  *b.xt   + a.yz  *b.zt   + a.xt  *b.xy   + a.yt  *b.q    - a.zt  *b.yz   ; 
c.zt   =  + a.q   *b.zt   - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   + a.zt  *b.q    ; 
c.xyz  =  + a.q   *b.xyz  - a.xt  *b.yzt  + a.yt  *b.xzt  - a.zt  *b.xyt  + a.xyz *b.q    + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   ; 
c.xyt  =  + a.q   *b.xyt  - a.xz  *b.yzt  + a.yz  *b.xzt  - a.zt  *b.xyz  + a.xyz *b.zt   + a.xyt *b.q    - a.xzt *b.yz   + a.yzt *b.xz   ; 
c.xzt  =  + a.q   *b.xzt  + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  - a.xyz *b.yt   + a.xyt *b.yz   + a.xzt *b.q    - a.yzt *b.xy   ; 
c.yzt  =  + a.q   *b.yzt  - a.xy  *b.xzt  + a.xz  *b.xyt  - a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.yzt *b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.xyzt*b.q    ; 


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);



	return c;
}

//////////////////////////////////////////////////////

Mink Shrinker(const Mink &a, const Mink &b)
{

	Mink c;

/*


Terms with reduced rank
 <    |  q     x     y     z     t     xy    xz    yz    xt    yt    zt    xyz   xyt   xzt   yzt   xyzt 
-------------------------------------------------------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     0     y     z     0     t     0     0     yz    yt    zt    0     yzt  
 y    |  0     0     q     0     0    -x     0     z     0     t     0    -xz   -xt    0     zt   -xzt  
 z    |  0     0     0     q     0     0    -x    -y     0     0     t     xy    0    -xt   -yt    xyt  
 t    |  0     0     0     0    -q     0     0     0     x     y     z     0    -xy   -xz   -yz    xyz  
 xy   |  0    -y     x     0     0    -q     0     0     0     0     0    -z    -t     0     0    -zt   
 xz   |  0    -z     0     x     0     0    -q     0     0     0     0     y     0    -t     0     yt   
 yz   |  0     0    -z     y     0     0     0    -q     0     0     0    -x     0     0    -t    -xt   
 xt   |  0    -t     0     0    -x     0     0     0     q     0     0     0    -y    -z     0     yz   
 yt   |  0     0    -t     0    -y     0     0     0     0     q     0     0     x     0    -z    -xz   
 zt   |  0     0     0    -t    -z     0     0     0     0     0     q     0     0     x     y     xy   
 xyz  |  0     yz   -xz    xy    0    -z     y    -x     0     0     0    -q    -zt    yt   -xt   -t    
 xyt  |  0     yt   -xt    0    -xy   -t     0     0    -y     x     0     zt    q     yz   -xz   -z    
 xzt  |  0     zt    0    -xt   -xz    0    -t     0    -z     0     x    -yt   -yz    q     xy    y    
 yzt  |  0     0     zt   -yt   -yz    0     0    -t     0    -z     y     xt    xz   -xy    q    -x    
 xyzt |  0    -yzt   xzt  -xyt  -xyz  -zt    yt   -xt    yz   -xz    xy    t     z    -y     x    -q    


*/


// Shrinker equation set for test purposes
c.q    =  + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    - a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   + a.xt  *b.xt   + a.yt  *b.yt   + a.zt  *b.zt   - a.xyz *b.xyz  + a.xyt *b.xyt  + a.xzt *b.xzt  + a.yzt *b.yzt  - a.xyzt*b.xyzt ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   + a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  - a.xt  *b.t    + a.yt  *b.xyt  + a.zt  *b.xzt  - a.xyz *b.yz   + a.xyt *b.yt   + a.xzt *b.zt   - a.yzt *b.xyzt + a.xyzt*b.yzt  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   + a.t   *b.yt   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    - a.xt  *b.xyt  - a.yt  *b.t    + a.zt  *b.yzt  + a.xyz *b.xz   - a.xyt *b.xt   + a.xzt *b.xyzt + a.yzt *b.zt   - a.xyzt*b.xzt  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   + a.t   *b.zt   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    - a.xt  *b.xzt  - a.yt  *b.yzt  - a.zt  *b.t    - a.xyz *b.xy   - a.xyt *b.xyzt - a.xzt *b.xt   - a.yzt *b.yt   + a.xyzt*b.xyt  ; 
c.t    =  + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   + a.xyzt*b.xyz  ; 
c.xy   =  + a.z   *b.xyz  - a.t   *b.xyt  + a.zt  *b.xyzt + a.xyz *b.z    - a.xyt *b.t    + a.xzt *b.yzt  - a.yzt *b.xzt  + a.xyzt*b.zt   ; 
c.xz   =  - a.y   *b.xyz  - a.t   *b.xzt  - a.yt  *b.xyzt - a.xyz *b.y    - a.xyt *b.yzt  - a.xzt *b.t    + a.yzt *b.xyt  - a.xyzt*b.yt   ; 
c.yz   =  + a.x   *b.xyz  - a.t   *b.yzt  + a.xt  *b.xyzt + a.xyz *b.x    + a.xyt *b.xzt  - a.xzt *b.xyt  - a.yzt *b.t    + a.xyzt*b.xt   ; 
c.xt   =  - a.y   *b.xyt  - a.z   *b.xzt  - a.yz  *b.xyzt - a.xyz *b.yzt  - a.xyt *b.y    - a.xzt *b.z    + a.yzt *b.xyz  - a.xyzt*b.yz   ; 
c.yt   =  + a.x   *b.xyt  - a.z   *b.yzt  + a.xz  *b.xyzt + a.xyz *b.xzt  + a.xyt *b.x    - a.xzt *b.xyz  - a.yzt *b.z    + a.xyzt*b.xz   ; 
c.zt   =  + a.x   *b.xzt  + a.y   *b.yzt  - a.xy  *b.xyzt - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy   ; 
c.xyz  =  + a.t   *b.xyzt - a.xyzt*b.t    ; 
c.xyt  =  + a.z   *b.xyzt - a.xyzt*b.z    ; 
c.xzt  =  - a.y   *b.xyzt + a.xyzt*b.y    ; 
c.yzt  =  + a.x   *b.xyzt - a.xyzt*b.x    ; 
c.xyzt = 0 ; 


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);


	return c;
}

///////////////////////////////////////////////////////

Mink Symmetric(const Mink &a, const Mink &b)
{

	Mink c;

//	c = (a*b + b*a)/2;

// Symmetric Product Equations

c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    - a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   + a.xt  *b.xt   + a.yt  *b.yt   + a.zt  *b.zt   - a.xyz *b.xyz  + a.xyt *b.xyt  + a.xzt *b.xzt  + a.yzt *b.yzt  - a.xyzt*b.xyzt ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.yz  *b.xyz  + a.yt  *b.xyt  + a.zt  *b.xzt  - a.xyz *b.yz   + a.xyt *b.yt   + a.xzt *b.zt   ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    + a.xz  *b.xyz  - a.xt  *b.xyt  + a.zt  *b.yzt  + a.xyz *b.xz   - a.xyt *b.xt   + a.yzt *b.zt   ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    - a.xy  *b.xyz  - a.xt  *b.xzt  - a.yt  *b.yzt  - a.xyz *b.xy   - a.xzt *b.xt   - a.yzt *b.yt   ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   ; 
c.xy   =  + a.q   *b.xy   + a.z   *b.xyz  - a.t   *b.xyt  + a.xy  *b.q    + a.zt  *b.xyzt + a.xyz *b.z    - a.xyt *b.t    + a.xyzt*b.zt   ; 
c.xz   =  + a.q   *b.xz   - a.y   *b.xyz  - a.t   *b.xzt  + a.xz  *b.q    - a.yt  *b.xyzt - a.xyz *b.y    - a.xzt *b.t    - a.xyzt*b.yt   ; 
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  - a.t   *b.yzt  + a.yz  *b.q    + a.xt  *b.xyzt + a.xyz *b.x    - a.yzt *b.t    + a.xyzt*b.xt   ; 
c.xt   =  + a.q   *b.xt   - a.y   *b.xyt  - a.z   *b.xzt  - a.yz  *b.xyzt + a.xt  *b.q    - a.xyt *b.y    - a.xzt *b.z    - a.xyzt*b.yz   ; 
c.yt   =  + a.q   *b.yt   + a.x   *b.xyt  - a.z   *b.yzt  + a.xz  *b.xyzt + a.yt  *b.q    + a.xyt *b.x    - a.yzt *b.z    + a.xyzt*b.xz   ; 
c.zt   =  + a.q   *b.zt   + a.x   *b.xzt  + a.y   *b.yzt  - a.xy  *b.xyzt + a.zt  *b.q    + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy   ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    + a.xyt *b.q    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    + a.xzt *b.q    ; 
c.yzt  =  + a.q   *b.yzt  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    + a.yzt *b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyzt*b.q    ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);


	return c;
}

//////////////////////////////////////////////////////

Mink AntiSymmetric(const Mink &a, const Mink &b)
{

	Mink c;

//	c = (a*b - b*a)/2;

// AntiSymmetric Product Equations

c.q    =  0 ; 
c.x    =  - a.y   *b.xy   - a.z   *b.xz   + a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    - a.xt  *b.t    - a.yzt *b.xyzt + a.xyzt*b.yzt  ; 
c.y    =  + a.x   *b.xy   - a.z   *b.yz   + a.t   *b.yt   - a.xy  *b.x    + a.yz  *b.z    - a.yt  *b.t    + a.xzt *b.xyzt - a.xyzt*b.xzt  ; 
c.z    =  + a.x   *b.xz   + a.y   *b.yz   + a.t   *b.zt   - a.xz  *b.x    - a.yz  *b.y    - a.zt  *b.t    - a.xyt *b.xyzt + a.xyzt*b.xyt  ; 
c.t    =  + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt + a.xyzt*b.xyz  ; 
c.xy   =  + a.x   *b.y    - a.y   *b.x    - a.xz  *b.yz   + a.yz  *b.xz   + a.xt  *b.yt   - a.yt  *b.xt   + a.xzt *b.yzt  - a.yzt *b.xzt  ; 
c.xz   =  + a.x   *b.z    - a.z   *b.x    + a.xy  *b.yz   - a.yz  *b.xy   + a.xt  *b.zt   - a.zt  *b.xt   - a.xyt *b.yzt  + a.yzt *b.xyt  ; 
c.yz   =  + a.y   *b.z    - a.z   *b.y    - a.xy  *b.xz   + a.xz  *b.xy   + a.yt  *b.zt   - a.zt  *b.yt   + a.xyt *b.xzt  - a.xzt *b.xyt  ; 
c.xt   =  + a.x   *b.t    - a.t   *b.x    + a.xy  *b.yt   + a.xz  *b.zt   - a.yt  *b.xy   - a.zt  *b.xz   - a.xyz *b.yzt  + a.yzt *b.xyz  ; 
c.yt   =  + a.y   *b.t    - a.t   *b.y    - a.xy  *b.xt   + a.yz  *b.zt   + a.xt  *b.xy   - a.zt  *b.yz   + a.xyz *b.xzt  - a.xzt *b.xyz  ; 
c.zt   =  + a.z   *b.t    - a.t   *b.z    - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   - a.xyz *b.xyt  + a.xyt *b.xyz  ; 
c.xyz  =  + a.t   *b.xyzt - a.xt  *b.yzt  + a.yt  *b.xzt  - a.zt  *b.xyt  + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   - a.xyzt*b.t    ; 
c.xyt  =  + a.z   *b.xyzt - a.xz  *b.yzt  + a.yz  *b.xzt  - a.zt  *b.xyz  + a.xyz *b.zt   - a.xzt *b.yz   + a.yzt *b.xz   - a.xyzt*b.z    ; 
c.xzt  =  - a.y   *b.xyzt + a.xy  *b.yzt  - a.yz  *b.xyt  + a.yt  *b.xyz  - a.xyz *b.yt   + a.xyt *b.yz   - a.yzt *b.xy   + a.xyzt*b.y    ; 
c.yzt  =  + a.x   *b.xyzt - a.xy  *b.xzt  + a.xz  *b.xyt  - a.xt  *b.xyz  + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   - a.xyzt*b.x    ; 
c.xyzt =  + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    ; 


	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);


	return c;
}

//////////////////////////////////////////////////////

Mink Inner(const Mink &a, const Mink &b)
{

	Mink c;

c.q    = -a.t*b.t + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt;

c.x    =  + (a.t*b.xt - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt);
c.y    =  + (a.t*b.yt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt);
c.z    =  + (a.t*b.zt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt - a.zt*b.t);
c.t    =  + (a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z);

c.xy   =  + (-a.t*b.xyt - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt + a.z*b.xyz + a.zt*b.xyzt);
c.xz   =  - ( a.t*b.xzt + a.xyz*b.y + a.xyzt*b.yt + a.xzt*b.t + a.y*b.xyz + a.yt*b.xyzt);
c.xt   =  - ( a.xyt*b.y + a.xyzt*b.yz + a.xzt*b.z + a.y*b.xyt + a.yz*b.xyzt + a.z*b.xzt);
c.yz   =  + (-a.t*b.yzt + a.x*b.xyz + a.xt*b.xyzt + a.xyz*b.x + a.xyzt*b.xt - a.yzt*b.t);
c.yt   =  + ( a.x*b.xyt + a.xyt*b.x + a.xyzt*b.xz + a.xz*b.xyzt - a.yzt*b.z - a.z*b.yzt);
c.zt   =  + ( a.x*b.xzt - a.xy*b.xyzt - a.xyzt*b.xy + a.xzt*b.x + a.y*b.yzt + a.yzt*b.y);

c.xyz  =  + ( a.t*b.xyzt - a.xyzt*b.t);
c.xyt  =  + (-a.xyzt*b.z + a.z*b.xyzt);
c.xzt  =  + ( a.xyzt*b.y - a.y*b.xyzt);
c.yzt  =  + ( a.x*b.xyzt - a.xyzt*b.x);

c.xyzt = 0;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);

	return c;
}

//////////////////////////////////////////////////////

Mink LeftContraction (const Mink &a, const Mink &b)
{

	Mink c;

c.q    = a.q*b.q - a.t*b.t + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt;

c.x    = + (a.q*b.x + a.t*b.xt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt);
c.y    = + (a.q*b.y + a.t*b.yt + a.x*b.xy - a.xt*b.xyt + a.xz*b.xyz + a.xzt*b.xyzt - a.z*b.yz + a.zt*b.yzt);
c.z    = + (a.q*b.z + a.t*b.zt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt + a.y*b.yz - a.yt*b.yzt);
c.t    = + (a.q*b.t + a.x*b.xt - a.xy*b.xyt - a.xyz*b.xyzt - a.xz*b.xzt + a.y*b.yt - a.yz*b.yzt + a.z*b.zt);

c.xy   = + (a.q*b.xy - a.t*b.xyt + a.z*b.xyz + a.zt*b.xyzt);
c.xz   = + (a.q*b.xz - a.t*b.xzt - a.y*b.xyz - a.yt*b.xyzt);
c.xt   = + (a.q*b.xt - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt);
c.yz   = + (a.q*b.yz - a.t*b.yzt + a.x*b.xyz + a.xt*b.xyzt);
c.yt   = + (a.q*b.yt + a.x*b.xyt + a.xz*b.xyzt - a.z*b.yzt);
c.zt   = + (a.q*b.zt + a.x*b.xzt - a.xy*b.xyzt + a.y*b.yzt);

c.xyz  = + (a.q*b.xyz + a.t*b.xyzt);
c.xyt  = + (a.q*b.xyt + a.z*b.xyzt);
c.xzt  = + (a.q*b.xzt - a.y*b.xyzt);
c.yzt  = + (a.q*b.yzt + a.x*b.xyzt);

c.xyzt = + a.q*b.xyzt;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);

	return c;
}

//////////////////////////////////////////////////////

Mink RightContraction (const Mink &a, const Mink &b)
{

	Mink c;

c.q   = a.q*b.q - a.t*b.t + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt;

c.x    = + (a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt);
c.y    = + (-a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt);
c.z    = + (-a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t);
c.t    = + (a.t*b.q - a.xt*b.x - a.xyt*b.xy + a.xyzt*b.xyz - a.xzt*b.xz - a.yt*b.y - a.yzt*b.yz - a.zt*b.z);

c.xy   = + (a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt);
c.xz   = + (-a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t);
c.xt   = + (a.xt*b.q - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z);
c.yz   = + (a.xyz*b.x + a.xyzt*b.xt + a.yz*b.q - a.yzt*b.t);
c.yt   = + (a.xyt*b.x + a.xyzt*b.xz + a.yt*b.q - a.yzt*b.z);
c.zt   = + (-a.xyzt*b.xy + a.xzt*b.x + a.yzt*b.y + a.zt*b.q);

c.xyz  = + (a.xyz*b.q - a.xyzt*b.t);
c.xyt  = + (a.xyt*b.q - a.xyzt*b.z);
c.xzt  = + (a.xyzt*b.y + a.xzt*b.q);
c.yzt  = + (-a.xyzt*b.x + a.yzt*b.q);

c.xyzt = + a.xyzt*b.q;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);
	c.xyzt = expand(c.xyzt);

	return c;
}

//////////////////////////////////////////////////////

ex Determinant(Mink A) {
/*

Using the reverse operator as an example, I want to find det(V).
\begin{verbatim}
V = Mink(A,  B, C, D, E,  F, G, H, J, K, L,  M, N, P, R,  S)
U = Mink(A,  B, C, D, E, -F,-G,-H,-J,-K,-L, -M,-N,-P,-R,  S)
W = U*V

a = W.q = A^2+B^2+C^2+D^2-E^2+F^2+G^2+H^2-J^2-K^2-L^2+M^2-N^2-P^2-R^2-S^2

b = W.x =  + 2*A*B + 2*C*F + 2*D*G - 2*E*J + 2*H*M - 2*K*N - 2*L*P - 2*R*S
c = W.y =  + 2*A*C - 2*B*F + 2*D*H - 2*E*K - 2*G*M + 2*J*N - 2*L*R + 2*P*S
d = W.z =  + 2*A*D - 2*B*G - 2*C*H - 2*E*L + 2*J*P + 2*K*R + 2*M*F - 2*N*S
e = W.t =  + 2*A*E - 2*B*J - 2*C*K - 2*D*L + 2*F*N + 2*G*P + 2*H*R - 2*M*S

s = W.xyzt =  + 2*A*S - 2*B*R + 2*C*P - 2*D*N + 2*E*M - 2*F*L + 2*G*K - 2*H*J

det(W) = det(V)*det(U) = (det(V))^2 = (a*a - b*b - c*c - d*d + e*e + s*s)^2
det(V) = a*a - b*b - c*c - d*d + e*e + s*s (sign verified)
\end{verbatim}

*/	
	Mink B, C;
	ex a, b, c, d, e, s, det;

	B = Reverse(A);
	C = A*B;

	a = C.q;
	b = C.x;
	c = C.y;
	d = C.z;
	e = C.t;

	s = C.xyzt;
	det = expand(a*a - b*b - c*c - d*d + e*e + s*s);

	return(det);
}

//////////////////////////////////////////////////////

Mink Adjugate(Mink V)
{
	Mink W;

//Product_Clifford(Reverse(r),Conjugate(Product_Clifford(r,Reverse(r)))) = Adjugate(r);

// V is the input multivector. Local variables used are

	ex a,  b,c,d,e,  f,g,h,j,k,l,  m,n,p,r,  s;

	a = V.q;
	b = V.x;   c = V.y;   d = V.z;   e = V.t;
	f = V.xy;  g = V.xz;  h = V.yz;  j = V.xt;  k = V.yt;  l = V.zt;
	m = V.xyz; n = V.xyt; p = V.xzt; r = V.yzt;
	s = V.xyzt;


W.q = 
 + a*( + a*a - b*b - c*c - d*d + e*e + f*f + g*g + h*h - j*j - k*k - l*l
                                          + m*m - n*n - p*p - r*r + s*s)
 + 2*( - b*h*m + b*k*n + b*l*p + c*g*m - c*j*n + c*l*r - d*f*m - d*j*p
                - d*k*r + e*f*n + e*g*p + e*h*r - f*l*s + g*k*s - h*j*s);
		
W.x = 
 + b*( - a*a + b*b + c*c + d*d - e*e - f*f - g*g + h*h + j*j - k*k - l*l
                                          + m*m - n*n - p*p + r*r - s*s)
 + 2*( - a*h*m + a*k*n + a*l*p - c*g*h + c*j*k - c*p*r + d*f*h + d*j*l
                + d*n*r - e*f*k - e*g*l - e*m*r + f*p*s - g*n*s + m*j*s);
		
W.y = 
 + c*( - a*a + b*b + c*c + d*d - e*e - f*f + g*g - h*h - j*j + k*k - l*l
                                          + m*m - n*n + p*p - r*r - s*s)
 + 2*( + a*g*m - a*j*n + a*r*l - b*g*h + b*j*k - b*p*r - d*f*g + d*k*l
               - d*n*p + e*f*j - e*h*l + e*m*p + f*r*s - h*n*s + m*k*s);
	       
W.z = 
 + d*( - a*a + b*b + c*c + d*d - e*e + f*f - g*g - h*h - j*j - k*k + l*l
                                          + m*m + n*n - p*p - r*r - s*s)
 + 2*( - a*f*m - a*j*p - a*k*r + b*f*h + b*j*l + b*n*r - c*f*g + c*k*l
               - c*n*p + e*g*j + e*h*k - e*m*n + g*r*s - h*p*s + l*m*s);
		
W.t = 
 + e*( - a*a + b*b + c*c + d*d - e*e + f*f + g*g + h*h + j*j + k*k + l*l
                                          - m*m - n*n - p*p - r*r - s*s)
 + 2*( - a*f*n - a*g*p - a*h*r + b*f*k + b*g*l + b*m*r - c*f*j + c*h*l
               - c*m*p - d*g*j - d*h*k + d*m*n + j*r*s - k*p*s + l*n*s);





W.xy = 
 + f*( - a*a + b*b + c*c - d*d + e*e - f*f - g*g - h*h + j*j + k*k - l*l
                                          + m*m - n*n + p*p + r*r + s*s)
 + 2*( + a*d*m - a*e*n + a*l*s - b*d*h + b*e*k - b*p*s + c*d*g - c*r*s
               - e*c*j + g*k*l - g*n*p - h*j*l - h*n*r + j*m*p + k*m*r);
	       
W.xz = 
 + g*( - a*a + b*b - c*c + d*d + e*e - f*f - g*g - h*h + j*j - k*k + l*l
                                          + m*m + n*n - p*p + r*r + s*s)
 + 2*( - a*c*m - a*e*p - a*k*s + b*c*h + b*e*l + b*n*s + c*d*f - d*e*j
               - d*r*s + f*k*l - f*n*p + h*j*k - h*p*r - j*m*n + l*m*r);
		
W.yz = 
 + h*( - a*a - b*b + c*c + d*d + e*e - f*f - g*g - h*h - j*j + k*k + l*l
                                          + m*m + n*n + p*p - r*r + s*s)
 + 2*( + a*b*m - a*e*r + a*j*s + b*c*g - b*d*f + c*e*l + c*n*s - d*e*k
               + d*p*s + g*j*k - g*p*r - f*j*l - f*n*r - k*m*n - m*l*p);
	       
W.xt = 
 + j*( - a*a + b*b - c*c - d*d - e*e - f*f - g*g + h*h + j*j + k*k + l*l
                                          - m*m - n*n - p*p + r*r + s*s) 
 + 2*( - a*c*n - a*d*p - a*h*s + b*c*k + b*d*l + b*m*s + c*e*f + d*e*g
               - e*r*s + f*h*l - f*m*p - g*h*k + g*m*n - k*p*r + l*n*r);
	       
W.yt = 
 + k*( - a*a - b*b + c*c - d*d - e*e - f*f + g*g - h*h + j*j + k*k + l*l
                                          - m*m - n*n + p*p - r*r + s*s)
 + 2*( + a*b*n - a*d*r + a*g*s + b*c*j - b*e*f + c*d*l + c*m*s + d*e*h
                + e*p*s - f*g*l - f*m*r - g*h*j + h*m*n - j*p*r - l*n*p);
		
W.zt = 
 + l*( - a*a - b*b - c*c + d*d - e*e + f*f - g*g - h*h + j*j + k*k + l*l
                                          - m*m + n*n - p*p - r*r + s*s)
 + 2*( + a*b*p + a*c*r - a*f*s + b*d*j - b*e*g + c*d*k - c*e*h + d*m*s
               - e*n*s - f*g*k + f*h*j - g*m*r + h*m*p + j*n*r - k*n*p);
	       
W.xyz = 
 + m*( - a*a - b*b - c*c - d*d - e*e + f*f + g*g + h*h + j*j + k*k + l*l
                                           - m*m + n*n + p*p + r*r - s*s)
 + 2*( + a*b*h - a*c*g + a*d*f + b*e*r - b*j*s - c*e*p - c*k*s + d*e*n
               - d*s*l + f*j*p + f*k*r - g*j*n + g*l*r - h*k*n - h*l*p);


	       
W.xyt = 
 + n*( - a*a - b*b - c*c + d*d + e*e + f*f - g*g - h*h - j*j - k*k + l*l
                                          - m*m + n*n + p*p + r*r - s*s)
 + 2*( + a*b*k - a*c*j + a*e*f + b*d*r - b*g*s - c*d*p - c*h*s - d*e*m
               - e*l*s + f*g*p + f*h*r + g*j*m + h*k*m + j*l*r - k*l*p);
	      
W.xzt = 
 + p*( - a*a - b*b + c*c - d*d + e*e - f*f + g*g - h*h - j*j + k*k - l*l
                                          - m*m + n*n + p*p + r*r - s*s)
 + 2*( + a*b*l - a*d*j + a*e*g - b*c*r + b*f*s - c*d*n + c*e*m - d*h*s
               + e*k*s + f*g*n - f*j*m + g*h*r + h*l*m - j*k*r - k*l*n);
	       
W.yzt = 
 + r*( - a*a + b*b - c*c - d*d + e*e - f*f - g*g + h*h + j*j - k*k - l*l
                                          - m*m + n*n + p*p + r*r - s*s)
 + 2*( + a*c*l - a*d*k + a*e*h - b*c*p + b*d*n - b*e*m + c*f*s + d*g*s
               - e*j*s + f*h*n - f*k*m + g*h*p - g*l*m - j*k*p + j*l*n);
	       
W.xyzt = 
 + s*( - a*a + b*b + c*c + d*d - e*e + f*f + g*g + h*h - j*j - k*k - l*l
                                          - m*m + n*n + p*p + r*r - s*s)
 + 2*( + a*f*l - a*g*k + a*h*j - b*f*p + b*g*n - b*j*m - c*f*r + c*h*n
              - c*k*m  - d*g*r + d*h*p - d*l*m + e*j*r - e*k*p + e*l*n);


	return W;
}

//////////////////////////////////////////////////////

Mink Reciprocal(Mink a)
{
	ex b;
	Mink c,d;

	b = Determinant(a);

	c = Adjugate(a);

	d.q = expand(c.q/b);
	d.x = expand(c.x/b);
	d.y = expand(c.y/b);
	d.z = expand(c.z/b);
	d.t = expand(c.t/b);
	d.xy = expand(c.xy/b);
	d.xz = expand(c.xz/b);
	d.yz = expand(c.yz/b);
	d.xt = expand(c.xt/b);
	d.yt = expand(c.yt/b);
	d.zt = expand(c.zt/b);
	d.xyz = expand(c.xyz/b);
	d.xyt = expand(c.xyt/b);
	d.xzt = expand(c.xzt/b);
	d.yzt = expand(c.yzt/b);
	d.xyzt = expand(c.xyzt/b);

	return d;
}


//////////////////////////////////////////////////////

Mink	Comp(Mink u, int i) {

	Mink MV;
	ex a, b,c,d,e, f,g,h,j,k,l, m,n,p,r, s;
	int Mask, Sign[16];

	Mask = 32;
	if((Mask & i) == 0) Sign[5] = +1; else Sign[5] = -1;  Mask >>= 1;  // q
	if((Mask & i) == 0) Sign[4] = +1; else Sign[4] = -1;  Mask >>= 1;  // x
	if((Mask & i) == 0) Sign[3] = +1; else Sign[3] = -1;  Mask >>= 1;  // y
	if((Mask & i) == 0) Sign[2] = +1; else Sign[2] = -1;  Mask >>= 1;  // z
	if((Mask & i) == 0) Sign[1] = +1; else Sign[1] = -1;  Mask >>= 1;  // t
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

	a = Sign[5]*u.q;   b = Sign[4]*u.x;   c = Sign[3]*u.y;   d = Sign[2]*u.z;   e = Sign[1]*u.t;
	f = Sign[0]*u.xy;  g = Sign[6]*u.xz;  j = Sign[7]*u.xt;  h = Sign[8]*u.yz;  k = Sign[9]*u.yt;  l = Sign[10]*u.zt;
	m = Sign[11]*u.xyz;  n = Sign[12]*u.xyt;  p = Sign[13]*u.xzt;   r = Sign[14]*u.yzt;  s = Sign[15]*u.xyzt;

	MV = Mink(a, b,c,d,e, f,g,h,j,k,l, m,n,p,r, s) ; 
	
	return MV;
}


//////////////////////////////////////////////////////

Mink Magic(const Mink &u, int i) {  // KN - sorted magic-square remappings

	ex  a,  b, c, d, e,  f, g, h, j, k, l,  m, n, p, r,  s;

	a = u.q;	b = u.x;	c = u.y;	d = u.z;	e = u.t;
	f = u.xy;	g = u.xz;	h = u.yz;	j = u.xt;	k = u.yt;	l = u.zt;
	m = u.xyz;	n = u.xyt;	p = u.xzt;	r = u.yzt;	s = u.xyzt;

	Mink w;
	w = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);

	if(i== 0) w = Mink( a,  b, c, d, e,  f, g, h, j, k, l,  m, n, p, r,  s) ; 
	if(i== 1) w = Mink( a,  b, c, d,-s,  f, g, h,-r, p,-n,  m, l,-k, j,  e) ; 
	if(i== 2) w = Mink( a,  b, d, c, e,  g, f,-h, j, l, k, -m, p, n,-r, -s) ; 
	if(i== 3) w = Mink( a,  b, d, c, s,  g, f,-h, r, n,-p, -m, k,-l, j,  e) ; 
	if(i== 4) w = Mink( a,  b, j, r, f,  e, s,-m, c, k,-p, -h, n,-l, d,  g) ; 
	if(i== 5) w = Mink( a,  b, j,-r, g,  e,-s, m, d, l,-n,  h, p,-k, c,  f) ; 
	if(i== 6) w = Mink( a,  b,-r, j, f, -s, e,-m, c, p, k, -h, l, n, d,  g) ; 
	if(i== 7) w = Mink( a,  b, r, j, g,  s, e, m, d, n, l,  h, k, p, c,  f) ; 
	if(i== 8) w = Mink( a,  c, b, d, e, -f, h, g, k, j, l, -m,-n, r, p, -s) ; 
	if(i== 9) w = Mink( a,  c, b, d, s, -f, h, g,-p, r, n, -m, l,-j, k,  e) ; 
	if(i==10) w = Mink( a,  c, d, b, e,  h,-f,-g, k, l, j,  m, r,-n,-p,  s) ; 
	if(i==11) w = Mink( a,  c, d, b,-s,  h,-f,-g, p,-n,-r,  m, j,-l, k,  e) ; 
	if(i==12) w = Mink( a,  c, k, p,-f,  e,-s, m, b, j,-r, -g,-n,-l, d,  h) ; 
	if(i==13) w = Mink( a,  c, k,-p, h,  e, s,-m, d, l, n,  g, r,-j, b, -f) ; 
	if(i==14) w = Mink( a,  c,-p, k,-f,  s, e, m, b, r, j, -g, l,-n, d,  h) ; 
	if(i==15) w = Mink( a,  c, p, k, h, -s, e,-m, d,-n, l,  g, j, r, b, -f) ; 
	if(i==16) w = Mink( a,  d, b, c, e, -g,-h, f, l, j, k,  m,-p,-r, n,  s) ; 
	if(i==17) w = Mink( a,  d, b, c,-s, -g,-h, f,-n,-r, p,  m, k,-j, l,  e) ; 
	if(i==18) w = Mink( a,  d, c, b, e, -h,-g,-f, l, k, j, -m,-r,-p,-n, -s) ; 
	if(i==19) w = Mink( a,  d, c, b, s, -h,-g,-f, n,-p, r, -m, j,-k, l,  e) ; 
	if(i==20) w = Mink( a,  d, l, n,-g,  e, s,-m, b, j, r, -f,-p,-k, c, -h) ; 
	if(i==21) w = Mink( a,  d, l,-n,-h,  e,-s, m, c, k, p,  f,-r,-j, b, -g) ; 
	if(i==22) w = Mink( a,  d,-n, l,-g, -s, e,-m, b,-r, j, -f, k,-p, c, -h) ; 
	if(i==23) w = Mink( a,  d, n, l,-h,  s, e, m, c,-p, k,  f, j,-r, b, -g) ; 
	if(i==24) w = Mink( a, -j, b,-r, f,  e,-m,-s,-k, c, p,  h, n, d, l, -g) ; 
	if(i==25) w = Mink( a, -j, b, r, g,  e, m, s,-l, d, n, -h, p, c, k, -f) ; 
	if(i==26) w = Mink( a, -j,-k,-l, e,  f, g, h, b, c, d, -s, n, p, r,  m) ; 
	if(i==27) w = Mink( a, -j, k,-l, m, -f, g,-h, r, p, n,  s, d, c, b,  e) ; 
	if(i==28) w = Mink( a, -j,-l,-k, e,  g, f,-h, b, d, c,  s, p, n,-r, -m) ; 
	if(i==29) w = Mink( a, -j, l,-k,-m, -g, f, h,-r, n, p, -s, c, d, b,  e) ; 
	if(i==30) w = Mink( a, -j,-r, b, f, -m, e, s,-k, p, c, -h, d, n,-l,  g) ; 
	if(i==31) w = Mink( a, -j, r, b, g,  m, e,-s,-l, n, d,  h, c, p,-k,  f) ; 
	if(i==32) w = Mink( a, -k, c,-p,-f,  e, m, s,-j, b, r,  g,-n, d, l, -h) ; 
	if(i==33) w = Mink( a, -k, c, p, h,  e,-m,-s,-l, d,-n, -g, r, b, j,  f) ; 
	if(i==34) w = Mink( a, -k,-j,-l, e, -f, h, g, c, b, d,  s,-n, r, p, -m) ; 
	if(i==35) w = Mink( a, -k, j,-l,-m,  f, h,-g, p, r,-n, -s, d, b, c,  e) ; 
	if(i==36) w = Mink( a, -k,-l,-j, e,  h,-f,-g, c, d, b, -s, r,-n,-p,  m) ; 
	if(i==37) w = Mink( a, -k, l,-j, m, -h,-f, g,-p,-n, r,  s, b, d, c,  e) ; 
	if(i==38) w = Mink( a, -k,-p, c,-f,  m, e,-s,-j, r, b, -g, d,-n,-l,  h) ; 
	if(i==39) w = Mink( a, -k, p, c, h, -m, e, s,-l,-n, d,  g, b, r,-j, -f) ; 
	if(i==40) w = Mink( a, -l, d,-n,-g,  e,-m,-s,-j, b,-r,  f,-p, c, k,  h) ; 
	if(i==41) w = Mink( a, -l, d, n,-h,  e, m, s,-k, c,-p, -f,-r, b, j,  g) ; 
	if(i==42) w = Mink( a, -l,-j,-k, e, -g,-h, f, d, b, c, -s,-p,-r, n,  m) ; 
	if(i==43) w = Mink( a, -l, j,-k, m,  g,-h,-f, n,-r,-p,  s, c, b, d,  e) ; 
	if(i==44) w = Mink( a, -l,-k,-j, e, -h,-g,-f, d, c, b,  s,-r,-p,-n, -m) ; 
	if(i==45) w = Mink( a, -l, k,-j,-m,  h,-g, f,-n,-p,-r, -s, b, c, d,  e) ; 
	if(i==46) w = Mink( a, -l,-n, d,-g, -m, e, s,-j,-r, b, -f, c,-p,-k, -h) ; 
	if(i==47) w = Mink( a, -l, n, d,-h,  m, e,-s,-k,-p, c,  f, b,-r,-j, -g) ; 
	if(i==48) w = Mink( a, -n, d, l,-g,  s,-m, e,-r, b, j,  f,-k, c,-p,  h) ; 
	if(i==49) w = Mink( a,  n, d, l,-h, -s, m, e,-p, c, k, -f,-j, b,-r,  g) ; 
	if(i==50) w = Mink( a,  n,-l, d,-g, -m,-s, e, r,-j, b, -f, c, k,-p, -h) ; 
	if(i==51) w = Mink( a, -n,-l, d,-h,  m, s, e, p,-k, c,  f, b, j,-r, -g) ; 
	if(i==52) w = Mink( a, -n,-p,-r,-m,  h,-g, f, l,-k, j,  e, b, c, d,  s) ; 
	if(i==53) w = Mink( a,  n,-p, r,-s, -h,-g,-f, d, c, b,  e,-j, k,-l, -m) ; 
	if(i==54) w = Mink( a,  n,-r,-p, m,  g,-h,-f, l,-j, k,  e, c, b, d, -s) ; 
	if(i==55) w = Mink( a, -n,-r, p, s, -g,-h, f, d, b, c,  e,-k, j,-l,  m) ; 
	if(i==56) w = Mink( a, -p, c, k,-f, -s, m, e, r, b, j,  g,-l, d,-n, -h) ; 
	if(i==57) w = Mink( a,  p, c, k, h,  s,-m, e,-n, d, l, -g,-j, b, r,  f) ; 
	if(i==58) w = Mink( a,  p,-k, c,-f,  m, s, e,-r,-j, b, -g, d, l,-n,  h) ; 
	if(i==59) w = Mink( a, -p,-k, c, h, -m,-s, e, n,-l, d,  g, b, j, r, -f) ; 
	if(i==60) w = Mink( a, -p,-n, r, m, -h,-f, g, k,-l, j,  e, b, d, c, -s) ; 
	if(i==61) w = Mink( a,  p,-n,-r, s,  h,-f,-g, c, d, b,  e,-j, l,-k,  m) ; 
	if(i==62) w = Mink( a,  p, r,-n,-m,  f, h,-g, k,-j, l,  e, d, b, c,  s) ; 
	if(i==63) w = Mink( a, -p, r, n,-s, -f, h, g, c, b, d,  e,-l, j,-k, -m) ; 
	if(i==64) w = Mink( a, -r, b, j, f,  s,-m, e, p, c, k,  h,-l, d, n, -g) ; 
	if(i==65) w = Mink( a,  r, b, j, g, -s, m, e, n, d, l, -h,-k, c, p, -f) ; 
	if(i==66) w = Mink( a,  r,-j, b, f, -m,-s, e,-p,-k, c, -h, d, l, n,  g) ; 
	if(i==67) w = Mink( a, -r,-j, b, g,  m, s, e,-n,-l, d,  h, c, k, p,  f) ; 
	if(i==68) w = Mink( a, -r, n, p,-m, -g, f, h, j,-l, k,  e, c, d, b,  s) ; 
	if(i==69) w = Mink( a,  r, n,-p,-s,  g, f,-h, b, d, c,  e,-k, l,-j, -m) ; 
	if(i==70) w = Mink( a,  r, p, n, m, -f, g,-h, j,-k, l,  e, d, c, b, -s) ; 
	if(i==71) w = Mink( a, -r, p,-n, s,  f, g, h, b, c, d,  e,-l, k,-j,  m) ; 

	return w;
}

////////////////////////////////////


