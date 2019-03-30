// Routines for Geometric Algebra in Five Dimensional metric(4,1) space
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

struct GA5_4_1{	
	ex q,   w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt,
		wxy,wxz,wxt,wyz,wyt,wzt,xyz,xyt,xzt,yzt,  
		wxyz,wxyt,wxzt,wyzt,xyzt,  wxyzt;
		
	GA5_4_1() { q = 0 ;    w = 0 ; x = 0 ; y = 0 ; z = 0 ; t = 0 ;  wx = 0 ; wy = 0 ; 
		wz = 0 ; wt = 0 ; xy = 0 ; xz = 0 ; xt = 0 ; yz = 0 ; yt = 0 ; zt = 0 ; 
		wxy = 0 ; wxz = 0 ; wxt = 0 ; wyz = 0 ; wyt = 0 ; wzt = 0 ; xyz = 0 ; xyt = 0 ; xzt = 0 ; yzt = 0 ;   
		wxyz = 0 ; wxyt = 0 ; wxzt = 0 ; wyzt = 0 ; xyzt = 0 ;   wxyzt = 0 ; }
	GA5_4_1(   ex q1, ex w1, ex x1, ex y1, ex z1, ex t1, ex wx1, ex wy1, ex wz1, ex wt1, ex xy1, ex xz1,
		 ex xt1, ex yz1, ex yt1, ex zt1, ex wxy1, ex wxz1, ex wxt1, ex wyz1, ex wyt1, ex wzt1, 
		 ex xyz1, ex xyt1, ex xzt1, ex yzt1, ex wxyz1, ex wxyt1, ex wxzt1, ex wyzt1, ex xyzt1, ex wxyzt1) 
		{
			q = q1;	
			w = w1;	x = x1;	y = y1;	z = z1;	t = t1;
			wx = wx1; wy = wy1; wz = wz1; wt = wt1; xy = xy1;
			xz = xz1; xt = xt1; yz = yz1; yt = yt1; zt = zt1;
			wxy = wxy1; wxz = wxz1; wxt = wxt1; wyz = wyz1; wyt = wyt1;
			wzt = wzt1; xyz = xyz1; xyt = xyt1; xzt = xzt1; yzt = yzt1;
			wxyz = wxyz1; wxyt = wxyt1; wxzt = wxzt1; wyzt = wyzt1; xyzt = xyzt1;
			wxyzt = wxyzt1;
		}
	};


//////////////////////////////////////////////////////

// Necessary forward declarations

GA5_4_1 LeftContraction (const GA5_4_1 &a, const GA5_4_1 &b) ;
GA5_4_1 Product(const GA5_4_1 &a, const GA5_4_1 &b) ;

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, GA5_4_1 &v) {
	return ff << "(" 
		<< v.q << ", " 

		<< v.w << "," 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.t << ", " 

		<< v.wx << "," 
		<< v.wy << "," 
		<< v.wz << "," 
		<< v.wt << "," 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.xt << "," 
		<< v.yz << "," 
		<< v.yt << "," 
		<< v.zt << ", " 

		<< v.wxy << ","
		<< v.wxz << ","
		<< v.wxt << ","
		<< v.wyz << ","
		<< v.wyt << ","
		<< v.wzt << ","
		<< v.xyz << ","
		<< v.xyt << ","
		<< v.xzt << ","
		<< v.yzt << ", "

		<< v.wxyz << ","
		<< v.wxyt << ","
		<< v.wxzt << ","
		<< v.wyzt << ","
		<< v.xyzt << ", "


		<< v.wxyzt << ")";
}


//////////////////////////////////////////////////////

void PrintMV(GA5_4_1 &v) {
	cout << "(" 
		<< v.q << ", \n" 

		<< v.w << "," 
		<< v.x << "," 
		<< v.y << "," 
		<< v.z << "," 
		<< v.t << ", \n" 

		<< v.wx << "," 
		<< v.wy << "," 
		<< v.wz << "," 
		<< v.wt << "," 
		<< v.xy << "," 
		<< v.xz << "," 
		<< v.xt << "," 
		<< v.yz << "," 
		<< v.yt << "," 
		<< v.zt << ", \n" 

		<< v.wxy << ","
		<< v.wxz << ","
		<< v.wxt << ","
		<< v.wyz << ","
		<< v.wyt << ","
		<< v.wzt << ","
		<< v.xyz << ","
		<< v.xyt << ","
		<< v.xzt << ","
		<< v.yzt << ", \n"

		<< v.wxyz << ","
		<< v.wxyt << ","
		<< v.wxzt << ","
		<< v.wyzt << ","
		<< v.xyzt << ", "


		<< v.wxyzt << ")\n";
}

//////////////////////////////////////////////////////

GA5_4_1 OverBar(GA5_4_1 a)
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA5_4_1 b;		// a blade wedge b blade = pseudovector blade wxyzt
// a^OverBar(a) = pseudovector blade wxyzt
// OverBar(Blade[i]) = Reciprocal(Blade[i])*wxyzt ;

	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t =  a.wxyz;

	b.wx =  a.yzt; 	b.wy = -a.xzt; 	b.wz =  a.xyt; 	b.wt = -a.xyz; 	b.xy =  a.wzt;
	b.xz = -a.wyt; 	b.xt =  a.wyz; 	b.yz =  a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy =  a.zt; 	b.wxz = -a.yt; 	b.wxt =  a.yz; 	b.wyz =  a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz = -a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt = -a.z; 	b.wxzt =  a.y; 	b.wyzt = -a.x; 	b.xyzt =  a.w;

	b.wxyzt = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 UnderBar(GA5_4_1 a)   
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA5_4_1 b;		// b blade wedge a blade = pseudovector blade

// UnderBar(Blade[i]) = wxyzt*Reciprocal(Blade[i]) ;

	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t =  a.wxyz;

	b.wx =  a.yzt; 	b.wy = -a.xzt; 	b.wz =  a.xyt; 	b.wt = -a.xyz; 	b.xy =  a.wzt;
	b.xz = -a.wyt; 	b.xt =  a.wyz; 	b.yz =  a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy =  a.zt; 	b.wxz = -a.yt; 	b.wxt =  a.yz; 	b.wyz =  a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz = -a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt = -a.z; 	b.wxzt =  a.y; 	b.wyzt = -a.x; 	b.xyzt =  a.w;

	b.wxyzt = a.q;

	return b;
}  


//////////////////////////////////////////////////////

GA5_4_1 Reverse(const GA5_4_1 &a)	// change sign based upon multivector product reversal.
// Preserves determinant faithfully
// Ex:  xy -> -yx, xyz -> -zyx, xyzt -> tzyx, wxyzt => tzyxw
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w =  a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = a.t;
	b.wx = - a.wx;		b.wy = - a.wy;		b.wz = - a.wz;		b.wt = - a.wt;		b.xy = - a.xy;
	b.xz = - a.xz;		b.yz = - a.yz;		b.xt = - a.xt;		b.yt = - a.yt;		b.zt = - a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = - a.xyz;	b.xyt = - a.xyt;	b.xzt = - a.xzt;	b.yzt = - a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt =  a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Reverse_Vector_Quad(const GA5_4_1 &a)	// change sign based upon multivector product reversal.
// only useful with 20 zero bivector and trivector elements. Otherwise, trashes determinant
// Ex:  xy -> -yx, xyz -> -zyx, xyzt -> tzyx, wxyzt => tzyxw
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = a.wxt;		b.wyz = a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt =  a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Transpose(const GA5_4_1 &a)	// verified
//  Determinant unchanged, as expected
//	+ a*q
//	- b*w + c*x + d*y + e*z - f*t
//	+ g*wx + h*wy + j*wz - k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	+ S*wxy + R*wxz - P*wxt + N*wyz - M*wyt - L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	- F*wxyz + E*wxyt + D*wxzt + C*wyzt - B*xyzt
//	+ A*wxyzt 
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = - a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = - a.wt;		b.xy = - a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = -a.wxt;		b.wyz = a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Conjugation(const GA5_4_1 &a)	// KN - destroys determinant. useless.
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t =  -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt =  -a.wt;		b.xy =  -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 CliffordConjugation(const GA5_4_1 &a)	
//(a, -b,-c,-d,-e,-f, -g,-h,-j,-k,-l, -m,-n,-p,-r,-s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, -A)
// Conjugates complex determinant. Magnitude matches.
{
	GA5_4_1 b;

	b.q =  a.q;

	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t =  -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = -a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;

	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = a.wxt;		b.wyz = a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;

	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Dual(const GA5_4_1 &a)	// dual =  a*I_inv 
{

//	GA5_4_1 I_inv,b;	// I_inv created initially zero
//	I_inv.wxyzt = -1;	// pure pseudoscalar
//	b = Product(a,I_inv);	// form dual

	GA5_4_1 b;

//Dual(MV1)  = (A, B,-C,D,-E,-F, -G,H,-J,-K,-L,  M,N,-P,-R,S, -s,r,p,-n,-m,  l,k,j,-h,g, f,e,-d,c,-b, -a)

	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t = -a.wxyz;

	b.wx = -a.yzt; 	b.wy =  a.xzt; 	b.wz = -a.xyt; 	b.wt = -a.xyz; 	b.xy = -a.wzt;
	b.xz =  a.wyt; 	b.xt =  a.wyz; 	b.yz = -a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy = -a.zt; 	b.wxz =  a.yt; 	b.wxt =  a.yz; 	b.wyz = -a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz =  a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt =  a.z; 	b.wxzt = -a.y; 	b.wyzt =  a.x; 	b.xyzt = -a.w;

	b.wxyzt =-a.q;


	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 DorstDual(GA5_4_1 a)
{
	GA5_4_1 b, I_inv;

	I_inv.wxyzt =-1;

	b = LeftContraction(a,I_inv);
	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 DorstUnDual(GA5_4_1 a)
{
	GA5_4_1 b, I;

	I.wxyzt = 1;

	b = LeftContraction(a,I);
	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 operator+(const GA5_4_1 &u, const GA5_4_1 &v)
{
	GA5_4_1 w;
	w.q =   u.q   + v.q  ;

	w.w   = u.w   + v.w  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.t   = u.t   + v.t  ;

	w.wx  = u.wx  + v.wx ;
	w.wy  = u.wy  + v.wy ;
	w.wz  = u.wz  + v.wz ;
	w.wt  = u.wt  + v.wt ;
	w.xy  = u.xy  + v.xy ;
	w.xz  = u.xz  + v.xz ;
	w.yz  = u.yz  + v.yz ;
	w.xt  = u.xt  + v.xt ;
	w.yt  = u.yt  + v.yt ;
	w.zt  = u.zt  + v.zt ;

	w.wxy = u.wxy + v.wxy;
	w.wxz = u.wxz + v.wxz;
	w.wxt = u.wxt + v.wxt;
	w.wyz = u.wyz + v.wyz;
	w.wyt = u.wyt + v.wyt;
	w.wzt = u.wzt + v.wzt;
	w.xyz = u.xyz + v.xyz;
	w.xyt = u.xyt + v.xyt;
	w.xzt = u.xzt + v.xzt;
	w.yzt = u.yzt + v.yzt;

	w.wxyz = u.wxyz + v.wxyz;
	w.wxyt = u.wxyt + v.wxyt;
	w.wxzt = u.wxzt + v.wxzt;
	w.wyzt = u.wyzt + v.wyzt;
	w.xyzt = u.xyzt + v.xyzt;

	w.wxyzt = u.wxyzt + v.wxyzt;

	return w;
}

//////////////////////////////////////////////////////

GA5_4_1 operator-(const GA5_4_1 &u, const GA5_4_1 &v)
{
	GA5_4_1 w;
	w.q =   u.q   - v.q  ;

	w.w   = u.w   - v.w  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.t   = u.t   - v.t  ;

	w.wx  = u.wx  - v.wx ;
	w.wy  = u.wy  - v.wy ;
	w.wz  = u.wz  - v.wz ;
	w.wt  = u.wt  - v.wt ;
	w.xy  = u.xy  - v.xy ;
	w.xz  = u.xz  - v.xz ;
	w.yz  = u.yz  - v.yz ;
	w.xt  = u.xt  - v.xt ;
	w.yt  = u.yt  - v.yt ;
	w.zt  = u.zt  - v.zt ;

	w.wxy = u.wxy - v.wxy;
	w.wxz = u.wxz - v.wxz;
	w.wxt = u.wxt - v.wxt;
	w.wyz = u.wyz - v.wyz;
	w.wyt = u.wyt - v.wyt;
	w.wzt = u.wzt - v.wzt;
	w.xyz = u.xyz - v.xyz;
	w.xyt = u.xyt - v.xyt;
	w.xzt = u.xzt - v.xzt;
	w.yzt = u.yzt - v.yzt;

	w.wxyz = u.wxyz - v.wxyz;
	w.wxyt = u.wxyt - v.wxyt;
	w.wxzt = u.wxzt - v.wxzt;
	w.wyzt = u.wyzt - v.wyzt;
	w.xyzt = u.xyzt - v.xyzt;

	w.wxyzt = u.wxyzt - v.wxyzt;

	return w;
}

//////////////////////////////////////////////////////

int operator==(const GA5_4_1 &u, const GA5_4_1 &v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.w ==v.w )&&
			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&
			(u.t ==v.t )&&

			(u.wx==v.wx)&&
			(u.wy==v.wy)&&
			(u.wz==v.wz)&&
			(u.wt==v.wt)&&
			(u.xy==v.xy)&&
			(u.xz==v.xz)&&
			(u.yz==v.yz)&&
			(u.xt==v.xt)&&
			(u.yt==v.yt)&&
			(u.zt==v.zt)&&

			(u.wxy==v.wxy)&&
			(u.wxz==v.wxz)&&
			(u.wxt==v.wxt)&&
			(u.wyz==v.wyz)&&
			(u.wyt==v.wyt)&&
			(u.wzt==v.wzt)&&
			(u.xyz==v.xyz)&&
			(u.xyt==v.xyt)&&
			(u.xzt==v.xzt)&&
			(u.yzt==v.yzt)&&

			(u.wxyz==v.wxyz)&&
			(u.wxyt==v.wxyt)&&
			(u.wxzt==v.wxzt)&&
			(u.wyzt==v.wyzt)&&
			(u.xyzt==v.xyzt)&&

			(u.wxyzt==v.wxyzt);
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const GA5_4_1 &u, const GA5_4_1 &v)
{
	int result;

	result = 	(u.q !=v.q )||

			(u.w !=v.w )||
			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||
			(u.t !=v.t )||

			(u.wx!=v.wx)||
			(u.wy!=v.wy)||
			(u.wz!=v.wz)||
			(u.wt!=v.wt)||
			(u.xy!=v.xy)||
			(u.xz!=v.xz)||
			(u.yz!=v.yz)||
			(u.xt!=v.xt)||
			(u.yt!=v.yt)||
			(u.zt!=v.zt)||

			(u.wxy!=v.wxy)||
			(u.wxz!=v.wxz)||
			(u.wxt!=v.wxt)||
			(u.wyz!=v.wyz)||
			(u.wyt!=v.wyt)||
			(u.wzt!=v.wzt)||
			(u.xyz!=v.xyz)||
			(u.xyt!=v.xyt)||
			(u.xzt!=v.xzt)||
			(u.yzt!=v.yzt)||

			(u.wxyz!=v.wxyz)||
			(u.wxyt!=v.wxyt)||
			(u.wxzt!=v.wxzt)||
			(u.wyzt!=v.wyzt)||
			(u.xyzt!=v.xyzt)||

			(u.wxyzt!=v.wxyzt);
	return result;
}

//////////////////////////////////////////////////////

GA5_4_1 operator*(const GA5_4_1 &a, const GA5_4_1 &b) {
	GA5_4_1 c;


c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.q*b.w + a.t*b.wt + a.w*b.q - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.q*b.x + a.t*b.xt + a.w*b.wx - a.wt*b.wxt - a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.q*b.y + a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt - a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt + a.wxz*b.wxyz - a.wxzt*b.wxyzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.q*b.z + a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt + a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wy*b.wyz - a.wyt*b.wyzt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t) ; 

c.t = + (a.q*b.t + a.t*b.q + a.w*b.wt - a.wt*b.w - a.wx*b.wxt - a.wxt*b.wx - a.wxy*b.wxyt + a.wxyt*b.wxy + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz - a.wxz*b.wxzt + a.wxzt*b.wxz - a.wy*b.wyt - a.wyt*b.wy - a.wyz*b.wyzt + a.wyzt*b.wyz - a.wz*b.wzt - a.wzt*b.wz + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z) ; 



c.wx = + (a.q*b.wx - a.t*b.wxt + a.w*b.x + a.wt*b.xt + a.wx*b.q - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt - a.wy*b.xy + a.wyt*b.xyt - a.wyz*b.xyz - a.wyzt*b.xyzt - a.wz*b.xz + a.wzt*b.xzt - a.x*b.w - a.xt*b.wt + a.xy*b.wy - a.xyt*b.wyt + a.xyz*b.wyz + a.xyzt*b.wyzt + a.xz*b.wz - a.xzt*b.wzt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 

c.wy = + (a.q*b.wy - a.t*b.wyt + a.w*b.y + a.wt*b.yt + a.wx*b.xy - a.wxt*b.xyt - a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt + a.wxz*b.xyz + a.wxzt*b.xyzt + a.wy*b.q - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt - a.wz*b.yz + a.wzt*b.yzt - a.x*b.wxy - a.xt*b.wxyt - a.xy*b.wx + a.xyt*b.wxt - a.xyz*b.wxz - a.xyzt*b.wxzt + a.xz*b.wxyz - a.xzt*b.wxyzt - a.y*b.w - a.yt*b.wt + a.yz*b.wz - a.yzt*b.wzt + a.z*b.wyz + a.zt*b.wyzt) ; 

c.wz = + (a.q*b.wz - a.t*b.wzt + a.w*b.z + a.wt*b.zt + a.wx*b.xz - a.wxt*b.xzt - a.wxy*b.xyz - a.wxyt*b.xyzt - a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt + a.wy*b.yz - a.wyt*b.yzt - a.wyz*b.y - a.wyzt*b.yt + a.wz*b.q - a.wzt*b.t - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt + a.xyz*b.wxy + a.xyzt*b.wxyt - a.xz*b.wx + a.xzt*b.wxt - a.y*b.wyz - a.yt*b.wyzt - a.yz*b.wy + a.yzt*b.wyt - a.z*b.w - a.zt*b.wt) ; 

c.wt = + (a.q*b.wt - a.t*b.w + a.w*b.t + a.wt*b.q + a.wx*b.xt - a.wxt*b.x - a.wxy*b.xyt - a.wxyt*b.xy - a.wxyz*b.xyzt + a.wxyzt*b.xyz - a.wxz*b.xzt - a.wxzt*b.xz + a.wy*b.yt - a.wyt*b.y - a.wyz*b.yzt - a.wyzt*b.yz + a.wz*b.zt - a.wzt*b.z - a.x*b.wxt - a.xt*b.wx - a.xy*b.wxyt + a.xyt*b.wxy + a.xyz*b.wxyzt + a.xyzt*b.wxyz - a.xz*b.wxzt + a.xzt*b.wxz - a.y*b.wyt - a.yt*b.wy - a.yz*b.wyzt + a.yzt*b.wyz - a.z*b.wzt - a.zt*b.wz) ; 

c.xy = + (a.q*b.xy - a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt - a.wx*b.wy + a.wxt*b.wyt + a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt - a.wxz*b.wyz - a.wxzt*b.wyzt + a.wy*b.wx - a.wyt*b.wxt + a.wyz*b.wxz + a.wyzt*b.wxzt - a.wz*b.wxyz + a.wzt*b.wxyzt + a.x*b.y + a.xt*b.yt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt - a.xz*b.yz + a.xzt*b.yzt - a.y*b.x - a.yt*b.xt + a.yz*b.xz - a.yzt*b.xzt + a.z*b.xyz + a.zt*b.xyzt) ; 

c.xz = + (a.q*b.xz - a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt - a.wx*b.wz + a.wxt*b.wzt + a.wxy*b.wyz + a.wxyt*b.wyzt + a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.wyz*b.wxy - a.wyzt*b.wxyt + a.wz*b.wx - a.wzt*b.wxt + a.x*b.z + a.xt*b.zt + a.xy*b.yz - a.xyt*b.yzt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt - a.yz*b.xy + a.yzt*b.xyt - a.z*b.x - a.zt*b.xt) ; 

c.xt = + (a.q*b.xt - a.t*b.x + a.w*b.wxt + a.wt*b.wx - a.wx*b.wt + a.wxt*b.w + a.wxy*b.wyt + a.wxyt*b.wy + a.wxyz*b.wyzt - a.wxyzt*b.wyz + a.wxz*b.wzt + a.wxzt*b.wz + a.wy*b.wxyt - a.wyt*b.wxy - a.wyz*b.wxyzt - a.wyzt*b.wxyz + a.wz*b.wxzt - a.wzt*b.wxz + a.x*b.t + a.xt*b.q + a.xy*b.yt - a.xyt*b.y - a.xyz*b.yzt - a.xyzt*b.yz + a.xz*b.zt - a.xzt*b.z - a.y*b.xyt - a.yt*b.xy - a.yz*b.xyzt + a.yzt*b.xyz - a.z*b.xzt - a.zt*b.xz) ; 

c.yz = + (a.q*b.yz - a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt - a.wxy*b.wxz - a.wxyt*b.wxzt - a.wxyz*b.wx + a.wxyzt*b.wxt + a.wxz*b.wxy + a.wxzt*b.wxyt - a.wy*b.wz + a.wyt*b.wzt + a.wyz*b.w + a.wyzt*b.wt + a.wz*b.wy - a.wzt*b.wyt + a.x*b.xyz + a.xt*b.xyzt - a.xy*b.xz + a.xyt*b.xzt + a.xyz*b.x + a.xyzt*b.xt + a.xz*b.xy - a.xzt*b.xyt + a.y*b.z + a.yt*b.zt + a.yz*b.q - a.yzt*b.t - a.z*b.y - a.zt*b.yt) ; 

c.yt = + (a.q*b.yt - a.t*b.y + a.w*b.wyt + a.wt*b.wy - a.wx*b.wxyt + a.wxt*b.wxy - a.wxy*b.wxt - a.wxyt*b.wx - a.wxyz*b.wxzt + a.wxyzt*b.wxz + a.wxz*b.wxyzt + a.wxzt*b.wxyz - a.wy*b.wt + a.wyt*b.w + a.wyz*b.wzt + a.wyzt*b.wz + a.wz*b.wyzt - a.wzt*b.wyz + a.x*b.xyt + a.xt*b.xy - a.xy*b.xt + a.xyt*b.x + a.xyz*b.xzt + a.xyzt*b.xz + a.xz*b.xyzt - a.xzt*b.xyz + a.y*b.t + a.yt*b.q + a.yz*b.zt - a.yzt*b.z - a.z*b.yzt - a.zt*b.yz) ; 

c.zt = + (a.q*b.zt - a.t*b.z + a.w*b.wzt + a.wt*b.wz - a.wx*b.wxzt + a.wxt*b.wxz - a.wxy*b.wxyzt - a.wxyt*b.wxyz + a.wxyz*b.wxyt - a.wxyzt*b.wxy - a.wxz*b.wxt - a.wxzt*b.wx - a.wy*b.wyzt + a.wyt*b.wyz - a.wyz*b.wyt - a.wyzt*b.wy - a.wz*b.wt + a.wzt*b.w + a.x*b.xzt + a.xt*b.xz - a.xy*b.xyzt + a.xyt*b.xyz - a.xyz*b.xyt - a.xyzt*b.xy - a.xz*b.xt + a.xzt*b.x + a.y*b.yzt + a.yt*b.yz - a.yz*b.yt + a.yzt*b.y + a.z*b.t + a.zt*b.q) ; 



c.wxy = + (a.q*b.wxy + a.t*b.wxyt + a.w*b.xy - a.wt*b.xyt + a.wx*b.y + a.wxt*b.yt + a.wxy*b.q - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt - a.wxz*b.yz + a.wxzt*b.yzt - a.wy*b.x - a.wyt*b.xt + a.wyz*b.xz - a.wyzt*b.xzt + a.wz*b.xyz + a.wzt*b.xyzt - a.x*b.wy + a.xt*b.wyt + a.xy*b.w + a.xyt*b.wt - a.xyz*b.wz + a.xyzt*b.wzt - a.xz*b.wyz - a.xzt*b.wyzt + a.y*b.wx - a.yt*b.wxt + a.yz*b.wxz + a.yzt*b.wxzt - a.z*b.wxyz + a.zt*b.wxyzt) ; 

c.wxz = + (a.q*b.wxz + a.t*b.wxzt + a.w*b.xz - a.wt*b.xzt + a.wx*b.z + a.wxt*b.zt + a.wxy*b.yz - a.wxyt*b.yzt - a.wxyz*b.y - a.wxyzt*b.yt + a.wxz*b.q - a.wxzt*b.t - a.wy*b.xyz - a.wyt*b.xyzt - a.wyz*b.xy + a.wyzt*b.xyt - a.wz*b.x - a.wzt*b.xt - a.x*b.wz + a.xt*b.wzt + a.xy*b.wyz + a.xyt*b.wyzt + a.xyz*b.wy - a.xyzt*b.wyt + a.xz*b.w + a.xzt*b.wt + a.y*b.wxyz - a.yt*b.wxyzt - a.yz*b.wxy - a.yzt*b.wxyt + a.z*b.wx - a.zt*b.wxt) ; 

c.wxt = + (a.q*b.wxt + a.t*b.wx + a.w*b.xt - a.wt*b.x + a.wx*b.t + a.wxt*b.q + a.wxy*b.yt - a.wxyt*b.y - a.wxyz*b.yzt - a.wxyzt*b.yz + a.wxz*b.zt - a.wxzt*b.z - a.wy*b.xyt - a.wyt*b.xy - a.wyz*b.xyzt + a.wyzt*b.xyz - a.wz*b.xzt - a.wzt*b.xz - a.x*b.wt + a.xt*b.w + a.xy*b.wyt + a.xyt*b.wy + a.xyz*b.wyzt - a.xyzt*b.wyz + a.xz*b.wzt + a.xzt*b.wz + a.y*b.wxyt - a.yt*b.wxy - a.yz*b.wxyzt - a.yzt*b.wxyz + a.z*b.wxzt - a.zt*b.wxz) ; 

c.wyz = + (a.q*b.wyz + a.t*b.wyzt + a.w*b.yz - a.wt*b.yzt + a.wx*b.xyz + a.wxt*b.xyzt - a.wxy*b.xz + a.wxyt*b.xzt + a.wxyz*b.x + a.wxyzt*b.xt + a.wxz*b.xy - a.wxzt*b.xyt + a.wy*b.z + a.wyt*b.zt + a.wyz*b.q - a.wyzt*b.t - a.wz*b.y - a.wzt*b.yt - a.x*b.wxyz + a.xt*b.wxyzt - a.xy*b.wxz - a.xyt*b.wxzt - a.xyz*b.wx + a.xyzt*b.wxt + a.xz*b.wxy + a.xzt*b.wxyt - a.y*b.wz + a.yt*b.wzt + a.yz*b.w + a.yzt*b.wt + a.z*b.wy - a.zt*b.wyt) ; 

c.wyt = + (a.q*b.wyt + a.t*b.wy + a.w*b.yt - a.wt*b.y + a.wx*b.xyt + a.wxt*b.xy - a.wxy*b.xt + a.wxyt*b.x + a.wxyz*b.xzt + a.wxyzt*b.xz + a.wxz*b.xyzt - a.wxzt*b.xyz + a.wy*b.t + a.wyt*b.q + a.wyz*b.zt - a.wyzt*b.z - a.wz*b.yzt - a.wzt*b.yz - a.x*b.wxyt + a.xt*b.wxy - a.xy*b.wxt - a.xyt*b.wx - a.xyz*b.wxzt + a.xyzt*b.wxz + a.xz*b.wxyzt + a.xzt*b.wxyz - a.y*b.wt + a.yt*b.w + a.yz*b.wzt + a.yzt*b.wz + a.z*b.wyzt - a.zt*b.wyz) ; 

c.wzt = + (a.q*b.wzt + a.t*b.wz + a.w*b.zt - a.wt*b.z + a.wx*b.xzt + a.wxt*b.xz - a.wxy*b.xyzt + a.wxyt*b.xyz - a.wxyz*b.xyt - a.wxyzt*b.xy - a.wxz*b.xt + a.wxzt*b.x + a.wy*b.yzt + a.wyt*b.yz - a.wyz*b.yt + a.wyzt*b.y + a.wz*b.t + a.wzt*b.q - a.x*b.wxzt + a.xt*b.wxz - a.xy*b.wxyzt - a.xyt*b.wxyz + a.xyz*b.wxyt - a.xyzt*b.wxy - a.xz*b.wxt - a.xzt*b.wx - a.y*b.wyzt + a.yt*b.wyz - a.yz*b.wyt - a.yzt*b.wy - a.z*b.wt + a.zt*b.w) ; 

c.xyz = + (a.q*b.xyz + a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt - a.wx*b.wyz - a.wxt*b.wyzt + a.wxy*b.wz - a.wxyt*b.wzt - a.wxyz*b.w - a.wxyzt*b.wt - a.wxz*b.wy + a.wxzt*b.wyt + a.wy*b.wxz + a.wyt*b.wxzt + a.wyz*b.wx - a.wyzt*b.wxt - a.wz*b.wxy - a.wzt*b.wxyt + a.x*b.yz - a.xt*b.yzt + a.xy*b.z + a.xyt*b.zt + a.xyz*b.q - a.xyzt*b.t - a.xz*b.y - a.xzt*b.yt - a.y*b.xz + a.yt*b.xzt + a.yz*b.x + a.yzt*b.xt + a.z*b.xy - a.zt*b.xyt) ; 

c.xyt = + (a.q*b.xyt + a.t*b.xy + a.w*b.wxyt - a.wt*b.wxy - a.wx*b.wyt - a.wxt*b.wy + a.wxy*b.wt - a.wxyt*b.w - a.wxyz*b.wzt - a.wxyzt*b.wz - a.wxz*b.wyzt + a.wxzt*b.wyz + a.wy*b.wxt + a.wyt*b.wx + a.wyz*b.wxzt - a.wyzt*b.wxz - a.wz*b.wxyzt - a.wzt*b.wxyz + a.x*b.yt - a.xt*b.y + a.xy*b.t + a.xyt*b.q + a.xyz*b.zt - a.xyzt*b.z - a.xz*b.yzt - a.xzt*b.yz - a.y*b.xt + a.yt*b.x + a.yz*b.xzt + a.yzt*b.xz + a.z*b.xyzt - a.zt*b.xyz) ; 

c.xzt = + (a.q*b.xzt + a.t*b.xz + a.w*b.wxzt - a.wt*b.wxz - a.wx*b.wzt - a.wxt*b.wz + a.wxy*b.wyzt - a.wxyt*b.wyz + a.wxyz*b.wyt + a.wxyzt*b.wy + a.wxz*b.wt - a.wxzt*b.w + a.wy*b.wxyzt + a.wyt*b.wxyz - a.wyz*b.wxyt + a.wyzt*b.wxy + a.wz*b.wxt + a.wzt*b.wx + a.x*b.zt - a.xt*b.z + a.xy*b.yzt + a.xyt*b.yz - a.xyz*b.yt + a.xyzt*b.y + a.xz*b.t + a.xzt*b.q - a.y*b.xyzt + a.yt*b.xyz - a.yz*b.xyt - a.yzt*b.xy - a.z*b.xt + a.zt*b.x) ; 

c.yzt = + (a.q*b.yzt + a.t*b.yz + a.w*b.wyzt - a.wt*b.wyz - a.wx*b.wxyzt - a.wxt*b.wxyz - a.wxy*b.wxzt + a.wxyt*b.wxz - a.wxyz*b.wxt - a.wxyzt*b.wx + a.wxz*b.wxyt - a.wxzt*b.wxy - a.wy*b.wzt - a.wyt*b.wz + a.wyz*b.wt - a.wyzt*b.w + a.wz*b.wyt + a.wzt*b.wy + a.x*b.xyzt - a.xt*b.xyz - a.xy*b.xzt - a.xyt*b.xz + a.xyz*b.xt - a.xyzt*b.x + a.xz*b.xyt + a.xzt*b.xy + a.y*b.zt - a.yt*b.z + a.yz*b.t + a.yzt*b.q - a.z*b.yt + a.zt*b.y) ; 



c.wxyz = + (a.q*b.wxyz - a.t*b.wxyzt + a.w*b.xyz + a.wt*b.xyzt + a.wx*b.yz - a.wxt*b.yzt + a.wxy*b.z + a.wxyt*b.zt + a.wxyz*b.q - a.wxyzt*b.t - a.wxz*b.y - a.wxzt*b.yt - a.wy*b.xz + a.wyt*b.xzt + a.wyz*b.x + a.wyzt*b.xt + a.wz*b.xy - a.wzt*b.xyt - a.x*b.wyz - a.xt*b.wyzt + a.xy*b.wz - a.xyt*b.wzt - a.xyz*b.w - a.xyzt*b.wt - a.xz*b.wy + a.xzt*b.wyt + a.y*b.wxz + a.yt*b.wxzt + a.yz*b.wx - a.yzt*b.wxt - a.z*b.wxy - a.zt*b.wxyt) ; 

c.wxyt = + (a.q*b.wxyt - a.t*b.wxy + a.w*b.xyt + a.wt*b.xy + a.wx*b.yt - a.wxt*b.y + a.wxy*b.t + a.wxyt*b.q + a.wxyz*b.zt - a.wxyzt*b.z - a.wxz*b.yzt - a.wxzt*b.yz - a.wy*b.xt + a.wyt*b.x + a.wyz*b.xzt + a.wyzt*b.xz + a.wz*b.xyzt - a.wzt*b.xyz - a.x*b.wyt - a.xt*b.wy + a.xy*b.wt - a.xyt*b.w - a.xyz*b.wzt - a.xyzt*b.wz - a.xz*b.wyzt + a.xzt*b.wyz + a.y*b.wxt + a.yt*b.wx + a.yz*b.wxzt - a.yzt*b.wxz - a.z*b.wxyzt - a.zt*b.wxyz) ; 

c.wxzt = + (a.q*b.wxzt - a.t*b.wxz + a.w*b.xzt + a.wt*b.xz + a.wx*b.zt - a.wxt*b.z + a.wxy*b.yzt + a.wxyt*b.yz - a.wxyz*b.yt + a.wxyzt*b.y + a.wxz*b.t + a.wxzt*b.q - a.wy*b.xyzt + a.wyt*b.xyz - a.wyz*b.xyt - a.wyzt*b.xy - a.wz*b.xt + a.wzt*b.x - a.x*b.wzt - a.xt*b.wz + a.xy*b.wyzt - a.xyt*b.wyz + a.xyz*b.wyt + a.xyzt*b.wy + a.xz*b.wt - a.xzt*b.w + a.y*b.wxyzt + a.yt*b.wxyz - a.yz*b.wxyt + a.yzt*b.wxy + a.z*b.wxt + a.zt*b.wx) ; 

c.wyzt = + (a.q*b.wyzt - a.t*b.wyz + a.w*b.yzt + a.wt*b.yz + a.wx*b.xyzt - a.wxt*b.xyz - a.wxy*b.xzt - a.wxyt*b.xz + a.wxyz*b.xt - a.wxyzt*b.x + a.wxz*b.xyt + a.wxzt*b.xy + a.wy*b.zt - a.wyt*b.z + a.wyz*b.t + a.wyzt*b.q - a.wz*b.yt + a.wzt*b.y - a.x*b.wxyzt - a.xt*b.wxyz - a.xy*b.wxzt + a.xyt*b.wxz - a.xyz*b.wxt - a.xyzt*b.wx + a.xz*b.wxyt - a.xzt*b.wxy - a.y*b.wzt - a.yt*b.wz + a.yz*b.wt - a.yzt*b.w + a.z*b.wyt + a.zt*b.wy) ; 

c.xyzt = + (a.q*b.xyzt - a.t*b.xyz + a.w*b.wxyzt + a.wt*b.wxyz - a.wx*b.wyzt + a.wxt*b.wyz + a.wxy*b.wzt + a.wxyt*b.wz - a.wxyz*b.wt + a.wxyzt*b.w - a.wxz*b.wyt - a.wxzt*b.wy + a.wy*b.wxzt - a.wyt*b.wxz + a.wyz*b.wxt + a.wyzt*b.wx - a.wz*b.wxyt + a.wzt*b.wxy + a.x*b.yzt + a.xt*b.yz + a.xy*b.zt - a.xyt*b.z + a.xyz*b.t + a.xyzt*b.q - a.xz*b.yt + a.xzt*b.y - a.y*b.xzt - a.yt*b.xz + a.yz*b.xt - a.yzt*b.x + a.z*b.xyt + a.zt*b.xy) ; 


c.wxyzt = + (a.q*b.wxyzt + a.t*b.wxyz + a.w*b.xyzt - a.wt*b.xyz + a.wx*b.yzt + a.wxt*b.yz + a.wxy*b.zt - a.wxyt*b.z + a.wxyz*b.t + a.wxyzt*b.q - a.wxz*b.yt + a.wxzt*b.y - a.wy*b.xzt - a.wyt*b.xz + a.wyz*b.xt - a.wyzt*b.x + a.wz*b.xyt + a.wzt*b.xy - a.x*b.wyzt + a.xt*b.wyz + a.xy*b.wzt + a.xyt*b.wz - a.xyz*b.wt + a.xyzt*b.w - a.xz*b.wyt - a.xzt*b.wy + a.y*b.wxzt - a.yt*b.wxz + a.yz*b.wxt + a.yzt*b.wx - a.z*b.wxyt + a.zt*b.wxy) ; 


	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 Product(const GA5_4_1 &a, const GA5_4_1 &b) {
	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.q*b.w + a.t*b.wt + a.w*b.q - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.q*b.x + a.t*b.xt + a.w*b.wx - a.wt*b.wxt - a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.q*b.y + a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt - a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt + a.wxz*b.wxyz - a.wxzt*b.wxyzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.q*b.z + a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt + a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wy*b.wyz - a.wyt*b.wyzt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t) ; 

c.t = + (a.q*b.t + a.t*b.q + a.w*b.wt - a.wt*b.w - a.wx*b.wxt - a.wxt*b.wx - a.wxy*b.wxyt + a.wxyt*b.wxy + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz - a.wxz*b.wxzt + a.wxzt*b.wxz - a.wy*b.wyt - a.wyt*b.wy - a.wyz*b.wyzt + a.wyzt*b.wyz - a.wz*b.wzt - a.wzt*b.wz + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z) ; 



c.wx = + (a.q*b.wx - a.t*b.wxt + a.w*b.x + a.wt*b.xt + a.wx*b.q - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt - a.wy*b.xy + a.wyt*b.xyt - a.wyz*b.xyz - a.wyzt*b.xyzt - a.wz*b.xz + a.wzt*b.xzt - a.x*b.w - a.xt*b.wt + a.xy*b.wy - a.xyt*b.wyt + a.xyz*b.wyz + a.xyzt*b.wyzt + a.xz*b.wz - a.xzt*b.wzt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 

c.wy = + (a.q*b.wy - a.t*b.wyt + a.w*b.y + a.wt*b.yt + a.wx*b.xy - a.wxt*b.xyt - a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt + a.wxz*b.xyz + a.wxzt*b.xyzt + a.wy*b.q - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt - a.wz*b.yz + a.wzt*b.yzt - a.x*b.wxy - a.xt*b.wxyt - a.xy*b.wx + a.xyt*b.wxt - a.xyz*b.wxz - a.xyzt*b.wxzt + a.xz*b.wxyz - a.xzt*b.wxyzt - a.y*b.w - a.yt*b.wt + a.yz*b.wz - a.yzt*b.wzt + a.z*b.wyz + a.zt*b.wyzt) ; 

c.wz = + (a.q*b.wz - a.t*b.wzt + a.w*b.z + a.wt*b.zt + a.wx*b.xz - a.wxt*b.xzt - a.wxy*b.xyz - a.wxyt*b.xyzt - a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt + a.wy*b.yz - a.wyt*b.yzt - a.wyz*b.y - a.wyzt*b.yt + a.wz*b.q - a.wzt*b.t - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt + a.xyz*b.wxy + a.xyzt*b.wxyt - a.xz*b.wx + a.xzt*b.wxt - a.y*b.wyz - a.yt*b.wyzt - a.yz*b.wy + a.yzt*b.wyt - a.z*b.w - a.zt*b.wt) ; 

c.wt = + (a.q*b.wt - a.t*b.w + a.w*b.t + a.wt*b.q + a.wx*b.xt - a.wxt*b.x - a.wxy*b.xyt - a.wxyt*b.xy - a.wxyz*b.xyzt + a.wxyzt*b.xyz - a.wxz*b.xzt - a.wxzt*b.xz + a.wy*b.yt - a.wyt*b.y - a.wyz*b.yzt - a.wyzt*b.yz + a.wz*b.zt - a.wzt*b.z - a.x*b.wxt - a.xt*b.wx - a.xy*b.wxyt + a.xyt*b.wxy + a.xyz*b.wxyzt + a.xyzt*b.wxyz - a.xz*b.wxzt + a.xzt*b.wxz - a.y*b.wyt - a.yt*b.wy - a.yz*b.wyzt + a.yzt*b.wyz - a.z*b.wzt - a.zt*b.wz) ; 

c.xy = + (a.q*b.xy - a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt - a.wx*b.wy + a.wxt*b.wyt + a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt - a.wxz*b.wyz - a.wxzt*b.wyzt + a.wy*b.wx - a.wyt*b.wxt + a.wyz*b.wxz + a.wyzt*b.wxzt - a.wz*b.wxyz + a.wzt*b.wxyzt + a.x*b.y + a.xt*b.yt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt - a.xz*b.yz + a.xzt*b.yzt - a.y*b.x - a.yt*b.xt + a.yz*b.xz - a.yzt*b.xzt + a.z*b.xyz + a.zt*b.xyzt) ; 

c.xz = + (a.q*b.xz - a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt - a.wx*b.wz + a.wxt*b.wzt + a.wxy*b.wyz + a.wxyt*b.wyzt + a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.wyz*b.wxy - a.wyzt*b.wxyt + a.wz*b.wx - a.wzt*b.wxt + a.x*b.z + a.xt*b.zt + a.xy*b.yz - a.xyt*b.yzt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt - a.yz*b.xy + a.yzt*b.xyt - a.z*b.x - a.zt*b.xt) ; 

c.xt = + (a.q*b.xt - a.t*b.x + a.w*b.wxt + a.wt*b.wx - a.wx*b.wt + a.wxt*b.w + a.wxy*b.wyt + a.wxyt*b.wy + a.wxyz*b.wyzt - a.wxyzt*b.wyz + a.wxz*b.wzt + a.wxzt*b.wz + a.wy*b.wxyt - a.wyt*b.wxy - a.wyz*b.wxyzt - a.wyzt*b.wxyz + a.wz*b.wxzt - a.wzt*b.wxz + a.x*b.t + a.xt*b.q + a.xy*b.yt - a.xyt*b.y - a.xyz*b.yzt - a.xyzt*b.yz + a.xz*b.zt - a.xzt*b.z - a.y*b.xyt - a.yt*b.xy - a.yz*b.xyzt + a.yzt*b.xyz - a.z*b.xzt - a.zt*b.xz) ; 

c.yz = + (a.q*b.yz - a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt - a.wxy*b.wxz - a.wxyt*b.wxzt - a.wxyz*b.wx + a.wxyzt*b.wxt + a.wxz*b.wxy + a.wxzt*b.wxyt - a.wy*b.wz + a.wyt*b.wzt + a.wyz*b.w + a.wyzt*b.wt + a.wz*b.wy - a.wzt*b.wyt + a.x*b.xyz + a.xt*b.xyzt - a.xy*b.xz + a.xyt*b.xzt + a.xyz*b.x + a.xyzt*b.xt + a.xz*b.xy - a.xzt*b.xyt + a.y*b.z + a.yt*b.zt + a.yz*b.q - a.yzt*b.t - a.z*b.y - a.zt*b.yt) ; 

c.yt = + (a.q*b.yt - a.t*b.y + a.w*b.wyt + a.wt*b.wy - a.wx*b.wxyt + a.wxt*b.wxy - a.wxy*b.wxt - a.wxyt*b.wx - a.wxyz*b.wxzt + a.wxyzt*b.wxz + a.wxz*b.wxyzt + a.wxzt*b.wxyz - a.wy*b.wt + a.wyt*b.w + a.wyz*b.wzt + a.wyzt*b.wz + a.wz*b.wyzt - a.wzt*b.wyz + a.x*b.xyt + a.xt*b.xy - a.xy*b.xt + a.xyt*b.x + a.xyz*b.xzt + a.xyzt*b.xz + a.xz*b.xyzt - a.xzt*b.xyz + a.y*b.t + a.yt*b.q + a.yz*b.zt - a.yzt*b.z - a.z*b.yzt - a.zt*b.yz) ; 

c.zt = + (a.q*b.zt - a.t*b.z + a.w*b.wzt + a.wt*b.wz - a.wx*b.wxzt + a.wxt*b.wxz - a.wxy*b.wxyzt - a.wxyt*b.wxyz + a.wxyz*b.wxyt - a.wxyzt*b.wxy - a.wxz*b.wxt - a.wxzt*b.wx - a.wy*b.wyzt + a.wyt*b.wyz - a.wyz*b.wyt - a.wyzt*b.wy - a.wz*b.wt + a.wzt*b.w + a.x*b.xzt + a.xt*b.xz - a.xy*b.xyzt + a.xyt*b.xyz - a.xyz*b.xyt - a.xyzt*b.xy - a.xz*b.xt + a.xzt*b.x + a.y*b.yzt + a.yt*b.yz - a.yz*b.yt + a.yzt*b.y + a.z*b.t + a.zt*b.q) ; 



c.wxy = + (a.q*b.wxy + a.t*b.wxyt + a.w*b.xy - a.wt*b.xyt + a.wx*b.y + a.wxt*b.yt + a.wxy*b.q - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt - a.wxz*b.yz + a.wxzt*b.yzt - a.wy*b.x - a.wyt*b.xt + a.wyz*b.xz - a.wyzt*b.xzt + a.wz*b.xyz + a.wzt*b.xyzt - a.x*b.wy + a.xt*b.wyt + a.xy*b.w + a.xyt*b.wt - a.xyz*b.wz + a.xyzt*b.wzt - a.xz*b.wyz - a.xzt*b.wyzt + a.y*b.wx - a.yt*b.wxt + a.yz*b.wxz + a.yzt*b.wxzt - a.z*b.wxyz + a.zt*b.wxyzt) ; 

c.wxz = + (a.q*b.wxz + a.t*b.wxzt + a.w*b.xz - a.wt*b.xzt + a.wx*b.z + a.wxt*b.zt + a.wxy*b.yz - a.wxyt*b.yzt - a.wxyz*b.y - a.wxyzt*b.yt + a.wxz*b.q - a.wxzt*b.t - a.wy*b.xyz - a.wyt*b.xyzt - a.wyz*b.xy + a.wyzt*b.xyt - a.wz*b.x - a.wzt*b.xt - a.x*b.wz + a.xt*b.wzt + a.xy*b.wyz + a.xyt*b.wyzt + a.xyz*b.wy - a.xyzt*b.wyt + a.xz*b.w + a.xzt*b.wt + a.y*b.wxyz - a.yt*b.wxyzt - a.yz*b.wxy - a.yzt*b.wxyt + a.z*b.wx - a.zt*b.wxt) ; 

c.wxt = + (a.q*b.wxt + a.t*b.wx + a.w*b.xt - a.wt*b.x + a.wx*b.t + a.wxt*b.q + a.wxy*b.yt - a.wxyt*b.y - a.wxyz*b.yzt - a.wxyzt*b.yz + a.wxz*b.zt - a.wxzt*b.z - a.wy*b.xyt - a.wyt*b.xy - a.wyz*b.xyzt + a.wyzt*b.xyz - a.wz*b.xzt - a.wzt*b.xz - a.x*b.wt + a.xt*b.w + a.xy*b.wyt + a.xyt*b.wy + a.xyz*b.wyzt - a.xyzt*b.wyz + a.xz*b.wzt + a.xzt*b.wz + a.y*b.wxyt - a.yt*b.wxy - a.yz*b.wxyzt - a.yzt*b.wxyz + a.z*b.wxzt - a.zt*b.wxz) ; 

c.wyz = + (a.q*b.wyz + a.t*b.wyzt + a.w*b.yz - a.wt*b.yzt + a.wx*b.xyz + a.wxt*b.xyzt - a.wxy*b.xz + a.wxyt*b.xzt + a.wxyz*b.x + a.wxyzt*b.xt + a.wxz*b.xy - a.wxzt*b.xyt + a.wy*b.z + a.wyt*b.zt + a.wyz*b.q - a.wyzt*b.t - a.wz*b.y - a.wzt*b.yt - a.x*b.wxyz + a.xt*b.wxyzt - a.xy*b.wxz - a.xyt*b.wxzt - a.xyz*b.wx + a.xyzt*b.wxt + a.xz*b.wxy + a.xzt*b.wxyt - a.y*b.wz + a.yt*b.wzt + a.yz*b.w + a.yzt*b.wt + a.z*b.wy - a.zt*b.wyt) ; 

c.wyt = + (a.q*b.wyt + a.t*b.wy + a.w*b.yt - a.wt*b.y + a.wx*b.xyt + a.wxt*b.xy - a.wxy*b.xt + a.wxyt*b.x + a.wxyz*b.xzt + a.wxyzt*b.xz + a.wxz*b.xyzt - a.wxzt*b.xyz + a.wy*b.t + a.wyt*b.q + a.wyz*b.zt - a.wyzt*b.z - a.wz*b.yzt - a.wzt*b.yz - a.x*b.wxyt + a.xt*b.wxy - a.xy*b.wxt - a.xyt*b.wx - a.xyz*b.wxzt + a.xyzt*b.wxz + a.xz*b.wxyzt + a.xzt*b.wxyz - a.y*b.wt + a.yt*b.w + a.yz*b.wzt + a.yzt*b.wz + a.z*b.wyzt - a.zt*b.wyz) ; 

c.wzt = + (a.q*b.wzt + a.t*b.wz + a.w*b.zt - a.wt*b.z + a.wx*b.xzt + a.wxt*b.xz - a.wxy*b.xyzt + a.wxyt*b.xyz - a.wxyz*b.xyt - a.wxyzt*b.xy - a.wxz*b.xt + a.wxzt*b.x + a.wy*b.yzt + a.wyt*b.yz - a.wyz*b.yt + a.wyzt*b.y + a.wz*b.t + a.wzt*b.q - a.x*b.wxzt + a.xt*b.wxz - a.xy*b.wxyzt - a.xyt*b.wxyz + a.xyz*b.wxyt - a.xyzt*b.wxy - a.xz*b.wxt - a.xzt*b.wx - a.y*b.wyzt + a.yt*b.wyz - a.yz*b.wyt - a.yzt*b.wy - a.z*b.wt + a.zt*b.w) ; 

c.xyz = + (a.q*b.xyz + a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt - a.wx*b.wyz - a.wxt*b.wyzt + a.wxy*b.wz - a.wxyt*b.wzt - a.wxyz*b.w - a.wxyzt*b.wt - a.wxz*b.wy + a.wxzt*b.wyt + a.wy*b.wxz + a.wyt*b.wxzt + a.wyz*b.wx - a.wyzt*b.wxt - a.wz*b.wxy - a.wzt*b.wxyt + a.x*b.yz - a.xt*b.yzt + a.xy*b.z + a.xyt*b.zt + a.xyz*b.q - a.xyzt*b.t - a.xz*b.y - a.xzt*b.yt - a.y*b.xz + a.yt*b.xzt + a.yz*b.x + a.yzt*b.xt + a.z*b.xy - a.zt*b.xyt) ; 

c.xyt = + (a.q*b.xyt + a.t*b.xy + a.w*b.wxyt - a.wt*b.wxy - a.wx*b.wyt - a.wxt*b.wy + a.wxy*b.wt - a.wxyt*b.w - a.wxyz*b.wzt - a.wxyzt*b.wz - a.wxz*b.wyzt + a.wxzt*b.wyz + a.wy*b.wxt + a.wyt*b.wx + a.wyz*b.wxzt - a.wyzt*b.wxz - a.wz*b.wxyzt - a.wzt*b.wxyz + a.x*b.yt - a.xt*b.y + a.xy*b.t + a.xyt*b.q + a.xyz*b.zt - a.xyzt*b.z - a.xz*b.yzt - a.xzt*b.yz - a.y*b.xt + a.yt*b.x + a.yz*b.xzt + a.yzt*b.xz + a.z*b.xyzt - a.zt*b.xyz) ; 

c.xzt = + (a.q*b.xzt + a.t*b.xz + a.w*b.wxzt - a.wt*b.wxz - a.wx*b.wzt - a.wxt*b.wz + a.wxy*b.wyzt - a.wxyt*b.wyz + a.wxyz*b.wyt + a.wxyzt*b.wy + a.wxz*b.wt - a.wxzt*b.w + a.wy*b.wxyzt + a.wyt*b.wxyz - a.wyz*b.wxyt + a.wyzt*b.wxy + a.wz*b.wxt + a.wzt*b.wx + a.x*b.zt - a.xt*b.z + a.xy*b.yzt + a.xyt*b.yz - a.xyz*b.yt + a.xyzt*b.y + a.xz*b.t + a.xzt*b.q - a.y*b.xyzt + a.yt*b.xyz - a.yz*b.xyt - a.yzt*b.xy - a.z*b.xt + a.zt*b.x) ; 

c.yzt = + (a.q*b.yzt + a.t*b.yz + a.w*b.wyzt - a.wt*b.wyz - a.wx*b.wxyzt - a.wxt*b.wxyz - a.wxy*b.wxzt + a.wxyt*b.wxz - a.wxyz*b.wxt - a.wxyzt*b.wx + a.wxz*b.wxyt - a.wxzt*b.wxy - a.wy*b.wzt - a.wyt*b.wz + a.wyz*b.wt - a.wyzt*b.w + a.wz*b.wyt + a.wzt*b.wy + a.x*b.xyzt - a.xt*b.xyz - a.xy*b.xzt - a.xyt*b.xz + a.xyz*b.xt - a.xyzt*b.x + a.xz*b.xyt + a.xzt*b.xy + a.y*b.zt - a.yt*b.z + a.yz*b.t + a.yzt*b.q - a.z*b.yt + a.zt*b.y) ; 



c.wxyz = + (a.q*b.wxyz - a.t*b.wxyzt + a.w*b.xyz + a.wt*b.xyzt + a.wx*b.yz - a.wxt*b.yzt + a.wxy*b.z + a.wxyt*b.zt + a.wxyz*b.q - a.wxyzt*b.t - a.wxz*b.y - a.wxzt*b.yt - a.wy*b.xz + a.wyt*b.xzt + a.wyz*b.x + a.wyzt*b.xt + a.wz*b.xy - a.wzt*b.xyt - a.x*b.wyz - a.xt*b.wyzt + a.xy*b.wz - a.xyt*b.wzt - a.xyz*b.w - a.xyzt*b.wt - a.xz*b.wy + a.xzt*b.wyt + a.y*b.wxz + a.yt*b.wxzt + a.yz*b.wx - a.yzt*b.wxt - a.z*b.wxy - a.zt*b.wxyt) ; 

c.wxyt = + (a.q*b.wxyt - a.t*b.wxy + a.w*b.xyt + a.wt*b.xy + a.wx*b.yt - a.wxt*b.y + a.wxy*b.t + a.wxyt*b.q + a.wxyz*b.zt - a.wxyzt*b.z - a.wxz*b.yzt - a.wxzt*b.yz - a.wy*b.xt + a.wyt*b.x + a.wyz*b.xzt + a.wyzt*b.xz + a.wz*b.xyzt - a.wzt*b.xyz - a.x*b.wyt - a.xt*b.wy + a.xy*b.wt - a.xyt*b.w - a.xyz*b.wzt - a.xyzt*b.wz - a.xz*b.wyzt + a.xzt*b.wyz + a.y*b.wxt + a.yt*b.wx + a.yz*b.wxzt - a.yzt*b.wxz - a.z*b.wxyzt - a.zt*b.wxyz) ; 

c.wxzt = + (a.q*b.wxzt - a.t*b.wxz + a.w*b.xzt + a.wt*b.xz + a.wx*b.zt - a.wxt*b.z + a.wxy*b.yzt + a.wxyt*b.yz - a.wxyz*b.yt + a.wxyzt*b.y + a.wxz*b.t + a.wxzt*b.q - a.wy*b.xyzt + a.wyt*b.xyz - a.wyz*b.xyt - a.wyzt*b.xy - a.wz*b.xt + a.wzt*b.x - a.x*b.wzt - a.xt*b.wz + a.xy*b.wyzt - a.xyt*b.wyz + a.xyz*b.wyt + a.xyzt*b.wy + a.xz*b.wt - a.xzt*b.w + a.y*b.wxyzt + a.yt*b.wxyz - a.yz*b.wxyt + a.yzt*b.wxy + a.z*b.wxt + a.zt*b.wx) ; 

c.wyzt = + (a.q*b.wyzt - a.t*b.wyz + a.w*b.yzt + a.wt*b.yz + a.wx*b.xyzt - a.wxt*b.xyz - a.wxy*b.xzt - a.wxyt*b.xz + a.wxyz*b.xt - a.wxyzt*b.x + a.wxz*b.xyt + a.wxzt*b.xy + a.wy*b.zt - a.wyt*b.z + a.wyz*b.t + a.wyzt*b.q - a.wz*b.yt + a.wzt*b.y - a.x*b.wxyzt - a.xt*b.wxyz - a.xy*b.wxzt + a.xyt*b.wxz - a.xyz*b.wxt - a.xyzt*b.wx + a.xz*b.wxyt - a.xzt*b.wxy - a.y*b.wzt - a.yt*b.wz + a.yz*b.wt - a.yzt*b.w + a.z*b.wyt + a.zt*b.wy) ; 

c.xyzt = + (a.q*b.xyzt - a.t*b.xyz + a.w*b.wxyzt + a.wt*b.wxyz - a.wx*b.wyzt + a.wxt*b.wyz + a.wxy*b.wzt + a.wxyt*b.wz - a.wxyz*b.wt + a.wxyzt*b.w - a.wxz*b.wyt - a.wxzt*b.wy + a.wy*b.wxzt - a.wyt*b.wxz + a.wyz*b.wxt + a.wyzt*b.wx - a.wz*b.wxyt + a.wzt*b.wxy + a.x*b.yzt + a.xt*b.yz + a.xy*b.zt - a.xyt*b.z + a.xyz*b.t + a.xyzt*b.q - a.xz*b.yt + a.xzt*b.y - a.y*b.xzt - a.yt*b.xz + a.yz*b.xt - a.yzt*b.x + a.z*b.xyt + a.zt*b.xy) ; 


c.wxyzt = + (a.q*b.wxyzt + a.t*b.wxyz + a.w*b.xyzt - a.wt*b.xyz + a.wx*b.yzt + a.wxt*b.yz + a.wxy*b.zt - a.wxyt*b.z + a.wxyz*b.t + a.wxyzt*b.q - a.wxz*b.yt + a.wxzt*b.y - a.wy*b.xzt - a.wyt*b.xz + a.wyz*b.xt - a.wyzt*b.x + a.wz*b.xyt + a.wzt*b.xy - a.x*b.wyzt + a.xt*b.wyz + a.xy*b.wzt + a.xyt*b.wz - a.xyz*b.wt + a.xyzt*b.w - a.xz*b.wyt - a.xzt*b.wy + a.y*b.wxzt - a.yt*b.wxz + a.yz*b.wxt + a.yzt*b.wx - a.z*b.wxyt + a.zt*b.wxy) ; 



	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 operator*(const GA5_4_1 &v, const ex &u) {
	GA5_4_1 w;

	w.q    = expand(u*v.q    ) ;

	w.w    = expand(u*v.w    ) ;
	w.x    = expand(u*v.x    ) ;
	w.y    = expand(u*v.y    ) ;
	w.z    = expand(u*v.z    ) ;
	w.t    = expand(u*v.t    ) ;

	w.wx   = expand(u*v.wx   ) ;
	w.wy   = expand(u*v.wy   ) ;
	w.wz   = expand(u*v.wz   ) ;
	w.wt   = expand(u*v.wt   ) ;
	w.xy   = expand(u*v.xy   ) ;
	w.xz   = expand(u*v.xz   ) ;
	w.xt   = expand(u*v.xt   ) ;
	w.yz   = expand(u*v.yz   ) ;
	w.yt   = expand(u*v.yt   ) ;
	w.zt   = expand(u*v.zt   ) ;

	w.wxy  = expand(u*v.wxy  ) ;
	w.wxz  = expand(u*v.wxz  ) ;
	w.wxt  = expand(u*v.wxt  ) ;
	w.wyz  = expand(u*v.wyz  ) ;
	w.wyt  = expand(u*v.wyt  ) ;
	w.wzt  = expand(u*v.wzt  ) ;
	w.xyz  = expand(u*v.xyz  ) ;
	w.xyt  = expand(u*v.xyt  ) ;
	w.xzt  = expand(u*v.xzt  ) ;
	w.yzt  = expand(u*v.yzt  ) ;

	w.wxyz = expand(u*v.wxyz ) ;
	w.wxyt = expand(u*v.wxyt ) ;
	w.wxzt = expand(u*v.wxzt ) ;
	w.wyzt = expand(u*v.wyzt ) ;
	w.xyzt = expand(u*v.xyzt ) ;

	w.wxyzt = expand(u*v.wxyzt ) ;

	return w;
}

//////////////////////////////////////////////////////

GA5_4_1 operator*(const ex &u, const GA5_4_1 &v) {
	GA5_4_1 w;

	w.q    = expand(u*v.q    ) ;

	w.w    = expand(u*v.w    ) ;
	w.x    = expand(u*v.x    ) ;
	w.y    = expand(u*v.y    ) ;
	w.z    = expand(u*v.z    ) ;
	w.t    = expand(u*v.t    ) ;

	w.wx   = expand(u*v.wx   ) ;
	w.wy   = expand(u*v.wy   ) ;
	w.wz   = expand(u*v.wz   ) ;
	w.wt   = expand(u*v.wt   ) ;
	w.xy   = expand(u*v.xy   ) ;
	w.xz   = expand(u*v.xz   ) ;
	w.xt   = expand(u*v.xt   ) ;
	w.yz   = expand(u*v.yz   ) ;
	w.yt   = expand(u*v.yt   ) ;
	w.zt   = expand(u*v.zt   ) ;

	w.wxy  = expand(u*v.wxy  ) ;
	w.wxz  = expand(u*v.wxz  ) ;
	w.wxt  = expand(u*v.wxt  ) ;
	w.wyz  = expand(u*v.wyz  ) ;
	w.wyt  = expand(u*v.wyt  ) ;
	w.wzt  = expand(u*v.wzt  ) ;
	w.xyz  = expand(u*v.xyz  ) ;
	w.xyt  = expand(u*v.xyt  ) ;
	w.xzt  = expand(u*v.xzt  ) ;
	w.yzt  = expand(u*v.yzt  ) ;

	w.wxyz = expand(u*v.wxyz ) ;
	w.wxyt = expand(u*v.wxyt ) ;
	w.wxzt = expand(u*v.wxzt ) ;
	w.wyzt = expand(u*v.wyzt ) ;
	w.xyzt = expand(u*v.xyzt ) ;

	w.wxyzt = expand(u*v.wxyzt ) ;

	return w;
}


//////////////////////////////////////////////////////

GA5_4_1 operator/(const GA5_4_1 &u, const int i)
{
	GA5_4_1 a;
	a.q = u.q/i;

	a.w = u.w/i;
	a.x = u.x/i;
	a.y = u.y/i;
	a.z = u.z/i;
	a.t = u.t/i;

	a.wx = u.wx/i;
	a.wy = u.wy/i;
	a.wz = u.wz/i;
	a.wt = u.wt/i;
	a.xy = u.xy/i;
	a.xz = u.xz/i;
	a.yz = u.yz/i;
	a.xt = u.xt/i;
	a.yt = u.yt/i;
	a.zt = u.zt/i;

	a.wxy = u.wxy/i;
	a.wxz = u.wxz/i;
	a.wxt = u.wxt/i;
	a.wyz = u.wyz/i;
	a.wyt = u.wyt/i;
	a.wzt = u.wzt/i;
	a.xyz = u.xyz/i;
	a.xyt = u.xyt/i;
	a.xzt = u.xzt/i;
	a.yzt = u.yzt/i;

	a.wxyz = u.wxyz/i;
	a.wxyt = u.wxyt/i;
	a.wxzt = u.wxzt/i;
	a.wyzt = u.wyzt/i;
	a.xyzt = u.xyzt/i;

	a.wxyzt = u.wxyzt/i;

	return a;
}

//////////////////////////////////////////////////////

GA5_4_1 operator^(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

c.q    =  + a.q   *b.q    ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    ; 
c.wx   =  + a.q   *b.wx   + a.w   *b.x    - a.x   *b.w    + a.wx  *b.q    ; 
c.wy   =  + a.q   *b.wy   + a.w   *b.y    - a.y   *b.w    + a.wy  *b.q    ; 
c.wz   =  + a.q   *b.wz   + a.w   *b.z    - a.z   *b.w    + a.wz  *b.q    ; 
c.wt   =  + a.q   *b.wt   + a.w   *b.t    - a.t   *b.w    + a.wt  *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.xy  *b.q    ; 
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.z   *b.x    + a.xz  *b.q    ; 
c.xt   =  + a.q   *b.xt   + a.x   *b.t    - a.t   *b.x    + a.xt  *b.q    ; 
c.yz   =  + a.q   *b.yz   + a.y   *b.z    - a.z   *b.y    + a.yz  *b.q    ; 
c.yt   =  + a.q   *b.yt   + a.y   *b.t    - a.t   *b.y    + a.yt  *b.q    ; 
c.zt   =  + a.q   *b.zt   + a.z   *b.t    - a.t   *b.z    + a.zt  *b.q    ; 
c.wxy  =  + a.q   *b.wxy  + a.w   *b.xy   - a.x   *b.wy   + a.y   *b.wx   + a.wx  *b.y    - a.wy  *b.x    + a.xy  *b.w    + a.wxy *b.q    ; 
c.wxz  =  + a.q   *b.wxz  + a.w   *b.xz   - a.x   *b.wz   + a.z   *b.wx   + a.wx  *b.z    - a.wz  *b.x    + a.xz  *b.w    + a.wxz *b.q    ; 
c.wxt  =  + a.q   *b.wxt  + a.w   *b.xt   - a.x   *b.wt   + a.t   *b.wx   + a.wx  *b.t    - a.wt  *b.x    + a.xt  *b.w    + a.wxt *b.q    ; 
c.wyz  =  + a.q   *b.wyz  + a.w   *b.yz   - a.y   *b.wz   + a.z   *b.wy   + a.wy  *b.z    - a.wz  *b.y    + a.yz  *b.w    + a.wyz *b.q    ; 
c.wyt  =  + a.q   *b.wyt  + a.w   *b.yt   - a.y   *b.wt   + a.t   *b.wy   + a.wy  *b.t    - a.wt  *b.y    + a.yt  *b.w    + a.wyt *b.q    ; 
c.wzt  =  + a.q   *b.wzt  + a.w   *b.zt   - a.z   *b.wt   + a.t   *b.wz   + a.wz  *b.t    - a.wt  *b.z    + a.zt  *b.w    + a.wzt *b.q    ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    + a.xyt *b.q    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    + a.xzt *b.q    ; 
c.yzt  =  + a.q   *b.yzt  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    + a.yzt *b.q    ; 
c.wxyz =  + a.q   *b.wxyz + a.w   *b.xyz  - a.x   *b.wyz  + a.y   *b.wxz  - a.z   *b.wxy  + a.wx  *b.yz   - a.wy  *b.xz   + a.wz  *b.xy   + a.xy  *b.wz   - a.xz  *b.wy   + a.yz  *b.wx   + a.wxy *b.z    - a.wxz *b.y    + a.wyz *b.x    - a.xyz *b.w    + a.wxyz*b.q    ; 
c.wxyt =  + a.q   *b.wxyt + a.w   *b.xyt  - a.x   *b.wyt  + a.y   *b.wxt  - a.t   *b.wxy  + a.wx  *b.yt   - a.wy  *b.xt   + a.wt  *b.xy   + a.xy  *b.wt   - a.xt  *b.wy   + a.yt  *b.wx   + a.wxy *b.t    - a.wxt *b.y    + a.wyt *b.x    - a.xyt *b.w    + a.wxyt*b.q    ; 
c.wxzt =  + a.q   *b.wxzt + a.w   *b.xzt  - a.x   *b.wzt  + a.z   *b.wxt  - a.t   *b.wxz  + a.wx  *b.zt   - a.wz  *b.xt   + a.wt  *b.xz   + a.xz  *b.wt   - a.xt  *b.wz   + a.zt  *b.wx   + a.wxz *b.t    - a.wxt *b.z    + a.wzt *b.x    - a.xzt *b.w    + a.wxzt*b.q    ; 
c.wyzt =  + a.q   *b.wyzt + a.w   *b.yzt  - a.y   *b.wzt  + a.z   *b.wyt  - a.t   *b.wyz  + a.wy  *b.zt   - a.wz  *b.yt   + a.wt  *b.yz   + a.yz  *b.wt   - a.yt  *b.wz   + a.zt  *b.wy   + a.wyz *b.t    - a.wyt *b.z    + a.wzt *b.y    - a.yzt *b.w    + a.wyzt*b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.xt  *b.yz   + a.yz  *b.xt   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q    ; 
c.wxyzt =  + a.q   *b.wxyzt + a.w   *b.xyzt - a.x   *b.wyzt + a.y   *b.wxzt - a.z   *b.wxyt + a.t   *b.wxyz + a.wx  *b.yzt  - a.wy  *b.xzt  + a.wz  *b.xyt  - a.wt  *b.xyz  + a.xy  *b.wzt  - a.xz  *b.wyt  + a.xt  *b.wyz  + a.yz  *b.wxt  - a.yt  *b.wxz  + a.zt  *b.wxy  + a.wxy *b.zt   - a.wxz *b.yt   + a.wxt *b.yz   + a.wyz *b.xt   - a.wyt *b.xz   + a.wzt *b.xy   - a.xyz *b.wt   + a.xyt *b.wz   - a.xzt *b.wy   + a.yzt *b.wx   + a.wxyz*b.t    - a.wxyt*b.z    + a.wxzt*b.y    - a.wyzt*b.x    + a.xyzt*b.w    + a.wxyzt*b.q    ; 


	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Wedge(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

c.q    =  + a.q   *b.q    ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    ; 
c.wx   =  + a.q   *b.wx   + a.w   *b.x    - a.x   *b.w    + a.wx  *b.q    ; 
c.wy   =  + a.q   *b.wy   + a.w   *b.y    - a.y   *b.w    + a.wy  *b.q    ; 
c.wz   =  + a.q   *b.wz   + a.w   *b.z    - a.z   *b.w    + a.wz  *b.q    ; 
c.wt   =  + a.q   *b.wt   + a.w   *b.t    - a.t   *b.w    + a.wt  *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.xy  *b.q    ; 
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.z   *b.x    + a.xz  *b.q    ; 
c.xt   =  + a.q   *b.xt   + a.x   *b.t    - a.t   *b.x    + a.xt  *b.q    ; 
c.yz   =  + a.q   *b.yz   + a.y   *b.z    - a.z   *b.y    + a.yz  *b.q    ; 
c.yt   =  + a.q   *b.yt   + a.y   *b.t    - a.t   *b.y    + a.yt  *b.q    ; 
c.zt   =  + a.q   *b.zt   + a.z   *b.t    - a.t   *b.z    + a.zt  *b.q    ; 
c.wxy  =  + a.q   *b.wxy  + a.w   *b.xy   - a.x   *b.wy   + a.y   *b.wx   + a.wx  *b.y    - a.wy  *b.x    + a.xy  *b.w    + a.wxy *b.q    ; 
c.wxz  =  + a.q   *b.wxz  + a.w   *b.xz   - a.x   *b.wz   + a.z   *b.wx   + a.wx  *b.z    - a.wz  *b.x    + a.xz  *b.w    + a.wxz *b.q    ; 
c.wxt  =  + a.q   *b.wxt  + a.w   *b.xt   - a.x   *b.wt   + a.t   *b.wx   + a.wx  *b.t    - a.wt  *b.x    + a.xt  *b.w    + a.wxt *b.q    ; 
c.wyz  =  + a.q   *b.wyz  + a.w   *b.yz   - a.y   *b.wz   + a.z   *b.wy   + a.wy  *b.z    - a.wz  *b.y    + a.yz  *b.w    + a.wyz *b.q    ; 
c.wyt  =  + a.q   *b.wyt  + a.w   *b.yt   - a.y   *b.wt   + a.t   *b.wy   + a.wy  *b.t    - a.wt  *b.y    + a.yt  *b.w    + a.wyt *b.q    ; 
c.wzt  =  + a.q   *b.wzt  + a.w   *b.zt   - a.z   *b.wt   + a.t   *b.wz   + a.wz  *b.t    - a.wt  *b.z    + a.zt  *b.w    + a.wzt *b.q    ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    + a.xyt *b.q    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    + a.xzt *b.q    ; 
c.yzt  =  + a.q   *b.yzt  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    + a.yzt *b.q    ; 
c.wxyz =  + a.q   *b.wxyz + a.w   *b.xyz  - a.x   *b.wyz  + a.y   *b.wxz  - a.z   *b.wxy  + a.wx  *b.yz   - a.wy  *b.xz   + a.wz  *b.xy   + a.xy  *b.wz   - a.xz  *b.wy   + a.yz  *b.wx   + a.wxy *b.z    - a.wxz *b.y    + a.wyz *b.x    - a.xyz *b.w    + a.wxyz*b.q    ; 
c.wxyt =  + a.q   *b.wxyt + a.w   *b.xyt  - a.x   *b.wyt  + a.y   *b.wxt  - a.t   *b.wxy  + a.wx  *b.yt   - a.wy  *b.xt   + a.wt  *b.xy   + a.xy  *b.wt   - a.xt  *b.wy   + a.yt  *b.wx   + a.wxy *b.t    - a.wxt *b.y    + a.wyt *b.x    - a.xyt *b.w    + a.wxyt*b.q    ; 
c.wxzt =  + a.q   *b.wxzt + a.w   *b.xzt  - a.x   *b.wzt  + a.z   *b.wxt  - a.t   *b.wxz  + a.wx  *b.zt   - a.wz  *b.xt   + a.wt  *b.xz   + a.xz  *b.wt   - a.xt  *b.wz   + a.zt  *b.wx   + a.wxz *b.t    - a.wxt *b.z    + a.wzt *b.x    - a.xzt *b.w    + a.wxzt*b.q    ; 
c.wyzt =  + a.q   *b.wyzt + a.w   *b.yzt  - a.y   *b.wzt  + a.z   *b.wyt  - a.t   *b.wyz  + a.wy  *b.zt   - a.wz  *b.yt   + a.wt  *b.yz   + a.yz  *b.wt   - a.yt  *b.wz   + a.zt  *b.wy   + a.wyz *b.t    - a.wyt *b.z    + a.wzt *b.y    - a.yzt *b.w    + a.wyzt*b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.xt  *b.yz   + a.yz  *b.xt   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q    ; 
c.wxyzt =  + a.q   *b.wxyzt + a.w   *b.xyzt - a.x   *b.wyzt + a.y   *b.wxzt - a.z   *b.wxyt + a.t   *b.wxyz + a.wx  *b.yzt  - a.wy  *b.xzt  + a.wz  *b.xyt  - a.wt  *b.xyz  + a.xy  *b.wzt  - a.xz  *b.wyt  + a.xt  *b.wyz  + a.yz  *b.wxt  - a.yt  *b.wxz  + a.zt  *b.wxy  + a.wxy *b.zt   - a.wxz *b.yt   + a.wxt *b.yz   + a.wyz *b.xt   - a.wyt *b.xz   + a.wzt *b.xy   - a.xyz *b.wt   + a.xyt *b.wz   - a.xzt *b.wy   + a.yzt *b.wx   + a.wxyz*b.t    - a.wxyt*b.z    + a.wxzt*b.y    - a.wyzt*b.x    + a.xyzt*b.w    + a.wxyzt*b.q    ; 


	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 AntiWedge(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

	c.q =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz
		 + a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt
		 + a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt
		 + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz
		 - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x
		 + a.xyzt*b.w + a.wxyzt*b.q ; 
	c.w =  + a.w*b.wxyzt + a.wx*b.wyzt - a.wy*b.wxzt + a.wz*b.wxyt - a.wt*b.wxyz + a.wxy*b.wzt
		 - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.wxyz*b.wt
		 - a.wxyt*b.wz + a.wxzt*b.wy - a.wyzt*b.wx + a.wxyzt*b.w ; 
	c.x =  + a.x*b.wxyzt + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz + a.wxy*b.xzt
		 - a.wxz*b.xyt + a.wxt*b.xyz + a.xyz*b.wxt - a.xyt*b.wxz + a.xzt*b.wxy + a.wxyz*b.xt
		 - a.wxyt*b.xz + a.wxzt*b.xy - a.xyzt*b.wx + a.wxyzt*b.x ; 
	c.y =  + a.y*b.wxyzt + a.wy*b.xyzt - a.xy*b.wyzt + a.yz*b.wxyt - a.yt*b.wxyz + a.wxy*b.yzt
		 - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz + a.yzt*b.wxy + a.wxyz*b.yt
		 - a.wxyt*b.yz + a.wyzt*b.xy - a.xyzt*b.wy + a.wxyzt*b.y ; 
	c.z =  + a.z*b.wxyzt + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz + a.wxz*b.yzt
		 - a.wyz*b.xzt + a.wzt*b.xyz + a.xyz*b.wzt - a.xzt*b.wyz + a.yzt*b.wxz + a.wxyz*b.zt
		 - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz + a.wxyzt*b.z ; 
	c.t =  + a.t*b.wxyzt + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt + a.wxt*b.yzt
		 - a.wyt*b.xzt + a.wzt*b.xyt + a.xyt*b.wzt - a.xzt*b.wyt + a.yzt*b.wxt + a.wxyt*b.zt
		 - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt + a.wxyzt*b.t ; 
	
	c.wx =  + a.wx*b.wxyzt + a.wxy*b.wxzt - a.wxz*b.wxyt + a.wxt*b.wxyz + a.wxyz*b.wxt - a.wxyt*b.wxz + a.wxzt*b.wxy + a.wxyzt*b.wx ; 
	c.wy =  + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
	c.wz =  + a.wz*b.wxyzt + a.wxz*b.wyzt - a.wyz*b.wxzt + a.wzt*b.wxyz + a.wxyz*b.wzt - a.wxzt*b.wyz + a.wyzt*b.wxz + a.wxyzt*b.wz ; 
	c.wt =  + a.wt*b.wxyzt + a.wxt*b.wyzt - a.wyt*b.wxzt + a.wzt*b.wxyt + a.wxyt*b.wzt - a.wxzt*b.wyt + a.wyzt*b.wxt + a.wxyzt*b.wt ; 
	c.xy =  + a.xy*b.wxyzt + a.wxy*b.xyzt - a.xyz*b.wxyt + a.xyt*b.wxyz + a.wxyz*b.xyt - a.wxyt*b.xyz + a.xyzt*b.wxy + a.wxyzt*b.xy ; 
	c.xz =  + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.xt =  + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.yz =  + a.yz*b.wxyzt + a.wyz*b.xyzt - a.xyz*b.wyzt + a.yzt*b.wxyz + a.wxyz*b.yzt - a.wyzt*b.xyz + a.xyzt*b.wyz + a.wxyzt*b.yz ; 
	c.yt =  + a.yt*b.wxyzt + a.wyt*b.xyzt - a.xyt*b.wyzt + a.yzt*b.wxyt + a.wxyt*b.yzt - a.wyzt*b.xyt + a.xyzt*b.wyt + a.wxyzt*b.yt ; 
	c.zt =  + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 

	c.wxy =  + a.wxy*b.wxyzt + a.wxyz*b.wxyt - a.wxyt*b.wxyz + a.wxyzt*b.wxy ; 
	c.wxz =  + a.wxz*b.wxyzt + a.wxyz*b.wxzt - a.wxzt*b.wxyz + a.wxyzt*b.wxz ; 
	c.wxt =  + a.wxt*b.wxyzt + a.wxyt*b.wxzt - a.wxzt*b.wxyt + a.wxyzt*b.wxt ; 
	c.wyz =  + a.wyz*b.wxyzt + a.wxyz*b.wyzt - a.wyzt*b.wxyz + a.wxyzt*b.wyz ; 
	c.wyt =  + a.wyt*b.wxyzt + a.wxyt*b.wyzt - a.wyzt*b.wxyt + a.wxyzt*b.wyt ; 
	c.wzt =  + a.wzt*b.wxyzt + a.wxzt*b.wyzt - a.wyzt*b.wxzt + a.wxyzt*b.wzt ; 
	c.xyz =  + a.xyz*b.wxyzt + a.wxyz*b.xyzt - a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xyt =  + a.xyt*b.wxyzt + a.wxyt*b.xyzt - a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.xzt =  + a.xzt*b.wxyzt + a.wxzt*b.xyzt - a.xyzt*b.wxzt + a.wxyzt*b.xzt ; 
	c.yzt =  + a.yzt*b.wxyzt + a.wyzt*b.xyzt - a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 

	c.wxyz =  + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
	c.wxyt =  + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
	c.wxzt =  + a.wxzt*b.wxyzt + a.wxyzt*b.wxzt ; 
	c.wyzt =  + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
	c.xyzt =  + a.xyzt*b.wxyzt + a.wxyzt*b.xyzt ; 

	c.wxyzt =  + a.wxyzt*b.wxyzt ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Regressive(GA5_4_1 a, GA5_4_1 b) 
{
	GA5_4_1 c;

//	GA5_4_1 I, I_inv;
//	I.wxyzt = 1;
//	I_inv.wxyzt = -1;
//	c = ((a*I_inv)^(b*I_inv))*I;


	c.q =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz 
		+ a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt 
		+ a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt
		 + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz
		 - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w + a.wxyzt*b.q ; 

	c.w =  + a.w*b.wxyzt - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz
		 + a.wxy*b.wzt - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz
		 + a.wzt*b.wxy - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx + a.wxyzt*b.w ; 
	c.x =  + a.x*b.wxyzt - a.wx*b.xyzt + a.xy*b.wxzt - a.xz*b.wxyt + a.xt*b.wxyz
		 + a.wxy*b.xzt - a.wxz*b.xyt + a.wxt*b.xyz + a.xyz*b.wxt - a.xyt*b.wxz
		 + a.xzt*b.wxy - a.wxyz*b.xt + a.wxyt*b.xz - a.wxzt*b.xy + a.xyzt*b.wx + a.wxyzt*b.x ; 
	c.y =  + a.y*b.wxyzt - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz
		 + a.wxy*b.yzt - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz
		 + a.yzt*b.wxy - a.wxyz*b.yt + a.wxyt*b.yz - a.wyzt*b.xy + a.xyzt*b.wy + a.wxyzt*b.y ; 
	c.z =  + a.z*b.wxyzt - a.wz*b.xyzt + a.xz*b.wyzt - a.yz*b.wxzt + a.zt*b.wxyz
		 + a.wxz*b.yzt - a.wyz*b.xzt + a.wzt*b.xyz + a.xyz*b.wzt - a.xzt*b.wyz
		 + a.yzt*b.wxz - a.wxyz*b.zt + a.wxzt*b.yz - a.wyzt*b.xz + a.xyzt*b.wz + a.wxyzt*b.z ; 
	c.t =  + a.t*b.wxyzt - a.wt*b.xyzt + a.xt*b.wyzt - a.yt*b.wxzt + a.zt*b.wxyt
		 + a.wxt*b.yzt - a.wyt*b.xzt + a.wzt*b.xyt + a.xyt*b.wzt - a.xzt*b.wyt
		 + a.yzt*b.wxt - a.wxyt*b.zt + a.wxzt*b.yt - a.wyzt*b.xt + a.xyzt*b.wt + a.wxyzt*b.t ; 

	c.wx =  + a.wx*b.wxyzt + a.wxy*b.wxzt - a.wxz*b.wxyt + a.wxt*b.wxyz + a.wxyz*b.wxt - a.wxyt*b.wxz + a.wxzt*b.wxy + a.wxyzt*b.wx ; 
	c.wy =  + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
	c.wz =  + a.wz*b.wxyzt + a.wxz*b.wyzt - a.wyz*b.wxzt + a.wzt*b.wxyz + a.wxyz*b.wzt - a.wxzt*b.wyz + a.wyzt*b.wxz + a.wxyzt*b.wz ; 
	c.wt =  + a.wt*b.wxyzt + a.wxt*b.wyzt - a.wyt*b.wxzt + a.wzt*b.wxyt + a.wxyt*b.wzt - a.wxzt*b.wyt + a.wyzt*b.wxt + a.wxyzt*b.wt ; 
	c.xy =  + a.xy*b.wxyzt + a.wxy*b.xyzt - a.xyz*b.wxyt + a.xyt*b.wxyz + a.wxyz*b.xyt - a.wxyt*b.xyz + a.xyzt*b.wxy + a.wxyzt*b.xy ; 
	c.xz =  + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.xt =  + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.yz =  + a.yz*b.wxyzt + a.wyz*b.xyzt - a.xyz*b.wyzt + a.yzt*b.wxyz + a.wxyz*b.yzt - a.wyzt*b.xyz + a.xyzt*b.wyz + a.wxyzt*b.yz ; 
	c.yt =  + a.yt*b.wxyzt + a.wyt*b.xyzt - a.xyt*b.wyzt + a.yzt*b.wxyt + a.wxyt*b.yzt - a.wyzt*b.xyt + a.xyzt*b.wyt + a.wxyzt*b.yt ; 
	c.zt =  + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 

	c.wxy =  + a.wxy*b.wxyzt - a.wxyz*b.wxyt + a.wxyt*b.wxyz + a.wxyzt*b.wxy ; 
	c.wxz =  + a.wxz*b.wxyzt - a.wxyz*b.wxzt + a.wxzt*b.wxyz + a.wxyzt*b.wxz ; 
	c.wxt =  + a.wxt*b.wxyzt - a.wxyt*b.wxzt + a.wxzt*b.wxyt + a.wxyzt*b.wxt ; 
	c.wyz =  + a.wyz*b.wxyzt - a.wxyz*b.wyzt + a.wyzt*b.wxyz + a.wxyzt*b.wyz ; 
	c.wyt =  + a.wyt*b.wxyzt - a.wxyt*b.wyzt + a.wyzt*b.wxyt + a.wxyzt*b.wyt ; 
	c.wzt =  + a.wzt*b.wxyzt - a.wxzt*b.wyzt + a.wyzt*b.wxzt + a.wxyzt*b.wzt ; 
	c.xyz =  + a.xyz*b.wxyzt - a.wxyz*b.xyzt + a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xyt =  + a.xyt*b.wxyzt - a.wxyt*b.xyzt + a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.xzt =  + a.xzt*b.wxyzt - a.wxzt*b.xyzt + a.xyzt*b.wxzt + a.wxyzt*b.xzt ; 
	c.yzt =  + a.yzt*b.wxyzt - a.wyzt*b.xyzt + a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 

	c.wxyz =  + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
	c.wxyt =  + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
	c.wxzt =  + a.wxzt*b.wxyzt + a.wxyzt*b.wxzt ; 
	c.wyzt =  + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
	c.xyzt =  + a.xyzt*b.wxyzt + a.wxyzt*b.xyzt ; 

	c.wxyzt =  + a.wxyzt*b.wxyzt ; 


	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 RegressiveViaFormula(GA5_4_1 a, GA5_4_1 b) 
{
	GA5_4_1 c;

	GA5_4_1 I, I_inv;
	I.wxyzt = 1;
	I_inv.wxyzt = -1;
	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 LowerRightViaFormula(GA5_4_1 a, GA5_4_1 b)  // Like a mirror of Wedge along rising diagonal
{

// -Wedge(Blade[i]*xyzt,xyzt*Blade[j])

	GA5_4_1 c;

	GA5_4_1 I,I_inv,d,e,f;
	
	I.wxyzt = 1;	I_inv.wxyzt = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Expander(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

//	Expander equation set
	c.q = 0 ; 

	c.w = 0 ; 
	c.x = 0 ; 
	c.y = 0 ; 
	c.z = 0 ; 
	c.t = 0 ; 

	c.wx =  + a.w*b.x - a.x*b.w ; 
	c.wy =  + a.w*b.y - a.y*b.w ; 
	c.wz =  + a.w*b.z - a.z*b.w ; 
	c.wt =  + a.w*b.t - a.t*b.w ; 
	c.xy =  + a.x*b.y - a.y*b.x ; 
	c.xz =  + a.x*b.z - a.z*b.x ; 
	c.xt =  + a.x*b.t - a.t*b.x ; 
	c.yz =  + a.y*b.z - a.z*b.y ; 
	c.yt =  + a.y*b.t - a.t*b.y ; 
	c.zt =  + a.z*b.t - a.t*b.z ; 

	c.wxy =  + a.w*b.xy - a.x*b.wy + a.y*b.wx + a.wx*b.y - a.wy*b.x + a.xy*b.w ; 
	c.wxz =  + a.w*b.xz - a.x*b.wz + a.z*b.wx + a.wx*b.z - a.wz*b.x + a.xz*b.w ; 
	c.wxt =  + a.w*b.xt - a.x*b.wt + a.t*b.wx + a.wx*b.t - a.wt*b.x + a.xt*b.w ; 
	c.wyz =  + a.w*b.yz - a.y*b.wz + a.z*b.wy + a.wy*b.z - a.wz*b.y + a.yz*b.w ; 
	c.wyt =  + a.w*b.yt - a.y*b.wt + a.t*b.wy + a.wy*b.t - a.wt*b.y + a.yt*b.w ; 
	c.wzt =  + a.w*b.zt - a.z*b.wt + a.t*b.wz + a.wz*b.t - a.wt*b.z + a.zt*b.w ; 
	c.xyz =  + a.x*b.yz - a.y*b.xz + a.z*b.xy + a.xy*b.z - a.xz*b.y + a.yz*b.x ; 
	c.xyt =  + a.x*b.yt - a.y*b.xt + a.t*b.xy + a.xy*b.t - a.xt*b.y + a.yt*b.x ; 
	c.xzt =  + a.x*b.zt - a.z*b.xt + a.t*b.xz + a.xz*b.t - a.xt*b.z + a.zt*b.x ; 
	c.yzt =  + a.y*b.zt - a.z*b.yt + a.t*b.yz + a.yz*b.t - a.yt*b.z + a.zt*b.y ; 

	c.wxyz =  + a.w*b.xyz - a.x*b.wyz + a.y*b.wxz - a.z*b.wxy + a.wx*b.yz - a.wy*b.xz 
		+ a.wz*b.xy + a.xy*b.wz - a.xz*b.wy + a.yz*b.wx + a.wxy*b.z - a.wxz*b.y 
		- a.wxt*b.yzt + a.wyz*b.x + a.wyt*b.xzt - a.wzt*b.xyt - a.xyz*b.w 
		- a.xyt*b.wzt + a.xzt*b.wyt - a.yzt*b.wxt ; 
	c.wxyt =  + a.w*b.xyt - a.x*b.wyt + a.y*b.wxt - a.t*b.wxy + a.wx*b.yt - a.wy*b.xt 
		+ a.wt*b.xy + a.xy*b.wt - a.xt*b.wy + a.yt*b.wx + a.wxy*b.t - a.wxz*b.yzt 
		- a.wxt*b.y + a.wyz*b.xzt + a.wyt*b.x - a.wzt*b.xyz - a.xyz*b.wzt 
		- a.xyt*b.w + a.xzt*b.wyz - a.yzt*b.wxz ; 
	c.wxzt =  + a.w*b.xzt - a.x*b.wzt + a.z*b.wxt - a.t*b.wxz + a.wx*b.zt - a.wz*b.xt 
		+ a.wt*b.xz + a.xz*b.wt - a.xt*b.wz + a.zt*b.wx + a.wxy*b.yzt + a.wxz*b.t 
		- a.wxt*b.z - a.wyz*b.xyt + a.wyt*b.xyz + a.wzt*b.x + a.xyz*b.wyt 
		- a.xyt*b.wyz - a.xzt*b.w + a.yzt*b.wxy ; 
	c.wyzt =  + a.w*b.yzt - a.y*b.wzt + a.z*b.wyt - a.t*b.wyz + a.wy*b.zt - a.wz*b.yt 
		+ a.wt*b.yz + a.yz*b.wt - a.yt*b.wz + a.zt*b.wy - a.wxy*b.xzt + a.wxz*b.xyt 
		- a.wxt*b.xyz + a.wyz*b.t - a.wyt*b.z + a.wzt*b.y - a.xyz*b.wxt + a.xyt*b.wxz 
		- a.xzt*b.wxy - a.yzt*b.w ; 
	c.xyzt =  + a.x*b.yzt - a.y*b.xzt + a.z*b.xyt - a.t*b.xyz + a.xy*b.zt - a.xz*b.yt 
		+ a.xt*b.yz + a.yz*b.xt - a.yt*b.xz + a.zt*b.xy + a.wxy*b.wzt - a.wxz*b.wyt 
		+ a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.xyz*b.t 
		- a.xyt*b.z + a.xzt*b.y - a.yzt*b.x ; 

	c.wxyzt =  + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz 
		+ a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt 
		+ a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt 
		+ a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz 
		- a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Conserver(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

// Conserver equation set

	c.q =  + a.q*b.q ; 

	c.w =  + a.q*b.w + a.w*b.q ; 
	c.x =  + a.q*b.x + a.x*b.q ; 
	c.y =  + a.q*b.y + a.y*b.q ; 
	c.z =  + a.q*b.z + a.z*b.q ; 
	c.t =  + a.q*b.t + a.t*b.q ; 

	c.wx =  + a.q*b.wx + a.wx*b.q - a.wy*b.xy - a.wz*b.xz + a.wt*b.xt + a.xy*b.wy + a.xz*b.wz - a.xt*b.wt ; 
	c.wy =  + a.q*b.wy + a.wx*b.xy + a.wy*b.q - a.wz*b.yz + a.wt*b.yt - a.xy*b.wx + a.yz*b.wz - a.yt*b.wt ; 
	c.wz =  + a.q*b.wz + a.wx*b.xz + a.wy*b.yz + a.wz*b.q + a.wt*b.zt - a.xz*b.wx - a.yz*b.wy - a.zt*b.wt ; 
	c.wt =  + a.q*b.wt + a.wx*b.xt + a.wy*b.yt + a.wz*b.zt + a.wt*b.q - a.xt*b.wx - a.yt*b.wy - a.zt*b.wz ; 
	c.xy =  + a.q*b.xy - a.wx*b.wy + a.wy*b.wx + a.xy*b.q - a.xz*b.yz + a.xt*b.yt + a.yz*b.xz - a.yt*b.xt ; 
	c.xz =  + a.q*b.xz - a.wx*b.wz + a.wz*b.wx + a.xy*b.yz + a.xz*b.q + a.xt*b.zt - a.yz*b.xy - a.zt*b.xt ; 
	c.xt =  + a.q*b.xt - a.wx*b.wt + a.wt*b.wx + a.xy*b.yt + a.xz*b.zt + a.xt*b.q - a.yt*b.xy - a.zt*b.xz ; 
	c.yz =  + a.q*b.yz - a.wy*b.wz + a.wz*b.wy - a.xy*b.xz + a.xz*b.xy + a.yz*b.q + a.yt*b.zt - a.zt*b.yt ; 
	c.yt =  + a.q*b.yt - a.wy*b.wt + a.wt*b.wy - a.xy*b.xt + a.xt*b.xy + a.yz*b.zt + a.yt*b.q - a.zt*b.yz ; 
	c.zt =  + a.q*b.zt - a.wz*b.wt + a.wt*b.wz - a.xz*b.xt + a.xt*b.xz - a.yz*b.yt + a.yt*b.yz + a.zt*b.q ; 

	c.wxy =  + a.q*b.wxy + a.wz*b.xyz - a.wt*b.xyt - a.xz*b.wyz + a.xt*b.wyt + a.yz*b.wxz - a.yt*b.wxt
		 + a.wxy*b.q - a.wxz*b.yz + a.wxt*b.yt + a.wyz*b.xz - a.wyt*b.xt - a.xyz*b.wz + a.xyt*b.wt ; 
	c.wxz =  + a.q*b.wxz - a.wy*b.xyz - a.wt*b.xzt + a.xy*b.wyz + a.xt*b.wzt - a.yz*b.wxy - a.zt*b.wxt
		 + a.wxy*b.yz + a.wxz*b.q + a.wxt*b.zt - a.wyz*b.xy - a.wzt*b.xt + a.xyz*b.wy + a.xzt*b.wt ; 
	c.wxt =  + a.q*b.wxt - a.wy*b.xyt - a.wz*b.xzt + a.xy*b.wyt + a.xz*b.wzt - a.yt*b.wxy - a.zt*b.wxz
		 + a.wxy*b.yt + a.wxz*b.zt + a.wxt*b.q - a.wyt*b.xy - a.wzt*b.xz + a.xyt*b.wy + a.xzt*b.wz ; 
	c.wyz =  + a.q*b.wyz + a.wx*b.xyz - a.wt*b.yzt - a.xy*b.wxz + a.xz*b.wxy + a.yt*b.wzt - a.zt*b.wyt
		 - a.wxy*b.xz + a.wxz*b.xy + a.wyz*b.q + a.wyt*b.zt - a.wzt*b.yt - a.xyz*b.wx + a.yzt*b.wt ; 
	c.wyt =  + a.q*b.wyt + a.wx*b.xyt - a.wz*b.yzt - a.xy*b.wxt + a.xt*b.wxy + a.yz*b.wzt - a.zt*b.wyz
		 - a.wxy*b.xt + a.wxt*b.xy + a.wyz*b.zt + a.wyt*b.q - a.wzt*b.yz - a.xyt*b.wx + a.yzt*b.wz ; 
	c.wzt =  + a.q*b.wzt + a.wx*b.xzt + a.wy*b.yzt - a.xz*b.wxt + a.xt*b.wxz - a.yz*b.wyt + a.yt*b.wyz
		 - a.wxz*b.xt + a.wxt*b.xz - a.wyz*b.yt + a.wyt*b.yz + a.wzt*b.q - a.xzt*b.wx - a.yzt*b.wy ; 
	c.xyz =  + a.q*b.xyz - a.wx*b.wyz + a.wy*b.wxz - a.wz*b.wxy - a.xt*b.yzt + a.yt*b.xzt - a.zt*b.xyt
		 + a.wxy*b.wz - a.wxz*b.wy + a.wyz*b.wx + a.xyz*b.q + a.xyt*b.zt - a.xzt*b.yt + a.yzt*b.xt ; 
	c.xyt =  + a.q*b.xyt - a.wx*b.wyt + a.wy*b.wxt - a.wt*b.wxy - a.xz*b.yzt + a.yz*b.xzt - a.zt*b.xyz
		 + a.wxy*b.wt - a.wxt*b.wy + a.wyt*b.wx + a.xyz*b.zt + a.xyt*b.q - a.xzt*b.yz + a.yzt*b.xz ; 
	c.xzt =  + a.q*b.xzt - a.wx*b.wzt + a.wz*b.wxt - a.wt*b.wxz + a.xy*b.yzt - a.yz*b.xyt + a.yt*b.xyz
		 + a.wxz*b.wt - a.wxt*b.wz + a.wzt*b.wx - a.xyz*b.yt + a.xyt*b.yz + a.xzt*b.q - a.yzt*b.xy ; 
	c.yzt =  + a.q*b.yzt - a.wy*b.wzt + a.wz*b.wyt - a.wt*b.wyz - a.xy*b.xzt + a.xz*b.xyt - a.xt*b.xyz
		 + a.wyz*b.wt - a.wyt*b.wz + a.wzt*b.wy + a.xyz*b.xt - a.xyt*b.xz + a.xzt*b.xy + a.yzt*b.q ; 

	c.wxyz =  + a.q*b.wxyz + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt
		 + a.wxyz*b.q + a.wxyt*b.zt - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt ; 
	c.wxyt =  + a.q*b.wxyt + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz
		 + a.wxyz*b.zt + a.wxyt*b.q - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz ; 
	c.wxzt =  + a.q*b.wxzt - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz
		 - a.wxyz*b.yt + a.wxyt*b.yz + a.wxzt*b.q - a.wyzt*b.xy + a.xyzt*b.wy ; 
	c.wyzt =  + a.q*b.wyzt + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz
		 + a.wxyz*b.xt - a.wxyt*b.xz + a.wxzt*b.xy + a.wyzt*b.q - a.xyzt*b.wx ; 
	c.xyzt =  + a.q*b.xyzt - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz
		 - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx + a.xyzt*b.q ; 

	c.wxyzt =  + a.q*b.wxyzt + a.wxyzt*b.q ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 Shrinker(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

// Shrinker equation set

	c.q =  + a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z - a.t*b.t - a.wx*b.wx - a.wy*b.wy
		 - a.wz*b.wz + a.wt*b.wt - a.xy*b.xy - a.xz*b.xz + a.xt*b.xt - a.yz*b.yz + a.yt*b.yt + a.zt*b.zt
		 - a.wxy*b.wxy - a.wxz*b.wxz + a.wxt*b.wxt - a.wyz*b.wyz + a.wyt*b.wyt + a.wzt*b.wzt - a.xyz*b.xyz
		 + a.xyt*b.xyt + a.xzt*b.xzt + a.yzt*b.yzt + a.wxyz*b.wxyz - a.wxyt*b.wxyt - a.wxzt*b.wxzt
		 - a.wyzt*b.wyzt - a.xyzt*b.xyzt - a.wxyzt*b.wxyzt ; 

	c.w =  - a.x*b.wx - a.y*b.wy - a.z*b.wz + a.t*b.wt + a.wx*b.x + a.wy*b.y + a.wz*b.z - a.wt*b.t
		 - a.xy*b.wxy - a.xz*b.wxz + a.xt*b.wxt - a.yz*b.wyz + a.yt*b.wyt + a.zt*b.wzt
		 - a.wxy*b.xy - a.wxz*b.xz + a.wxt*b.xt - a.wyz*b.yz + a.wyt*b.yt + a.wzt*b.zt
		 + a.xyz*b.wxyz - a.xyt*b.wxyt - a.xzt*b.wxzt - a.yzt*b.wyzt - a.wxyz*b.xyz
		 + a.wxyt*b.xyt + a.wxzt*b.xzt + a.wyzt*b.yzt - a.xyzt*b.wxyzt - a.wxyzt*b.xyzt ; 
	c.x =  + a.w*b.wx - a.y*b.xy - a.z*b.xz + a.t*b.xt - a.wx*b.w + a.wy*b.wxy + a.wz*b.wxz
		 - a.wt*b.wxt + a.xy*b.y + a.xz*b.z - a.xt*b.t - a.yz*b.xyz + a.yt*b.xyt + a.zt*b.xzt
		 + a.wxy*b.wy + a.wxz*b.wz - a.wxt*b.wt - a.wyz*b.wxyz + a.wyt*b.wxyt + a.wzt*b.wxzt
		 - a.xyz*b.yz + a.xyt*b.yt + a.xzt*b.zt - a.yzt*b.xyzt + a.wxyz*b.wyz - a.wxyt*b.wyt
		 - a.wxzt*b.wzt + a.wyzt*b.wxyzt + a.xyzt*b.yzt + a.wxyzt*b.wyzt ; 
	c.y =  + a.w*b.wy + a.x*b.xy - a.z*b.yz + a.t*b.yt - a.wx*b.wxy - a.wy*b.w + a.wz*b.wyz
		 - a.wt*b.wyt - a.xy*b.x + a.xz*b.xyz - a.xt*b.xyt + a.yz*b.z - a.yt*b.t + a.zt*b.yzt
		 - a.wxy*b.wx + a.wxz*b.wxyz - a.wxt*b.wxyt + a.wyz*b.wz - a.wyt*b.wt + a.wzt*b.wyzt
		 + a.xyz*b.xz - a.xyt*b.xt + a.xzt*b.xyzt + a.yzt*b.zt - a.wxyz*b.wxz + a.wxyt*b.wxt
		 - a.wxzt*b.wxyzt - a.wyzt*b.wzt - a.xyzt*b.xzt - a.wxyzt*b.wxzt ; 
	c.z =  + a.w*b.wz + a.x*b.xz + a.y*b.yz + a.t*b.zt - a.wx*b.wxz - a.wy*b.wyz - a.wz*b.w
		 - a.wt*b.wzt - a.xy*b.xyz - a.xz*b.x - a.xt*b.xzt - a.yz*b.y - a.yt*b.yzt - a.zt*b.t
		 - a.wxy*b.wxyz - a.wxz*b.wx - a.wxt*b.wxzt - a.wyz*b.wy - a.wyt*b.wyzt - a.wzt*b.wt
		 - a.xyz*b.xy - a.xyt*b.xyzt - a.xzt*b.xt - a.yzt*b.yt + a.wxyz*b.wxy + a.wxyt*b.wxyzt
		 + a.wxzt*b.wxt + a.wyzt*b.wyt + a.xyzt*b.xyt + a.wxyzt*b.wxyt ; 
	c.t =  + a.w*b.wt + a.x*b.xt + a.y*b.yt + a.z*b.zt - a.wx*b.wxt - a.wy*b.wyt - a.wz*b.wzt
		 - a.wt*b.w - a.xy*b.xyt - a.xz*b.xzt - a.xt*b.x - a.yz*b.yzt - a.yt*b.y - a.zt*b.z
		 - a.wxy*b.wxyt - a.wxz*b.wxzt - a.wxt*b.wx - a.wyz*b.wyzt - a.wyt*b.wy - a.wzt*b.wz
		 - a.xyz*b.xyzt - a.xyt*b.xy - a.xzt*b.xz - a.yzt*b.yz + a.wxyz*b.wxyzt + a.wxyt*b.wxy
		 + a.wxzt*b.wxz + a.wyzt*b.wyz + a.xyzt*b.xyz + a.wxyzt*b.wxyz ; 

	c.wx =  + a.y*b.wxy + a.z*b.wxz - a.t*b.wxt - a.yz*b.wxyz + a.yt*b.wxyt + a.zt*b.wxzt + a.wxy*b.y
		 + a.wxz*b.z - a.wxt*b.t - a.wyz*b.xyz + a.wyt*b.xyt + a.wzt*b.xzt + a.xyz*b.wyz - a.xyt*b.wyt
		 - a.xzt*b.wzt + a.yzt*b.wxyzt - a.wxyz*b.yz + a.wxyt*b.yt + a.wxzt*b.zt - a.wyzt*b.xyzt
		 + a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 
	c.wy =  - a.x*b.wxy + a.z*b.wyz - a.t*b.wyt + a.xz*b.wxyz - a.xt*b.wxyt + a.zt*b.wyzt - a.wxy*b.x
		 + a.wxz*b.xyz - a.wxt*b.xyt + a.wyz*b.z - a.wyt*b.t + a.wzt*b.yzt - a.xyz*b.wxz + a.xyt*b.wxt
		 - a.xzt*b.wxyzt - a.yzt*b.wzt + a.wxyz*b.xz - a.wxyt*b.xt + a.wxzt*b.xyzt + a.wyzt*b.zt
		 - a.xyzt*b.wxzt - a.wxyzt*b.xzt ; 
	c.wz =  - a.x*b.wxz - a.y*b.wyz - a.t*b.wzt - a.xy*b.wxyz - a.xt*b.wxzt - a.yt*b.wyzt - a.wxy*b.xyz
		 - a.wxz*b.x - a.wxt*b.xzt - a.wyz*b.y - a.wyt*b.yzt - a.wzt*b.t + a.xyz*b.wxy + a.xyt*b.wxyzt
		 + a.xzt*b.wxt + a.yzt*b.wyt - a.wxyz*b.xy - a.wxyt*b.xyzt - a.wxzt*b.xt - a.wyzt*b.yt
		 + a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.wt =  - a.x*b.wxt - a.y*b.wyt - a.z*b.wzt - a.xy*b.wxyt - a.xz*b.wxzt - a.yz*b.wyzt - a.wxy*b.xyt
		 - a.wxz*b.xzt - a.wxt*b.x - a.wyz*b.yzt - a.wyt*b.y - a.wzt*b.z + a.xyz*b.wxyzt + a.xyt*b.wxy
		 + a.xzt*b.wxz + a.yzt*b.wyz - a.wxyz*b.xyzt - a.wxyt*b.xy - a.wxzt*b.xz - a.wyzt*b.yz
		 + a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xy =  + a.w*b.wxy + a.z*b.xyz - a.t*b.xyt - a.wz*b.wxyz + a.wt*b.wxyt + a.zt*b.xyzt + a.wxy*b.w
		 - a.wxz*b.wyz + a.wxt*b.wyt + a.wyz*b.wxz - a.wyt*b.wxt + a.wzt*b.wxyzt + a.xyz*b.z - a.xyt*b.t
		 + a.xzt*b.yzt - a.yzt*b.xzt - a.wxyz*b.wz + a.wxyt*b.wt - a.wxzt*b.wyzt + a.wyzt*b.wxzt
		 + a.xyzt*b.zt + a.wxyzt*b.wzt ; 
	c.xz =  + a.w*b.wxz - a.y*b.xyz - a.t*b.xzt + a.wy*b.wxyz + a.wt*b.wxzt - a.yt*b.xyzt + a.wxy*b.wyz
		 + a.wxz*b.w + a.wxt*b.wzt - a.wyz*b.wxy - a.wyt*b.wxyzt - a.wzt*b.wxt - a.xyz*b.y - a.xyt*b.yzt
		 - a.xzt*b.t + a.yzt*b.xyt + a.wxyz*b.wy + a.wxyt*b.wyzt + a.wxzt*b.wt - a.wyzt*b.wxyt
		 - a.xyzt*b.yt - a.wxyzt*b.wyt ; 
	c.xt =  + a.w*b.wxt - a.y*b.xyt - a.z*b.xzt + a.wy*b.wxyt + a.wz*b.wxzt - a.yz*b.xyzt + a.wxy*b.wyt
		 + a.wxz*b.wzt + a.wxt*b.w - a.wyz*b.wxyzt - a.wyt*b.wxy - a.wzt*b.wxz - a.xyz*b.yzt - a.xyt*b.y
		 - a.xzt*b.z + a.yzt*b.xyz + a.wxyz*b.wyzt + a.wxyt*b.wy + a.wxzt*b.wz - a.wyzt*b.wxyz
		 - a.xyzt*b.yz - a.wxyzt*b.wyz ; 
	c.yz =  + a.w*b.wyz + a.x*b.xyz - a.t*b.yzt - a.wx*b.wxyz + a.wt*b.wyzt + a.xt*b.xyzt - a.wxy*b.wxz
		 + a.wxz*b.wxy + a.wxt*b.wxyzt + a.wyz*b.w + a.wyt*b.wzt - a.wzt*b.wyt + a.xyz*b.x + a.xyt*b.xzt
		 - a.xzt*b.xyt - a.yzt*b.t - a.wxyz*b.wx - a.wxyt*b.wxzt + a.wxzt*b.wxyt + a.wyzt*b.wt
		 + a.xyzt*b.xt + a.wxyzt*b.wxt ; 
	c.yt =  + a.w*b.wyt + a.x*b.xyt - a.z*b.yzt - a.wx*b.wxyt + a.wz*b.wyzt + a.xz*b.xyzt - a.wxy*b.wxt
		 + a.wxz*b.wxyzt + a.wxt*b.wxy + a.wyz*b.wzt + a.wyt*b.w - a.wzt*b.wyz + a.xyz*b.xzt + a.xyt*b.x
		 - a.xzt*b.xyz - a.yzt*b.z - a.wxyz*b.wxzt - a.wxyt*b.wx + a.wxzt*b.wxyz + a.wyzt*b.wz
		 + a.xyzt*b.xz + a.wxyzt*b.wxz ; 
	c.zt =  + a.w*b.wzt + a.x*b.xzt + a.y*b.yzt - a.wx*b.wxzt - a.wy*b.wyzt - a.xy*b.xyzt - a.wxy*b.wxyzt
		 - a.wxz*b.wxt + a.wxt*b.wxz - a.wyz*b.wyt + a.wyt*b.wyz + a.wzt*b.w - a.xyz*b.xyt + a.xyt*b.xyz
		 + a.xzt*b.x + a.yzt*b.y + a.wxyz*b.wxyt - a.wxyt*b.wxyz - a.wxzt*b.wx - a.wyzt*b.wy
		 - a.xyzt*b.xy - a.wxyzt*b.wxy ; 

	c.wxy =  - a.z*b.wxyz + a.t*b.wxyt + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt
		 + a.wxyz*b.z - a.wxyt*b.t + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 
	c.wxz =  + a.y*b.wxyz + a.t*b.wxzt - a.yt*b.wxyzt - a.wyt*b.xyzt + a.xyt*b.wyzt - a.yzt*b.wxyt
		 - a.wxyz*b.y - a.wxyt*b.yzt - a.wxzt*b.t + a.wyzt*b.xyt - a.xyzt*b.wyt - a.wxyzt*b.yt ; 
	c.wxt =  + a.y*b.wxyt + a.z*b.wxzt - a.yz*b.wxyzt - a.wyz*b.xyzt + a.xyz*b.wyzt - a.yzt*b.wxyz
		 - a.wxyz*b.yzt - a.wxyt*b.y - a.wxzt*b.z + a.wyzt*b.xyz - a.xyzt*b.wyz - a.wxyzt*b.yz ; 
	c.wyz =  - a.x*b.wxyz + a.t*b.wyzt + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt
		 + a.wxyz*b.x + a.wxyt*b.xzt - a.wxzt*b.xyt - a.wyzt*b.t + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.wyt =  - a.x*b.wxyt + a.z*b.wyzt + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz
		 + a.wxyz*b.xzt + a.wxyt*b.x - a.wxzt*b.xyz - a.wyzt*b.z + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.wzt =  - a.x*b.wxzt - a.y*b.wyzt - a.xy*b.wxyzt - a.wxy*b.xyzt + a.xyz*b.wxyt - a.xyt*b.wxyz
		 - a.wxyz*b.xyt + a.wxyt*b.xyz + a.wxzt*b.x + a.wyzt*b.y - a.xyzt*b.wxy - a.wxyzt*b.xy ; 
	c.xyz =  + a.w*b.wxyz + a.t*b.xyzt - a.wt*b.wxyzt - a.wxt*b.wyzt + a.wyt*b.wxzt - a.wzt*b.wxyt
		 - a.wxyz*b.w - a.wxyt*b.wzt + a.wxzt*b.wyt - a.wyzt*b.wxt - a.xyzt*b.t - a.wxyzt*b.wt ; 
	c.xyt =  + a.w*b.wxyt + a.z*b.xyzt - a.wz*b.wxyzt - a.wxz*b.wyzt + a.wyz*b.wxzt - a.wzt*b.wxyz
		 - a.wxyz*b.wzt - a.wxyt*b.w + a.wxzt*b.wyz - a.wyzt*b.wxz - a.xyzt*b.z - a.wxyzt*b.wz ; 
	c.xzt =  + a.w*b.wxzt - a.y*b.xyzt + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz
		 + a.wxyz*b.wyt - a.wxyt*b.wyz - a.wxzt*b.w + a.wyzt*b.wxy + a.xyzt*b.y + a.wxyzt*b.wy ; 
	c.yzt =  + a.w*b.wyzt + a.x*b.xyzt - a.wx*b.wxyzt - a.wxy*b.wxzt + a.wxz*b.wxyt - a.wxt*b.wxyz
		 - a.wxyz*b.wxt + a.wxyt*b.wxz - a.wxzt*b.wxy - a.wyzt*b.w - a.xyzt*b.x - a.wxyzt*b.wx ; 

	c.wxyz =  - a.t*b.wxyzt - a.wxyzt*b.t ; 
	c.wxyt =  - a.z*b.wxyzt - a.wxyzt*b.z ; 
	c.wxzt =  + a.y*b.wxyzt + a.wxyzt*b.y ; 
	c.wyzt =  - a.x*b.wxyzt - a.wxyzt*b.x ; 
	c.xyzt =  + a.w*b.wxyzt + a.wxyzt*b.w ; 

	c.wxyzt = 0 ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Symmetric(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

//	c = a*b/2 + b*a/2;


// Symmetric Product
c.q =  + a.q*b.q + a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z - a.t*b.t - a.wx*b.wx - a.wy*b.wy - a.wz*b.wz + a.wt*b.wt - a.xy*b.xy - a.xz*b.xz + a.xt*b.xt - a.yz*b.yz + a.yt*b.yt + a.zt*b.zt - a.wxy*b.wxy - a.wxz*b.wxz + a.wxt*b.wxt - a.wyz*b.wyz + a.wyt*b.wyt + a.wzt*b.wzt - a.xyz*b.xyz + a.xyt*b.xyt + a.xzt*b.xzt + a.yzt*b.yzt + a.wxyz*b.wxyz - a.wxyt*b.wxyt - a.wxzt*b.wxzt - a.wyzt*b.wyzt - a.xyzt*b.xyzt - a.wxyzt*b.wxyzt ; 
c.w =  + a.q*b.w + a.w*b.q - a.xy*b.wxy - a.xz*b.wxz + a.xt*b.wxt - a.yz*b.wyz + a.yt*b.wyt + a.zt*b.wzt - a.wxy*b.xy - a.wxz*b.xz + a.wxt*b.xt - a.wyz*b.yz + a.wyt*b.yt + a.wzt*b.zt - a.xyzt*b.wxyzt - a.wxyzt*b.xyzt ; 
c.x =  + a.q*b.x + a.x*b.q + a.wy*b.wxy + a.wz*b.wxz - a.wt*b.wxt - a.yz*b.xyz + a.yt*b.xyt + a.zt*b.xzt + a.wxy*b.wy + a.wxz*b.wz - a.wxt*b.wt - a.xyz*b.yz + a.xyt*b.yt + a.xzt*b.zt + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
c.y =  + a.q*b.y + a.y*b.q - a.wx*b.wxy + a.wz*b.wyz - a.wt*b.wyt + a.xz*b.xyz - a.xt*b.xyt + a.zt*b.yzt - a.wxy*b.wx + a.wyz*b.wz - a.wyt*b.wt + a.xyz*b.xz - a.xyt*b.xt + a.yzt*b.zt - a.wxzt*b.wxyzt - a.wxyzt*b.wxzt ; 
c.z =  + a.q*b.z + a.z*b.q - a.wx*b.wxz - a.wy*b.wyz - a.wt*b.wzt - a.xy*b.xyz - a.xt*b.xzt - a.yt*b.yzt - a.wxz*b.wx - a.wyz*b.wy - a.wzt*b.wt - a.xyz*b.xy - a.xzt*b.xt - a.yzt*b.yt + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
c.t =  + a.q*b.t + a.t*b.q - a.wx*b.wxt - a.wy*b.wyt - a.wz*b.wzt - a.xy*b.xyt - a.xz*b.xzt - a.yz*b.yzt - a.wxt*b.wx - a.wyt*b.wy - a.wzt*b.wz - a.xyt*b.xy - a.xzt*b.xz - a.yzt*b.yz + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
c.wx =  + a.q*b.wx + a.y*b.wxy + a.z*b.wxz - a.t*b.wxt + a.wx*b.q - a.yz*b.wxyz + a.yt*b.wxyt + a.zt*b.wxzt + a.wxy*b.y + a.wxz*b.z - a.wxt*b.t + a.yzt*b.wxyzt - a.wxyz*b.yz + a.wxyt*b.yt + a.wxzt*b.zt + a.wxyzt*b.yzt ; 
c.wy =  + a.q*b.wy - a.x*b.wxy + a.z*b.wyz - a.t*b.wyt + a.wy*b.q + a.xz*b.wxyz - a.xt*b.wxyt + a.zt*b.wyzt - a.wxy*b.x + a.wyz*b.z - a.wyt*b.t - a.xzt*b.wxyzt + a.wxyz*b.xz - a.wxyt*b.xt + a.wyzt*b.zt - a.wxyzt*b.xzt ; 
c.wz =  + a.q*b.wz - a.x*b.wxz - a.y*b.wyz - a.t*b.wzt + a.wz*b.q - a.xy*b.wxyz - a.xt*b.wxzt - a.yt*b.wyzt - a.wxz*b.x - a.wyz*b.y - a.wzt*b.t + a.xyt*b.wxyzt - a.wxyz*b.xy - a.wxzt*b.xt - a.wyzt*b.yt + a.wxyzt*b.xyt ; 
c.wt =  + a.q*b.wt - a.x*b.wxt - a.y*b.wyt - a.z*b.wzt + a.wt*b.q - a.xy*b.wxyt - a.xz*b.wxzt - a.yz*b.wyzt - a.wxt*b.x - a.wyt*b.y - a.wzt*b.z + a.xyz*b.wxyzt - a.wxyt*b.xy - a.wxzt*b.xz - a.wyzt*b.yz + a.wxyzt*b.xyz ; 
c.xy =  + a.q*b.xy + a.w*b.wxy + a.z*b.xyz - a.t*b.xyt - a.wz*b.wxyz + a.wt*b.wxyt + a.xy*b.q + a.zt*b.xyzt + a.wxy*b.w + a.wzt*b.wxyzt + a.xyz*b.z - a.xyt*b.t - a.wxyz*b.wz + a.wxyt*b.wt + a.xyzt*b.zt + a.wxyzt*b.wzt ; 
c.xz =  + a.q*b.xz + a.w*b.wxz - a.y*b.xyz - a.t*b.xzt + a.wy*b.wxyz + a.wt*b.wxzt + a.xz*b.q - a.yt*b.xyzt + a.wxz*b.w - a.wyt*b.wxyzt - a.xyz*b.y - a.xzt*b.t + a.wxyz*b.wy + a.wxzt*b.wt - a.xyzt*b.yt - a.wxyzt*b.wyt ; 
c.xt =  + a.q*b.xt + a.w*b.wxt - a.y*b.xyt - a.z*b.xzt + a.wy*b.wxyt + a.wz*b.wxzt + a.xt*b.q - a.yz*b.xyzt + a.wxt*b.w - a.wyz*b.wxyzt - a.xyt*b.y - a.xzt*b.z + a.wxyt*b.wy + a.wxzt*b.wz - a.xyzt*b.yz - a.wxyzt*b.wyz ; 
c.yz =  + a.q*b.yz + a.w*b.wyz + a.x*b.xyz - a.t*b.yzt - a.wx*b.wxyz + a.wt*b.wyzt + a.xt*b.xyzt + a.yz*b.q + a.wxt*b.wxyzt + a.wyz*b.w + a.xyz*b.x - a.yzt*b.t - a.wxyz*b.wx + a.wyzt*b.wt + a.xyzt*b.xt + a.wxyzt*b.wxt ; 
c.yt =  + a.q*b.yt + a.w*b.wyt + a.x*b.xyt - a.z*b.yzt - a.wx*b.wxyt + a.wz*b.wyzt + a.xz*b.xyzt + a.yt*b.q + a.wxz*b.wxyzt + a.wyt*b.w + a.xyt*b.x - a.yzt*b.z - a.wxyt*b.wx + a.wyzt*b.wz + a.xyzt*b.xz + a.wxyzt*b.wxz ; 
c.zt =  + a.q*b.zt + a.w*b.wzt + a.x*b.xzt + a.y*b.yzt - a.wx*b.wxzt - a.wy*b.wyzt - a.xy*b.xyzt + a.zt*b.q - a.wxy*b.wxyzt + a.wzt*b.w + a.xzt*b.x + a.yzt*b.y - a.wxzt*b.wx - a.wyzt*b.wy - a.xyzt*b.xy - a.wxyzt*b.wxy ; 
c.wxy =  + a.q*b.wxy + a.w*b.xy - a.x*b.wy + a.y*b.wx + a.wx*b.y - a.wy*b.x + a.xy*b.w + a.zt*b.wxyzt + a.wxy*b.q + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 
c.wxz =  + a.q*b.wxz + a.w*b.xz - a.x*b.wz + a.z*b.wx + a.wx*b.z - a.wz*b.x + a.xz*b.w - a.yt*b.wxyzt + a.wxz*b.q - a.wyt*b.xyzt + a.xyt*b.wyzt - a.yzt*b.wxyt - a.wxyt*b.yzt + a.wyzt*b.xyt - a.xyzt*b.wyt - a.wxyzt*b.yt ; 
c.wxt =  + a.q*b.wxt + a.w*b.xt - a.x*b.wt + a.t*b.wx + a.wx*b.t - a.wt*b.x + a.xt*b.w - a.yz*b.wxyzt + a.wxt*b.q - a.wyz*b.xyzt + a.xyz*b.wyzt - a.yzt*b.wxyz - a.wxyz*b.yzt + a.wyzt*b.xyz - a.xyzt*b.wyz - a.wxyzt*b.yz ; 
c.wyz =  + a.q*b.wyz + a.w*b.yz - a.y*b.wz + a.z*b.wy + a.wy*b.z - a.wz*b.y + a.xt*b.wxyzt + a.yz*b.w + a.wxt*b.xyzt + a.wyz*b.q - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
c.wyt =  + a.q*b.wyt + a.w*b.yt - a.y*b.wt + a.t*b.wy + a.wy*b.t - a.wt*b.y + a.xz*b.wxyzt + a.yt*b.w + a.wxz*b.xyzt + a.wyt*b.q - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
c.wzt =  + a.q*b.wzt + a.w*b.zt - a.z*b.wt + a.t*b.wz + a.wz*b.t - a.wt*b.z - a.xy*b.wxyzt + a.zt*b.w - a.wxy*b.xyzt + a.wzt*b.q + a.xyz*b.wxyt - a.xyt*b.wxyz - a.wxyz*b.xyt + a.wxyt*b.xyz - a.xyzt*b.wxy - a.wxyzt*b.xy ; 
c.xyz =  + a.q*b.xyz + a.x*b.yz - a.y*b.xz + a.z*b.xy - a.wt*b.wxyzt + a.xy*b.z - a.xz*b.y + a.yz*b.x - a.wxt*b.wyzt + a.wyt*b.wxzt - a.wzt*b.wxyt + a.xyz*b.q - a.wxyt*b.wzt + a.wxzt*b.wyt - a.wyzt*b.wxt - a.wxyzt*b.wt ; 
c.xyt =  + a.q*b.xyt + a.x*b.yt - a.y*b.xt + a.t*b.xy - a.wz*b.wxyzt + a.xy*b.t - a.xt*b.y + a.yt*b.x - a.wxz*b.wyzt + a.wyz*b.wxzt - a.wzt*b.wxyz + a.xyt*b.q - a.wxyz*b.wzt + a.wxzt*b.wyz - a.wyzt*b.wxz - a.wxyzt*b.wz ; 
c.xzt =  + a.q*b.xzt + a.x*b.zt - a.z*b.xt + a.t*b.xz + a.wy*b.wxyzt + a.xz*b.t - a.xt*b.z + a.zt*b.x + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.xzt*b.q + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
c.yzt =  + a.q*b.yzt + a.y*b.zt - a.z*b.yt + a.t*b.yz - a.wx*b.wxyzt + a.yz*b.t - a.yt*b.z + a.zt*b.y - a.wxy*b.wxzt + a.wxz*b.wxyt - a.wxt*b.wxyz + a.yzt*b.q - a.wxyz*b.wxt + a.wxyt*b.wxz - a.wxzt*b.wxy - a.wxyzt*b.wx ; 
c.wxyz =  + a.q*b.wxyz - a.t*b.wxyzt + a.wx*b.yz - a.wy*b.xz + a.wz*b.xy + a.xy*b.wz - a.xz*b.wy + a.yz*b.wx - a.wxt*b.yzt + a.wyt*b.xzt - a.wzt*b.xyt - a.xyt*b.wzt + a.xzt*b.wyt - a.yzt*b.wxt + a.wxyz*b.q - a.wxyzt*b.t ; 
c.wxyt =  + a.q*b.wxyt - a.z*b.wxyzt + a.wx*b.yt - a.wy*b.xt + a.wt*b.xy + a.xy*b.wt - a.xt*b.wy + a.yt*b.wx - a.wxz*b.yzt + a.wyz*b.xzt - a.wzt*b.xyz - a.xyz*b.wzt + a.xzt*b.wyz - a.yzt*b.wxz + a.wxyt*b.q - a.wxyzt*b.z ; 
c.wxzt =  + a.q*b.wxzt + a.y*b.wxyzt + a.wx*b.zt - a.wz*b.xt + a.wt*b.xz + a.xz*b.wt - a.xt*b.wz + a.zt*b.wx + a.wxy*b.yzt - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz + a.yzt*b.wxy + a.wxzt*b.q + a.wxyzt*b.y ; 
c.wyzt =  + a.q*b.wyzt - a.x*b.wxyzt + a.wy*b.zt - a.wz*b.yt + a.wt*b.yz + a.yz*b.wt - a.yt*b.wz + a.zt*b.wy - a.wxy*b.xzt + a.wxz*b.xyt - a.wxt*b.xyz - a.xyz*b.wxt + a.xyt*b.wxz - a.xzt*b.wxy + a.wyzt*b.q - a.wxyzt*b.x ; 
c.xyzt =  + a.q*b.xyzt + a.w*b.wxyzt + a.xy*b.zt - a.xz*b.yt + a.xt*b.yz + a.yz*b.xt - a.yt*b.xz + a.zt*b.xy + a.wxy*b.wzt - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.xyzt*b.q + a.wxyzt*b.w ; 
c.wxyzt =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz + a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt + a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w + a.wxyzt*b.q ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 AntiSymmetric(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

//	c = a*b/2 - b*a/2;


// AntiSymmetric Product
c.q =  0; 
c.w =  - a.x*b.wx - a.y*b.wy - a.z*b.wz + a.t*b.wt + a.wx*b.x + a.wy*b.y + a.wz*b.z - a.wt*b.t + a.xyz*b.wxyz - a.xyt*b.wxyt - a.xzt*b.wxzt - a.yzt*b.wyzt - a.wxyz*b.xyz + a.wxyt*b.xyt + a.wxzt*b.xzt + a.wyzt*b.yzt ; 
c.x =  + a.w*b.wx - a.y*b.xy - a.z*b.xz + a.t*b.xt - a.wx*b.w + a.xy*b.y + a.xz*b.z - a.xt*b.t - a.wyz*b.wxyz + a.wyt*b.wxyt + a.wzt*b.wxzt - a.yzt*b.xyzt + a.wxyz*b.wyz - a.wxyt*b.wyt - a.wxzt*b.wzt + a.xyzt*b.yzt ; 
c.y =  + a.w*b.wy + a.x*b.xy - a.z*b.yz + a.t*b.yt - a.wy*b.w - a.xy*b.x + a.yz*b.z - a.yt*b.t + a.wxz*b.wxyz - a.wxt*b.wxyt + a.wzt*b.wyzt + a.xzt*b.xyzt - a.wxyz*b.wxz + a.wxyt*b.wxt - a.wyzt*b.wzt - a.xyzt*b.xzt ; 
c.z =  + a.w*b.wz + a.x*b.xz + a.y*b.yz + a.t*b.zt - a.wz*b.w - a.xz*b.x - a.yz*b.y - a.zt*b.t - a.wxy*b.wxyz - a.wxt*b.wxzt - a.wyt*b.wyzt - a.xyt*b.xyzt + a.wxyz*b.wxy + a.wxzt*b.wxt + a.wyzt*b.wyt + a.xyzt*b.xyt ; 
c.t =  + a.w*b.wt + a.x*b.xt + a.y*b.yt + a.z*b.zt - a.wt*b.w - a.xt*b.x - a.yt*b.y - a.zt*b.z - a.wxy*b.wxyt - a.wxz*b.wxzt - a.wyz*b.wyzt - a.xyz*b.xyzt + a.wxyt*b.wxy + a.wxzt*b.wxz + a.wyzt*b.wyz + a.xyzt*b.xyz ; 
c.wx =  + a.w*b.x - a.x*b.w - a.wy*b.xy - a.wz*b.xz + a.wt*b.xt + a.xy*b.wy + a.xz*b.wz - a.xt*b.wt - a.wyz*b.xyz + a.wyt*b.xyt + a.wzt*b.xzt + a.xyz*b.wyz - a.xyt*b.wyt - a.xzt*b.wzt - a.wyzt*b.xyzt + a.xyzt*b.wyzt ; 
c.wy =  + a.w*b.y - a.y*b.w + a.wx*b.xy - a.wz*b.yz + a.wt*b.yt - a.xy*b.wx + a.yz*b.wz - a.yt*b.wt + a.wxz*b.xyz - a.wxt*b.xyt + a.wzt*b.yzt - a.xyz*b.wxz + a.xyt*b.wxt - a.yzt*b.wzt + a.wxzt*b.xyzt - a.xyzt*b.wxzt ; 
c.wz =  + a.w*b.z - a.z*b.w + a.wx*b.xz + a.wy*b.yz + a.wt*b.zt - a.xz*b.wx - a.yz*b.wy - a.zt*b.wt - a.wxy*b.xyz - a.wxt*b.xzt - a.wyt*b.yzt + a.xyz*b.wxy + a.xzt*b.wxt + a.yzt*b.wyt - a.wxyt*b.xyzt + a.xyzt*b.wxyt ; 
c.wt =  + a.w*b.t - a.t*b.w + a.wx*b.xt + a.wy*b.yt + a.wz*b.zt - a.xt*b.wx - a.yt*b.wy - a.zt*b.wz - a.wxy*b.xyt - a.wxz*b.xzt - a.wyz*b.yzt + a.xyt*b.wxy + a.xzt*b.wxz + a.yzt*b.wyz - a.wxyz*b.xyzt + a.xyzt*b.wxyz ; 
c.xy =  + a.x*b.y - a.y*b.x - a.wx*b.wy + a.wy*b.wx - a.xz*b.yz + a.xt*b.yt + a.yz*b.xz - a.yt*b.xt - a.wxz*b.wyz + a.wxt*b.wyt + a.wyz*b.wxz - a.wyt*b.wxt + a.xzt*b.yzt - a.yzt*b.xzt - a.wxzt*b.wyzt + a.wyzt*b.wxzt ; 
c.xz =  + a.x*b.z - a.z*b.x - a.wx*b.wz + a.wz*b.wx + a.xy*b.yz + a.xt*b.zt - a.yz*b.xy - a.zt*b.xt + a.wxy*b.wyz + a.wxt*b.wzt - a.wyz*b.wxy - a.wzt*b.wxt - a.xyt*b.yzt + a.yzt*b.xyt + a.wxyt*b.wyzt - a.wyzt*b.wxyt ; 
c.xt =  + a.x*b.t - a.t*b.x - a.wx*b.wt + a.wt*b.wx + a.xy*b.yt + a.xz*b.zt - a.yt*b.xy - a.zt*b.xz + a.wxy*b.wyt + a.wxz*b.wzt - a.wyt*b.wxy - a.wzt*b.wxz - a.xyz*b.yzt + a.yzt*b.xyz + a.wxyz*b.wyzt - a.wyzt*b.wxyz ; 
c.yz =  + a.y*b.z - a.z*b.y - a.wy*b.wz + a.wz*b.wy - a.xy*b.xz + a.xz*b.xy + a.yt*b.zt - a.zt*b.yt - a.wxy*b.wxz + a.wxz*b.wxy + a.wyt*b.wzt - a.wzt*b.wyt + a.xyt*b.xzt - a.xzt*b.xyt - a.wxyt*b.wxzt + a.wxzt*b.wxyt ; 
c.yt =  + a.y*b.t - a.t*b.y - a.wy*b.wt + a.wt*b.wy - a.xy*b.xt + a.xt*b.xy + a.yz*b.zt - a.zt*b.yz - a.wxy*b.wxt + a.wxt*b.wxy + a.wyz*b.wzt - a.wzt*b.wyz + a.xyz*b.xzt - a.xzt*b.xyz - a.wxyz*b.wxzt + a.wxzt*b.wxyz ; 
c.zt =  + a.z*b.t - a.t*b.z - a.wz*b.wt + a.wt*b.wz - a.xz*b.xt + a.xt*b.xz - a.yz*b.yt + a.yt*b.yz - a.wxz*b.wxt + a.wxt*b.wxz - a.wyz*b.wyt + a.wyt*b.wyz - a.xyz*b.xyt + a.xyt*b.xyz + a.wxyz*b.wxyt - a.wxyt*b.wxyz ; 
c.wxy =  - a.z*b.wxyz + a.t*b.wxyt + a.wz*b.xyz - a.wt*b.xyt - a.xz*b.wyz + a.xt*b.wyt + a.yz*b.wxz - a.yt*b.wxt - a.wxz*b.yz + a.wxt*b.yt + a.wyz*b.xz - a.wyt*b.xt - a.xyz*b.wz + a.xyt*b.wt + a.wxyz*b.z - a.wxyt*b.t ; 
c.wxz =  + a.y*b.wxyz + a.t*b.wxzt - a.wy*b.xyz - a.wt*b.xzt + a.xy*b.wyz + a.xt*b.wzt - a.yz*b.wxy - a.zt*b.wxt + a.wxy*b.yz + a.wxt*b.zt - a.wyz*b.xy - a.wzt*b.xt + a.xyz*b.wy + a.xzt*b.wt - a.wxyz*b.y - a.wxzt*b.t ; 
c.wxt =  + a.y*b.wxyt + a.z*b.wxzt - a.wy*b.xyt - a.wz*b.xzt + a.xy*b.wyt + a.xz*b.wzt - a.yt*b.wxy - a.zt*b.wxz + a.wxy*b.yt + a.wxz*b.zt - a.wyt*b.xy - a.wzt*b.xz + a.xyt*b.wy + a.xzt*b.wz - a.wxyt*b.y - a.wxzt*b.z ; 
c.wyz =  - a.x*b.wxyz + a.t*b.wyzt + a.wx*b.xyz - a.wt*b.yzt - a.xy*b.wxz + a.xz*b.wxy + a.yt*b.wzt - a.zt*b.wyt - a.wxy*b.xz + a.wxz*b.xy + a.wyt*b.zt - a.wzt*b.yt - a.xyz*b.wx + a.yzt*b.wt + a.wxyz*b.x - a.wyzt*b.t ; 
c.wyt =  - a.x*b.wxyt + a.z*b.wyzt + a.wx*b.xyt - a.wz*b.yzt - a.xy*b.wxt + a.xt*b.wxy + a.yz*b.wzt - a.zt*b.wyz - a.wxy*b.xt + a.wxt*b.xy + a.wyz*b.zt - a.wzt*b.yz - a.xyt*b.wx + a.yzt*b.wz + a.wxyt*b.x - a.wyzt*b.z ; 
c.wzt =  - a.x*b.wxzt - a.y*b.wyzt + a.wx*b.xzt + a.wy*b.yzt - a.xz*b.wxt + a.xt*b.wxz - a.yz*b.wyt + a.yt*b.wyz - a.wxz*b.xt + a.wxt*b.xz - a.wyz*b.yt + a.wyt*b.yz - a.xzt*b.wx - a.yzt*b.wy + a.wxzt*b.x + a.wyzt*b.y ; 
c.xyz =  + a.w*b.wxyz + a.t*b.xyzt - a.wx*b.wyz + a.wy*b.wxz - a.wz*b.wxy - a.xt*b.yzt + a.yt*b.xzt - a.zt*b.xyt + a.wxy*b.wz - a.wxz*b.wy + a.wyz*b.wx + a.xyt*b.zt - a.xzt*b.yt + a.yzt*b.xt - a.wxyz*b.w - a.xyzt*b.t ; 
c.xyt =  + a.w*b.wxyt + a.z*b.xyzt - a.wx*b.wyt + a.wy*b.wxt - a.wt*b.wxy - a.xz*b.yzt + a.yz*b.xzt - a.zt*b.xyz + a.wxy*b.wt - a.wxt*b.wy + a.wyt*b.wx + a.xyz*b.zt - a.xzt*b.yz + a.yzt*b.xz - a.wxyt*b.w - a.xyzt*b.z ; 
c.xzt =  + a.w*b.wxzt - a.y*b.xyzt - a.wx*b.wzt + a.wz*b.wxt - a.wt*b.wxz + a.xy*b.yzt - a.yz*b.xyt + a.yt*b.xyz + a.wxz*b.wt - a.wxt*b.wz + a.wzt*b.wx - a.xyz*b.yt + a.xyt*b.yz - a.yzt*b.xy - a.wxzt*b.w + a.xyzt*b.y ; 
c.yzt =  + a.w*b.wyzt + a.x*b.xyzt - a.wy*b.wzt + a.wz*b.wyt - a.wt*b.wyz - a.xy*b.xzt + a.xz*b.xyt - a.xt*b.xyz + a.wyz*b.wt - a.wyt*b.wz + a.wzt*b.wy + a.xyz*b.xt - a.xyt*b.xz + a.xzt*b.xy - a.wyzt*b.w - a.xyzt*b.x ; 
c.wxyz =  + a.w*b.xyz - a.x*b.wyz + a.y*b.wxz - a.z*b.wxy + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt + a.wxy*b.z - a.wxz*b.y + a.wyz*b.x - a.xyz*b.w + a.wxyt*b.zt - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt ; 
c.wxyt =  + a.w*b.xyt - a.x*b.wyt + a.y*b.wxt - a.t*b.wxy + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz + a.wxy*b.t - a.wxt*b.y + a.wyt*b.x - a.xyt*b.w + a.wxyz*b.zt - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz ; 
c.wxzt =  + a.w*b.xzt - a.x*b.wzt + a.z*b.wxt - a.t*b.wxz - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz + a.wxz*b.t - a.wxt*b.z + a.wzt*b.x - a.xzt*b.w - a.wxyz*b.yt + a.wxyt*b.yz - a.wyzt*b.xy + a.xyzt*b.wy ; 
c.wyzt =  + a.w*b.yzt - a.y*b.wzt + a.z*b.wyt - a.t*b.wyz + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz + a.wyz*b.t - a.wyt*b.z + a.wzt*b.y - a.yzt*b.w + a.wxyz*b.xt - a.wxyt*b.xz + a.wxzt*b.xy - a.xyzt*b.wx ; 
c.xyzt =  + a.x*b.yzt - a.y*b.xzt + a.z*b.xyt - a.t*b.xyz - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz + a.xyz*b.t - a.xyt*b.z + a.xzt*b.y - a.yzt*b.x - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx ; 
c.wxyzt =  0; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 Inner(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

c.q = -a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.t*b.wt - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.t*b.xt + a.w*b.wx - a.wt*b.wxt - a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt - a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt + a.wxz*b.wxyz - a.wxzt*b.wxyzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt + a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wy*b.wyz - a.wyt*b.wyzt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt - a.zt*b.t) ; 

c.t = + (a.w*b.wt - a.wt*b.w - a.wx*b.wxt - a.wxt*b.wx - a.wxy*b.wxyt + a.wxyt*b.wxy + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz - a.wxz*b.wxzt + a.wxzt*b.wxz - a.wy*b.wyt - a.wyt*b.wy - a.wyz*b.wyzt + a.wyzt*b.wyz - a.wz*b.wzt - a.wzt*b.wz + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z) ; 



c.wx = + (-a.t*b.wxt - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 

c.wy = + (-a.t*b.wyt - a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt - a.x*b.wxy - a.xt*b.wxyt + a.xz*b.wxyz - a.xzt*b.wxyzt + a.z*b.wyz + a.zt*b.wyzt) ; 

c.wz = + (-a.t*b.wzt - a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt - a.wyz*b.y - a.wyzt*b.yt - a.wzt*b.t - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt - a.y*b.wyz - a.yt*b.wyzt) ; 

c.wt = + (-a.wxt*b.x - a.wxyt*b.xy + a.wxyzt*b.xyz - a.wxzt*b.xz - a.wyt*b.y - a.wyzt*b.yz - a.wzt*b.z - a.x*b.wxt - a.xy*b.wxyt + a.xyz*b.wxyzt - a.xz*b.wxzt - a.y*b.wyt - a.yz*b.wyzt - a.z*b.wzt) ; 

c.xy = + (-a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt + a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt - a.wz*b.wxyz + a.wzt*b.wxyzt - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt + a.z*b.xyz + a.zt*b.xyzt) ; 

c.xz = + (-a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt + a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.xyz*b.y - a.xyzt*b.yt - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt) ; 

c.xt = + (a.w*b.wxt + a.wxt*b.w + a.wxyt*b.wy - a.wxyzt*b.wyz + a.wxzt*b.wz + a.wy*b.wxyt - a.wyz*b.wxyzt + a.wz*b.wxzt - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt) ; 

c.yz = + (-a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt - a.wxyz*b.wx + a.wxyzt*b.wxt + a.wyz*b.w + a.wyzt*b.wt + a.x*b.xyz + a.xt*b.xyzt + a.xyz*b.x + a.xyzt*b.xt - a.yzt*b.t) ; 

c.yt = + (a.w*b.wyt - a.wx*b.wxyt - a.wxyt*b.wx + a.wxyzt*b.wxz + a.wxz*b.wxyzt + a.wyt*b.w + a.wyzt*b.wz + a.wz*b.wyzt + a.x*b.xyt + a.xyt*b.x + a.xyzt*b.xz + a.xz*b.xyzt - a.yzt*b.z - a.z*b.yzt) ; 

c.zt = + (a.w*b.wzt - a.wx*b.wxzt - a.wxy*b.wxyzt - a.wxyzt*b.wxy - a.wxzt*b.wx - a.wy*b.wyzt - a.wyzt*b.wy + a.wzt*b.w + a.x*b.xzt - a.xy*b.xyzt - a.xyzt*b.xy + a.xzt*b.x + a.y*b.yzt + a.yzt*b.y) ; 


	c.wxy =  + (a.t*b.wxyt - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt - a.z*b.wxyz + a.zt*b.wxyzt) ; 
	c.wxz = + (a.t*b.wxzt - a.wxyz*b.y - a.wxyzt*b.yt - a.wxzt*b.t + a.y*b.wxyz - a.yt*b.wxyzt) ; 
	c.wxt = + (-a.wxyt*b.y - a.wxyzt*b.yz - a.wxzt*b.z + a.y*b.wxyt - a.yz*b.wxyzt + a.z*b.wxzt) ; 
	c.wyz = + (a.t*b.wyzt + a.wxyz*b.x + a.wxyzt*b.xt - a.wyzt*b.t - a.x*b.wxyz + a.xt*b.wxyzt) ; 
	c.wyt = + (a.wxyt*b.x + a.wxyzt*b.xz - a.wyzt*b.z - a.x*b.wxyt + a.xz*b.wxyzt + a.z*b.wyzt) ; 
	c.wzt = + (-a.wxyzt*b.xy + a.wxzt*b.x + a.wyzt*b.y - a.x*b.wxzt - a.xy*b.wxyzt - a.y*b.wyzt) ; 
	c.xyz = + (a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt - a.wxyz*b.w - a.wxyzt*b.wt - a.xyzt*b.t) ; 
	c.xyt = + (a.w*b.wxyt - a.wxyt*b.w - a.wxyzt*b.wz - a.wz*b.wxyzt - a.xyzt*b.z + a.z*b.xyzt) ; 
	c.xzt = + (a.w*b.wxzt + a.wxyzt*b.wy - a.wxzt*b.w + a.wy*b.wxyzt + a.xyzt*b.y - a.y*b.xyzt) ; 
	c.yzt = + (a.w*b.wyzt - a.wx*b.wxyzt - a.wxyzt*b.wx - a.wyzt*b.w + a.x*b.xyzt - a.xyzt*b.x) ; 

	c.wxyz = + (-a.t*b.wxyzt - a.wxyzt*b.t) ; 
	c.wxyt = + (-a.wxyzt*b.z - a.z*b.wxyzt) ; 
	c.wxzt = + (a.wxyzt*b.y + a.y*b.wxyzt) ; 
	c.wyzt = + (-a.wxyzt*b.x - a.x*b.wxyzt) ; 
	c.xyzt = + (a.w*b.wxyzt + a.wxyzt*b.w) ; 
	
	c.wxyzt = 0 ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 LeftContraction(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.q*b.w + a.t*b.wt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.q*b.x + a.t*b.xt + a.w*b.wx - a.wt*b.wxt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.q*b.y + a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt + a.wxz*b.wxyz - a.wxzt*b.wxyzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt + a.xz*b.xyz + a.xzt*b.xyzt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.q*b.z + a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt - a.wy*b.wyz - a.wyt*b.wyzt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt + a.y*b.yz - a.yt*b.yzt) ; 

c.t = + (a.q*b.t + a.w*b.wt - a.wx*b.wxt - a.wxy*b.wxyt + a.wxyz*b.wxyzt - a.wxz*b.wxzt - a.wy*b.wyt - a.wyz*b.wyzt - a.wz*b.wzt + a.x*b.xt - a.xy*b.xyt - a.xyz*b.xyzt - a.xz*b.xzt + a.y*b.yt - a.yz*b.yzt + a.z*b.zt) ; 

	c.wx = + (a.q*b.wx - a.t*b.wxt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 
	c.wy = + (a.q*b.wy - a.t*b.wyt - a.x*b.wxy - a.xt*b.wxyt + a.xz*b.wxyz - a.xzt*b.wxyzt + a.z*b.wyz + a.zt*b.wyzt) ; 
	c.wz = + (a.q*b.wz - a.t*b.wzt - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt - a.y*b.wyz - a.yt*b.wyzt) ; 
	c.wt = + (a.q*b.wt - a.x*b.wxt - a.xy*b.wxyt + a.xyz*b.wxyzt - a.xz*b.wxzt - a.y*b.wyt - a.yz*b.wyzt - a.z*b.wzt) ; 
	c.xy = + (a.q*b.xy - a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt - a.wz*b.wxyz + a.wzt*b.wxyzt + a.z*b.xyz + a.zt*b.xyzt) ; 
	c.xz = + (a.q*b.xz - a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.y*b.xyz - a.yt*b.xyzt) ; 
	c.xt = + (a.q*b.xt + a.w*b.wxt + a.wy*b.wxyt - a.wyz*b.wxyzt + a.wz*b.wxzt - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt) ; 
	c.yz = + (a.q*b.yz - a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt + a.x*b.xyz + a.xt*b.xyzt) ; 
	c.yt = + (a.q*b.yt + a.w*b.wyt - a.wx*b.wxyt + a.wxz*b.wxyzt + a.wz*b.wyzt + a.x*b.xyt + a.xz*b.xyzt - a.z*b.yzt) ; 
	c.zt = + (a.q*b.zt + a.w*b.wzt - a.wx*b.wxzt - a.wxy*b.wxyzt - a.wy*b.wyzt + a.x*b.xzt - a.xy*b.xyzt + a.y*b.yzt) ; 

	c.wxy = + (a.q*b.wxy + a.t*b.wxyt - a.z*b.wxyz + a.zt*b.wxyzt) ; 
	c.wxz = + (a.q*b.wxz + a.t*b.wxzt + a.y*b.wxyz - a.yt*b.wxyzt) ; 
	c.wxt = + (a.q*b.wxt + a.y*b.wxyt - a.yz*b.wxyzt + a.z*b.wxzt) ; 
	c.wyz = + (a.q*b.wyz + a.t*b.wyzt - a.x*b.wxyz + a.xt*b.wxyzt) ; 
	c.wyt = + (a.q*b.wyt - a.x*b.wxyt + a.xz*b.wxyzt + a.z*b.wyzt) ; 
	c.wzt = + (a.q*b.wzt - a.x*b.wxzt - a.xy*b.wxyzt - a.y*b.wyzt) ; 
	c.xyz = + (a.q*b.xyz + a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt) ; 
	c.xyt = + (a.q*b.xyt + a.w*b.wxyt - a.wz*b.wxyzt + a.z*b.xyzt) ; 
	c.xzt = + (a.q*b.xzt + a.w*b.wxzt + a.wy*b.wxyzt - a.y*b.xyzt) ; 
	c.yzt = + (a.q*b.yzt + a.w*b.wyzt - a.wx*b.wxyzt + a.x*b.xyzt) ; 

	c.wxyz = + (a.q*b.wxyz - a.t*b.wxyzt) ; 
	c.wxyt = + (a.q*b.wxyt - a.z*b.wxyzt) ; 
	c.wxzt = + (a.q*b.wxzt + a.y*b.wxyzt) ; 
	c.wyzt = + (a.q*b.wyzt - a.x*b.wxyzt) ; 
	c.xyzt = + (a.q*b.xyzt + a.w*b.wxyzt) ; 

	c.wxyzt = + a.q*b.wxyzt ; 

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 RightContraction(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.w*b.q - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt) ; 

c.x = + (-a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt) ; 

c.y = + (-a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt) ; 

c.z = + (a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t) ; 

c.t = + (a.t*b.q - a.wt*b.w - a.wxt*b.wx + a.wxyt*b.wxy + a.wxyzt*b.wxyz + a.wxzt*b.wxz - a.wyt*b.wy + a.wyzt*b.wyz - a.wzt*b.wz - a.xt*b.x - a.xyt*b.xy + a.xyzt*b.xyz - a.xzt*b.xz - a.yt*b.y - a.yzt*b.yz - a.zt*b.z) ; 

	c.wx = + (a.wx*b.q - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt) ; 
	c.wy = + (-a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt + a.wy*b.q - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt) ; 
	c.wz = + (-a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt - a.wyz*b.y - a.wyzt*b.yt + a.wz*b.q - a.wzt*b.t) ; 
	c.wt = + (a.wt*b.q - a.wxt*b.x - a.wxyt*b.xy + a.wxyzt*b.xyz - a.wxzt*b.xz - a.wyt*b.y - a.wyzt*b.yz - a.wzt*b.z) ; 
	c.xy = + (a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt) ; 
	c.xz = + (a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t) ; 
	c.xt = + (a.wxt*b.w + a.wxyt*b.wy - a.wxyzt*b.wyz + a.wxzt*b.wz + a.xt*b.q - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z) ; 
	c.yz = + (-a.wxyz*b.wx + a.wxyzt*b.wxt + a.wyz*b.w + a.wyzt*b.wt + a.xyz*b.x + a.xyzt*b.xt + a.yz*b.q - a.yzt*b.t) ; 
	c.yt = + (-a.wxyt*b.wx + a.wxyzt*b.wxz + a.wyt*b.w + a.wyzt*b.wz + a.xyt*b.x + a.xyzt*b.xz + a.yt*b.q - a.yzt*b.z) ; 
	c.zt = + (-a.wxyzt*b.wxy - a.wxzt*b.wx - a.wyzt*b.wy + a.wzt*b.w - a.xyzt*b.xy + a.xzt*b.x + a.yzt*b.y + a.zt*b.q) ; 

	c.wxy = + (a.wxy*b.q - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt) ; 
	c.wxz = + (-a.wxyz*b.y - a.wxyzt*b.yt + a.wxz*b.q - a.wxzt*b.t) ; 
	c.wxt = + (a.wxt*b.q - a.wxyt*b.y - a.wxyzt*b.yz - a.wxzt*b.z) ; 
	c.wyz = + (a.wxyz*b.x + a.wxyzt*b.xt + a.wyz*b.q - a.wyzt*b.t) ; 
	c.wyt = + (a.wxyt*b.x + a.wxyzt*b.xz + a.wyt*b.q - a.wyzt*b.z) ; 
	c.wzt = + (-a.wxyzt*b.xy + a.wxzt*b.x + a.wyzt*b.y + a.wzt*b.q) ; 
	c.xyz = + (-a.wxyz*b.w - a.wxyzt*b.wt + a.xyz*b.q - a.xyzt*b.t) ; 
	c.xyt = + (-a.wxyt*b.w - a.wxyzt*b.wz + a.xyt*b.q - a.xyzt*b.z) ; 
	c.xzt = + (a.wxyzt*b.wy - a.wxzt*b.w + a.xyzt*b.y + a.xzt*b.q) ; 
	c.yzt = + (-a.wxyzt*b.wx - a.wyzt*b.w - a.xyzt*b.x + a.yzt*b.q) ; 

	c.wxyz = + (a.wxyz*b.q - a.wxyzt*b.t) ; 
	c.wxyt = + (a.wxyt*b.q - a.wxyzt*b.z) ; 
	c.wxzt = + (a.wxyzt*b.y + a.wxzt*b.q) ; 
	c.wyzt = + (-a.wxyzt*b.x + a.wyzt*b.q) ; 
	c.xyzt = + (a.wxyzt*b.w + a.xyzt*b.q) ; 

	c.wxyzt = + a.wxyzt*b.q ;

	c.q = expand(c.q);
	c.w = expand(c.w);	c.x = expand(c.x);	c.y = expand(c.y);	c.z = expand(c.z);	c.t = expand(c.t);
	c.wx = expand(c.wx);	c.wy = expand(c.wy);	c.wz = expand(c.wz);	c.wt = expand(c.wt);	c.xy = expand(c.xy);
	c.xz = expand(c.xz);	c.yz = expand(c.yz);	c.xt = expand(c.xt);	c.yt = expand(c.yt);	c.zt = expand(c.zt);
	c.wxy = expand(c.wxy);	c.wxz = expand(c.wxz);	c.wxt = expand(c.wxt);	c.wyz = expand(c.wyz);	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);	c.xyz = expand(c.xyz);	c.xyt = expand(c.xyt);	c.xzt = expand(c.xzt);	c.yzt = expand(c.yzt);
	c.wxyz = expand(c.wxyz); c.wxyt = expand(c.wxyt); c.wxzt = expand(c.wxzt); c.wyzt = expand(c.wyzt); c.xyzt = expand(c.xyzt);
	c.wxyzt = expand(c.wxyzt);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Parity(const GA5_4_1 &a)	// Corresponds to parity transform: w -> -w, x -> -x, y -> -y, z ->-z, t -> -t
// Conjugates complex determinant. Magnitude matches.
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = - a.xyz;	b.xyt = - a.xyt;	b.xzt = - a.xzt;	b.yzt = - a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = - a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 ComplexConjugate(const GA5_4_1 &a)	// negate every term including w
// Conjugates complex determinant. Magnitude matches.
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = - a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = a.t;
	b.wx = - a.wx;		b.wy = - a.wy;		b.wz = - a.wz;		b.wt = - a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz = - a.wxyz;	b.wxyt = - a.wxyt;	b.wxzt = - a.wxzt;	b.wyzt = - a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = - a.wxyzt;

	return b;
}
//////////////////////////////////////////////////////

GA5_4_1 Hermitian(const GA5_4_1 &a)	// InverseBasis = reverse and time conjugation. Each basis times its InverseBasis = 1
// Conjugates complex determinant. Magnitude matches.
//	+ a*q
//	+ b*w + c*x + d*y + e*z - f*t
//	- g*wx - h*wy - j*wz + k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	- S*wxy - R*wxz + P*wxt - N*wyz + M*wyt + L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	+ F*wxyz - E*wxyt - D*wxzt - C*wyzt - B*xyzt
//	- A*wxyzt 
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = a.wxt;		b.wyz = -a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Hermitian_Component(const GA5_4_1 &a)

//	+ a*q
//	+ b*w + c*x + d*y + e*z - f*t
//	- g*wx - h*wy - j*wz + k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	- S*wxy - R*wxz + P*wxt - N*wyz + M*wyt + L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	+ F*wxyz - E*wxyt - D*wxzt - C*wyzt - B*xyzt
//	- A*wxyzt 
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = 0;
	b.wx = 0;		b.wy = 0;		b.wz = 0;		b.wt = a.wt;		b.xy = 0;
	b.xz = 0;		b.yz = 0;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = 0;		b.wxz = 0;		b.wxt = a.wxt;		b.wyz = 0;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = 0;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt = 0;		b.wxzt = 0;		b.wyzt = 0;		b.xyzt =  0;
	b.wxyzt = 0;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Anti_Hermitian_Component(const GA5_4_1 &a)	// InverseBasis = reverse and time conjugation. Each basis times its InverseBasis = 1
// Conjugates complex determinant. Magnitude matches.
//	+ a*q
//	+ b*w + c*x + d*y + e*z - f*t
//	- g*wx - h*wy - j*wz + k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	- S*wxy - R*wxz + P*wxt - N*wyz + M*wyt + L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	+ F*wxyz - E*wxyt - D*wxzt - C*wyzt - B*xyzt
//	- A*wxyzt 
{
	GA5_4_1 b;

	b.q = 0;
	b.w = 0;		b.x = 0;		b.y = 0;		b.z = 0;		b.t = a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = 0;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = 0;		b.yt = 0;		b.zt = 0;
	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = 0;		b.wyz = a.wyz;		b.wyt = 0;
	b.wzt = 0;		b.xyz = a.xyz;		b.xyt = 0;		b.xzt = 0;		b.yzt = 0;
	b.wxyz = 0;		b.wxyt = a.wxyt;	b.wxzt = a.wxzt;	b.wyzt = a.wyzt;	b.xyzt = a.xyzt;
	b.wxyzt = a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Conjugate(const GA5_4_1 &a)	// This is not the ComplexConjugate - still useless. destroys determinant
// This changed all non-scalar terms, then flipped t terms
{
	GA5_4_1 b;

/*	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = -a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;
*/
	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = a.wxt;		b.wyz = -a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 SymmetricProductViaFormula(const GA5_4_1 &u, const GA5_4_1 &v)	// symmetric product
{
	GA5_4_1 a;
	a = (u*v + v*u)/2;
	return a;
}

//////////////////////////////////////////////////////

GA5_4_1 AntiSymmetricProductViaFormula(const GA5_4_1 &u, const GA5_4_1 &v)	// antisymmetric product
{
	GA5_4_1 a;
	a = (u*v - v*u)/2;
	return a;
}

//////////////////////////////////////////////////////

ex Determinant(const GA5_4_1 &MV1)
{
	ex det;
	ex a, b,c,d,e,f, g,G, F,E,D,C,B, A;
//	GA5_4_1 MV2, MV3, MV4, MV5;

//	MV2 = Reverse(MV1);
//	MV3 = MV1*MV2;	// result has zero in bivector and trivector components

//	MV4 = Reverse_Vector_Quad(MV3);
//	MV5 = MV3*MV4;
//	det = MV5.q + I*MV5.wxyzt;

//	deltat for 100X MV determinant = 10 => 100 msec/determinant using above
//			3.8 msec using code below..   25X improvement


	a = + MV1.q*MV1.q + MV1.w*MV1.w + MV1.x*MV1.x + MV1.y*MV1.y + MV1.z*MV1.z - MV1.t*MV1.t
		 + MV1.wx*MV1.wx + MV1.wy*MV1.wy + MV1.wz*MV1.wz - MV1.wt*MV1.wt + MV1.xy*MV1.xy
		 + MV1.xz*MV1.xz - MV1.xt*MV1.xt + MV1.yz*MV1.yz - MV1.yt*MV1.yt - MV1.zt*MV1.zt
		 + MV1.wxy*MV1.wxy + MV1.wxz*MV1.wxz - MV1.wxt*MV1.wxt + MV1.wyz*MV1.wyz - MV1.wyt*MV1.wyt 
		 - MV1.wzt*MV1.wzt + MV1.xyz*MV1.xyz - MV1.xyt*MV1.xyt - MV1.xzt*MV1.xzt - MV1.yzt*MV1.yzt
		 + MV1.wxyz*MV1.wxyz - MV1.wxyt*MV1.wxyt - MV1.wxzt*MV1.wxzt - MV1.wyzt*MV1.wyzt - MV1.xyzt*MV1.xyzt
		 - MV1.wxyzt*MV1.wxyzt ; 

	b = 2*MV1.y*MV1.wy-2*MV1.xzt*MV1.wxzt+2*MV1.wyz*MV1.yz+2*MV1.wxz*MV1.xz-2*MV1.yzt*MV1.wyzt
		-2*MV1.t*MV1.wt+2*MV1.wxy*MV1.xy-2*MV1.wxyt*MV1.xyt-2*MV1.yt*MV1.wyt-2*MV1.wxt*MV1.xt+2*MV1.z*MV1.wz
		-2*MV1.xyzt*MV1.wxyzt+2*MV1.xyz*MV1.wxyz+2*MV1.x*MV1.wx+2*MV1.w*MV1.q-2*MV1.wzt*MV1.zt ; 

	c = -2*MV1.wyz*MV1.wxyz-2*MV1.xyzt*MV1.yzt-2*MV1.xzt*MV1.zt-2*MV1.w*MV1.wx+2*MV1.wxyzt*MV1.wyzt
		-2*MV1.wxy*MV1.wy+2*MV1.x*MV1.q+2*MV1.wxyt*MV1.wyt+2*MV1.y*MV1.xy+2*MV1.z*MV1.xz+2*MV1.wxt*MV1.wt
		-2*MV1.yt*MV1.xyt-2*MV1.t*MV1.xt-2*MV1.wxz*MV1.wz+2*MV1.xyz*MV1.yz+2*MV1.wzt*MV1.wxzt ; 

	d = 2*MV1.xt*MV1.xyt-2*MV1.wxyt*MV1.wxt-2*MV1.x*MV1.xy+2*MV1.wzt*MV1.wyzt+2*MV1.xyzt*MV1.xzt
		-2*MV1.wz*MV1.wyz+2*MV1.y*MV1.q-2*MV1.w*MV1.wy+2*MV1.wt*MV1.wyt-2*MV1.wxyzt*MV1.wxzt-2*MV1.zt*MV1.yzt
		-2*MV1.xz*MV1.xyz+2*MV1.z*MV1.yz+2*MV1.wxz*MV1.wxyz+2*MV1.wx*MV1.wxy-2*MV1.yt*MV1.t ; 

	e = 2*MV1.wxyt*MV1.wxyzt+2*MV1.xzt*MV1.xt+2*MV1.wzt*MV1.wt-2*MV1.wyzt*MV1.wyt
		+2*MV1.wyz*MV1.wy-2*MV1.xz*MV1.x+2*MV1.yt*MV1.yzt+2*MV1.wxz*MV1.wx-2*MV1.y*MV1.yz-2*MV1.wxy*MV1.wxyz
		-2*MV1.wxt*MV1.wxzt-2*MV1.t*MV1.zt+2*MV1.z*MV1.q-2*MV1.xyzt*MV1.xyt+2*MV1.xyz*MV1.xy-2*MV1.w*MV1.wz ; 

	f = 2*MV1.wy*MV1.wyt+2*MV1.wxt*MV1.wx-2*MV1.xyzt*MV1.xyz-2*MV1.yt*MV1.y-2*MV1.w*MV1.wt
		-2*MV1.x*MV1.xt+2*MV1.xy*MV1.xyt+2*MV1.yzt*MV1.yz+2*MV1.wz*MV1.wzt-2*MV1.wxz*MV1.wxzt-2*MV1.wxyt*MV1.wxy
		-2*MV1.wyz*MV1.wyzt+2*MV1.xz*MV1.xzt+2*MV1.wxyzt*MV1.wxyz-2*MV1.z*MV1.zt+2*MV1.t*MV1.q ; 

	F = -2*MV1.xzt*MV1.wyt+2*MV1.xyzt*MV1.wt-2*MV1.wx*MV1.yz-2*MV1.w*MV1.xyz+2*MV1.wxt*MV1.yzt
		+2*MV1.z*MV1.wxy+2*MV1.xz*MV1.wy-2*MV1.wxyt*MV1.zt-2*MV1.wxz*MV1.y-2*MV1.wz*MV1.xy+2*MV1.yt*MV1.wxzt
		-2*MV1.xt*MV1.wyzt+2*MV1.q*MV1.wxyz-2*MV1.t*MV1.wxyzt+2*MV1.wzt*MV1.xyt+2*MV1.x*MV1.wyz ; 

	E = 2*MV1.wzt*MV1.xyz-2*MV1.yt*MV1.wx+2*MV1.wxzt*MV1.yz+2*MV1.t*MV1.wxy+2*MV1.wxz*MV1.yzt
		+2*MV1.xt*MV1.wy+2*MV1.wxyt*MV1.q-2*MV1.zt*MV1.wxyz-2*MV1.wxt*MV1.y+2*MV1.x*MV1.wyt+2*MV1.xyzt*MV1.wz
		-2*MV1.z*MV1.wxyzt-2*MV1.w*MV1.xyt-2*MV1.xzt*MV1.wyz-2*MV1.wt*MV1.xy-2*MV1.xz*MV1.wyzt ; 

	D = -2*MV1.w*MV1.xzt+2*MV1.wxzt*MV1.q-2*MV1.wxyt*MV1.yz+2*MV1.wxyzt*MV1.y-2*MV1.yzt*MV1.wxy
		+2*MV1.wyz*MV1.xyt+2*MV1.wyzt*MV1.xy-2*MV1.xyz*MV1.wyt-2*MV1.z*MV1.wxt-2*MV1.xz*MV1.wt+2*MV1.wxz*MV1.t
		+2*MV1.wz*MV1.xt-2*MV1.zt*MV1.wx+2*MV1.wzt*MV1.x-2*MV1.xyzt*MV1.wy+2*MV1.yt*MV1.wxyz ; 

	C = 2*MV1.t*MV1.wyz-2*MV1.wxzt*MV1.xy+2*MV1.wzt*MV1.y+2*MV1.xyzt*MV1.wx+2*MV1.wyzt*MV1.q
		-2*MV1.w*MV1.yzt-2*MV1.xt*MV1.wxyz+2*MV1.wxt*MV1.xyz-2*MV1.zt*MV1.wy-2*MV1.wxyzt*MV1.x+2*MV1.xzt*MV1.wxy
		-2*MV1.z*MV1.wyt+2*MV1.yt*MV1.wz-2*MV1.wxz*MV1.xyt+2*MV1.wxyt*MV1.xz-2*MV1.wt*MV1.yz ; 

	B = 2*MV1.yt*MV1.xz-2*MV1.wxt*MV1.wyz-2*MV1.zt*MV1.xy-2*MV1.wxyt*MV1.wz+2*MV1.t*MV1.xyz
		-2*MV1.xt*MV1.yz-2*MV1.wzt*MV1.wxy+2*MV1.wxzt*MV1.wy+2*MV1.xzt*MV1.y+2*MV1.xyzt*MV1.q-2*MV1.x*MV1.yzt
		-2*MV1.wx*MV1.wyzt+2*MV1.w*MV1.wxyzt+2*MV1.wt*MV1.wxyz-2*MV1.z*MV1.xyt+2*MV1.wxz*MV1.wyt ; 

	A = 2*MV1.xzt*MV1.wy+2*MV1.wxz*MV1.yt+2*MV1.xyzt*MV1.w+2*MV1.wt*MV1.xyz+2*MV1.xz*MV1.wyt
		-2*MV1.z*MV1.wxyt-2*MV1.yzt*MV1.wx-2*MV1.xt*MV1.wyz-2*MV1.wzt*MV1.xy-2*MV1.wxt*MV1.yz+2*MV1.t*MV1.wxyz
		-2*MV1.x*MV1.wyzt+2*MV1.wxyzt*MV1.q-2*MV1.zt*MV1.wxy+2*MV1.wxzt*MV1.y-2*MV1.wz*MV1.xyt ; 


//	g = C^2-c^2-F^2-A^2+f^2+a^2+D^2-d^2+B^2-b^2+E^2-e^2 ; 

	g = a*a - b*b -c*c - d*d - e*e + f*f - F*F + E*E + D*D + C*C + B*B - A*A;

	G = + 2*a*A - 2*b*B + 2*c*C - 2*d*D + 2*e*E - 2*f*F ; 

	det = expand(g + G*I);
	
	return det;
}

/////////////////////////////////////////////////////////////////////

GA5_4_1 Adjugate(const GA5_4_1 &r)
{

// Adjugate = Reverse(r)*Reverse_Vector_Quad(r*Reverse(r));

	ex a, b,c,d,e,f, F,E,D,C,B, A;
	GA5_4_1 MV1, MV2, Adjugate;
	
	MV1 = Reverse(r);
//	MV2 = Reverse_Vector_Quad(r*Reverse(r););
//	Adjugate = MV1*MV2;
	
	a = + r.q*r.q + r.w*r.w + r.x*r.x + r.y*r.y + r.z*r.z - r.t*r.t
		 + r.wx*r.wx + r.wy*r.wy + r.wz*r.wz - r.wt*r.wt + r.xy*r.xy
		 + r.xz*r.xz - r.xt*r.xt + r.yz*r.yz - r.yt*r.yt - r.zt*r.zt
		 + r.wxy*r.wxy + r.wxz*r.wxz - r.wxt*r.wxt + r.wyz*r.wyz - r.wyt*r.wyt 
		 - r.wzt*r.wzt + r.xyz*r.xyz - r.xyt*r.xyt - r.xzt*r.xzt - r.yzt*r.yzt
		 + r.wxyz*r.wxyz - r.wxyt*r.wxyt - r.wxzt*r.wxzt - r.wyzt*r.wyzt - r.xyzt*r.xyzt
		 - r.wxyzt*r.wxyzt ; 

	b = 2*r.y*r.wy-2*r.xzt*r.wxzt+2*r.wyz*r.yz+2*r.wxz*r.xz-2*r.yzt*r.wyzt
		-2*r.t*r.wt+2*r.wxy*r.xy-2*r.wxyt*r.xyt-2*r.yt*r.wyt-2*r.wxt*r.xt+2*r.z*r.wz
		-2*r.xyzt*r.wxyzt+2*r.xyz*r.wxyz+2*r.x*r.wx+2*r.w*r.q-2*r.wzt*r.zt ; 

	c = -2*r.wyz*r.wxyz-2*r.xyzt*r.yzt-2*r.xzt*r.zt-2*r.w*r.wx+2*r.wxyzt*r.wyzt
		-2*r.wxy*r.wy+2*r.x*r.q+2*r.wxyt*r.wyt+2*r.y*r.xy+2*r.z*r.xz+2*r.wxt*r.wt
		-2*r.yt*r.xyt-2*r.t*r.xt-2*r.wxz*r.wz+2*r.xyz*r.yz+2*r.wzt*r.wxzt ; 

	d = 2*r.xt*r.xyt-2*r.wxyt*r.wxt-2*r.x*r.xy+2*r.wzt*r.wyzt+2*r.xyzt*r.xzt
		-2*r.wz*r.wyz+2*r.y*r.q-2*r.w*r.wy+2*r.wt*r.wyt-2*r.wxyzt*r.wxzt-2*r.zt*r.yzt
		-2*r.xz*r.xyz+2*r.z*r.yz+2*r.wxz*r.wxyz+2*r.wx*r.wxy-2*r.yt*r.t ; 

	e = 2*r.wxyt*r.wxyzt+2*r.xzt*r.xt+2*r.wzt*r.wt-2*r.wyzt*r.wyt
		+2*r.wyz*r.wy-2*r.xz*r.x+2*r.yt*r.yzt+2*r.wxz*r.wx-2*r.y*r.yz-2*r.wxy*r.wxyz
		-2*r.wxt*r.wxzt-2*r.t*r.zt+2*r.z*r.q-2*r.xyzt*r.xyt+2*r.xyz*r.xy-2*r.w*r.wz ; 

	f = 2*r.wy*r.wyt+2*r.wxt*r.wx-2*r.xyzt*r.xyz-2*r.yt*r.y-2*r.w*r.wt
		-2*r.x*r.xt+2*r.xy*r.xyt+2*r.yzt*r.yz+2*r.wz*r.wzt-2*r.wxz*r.wxzt-2*r.wxyt*r.wxy
		-2*r.wyz*r.wyzt+2*r.xz*r.xzt+2*r.wxyzt*r.wxyz-2*r.z*r.zt+2*r.t*r.q ; 

	F = -2*r.xzt*r.wyt+2*r.xyzt*r.wt-2*r.wx*r.yz-2*r.w*r.xyz+2*r.wxt*r.yzt
		+2*r.z*r.wxy+2*r.xz*r.wy-2*r.wxyt*r.zt-2*r.wxz*r.y-2*r.wz*r.xy+2*r.yt*r.wxzt
		-2*r.xt*r.wyzt+2*r.q*r.wxyz-2*r.t*r.wxyzt+2*r.wzt*r.xyt+2*r.x*r.wyz ; 

	E = 2*r.wzt*r.xyz-2*r.yt*r.wx+2*r.wxzt*r.yz+2*r.t*r.wxy+2*r.wxz*r.yzt
		+2*r.xt*r.wy+2*r.wxyt*r.q-2*r.zt*r.wxyz-2*r.wxt*r.y+2*r.x*r.wyt+2*r.xyzt*r.wz
		-2*r.z*r.wxyzt-2*r.w*r.xyt-2*r.xzt*r.wyz-2*r.wt*r.xy-2*r.xz*r.wyzt ; 

	D = -2*r.w*r.xzt+2*r.wxzt*r.q-2*r.wxyt*r.yz+2*r.wxyzt*r.y-2*r.yzt*r.wxy
		+2*r.wyz*r.xyt+2*r.wyzt*r.xy-2*r.xyz*r.wyt-2*r.z*r.wxt-2*r.xz*r.wt+2*r.wxz*r.t
		+2*r.wz*r.xt-2*r.zt*r.wx+2*r.wzt*r.x-2*r.xyzt*r.wy+2*r.yt*r.wxyz ; 

	C = 2*r.t*r.wyz-2*r.wxzt*r.xy+2*r.wzt*r.y+2*r.xyzt*r.wx+2*r.wyzt*r.q
		-2*r.w*r.yzt-2*r.xt*r.wxyz+2*r.wxt*r.xyz-2*r.zt*r.wy-2*r.wxyzt*r.x+2*r.xzt*r.wxy
		-2*r.z*r.wyt+2*r.yt*r.wz-2*r.wxz*r.xyt+2*r.wxyt*r.xz-2*r.wt*r.yz ; 

	B = 2*r.yt*r.xz-2*r.wxt*r.wyz-2*r.zt*r.xy-2*r.wxyt*r.wz+2*r.t*r.xyz
		-2*r.xt*r.yz-2*r.wzt*r.wxy+2*r.wxzt*r.wy+2*r.xzt*r.y+2*r.xyzt*r.q-2*r.x*r.yzt
		-2*r.wx*r.wyzt+2*r.w*r.wxyzt+2*r.wt*r.wxyz-2*r.z*r.xyt+2*r.wxz*r.wyt ; 

	A = 2*r.xzt*r.wy+2*r.wxz*r.yt+2*r.xyzt*r.w+2*r.wt*r.xyz+2*r.xz*r.wyt
		-2*r.z*r.wxyt-2*r.yzt*r.wx-2*r.xt*r.wyz-2*r.wzt*r.xy-2*r.wxt*r.yz+2*r.t*r.wxyz
		-2*r.x*r.wyzt+2*r.wxyzt*r.q-2*r.zt*r.wxy+2*r.wxzt*r.y-2*r.wz*r.xyt ; 

	// MV2 = Reverse_Vector_Quad(r*Reverse(r));
	MV2 = GA5_4_1(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);
	
	Adjugate = MV1*MV2;  // Adjugate = Reverse(r)*Reverse_Vector_Quad(r*Reverse(r));	
	Adjugate.q = expand(Adjugate.q);

	Adjugate.w = expand(Adjugate.w);
	Adjugate.x = expand(Adjugate.x);
	Adjugate.y = expand(Adjugate.y);
	Adjugate.z = expand(Adjugate.z);
	Adjugate.t = expand(Adjugate.t);

	Adjugate.wx = expand(Adjugate.wx);
	Adjugate.wy = expand(Adjugate.wy);
	Adjugate.wz = expand(Adjugate.wz);
	Adjugate.wt = expand(Adjugate.wt);
	Adjugate.xy = expand(Adjugate.xy);
	Adjugate.xz = expand(Adjugate.xz);
	Adjugate.yz = expand(Adjugate.yz);
	Adjugate.xt = expand(Adjugate.xt);
	Adjugate.yt = expand(Adjugate.yt);
	Adjugate.zt = expand(Adjugate.zt);

	Adjugate.wxy = expand(Adjugate.wxy);
	Adjugate.wxz = expand(Adjugate.wxz);
	Adjugate.wxt = expand(Adjugate.wxt);
	Adjugate.wyz = expand(Adjugate.wyz);
	Adjugate.wyt = expand(Adjugate.wyt);
	Adjugate.wzt = expand(Adjugate.wzt);
	Adjugate.xyz = expand(Adjugate.xyz);
	Adjugate.xyt = expand(Adjugate.xyt);
	Adjugate.xzt = expand(Adjugate.xzt);
	Adjugate.yzt = expand(Adjugate.yzt);

	Adjugate.wxyz = expand(Adjugate.wxyz);
	Adjugate.wxyt = expand(Adjugate.wxyt);
	Adjugate.wxzt = expand(Adjugate.wxzt);
	Adjugate.wyzt = expand(Adjugate.wyzt);
	Adjugate.xyzt = expand(Adjugate.xyzt);

	Adjugate.wxyzt = expand(Adjugate.wxyzt);

	return Adjugate;
}

/////////////////////////////////////////////////////////////////////

GA5_4_1 Reciprocal(const GA5_4_1 &r)
{
	ex det;
	ex a, b,c,d,e,f, g,G, F,E,D,C,B, A;
	GA5_4_1 MV1, MV2, MV3, Result;
	
	MV1 = Reverse(r);
//	MV2 = r*Reverse(r);
//	MV3 = Reverse_Vector_Quad(MV2);
	
	a = + r.q*r.q + r.w*r.w + r.x*r.x + r.y*r.y + r.z*r.z - r.t*r.t
		 + r.wx*r.wx + r.wy*r.wy + r.wz*r.wz - r.wt*r.wt + r.xy*r.xy
		 + r.xz*r.xz - r.xt*r.xt + r.yz*r.yz - r.yt*r.yt - r.zt*r.zt
		 + r.wxy*r.wxy + r.wxz*r.wxz - r.wxt*r.wxt + r.wyz*r.wyz - r.wyt*r.wyt 
		 - r.wzt*r.wzt + r.xyz*r.xyz - r.xyt*r.xyt - r.xzt*r.xzt - r.yzt*r.yzt
		 + r.wxyz*r.wxyz - r.wxyt*r.wxyt - r.wxzt*r.wxzt - r.wyzt*r.wyzt - r.xyzt*r.xyzt
		 - r.wxyzt*r.wxyzt ; 

	b = 2*r.y*r.wy-2*r.xzt*r.wxzt+2*r.wyz*r.yz+2*r.wxz*r.xz-2*r.yzt*r.wyzt
		-2*r.t*r.wt+2*r.wxy*r.xy-2*r.wxyt*r.xyt-2*r.yt*r.wyt-2*r.wxt*r.xt+2*r.z*r.wz
		-2*r.xyzt*r.wxyzt+2*r.xyz*r.wxyz+2*r.x*r.wx+2*r.w*r.q-2*r.wzt*r.zt ; 

	c = -2*r.wyz*r.wxyz-2*r.xyzt*r.yzt-2*r.xzt*r.zt-2*r.w*r.wx+2*r.wxyzt*r.wyzt
		-2*r.wxy*r.wy+2*r.x*r.q+2*r.wxyt*r.wyt+2*r.y*r.xy+2*r.z*r.xz+2*r.wxt*r.wt
		-2*r.yt*r.xyt-2*r.t*r.xt-2*r.wxz*r.wz+2*r.xyz*r.yz+2*r.wzt*r.wxzt ; 

	d = 2*r.xt*r.xyt-2*r.wxyt*r.wxt-2*r.x*r.xy+2*r.wzt*r.wyzt+2*r.xyzt*r.xzt
		-2*r.wz*r.wyz+2*r.y*r.q-2*r.w*r.wy+2*r.wt*r.wyt-2*r.wxyzt*r.wxzt-2*r.zt*r.yzt
		-2*r.xz*r.xyz+2*r.z*r.yz+2*r.wxz*r.wxyz+2*r.wx*r.wxy-2*r.yt*r.t ; 

	e = 2*r.wxyt*r.wxyzt+2*r.xzt*r.xt+2*r.wzt*r.wt-2*r.wyzt*r.wyt
		+2*r.wyz*r.wy-2*r.xz*r.x+2*r.yt*r.yzt+2*r.wxz*r.wx-2*r.y*r.yz-2*r.wxy*r.wxyz
		-2*r.wxt*r.wxzt-2*r.t*r.zt+2*r.z*r.q-2*r.xyzt*r.xyt+2*r.xyz*r.xy-2*r.w*r.wz ; 

	f = 2*r.wy*r.wyt+2*r.wxt*r.wx-2*r.xyzt*r.xyz-2*r.yt*r.y-2*r.w*r.wt
		-2*r.x*r.xt+2*r.xy*r.xyt+2*r.yzt*r.yz+2*r.wz*r.wzt-2*r.wxz*r.wxzt-2*r.wxyt*r.wxy
		-2*r.wyz*r.wyzt+2*r.xz*r.xzt+2*r.wxyzt*r.wxyz-2*r.z*r.zt+2*r.t*r.q ; 

	F = -2*r.xzt*r.wyt+2*r.xyzt*r.wt-2*r.wx*r.yz-2*r.w*r.xyz+2*r.wxt*r.yzt
		+2*r.z*r.wxy+2*r.xz*r.wy-2*r.wxyt*r.zt-2*r.wxz*r.y-2*r.wz*r.xy+2*r.yt*r.wxzt
		-2*r.xt*r.wyzt+2*r.q*r.wxyz-2*r.t*r.wxyzt+2*r.wzt*r.xyt+2*r.x*r.wyz ; 

	E = 2*r.wzt*r.xyz-2*r.yt*r.wx+2*r.wxzt*r.yz+2*r.t*r.wxy+2*r.wxz*r.yzt
		+2*r.xt*r.wy+2*r.wxyt*r.q-2*r.zt*r.wxyz-2*r.wxt*r.y+2*r.x*r.wyt+2*r.xyzt*r.wz
		-2*r.z*r.wxyzt-2*r.w*r.xyt-2*r.xzt*r.wyz-2*r.wt*r.xy-2*r.xz*r.wyzt ; 

	D = -2*r.w*r.xzt+2*r.wxzt*r.q-2*r.wxyt*r.yz+2*r.wxyzt*r.y-2*r.yzt*r.wxy
		+2*r.wyz*r.xyt+2*r.wyzt*r.xy-2*r.xyz*r.wyt-2*r.z*r.wxt-2*r.xz*r.wt+2*r.wxz*r.t
		+2*r.wz*r.xt-2*r.zt*r.wx+2*r.wzt*r.x-2*r.xyzt*r.wy+2*r.yt*r.wxyz ; 

	C = 2*r.t*r.wyz-2*r.wxzt*r.xy+2*r.wzt*r.y+2*r.xyzt*r.wx+2*r.wyzt*r.q
		-2*r.w*r.yzt-2*r.xt*r.wxyz+2*r.wxt*r.xyz-2*r.zt*r.wy-2*r.wxyzt*r.x+2*r.xzt*r.wxy
		-2*r.z*r.wyt+2*r.yt*r.wz-2*r.wxz*r.xyt+2*r.wxyt*r.xz-2*r.wt*r.yz ; 

	B = 2*r.yt*r.xz-2*r.wxt*r.wyz-2*r.zt*r.xy-2*r.wxyt*r.wz+2*r.t*r.xyz
		-2*r.xt*r.yz-2*r.wzt*r.wxy+2*r.wxzt*r.wy+2*r.xzt*r.y+2*r.xyzt*r.q-2*r.x*r.yzt
		-2*r.wx*r.wyzt+2*r.w*r.wxyzt+2*r.wt*r.wxyz-2*r.z*r.xyt+2*r.wxz*r.wyt ; 

	A = 2*r.xzt*r.wy+2*r.wxz*r.yt+2*r.xyzt*r.w+2*r.wt*r.xyz+2*r.xz*r.wyt
		-2*r.z*r.wxyt-2*r.yzt*r.wx-2*r.xt*r.wyz-2*r.wzt*r.xy-2*r.wxt*r.yz+2*r.t*r.wxyz
		-2*r.x*r.wyzt+2*r.wxyzt*r.q-2*r.zt*r.wxy+2*r.wxzt*r.y-2*r.wz*r.xyt ; 

	// MV2 = Reverse_Vector_Quad(r*Reverse(r));
	MV2 = GA5_4_1(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);

	g = a*a - b*b -c*c - d*d - e*e + f*f - F*F + E*E + D*D + C*C + B*B - A*A;

	G = + 2*a*A - 2*b*B + 2*c*C - 2*d*D + 2*e*E - 2*f*F ; 

	MV3 = GA5_4_1(g, 0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0, -G);

	det = g*g + G*G;
	
	Result = MV1*MV2*MV3;	// need to divide by det. (Result is the Adjugate.)
	Result.q = expand(Result.q)/det;

	Result.w = expand(Result.w)/det;
	Result.x = expand(Result.x)/det;
	Result.y = expand(Result.y)/det;
	Result.z = expand(Result.z)/det;
	Result.t = expand(Result.t)/det;

	Result.wx = expand(Result.wx)/det;
	Result.wy = expand(Result.wy)/det;
	Result.wz = expand(Result.wz)/det;
	Result.wt = expand(Result.wt)/det;
	Result.xy = expand(Result.xy)/det;
	Result.xz = expand(Result.xz)/det;
	Result.yz = expand(Result.yz)/det;
	Result.xt = expand(Result.xt)/det;
	Result.yt = expand(Result.yt)/det;
	Result.zt = expand(Result.zt)/det;

	Result.wxy = expand(Result.wxy)/det;
	Result.wxz = expand(Result.wxz)/det;
	Result.wxt = expand(Result.wxt)/det;
	Result.wyz = expand(Result.wyz)/det;
	Result.wyt = expand(Result.wyt)/det;
	Result.wzt = expand(Result.wzt)/det;
	Result.xyz = expand(Result.xyz)/det;
	Result.xyt = expand(Result.xyt)/det;
	Result.xzt = expand(Result.xzt)/det;
	Result.yzt = expand(Result.yzt)/det;

	Result.wxyz = expand(Result.wxyz)/det;
	Result.wxyt = expand(Result.wxyt)/det;
	Result.wxzt = expand(Result.wxzt)/det;
	Result.wyzt = expand(Result.wyzt)/det;
	Result.xyzt = expand(Result.xyzt)/det;

	Result.wxyzt = expand(Result.wxyzt)/det;

	return Result;
}

///////////////////////////////////

GA5_4_1 Magic(GA5_4_1 V, int i) {

	GA5_4_1 MV;

	ex a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A;

	a = V.q;	b = V.w;	c = V.x;	d = V.y;	e = V.z;	f = V.t;
	g = V.wx;	h = V.wy;	j = V.wz;	k = V.wt;	l = V.xy;
	m = V.xz;	n = V.xt;	p = V.yz;	r = V.yt;	s = V.zt;
	S = V.wxy;	R = V.wxz;	P = V.wxt;	N = V.wyz;	M = V.wyt;
	L = V.wzt;	K = V.xyz;	J = V.xyt;	H = V.xzt;	G = V.yzt;
	F = V.wxyz;	E = V.wxyt;	D = V.wxzt;	C = V.wyzt;	B = V.xyzt;	A = V.wxyzt;

	MV = V;	// default
	if (i<  0) printf("Magic - error: index below zero\n");
	if (i>359) printf("Magic - error: index above 359\n");

	if(i==  0) MV = GA5_4_1(+a,+b,+c,+d,+e,+f,+g,+h,+j,+k,+l,+m,+n,+p,+r,+s,+S,+R,+P,+N,+M,+L,+K,+J,+H,+G,+F,+E,+D,+C,+B,+A) ;
	if(i==  1) MV = GA5_4_1(+a,+b,+c,+d,-F,-E,+g,+h,-K,-J,+l,+N,+M,-R,-P,+s,+S,+p,+r,-m,-n,+L,+j,+k,+H,+G,+e,+f,+D,+C,+B,+A) ;
	if(i==  2) MV = GA5_4_1(+a,+b,+c,+F,+d,-E,+g,+K,+h,-J,-N,+l,+M,-R,-s,-P,-p,+S,+r,-m,-L,-n,+j,-H,+k,+G,+e,-D,+f,+C,+B,+A) ;
	if(i==  3) MV = GA5_4_1(+a,+b,+e,+c,+d,+f,+j,+g,+h,+k,-m,-p,+s,+l,+n,+r,-R,-N,+L,+S,+P,+M,+K,-H,-G,+J,+F,-D,-C,+E,+B,+A) ;
	if(i==  4) MV = GA5_4_1(+a,+b,-F,+c,+d,-E,-K,+g,+h,-J,-N,+R,+s,+l,+M,-P,-p,+m,+L,+S,+r,-n,+j,-H,-G,+k,+e,-D,-C,+f,+B,+A) ;
	if(i==  5) MV = GA5_4_1(+a,+b,+e,+c,-F,+D,+j,+g,-K,+H,-m,+S,+P,+N,-L,+r,-R,+l,+n,+p,-s,+M,+h,+k,-G,+J,+d,+f,-C,+E,+B,+A) ;
	if(i==  6) MV = GA5_4_1(+a,+b,+d,+e,+c,+f,+h,+j,+g,+k,+p,-l,+r,-m,+s,+n,+N,-S,+M,-R,+L,+P,+K,+G,-J,-H,+F,+C,-E,-D,+B,+A) ;
	if(i==  7) MV = GA5_4_1(+a,+b,+e,+F,+c,+D,+j,+K,+g,+H,-S,-m,+P,+N,-r,-L,-l,-R,+n,+p,-M,-s,+h,+G,+k,+J,+d,+C,+f,+E,+B,+A) ;
	if(i==  8) MV = GA5_4_1(+a,+b,-F,+e,+c,+D,-K,+j,+g,+H,-S,-N,+r,-m,+P,-L,-l,-p,+M,-R,+n,-s,+h,+G,-J,+k,+d,+C,-E,+f,+B,+A) ;
	if(i==  9) MV = GA5_4_1(+a,+b,+k,+J,+H,+g,+f,+E,+D,+c,-S,-R,+n,+p,-M,-L,-l,-m,+P,+N,-r,-s,+C,+d,+e,+F,+G,+h,+j,+K,+B,+A) ;
	if(i== 10) MV = GA5_4_1(+a,+b,-J,+k,+H,+g,-E,+f,+D,+c,-S,-p,+M,-R,+n,-L,-l,-N,+r,-m,+P,-s,+C,+d,-F,+e,+G,+h,-K,+j,+B,+A) ;
	if(i== 11) MV = GA5_4_1(+a,+b,-J,-H,+k,+g,-E,-D,+f,+c,+p,-S,+M,-R,+L,+n,+N,-l,+r,-m,+s,+P,+C,+F,+d,+e,+G,+K,+h,+j,+B,+A) ;
	if(i== 12) MV = GA5_4_1(+a,+b,+k,-H,-G,+j,+f,-D,-C,+e,+R,+N,+s,+l,-P,-M,+m,+p,+L,+S,-n,-r,+E,+c,+d,+F,+J,+g,+h,+K,+B,+A) ;
	if(i== 13) MV = GA5_4_1(+a,+b,-J,-H,-G,-K,-E,-D,-C,-F,+p,-m,+s,+l,-r,+n,+N,-R,+L,+S,-M,+P,+f,+c,+d,+e,+k,+g,+h,+j,+B,+A) ;
	if(i== 14) MV = GA5_4_1(+a,+b,+H,+k,-G,+j,+D,+f,-C,+e,+R,-l,+P,+N,+s,-M,+m,-S,+n,+p,+L,-r,+E,+c,-F,+d,+J,+g,-K,+h,+B,+A) ;
	if(i== 15) MV = GA5_4_1(+a,+b,+k,+G,-J,+h,+f,+C,-E,+d,-N,+S,+r,-m,-L,-P,-p,+l,+M,-R,-s,-n,-D,+e,+c,+F,-H,+j,+g,+K,+B,+A) ;
	if(i== 16) MV = GA5_4_1(+a,+b,+H,+G,+k,+j,+D,+C,+f,+e,+l,+R,+P,+N,+M,+s,+S,+m,+n,+p,+r,+L,+E,+F,+c,+d,+J,+K,+g,+h,+B,+A) ;
	if(i== 17) MV = GA5_4_1(+a,+b,+H,+G,-J,-K,+D,+C,-E,-F,+l,+p,+r,-m,-n,+s,+S,+N,+M,-R,-P,+L,+f,+e,+c,+d,+k,+j,+g,+h,+B,+A) ;
	if(i== 18) MV = GA5_4_1(+a,+b,-G,+k,-J,+h,-C,+f,-E,+d,-N,+m,+L,+S,+r,-P,-p,+R,+s,+l,+M,-n,-D,+e,-F,+c,-H,+j,-K,+g,+B,+A) ;
	if(i== 19) MV = GA5_4_1(+a,+b,-G,+J,+k,+h,-C,+E,+f,+d,-m,-N,+L,+S,+P,+r,-R,-p,+s,+l,+n,+M,-D,+F,+e,+c,-H,+K,+j,+g,+B,+A) ;
	if(i== 20) MV = GA5_4_1(+a,+b,-G,+J,+H,-K,-C,+E,+D,-F,-m,+l,+n,+p,-s,+r,-R,+S,+P,+N,-L,+M,+f,+d,+e,+c,+k,+h,+j,+g,+B,+A) ;
	if(i== 21) MV = GA5_4_1(+a,+b,+d,+e,-F,-C,+h,+j,-K,-G,+p,-R,+L,+S,-M,+n,+N,-m,+s,+l,-r,+P,+g,+k,-J,-H,+c,+f,-E,-D,+B,+A) ;
	if(i== 22) MV = GA5_4_1(+a,+b,+d,+F,+e,-C,+h,+K,+j,-G,+R,+p,+L,+S,-n,-M,+m,+N,+s,+l,-P,-r,+g,+J,+k,-H,+c,+E,+f,-D,+B,+A) ;
	if(i== 23) MV = GA5_4_1(+a,+b,-F,+d,+e,-C,-K,+h,+j,-G,+R,-S,+n,+p,+L,-M,+m,-l,+P,+N,+s,-r,+g,+J,+H,+k,+c,+E,+D,+f,+B,+A) ;
	if(i== 24) MV = GA5_4_1(+a,+c,+b,+e,+d,+f,-g,+m,+l,+n,+j,+h,+k,-p,+s,+r,-R,-S,-P,-K,+H,+J,-N,+L,+M,-G,+F,-D,-E,-B,-C,+A) ;
	if(i== 25) MV = GA5_4_1(+a,+c,+b,+e,-F,+D,-g,+m,+N,-L,+j,-K,+H,+S,+P,+r,-R,-p,+s,-h,-k,+J,+l,+n,+M,-G,+d,+f,-E,-B,-C,+A) ;
	if(i== 26) MV = GA5_4_1(+a,+c,+b,+F,+e,+D,-g,-N,+m,-L,+K,+j,+H,+S,-r,+P,+p,-R,+s,-h,-J,-k,+l,-M,+n,-G,+d,+E,+f,-B,-C,+A) ;
	if(i== 27) MV = GA5_4_1(+a,+d,+b,+c,+e,+f,-h,-l,+p,+r,+g,+j,+k,+m,+n,+s,+S,-N,-M,-K,-J,+G,+R,+P,+L,+H,+F,+E,-C,-B,+D,+A) ;
	if(i== 28) MV = GA5_4_1(+a,+d,+b,+c,-F,-E,-h,-l,-R,-P,+g,-K,-J,+N,+M,+s,+S,+m,+n,-j,-k,+G,+p,+r,+L,+H,+e,+f,-C,-B,+D,+A) ;
	if(i== 29) MV = GA5_4_1(+a,+F,+b,+c,+d,-E,-K,+N,-R,-s,+g,+h,-J,+l,+M,-P,-p,+m,+L,-j,+H,+G,+S,+r,-n,+k,+e,-D,-C,-B,+f,+A) ;
	if(i== 30) MV = GA5_4_1(+a,+d,+b,+F,+c,-E,-h,+R,-l,-P,+K,+g,-J,+N,-s,+M,-m,+S,+n,-j,-G,-k,+p,-L,+r,+H,+e,+C,+f,-B,+D,+A) ;
	if(i== 31) MV = GA5_4_1(+a,+e,+b,+d,+c,+f,-j,-p,-m,+s,+h,+g,+k,-l,+r,+n,+N,+R,-L,-K,-G,-H,-S,+M,+P,-J,+F,+C,+D,-B,-E,+A) ;
	if(i== 32) MV = GA5_4_1(+a,+F,+b,+e,+c,+D,-K,+S,+N,-r,+j,+g,+H,-m,+P,-L,-l,-p,+M,-h,-G,+J,-R,+n,-s,+k,+d,+C,-E,-B,+f,+A) ;
	if(i== 33) MV = GA5_4_1(+a,-k,+b,-J,-H,+g,+f,-S,-R,-n,-E,-D,+c,+p,+M,+L,+l,+m,+P,-C,+d,+e,+N,+r,+s,+F,+G,-h,-j,-B,+K,+A) ;
	if(i== 34) MV = GA5_4_1(+a,-H,+b,+k,+J,+g,+D,-R,+p,+L,+f,+E,+c,-S,+n,-M,+m,-N,-s,-C,+e,+F,-l,+P,-r,+d,+G,-j,-K,-B,+h,+A) ;
	if(i== 35) MV = GA5_4_1(+a,-H,+b,-J,+k,+g,+D,-p,-R,+L,-E,+f,+c,-S,+M,+n,+N,+m,-s,-C,-F,+e,-l,+r,+P,+d,+G,+K,-j,-B,+h,+A) ;
	if(i== 36) MV = GA5_4_1(+a,+J,+b,+k,+G,+h,-E,+S,-m,+P,+f,+C,+d,-N,+r,-L,-l,+R,-n,+D,+c,+F,-p,+M,-s,+e,-H,-g,-K,-B,+j,+A) ;
	if(i== 37) MV = GA5_4_1(+a,-k,+b,+H,+G,+j,+f,+R,+N,-s,+D,+C,+e,+l,+P,+M,-m,-p,+L,-E,+c,+d,+S,+n,+r,+F,+J,-g,-h,-B,+K,+A) ;
	if(i== 38) MV = GA5_4_1(+a,-H,+b,-J,+G,+K,+D,-p,-l,+r,-E,+C,+F,+m,-s,+n,+N,+S,-M,+f,+c,+e,+R,-L,+P,+d,-k,-g,-j,-B,+h,+A) ;
	if(i== 39) MV = GA5_4_1(+a,-k,+b,-G,+J,+h,+f,-N,+S,-r,-C,+E,+d,-m,+L,+P,+p,-l,+M,+D,+e,+c,-R,+s,+n,+F,-H,-j,-g,-B,+K,+A) ;
	if(i== 40) MV = GA5_4_1(+a,+J,+b,-G,+k,+h,-E,+m,+S,+P,-C,+f,+d,-N,+L,+r,-R,-l,-n,+D,-F,+c,-p,+s,+M,+e,-H,+K,-g,-B,+j,+A) ;
	if(i== 41) MV = GA5_4_1(+a,+J,+b,-G,-H,+K,-E,+m,-p,+s,-C,-D,+F,-l,-n,+r,-R,+N,-L,+f,+d,+c,-S,-P,+M,+e,-k,-h,-g,-B,+j,+A) ;
	if(i== 42) MV = GA5_4_1(+a,+G,+b,+k,-H,+j,-C,+N,+l,+M,+f,-D,+e,+R,+s,-P,-p,-S,-r,-E,+d,+F,+m,+L,-n,+c,+J,-h,-K,-B,+g,+A) ;
	if(i== 43) MV = GA5_4_1(+a,+G,+b,+H,+k,+j,-C,-l,+N,+M,+D,+f,+e,+R,+P,+s,+S,-p,-r,-E,-F,+d,+m,+n,+L,+c,+J,+K,-h,-B,+g,+A) ;
	if(i== 44) MV = GA5_4_1(+a,+G,+b,+H,+J,+K,-C,-l,+m,+n,+D,+E,+F,-p,-r,+s,+S,-R,-P,+f,+e,+d,-N,-M,+L,+c,-k,-j,-h,-B,+g,+A) ;
	if(i== 45) MV = GA5_4_1(+a,+e,+b,+d,-F,-C,-j,-p,+S,-M,+h,-K,-G,-R,+L,+n,+N,-l,+r,-g,-k,-H,-m,+s,+P,-J,+c,+f,+D,-B,-E,+A) ;
	if(i== 46) MV = GA5_4_1(+a,+F,+b,+d,+e,-C,-K,-R,+S,-n,+h,+j,-G,+p,+L,-M,+m,-l,+P,-g,-J,-H,+N,+s,-r,+k,+c,+E,+D,-B,+f,+A) ;
	if(i== 47) MV = GA5_4_1(+a,+e,+b,+F,+d,-C,-j,-S,-p,-M,+K,+h,-G,-R,-n,+L,+l,+N,+r,-g,+H,-k,-m,-P,+s,-J,+c,-D,+f,-B,-E,+A) ;
	if(i== 48) MV = GA5_4_1(+a,+c,+d,+b,+e,+f,+l,-g,+m,+n,-h,+p,+r,+j,+k,+s,+S,+K,+J,-R,-P,+H,-N,-M,+G,+L,+F,+E,+B,-D,-C,+A) ;
	if(i== 49) MV = GA5_4_1(+a,+c,+d,+b,-F,-E,+l,-g,+N,+M,-h,-R,-P,-K,-J,+s,+S,+j,+k,-p,-r,+H,+m,+n,+G,+L,+e,+f,+B,-D,-C,+A) ;
	if(i== 50) MV = GA5_4_1(+a,+c,-F,+b,+e,+D,+N,-g,+m,-L,+K,-S,+r,+j,+H,+P,+p,+h,+J,-R,+s,-k,+l,-M,+G,+n,+d,+E,+B,+f,-C,+A) ;
	if(i== 51) MV = GA5_4_1(+a,+e,+c,+b,+d,+f,-m,-j,-p,+s,-g,+l,+n,+h,+k,+r,-R,+K,-H,+N,-L,-G,-S,-P,+J,+M,+F,-D,+B,+C,-E,+A) ;
	if(i== 52) MV = GA5_4_1(+a,+e,+c,+b,-F,+D,-m,-j,+S,+P,-g,+N,-L,-K,+H,+r,-R,+h,+k,-l,-n,-G,-p,+s,+J,+M,+d,+f,+B,+C,-E,+A) ;
	if(i== 53) MV = GA5_4_1(+a,+F,+c,+b,+e,+D,+N,-K,+S,-r,-g,+m,-L,+j,+H,+P,+p,+h,+J,-l,+M,-G,-R,+s,-k,+n,+d,+E,+B,+C,+f,+A) ;
	if(i== 54) MV = GA5_4_1(+a,+d,+e,+b,+c,+f,+p,-h,-l,+r,-j,-m,+s,+g,+k,+n,+N,+K,+G,+S,-M,-J,+R,-L,-H,+P,+F,+C,+B,+E,+D,+A) ;
	if(i== 55) MV = GA5_4_1(+a,+d,-F,+b,+c,-E,-R,-h,-l,-P,+K,-N,+s,+g,-J,+M,-m,+j,+G,+S,+n,-k,+p,-L,-H,+r,+e,+C,+B,+f,+D,+A) ;
	if(i== 56) MV = GA5_4_1(+a,+F,+d,+b,+c,-E,-R,-K,+N,-s,-h,-l,-P,+g,-J,+M,-m,+j,+G,-p,+L,+H,+S,+n,-k,+r,+e,+C,+B,-D,+f,+A) ;
	if(i== 57) MV = GA5_4_1(+a,-k,-H,+b,-J,+g,-R,+f,-S,-n,+D,-p,+L,-E,+c,+M,-m,+C,+e,+l,+P,+d,+N,-s,-F,+r,+G,+j,+B,-h,+K,+A) ;
	if(i== 58) MV = GA5_4_1(+a,-J,-H,+b,+k,+g,+p,+E,-S,+M,+D,-R,+L,+f,+c,+n,+N,+C,+F,+l,-r,+d,+m,-s,+e,+P,+G,+K,+B,-h,-j,+A) ;
	if(i== 59) MV = GA5_4_1(+a,+H,-k,+b,-J,+g,-R,-D,+p,-L,+f,-S,-n,-E,+c,+M,-m,+C,+e,-N,+s,+F,+l,+P,+d,+r,+G,+j,+B,-K,-h,+A) ;
	if(i== 60) MV = GA5_4_1(+a,-k,+J,+b,-G,+h,+S,+f,-N,-r,-E,+m,+P,-C,+d,+L,+l,-D,+c,+p,+M,+e,-R,-n,-F,+s,-H,+g,+B,-j,+K,+A) ;
	if(i== 61) MV = GA5_4_1(+a,-J,-k,+b,-G,+h,+S,+E,-m,-P,+f,-N,-r,-C,+d,+L,+l,-D,+c,+R,+n,+F,+p,+M,+e,+s,-H,+g,+B,-K,-j,+A) ;
	if(i== 62) MV = GA5_4_1(+a,-J,-H,+b,-G,-K,+p,+E,-m,+s,+D,+l,-r,-C,-F,+n,+N,+f,+c,+R,-L,+d,-S,+M,+e,+P,+k,+g,+B,-h,-j,+A) ;
	if(i== 63) MV = GA5_4_1(+a,-k,+G,+b,+H,+j,+N,+f,+R,-s,-C,-l,+M,+D,+e,+P,+p,+E,+d,-m,+L,+c,+S,-r,-F,+n,+J,+h,+B,-g,+K,+A) ;
	if(i== 64) MV = GA5_4_1(+a,+H,+G,+b,+k,+j,+l,-D,+R,+P,-C,+N,+M,+f,+e,+s,+S,+E,+F,-m,-n,+c,-p,-r,+d,+L,+J,+K,+B,-g,-h,+A) ;
	if(i== 65) MV = GA5_4_1(+a,+H,+G,+b,-J,-K,+l,-D,+p,+r,-C,-m,-n,-E,-F,+s,+S,+f,+e,-N,-M,+c,+R,+P,+d,+L,+k,+j,+B,-g,-h,+A) ;
	if(i== 66) MV = GA5_4_1(+a,-G,+J,+b,+k,+h,-m,+C,-N,+L,-E,+S,+P,+f,+d,+r,-R,-D,+F,+p,-s,+e,-l,-n,+c,+M,-H,+K,+B,-j,-g,+A) ;
	if(i== 67) MV = GA5_4_1(+a,-G,-k,+b,+H,+j,+N,+C,+l,-M,+f,+R,-s,+D,+e,+P,+p,+E,+d,-S,+r,+F,-m,+L,+c,+n,+J,+h,+B,-K,-g,+A) ;
	if(i== 68) MV = GA5_4_1(+a,-G,+J,+b,+H,-K,-m,+C,+l,+n,-E,+p,-s,+D,-F,+r,-R,+f,+d,-S,-P,+e,-N,+L,+c,+M,+k,+h,+B,-j,-g,+A) ;
	if(i== 69) MV = GA5_4_1(+a,+d,+e,+b,-F,-C,+p,-h,-R,+L,-j,+S,-M,-K,-G,+n,+N,+g,+k,+m,-s,-J,-l,+r,-H,+P,+c,+f,+B,+E,+D,+A) ;
	if(i== 70) MV = GA5_4_1(+a,+e,-F,+b,+d,-C,+S,-j,-p,-M,+K,+R,+n,+h,-G,+L,+l,+g,-H,+N,+r,-k,-m,-P,+J,+s,+c,-D,+B,+f,-E,+A) ;
	if(i== 71) MV = GA5_4_1(+a,+F,+e,+b,+d,-C,+S,-K,-R,-n,-j,-p,-M,+h,-G,+L,+l,+g,-H,+m,+P,-J,+N,+r,-k,+s,+c,-D,+B,+E,+f,+A) ;
	if(i== 72) MV = GA5_4_1(+a,+c,+d,+F,+b,-E,+l,-N,-g,+M,+R,-h,-P,-K,-s,-J,-j,+S,+k,-p,-H,-r,+m,-G,+n,+L,+e,-B,+f,-D,-C,+A) ;
	if(i== 73) MV = GA5_4_1(+a,+c,+e,+d,+b,+f,+m,+l,-g,+n,-p,-j,+s,-h,+r,+k,-K,+R,+H,+S,+J,-P,-N,-G,-L,-M,+F,-B,+D,+E,-C,+A) ;
	if(i== 74) MV = GA5_4_1(+a,+c,-F,+d,+b,-E,+N,+l,-g,+M,+R,+K,+s,-h,-P,-J,-j,+p,+H,+S,+k,-r,+m,-G,-L,+n,+e,-B,+D,+f,-C,+A) ;
	if(i== 75) MV = GA5_4_1(+a,+d,+c,+e,+b,+f,-l,+p,-h,+r,+m,-g,+n,-j,+s,+k,-K,-S,-J,+N,+G,-M,+R,+H,-P,-L,+F,-B,-E,+C,+D,+A) ;
	if(i== 76) MV = GA5_4_1(+a,+F,+c,+d,+b,-E,+N,-R,-K,-s,+l,-g,+M,-h,-P,-J,-j,+p,+H,-m,+G,+L,+S,+k,-r,+n,+e,-B,+D,+C,+f,+A) ;
	if(i== 77) MV = GA5_4_1(+a,+e,+c,+F,+b,+D,-m,-S,-j,+P,-N,-g,-L,-K,-r,+H,-h,-R,+k,-l,+G,-n,-p,-J,+s,+M,+d,-B,+f,+C,-E,+A) ;
	if(i== 78) MV = GA5_4_1(+a,+e,+d,+c,+b,+f,-p,-m,-j,+s,-l,-h,+r,-g,+n,+k,-K,-N,-G,-R,-H,-L,-S,-J,-M,-P,+F,-B,-C,-D,-E,+A) ;
	if(i== 79) MV = GA5_4_1(+a,+e,-F,+c,+b,+D,+S,-m,-j,+P,-N,+K,+r,-g,-L,+H,-h,+l,-G,-R,+k,-n,-p,-J,-M,+s,+d,-B,-C,+f,-E,+A) ;
	if(i== 80) MV = GA5_4_1(+a,+F,+e,+c,+b,+D,+S,+N,-K,-r,-m,-j,+P,-g,-L,+H,-h,+l,-G,+p,+J,+M,-R,+k,-n,+s,+d,-B,-C,+E,+f,+A) ;
	if(i== 81) MV = GA5_4_1(+a,-k,-J,-H,+b,+g,-S,-R,+f,-n,+p,+E,+M,+D,+L,+c,-C,-l,+d,-m,+e,+P,+N,+F,-r,-s,+G,-B,+h,+j,+K,+A) ;
	if(i== 82) MV = GA5_4_1(+a,+J,-k,-H,+b,+g,-S,-p,-E,-M,-R,+f,-n,+D,+L,+c,-C,-l,+d,-N,-F,+r,-m,+e,+P,-s,+G,-B,+h,-K,+j,+A) ;
	if(i== 83) MV = GA5_4_1(+a,+J,+H,-k,+b,+g,+p,-S,-E,-M,-R,-D,-L,+f,-n,+c,-C,+N,+F,-l,+d,+r,-m,+e,+s,+P,+G,-B,+K,+h,+j,+A) ;
	if(i== 84) MV = GA5_4_1(+a,-k,+H,+G,+b,+j,+R,+N,+f,-s,+l,-D,+P,-C,+M,+e,-E,+m,+c,+p,+d,+L,+S,+F,-n,-r,+J,-B,+g,+h,+K,+A) ;
	if(i== 85) MV = GA5_4_1(+a,-H,-k,+G,+b,+j,+R,-l,+D,-P,+N,+f,-s,-C,+M,+e,-E,+m,+c,-S,-F,+n,+p,+d,+L,-r,+J,-B,+g,-K,+h,+A) ;
	if(i== 86) MV = GA5_4_1(+a,-H,-J,+G,+b,+K,-p,-l,+D,+r,+m,+E,-s,-C,+n,+F,+f,-N,+c,-S,+e,-M,+R,+d,+L,-P,-k,-B,+g,+j,+h,+A) ;
	if(i== 87) MV = GA5_4_1(+a,-k,-G,+J,+b,+h,-N,+S,+f,-r,-m,+C,+L,-E,+P,+d,+D,-p,+e,+l,+c,+M,-R,+F,-s,-n,-H,-B,+j,+g,+K,+A) ;
	if(i== 88) MV = GA5_4_1(+a,+J,-G,-H,+b,+K,+m,-p,-E,+s,-l,+C,-n,+D,+r,+F,+f,+R,+d,-N,+c,-L,-S,+e,+P,-M,-k,-B,+h,+g,+j,+A) ;
	if(i== 89) MV = GA5_4_1(+a,-H,-G,-k,+b,+j,+l,+R,+D,-P,+N,+C,-M,+f,-s,+e,-E,+S,+F,+m,+c,+n,+p,+d,+r,+L,+J,-B,+K,+g,+h,+A) ;
	if(i== 90) MV = GA5_4_1(+a,+G,-k,+J,+b,+h,-N,+m,-C,-L,+S,+f,-r,-E,+P,+d,+D,-p,+e,+R,-F,+s,+l,+c,+M,-n,-H,-B,+j,-K,+g,+A) ;
	if(i== 91) MV = GA5_4_1(+a,+G,-J,-k,+b,+h,-m,-N,-C,-L,+S,+E,-P,+f,-r,+d,+D,-R,+F,-p,+e,+s,+l,+c,+n,+M,-H,-B,+K,+j,+g,+A) ;
	if(i== 92) MV = GA5_4_1(+a,+G,+H,+J,+b,+K,-l,+m,-C,+n,-p,-D,-r,-E,+s,+F,+f,-S,+e,+R,+d,-P,-N,+c,+M,-L,-k,-B,+j,+h,+g,+A) ;
	if(i== 93) MV = GA5_4_1(+a,+d,+e,+F,+b,-C,+p,+R,-h,+L,-S,-j,-M,-K,-n,-G,-g,+N,+k,+m,+J,-s,-l,+H,+r,+P,+c,-B,+f,+E,+D,+A) ;
	if(i== 94) MV = GA5_4_1(+a,+d,-F,+e,+b,-C,-R,+p,-h,+L,-S,+K,+n,-j,-M,-G,-g,-m,-J,+N,+k,-s,-l,+H,-P,+r,+c,-B,-E,+f,+D,+A) ;
	if(i== 95) MV = GA5_4_1(+a,+F,+d,+e,+b,-C,-R,+S,-K,-n,+p,-h,+L,-j,-M,-G,-g,-m,-J,+l,-H,+P,+N,+k,-s,+r,+c,-B,-E,-D,+f,+A) ;
	if(i== 96) MV = GA5_4_1(+a,+c,+n,+L,+M,-g,+f,-D,-E,+b,+R,+S,+k,-p,-H,-J,-j,-h,-P,-K,-s,-r,-B,+e,+d,+F,-G,+m,+l,-N,-C,+A) ;
	if(i== 97) MV = GA5_4_1(+a,+c,-L,+n,+M,-g,+D,+f,-E,+b,+R,+p,+H,+S,+k,-J,-j,+K,+s,-h,-P,-r,-B,+e,-F,+d,-G,+m,+N,+l,-C,+A) ;
	if(i== 98) MV = GA5_4_1(+a,+c,-L,-M,+n,-g,+D,+E,+f,+b,-p,+R,+H,+S,+J,+k,-K,-j,+s,-h,+r,-P,-B,+F,+e,+d,-G,-N,+m,+l,-C,+A) ;
	if(i== 99) MV = GA5_4_1(+a,-k,-n,-r,-s,+f,+g,+h,+j,+b,+l,+m,+c,+p,+d,+e,-E,-D,+P,-C,+M,+L,-B,+J,+H,+G,+F,+S,+R,+N,+K,+A) ;
	if(i==100) MV = GA5_4_1(+a,-H,+L,-P,+r,+D,+g,+j,+K,+b,+m,-N,+c,-S,+e,+F,+f,+E,-s,-C,+n,-M,-B,-k,-J,+G,+d,+R,-p,-l,+h,+A) ;
	if(i==101) MV = GA5_4_1(+a,-H,+L,-r,-P,+D,+g,-K,+j,+b,+N,+m,+c,-S,-F,+e,-E,+f,-s,-C,+M,+n,-B,+J,-k,+G,+d,+p,+R,-l,+h,+A) ;
	if(i==102) MV = GA5_4_1(+a,+J,+P,-M,+s,-E,+h,+g,+K,+b,-l,+R,+d,-N,+c,+F,+f,+C,-n,+D,+r,-L,-B,-k,-G,-H,+e,-S,+m,-p,+j,+A) ;
	if(i==103) MV = GA5_4_1(+a,-k,-s,-n,-r,+f,+j,+g,+h,+b,-m,-p,+e,+l,+c,+d,+D,+C,+L,-E,+P,+M,-B,-H,-G,+J,+F,-R,-N,+S,+K,+A) ;
	if(i==104) MV = GA5_4_1(+a,-H,+r,+L,-P,+D,+K,+g,+j,+b,+N,+S,+F,+m,+c,+e,-E,+C,-M,+f,-s,+n,-B,+J,-G,-k,+d,+p,+l,+R,+h,+A) ;
	if(i==105) MV = GA5_4_1(+a,-k,-r,-s,-n,+f,+h,+j,+g,+b,+p,-l,+d,-m,+e,+c,-C,+E,+M,+D,+L,+P,-B,+G,-J,-H,+F,+N,-S,-R,+K,+A) ;
	if(i==106) MV = GA5_4_1(+a,+J,+P,-s,-M,-E,+h,-K,+g,+b,-R,-l,+d,-N,-F,+c,-C,+f,-n,+D,+L,+r,-B,+G,-k,-H,+e,-m,-S,-p,+j,+A) ;
	if(i==107) MV = GA5_4_1(+a,+J,+s,+P,-M,-E,+K,+h,+g,+b,-R,+N,+F,-l,+d,+c,-C,-D,-L,+f,-n,+r,-B,+G,+H,-k,+e,-m,+p,-S,+j,+A) ;
	if(i==108) MV = GA5_4_1(+a,+d,+r,+P,+L,-h,+f,+E,-C,+b,-S,+N,+k,+m,+J,-G,-g,-j,-M,-K,-n,-s,-B,+c,+e,+F,+H,-l,+p,+R,+D,+A) ;
	if(i==109) MV = GA5_4_1(+a,+d,-P,+r,+L,-h,-E,+f,-C,+b,-S,-m,-J,+N,+k,-G,-g,+K,+n,-j,-M,-s,-B,+c,-F,+e,+H,-l,-R,+p,+D,+A) ;
	if(i==110) MV = GA5_4_1(+a,+F,-s,+r,-n,-K,-E,-D,-C,+b,+p,-m,-J,+l,-H,-G,-g,-h,+L,-j,-M,+P,-B,+c,+d,+e,+k,+N,-R,+S,+f,+A) ;
	if(i==111) MV = GA5_4_1(+a,+d,-P,-L,+r,-h,-E,+C,+f,+b,+m,-S,-J,+N,+G,+k,-K,-g,+n,-j,+s,-M,-B,+F,+c,+e,+H,+R,-l,+p,+D,+A) ;
	if(i==112) MV = GA5_4_1(+a,+e,+s,+M,+P,-j,+f,+C,+D,+b,-N,-R,+k,-l,+G,+H,-h,-g,-L,-K,-r,-n,-B,+d,+c,+F,-J,-p,-m,-S,-E,+A) ;
	if(i==113) MV = GA5_4_1(+a,+F,-r,+n,-s,-K,+D,+C,-E,+b,+l,+p,+H,-m,+G,-J,-j,-g,+M,-h,-P,+L,-B,+e,+c,+d,+k,+S,+N,-R,+f,+A) ;
	if(i==114) MV = GA5_4_1(+a,+e,-M,+s,+P,-j,-C,+f,+D,+b,-N,+l,-G,-R,+k,+H,-h,+K,+r,-g,-L,-n,-B,+d,-F,+c,-J,-p,+S,-m,-E,+A) ;
	if(i==115) MV = GA5_4_1(+a,+F,-n,+s,-r,-K,-C,+E,+D,+b,-m,+l,-G,+p,+J,+H,-h,-j,+P,-g,-L,+M,-B,+d,+e,+c,+k,-R,+S,+N,+f,+A) ;
	if(i==116) MV = GA5_4_1(+a,+e,-M,-P,+s,-j,-C,-D,+f,+b,-l,-N,-G,-R,-H,+k,-K,-h,+r,-g,+n,-L,-B,+F,+d,+c,-J,-S,-p,-m,-E,+A) ;
	if(i==117) MV = GA5_4_1(+a,+G,+M,-L,+n,-C,+j,+h,+K,+b,-p,-S,+e,+R,+d,+F,+f,-D,-r,-E,+s,-P,-B,-k,+H,+J,+c,-N,-l,+m,+g,+A) ;
	if(i==118) MV = GA5_4_1(+a,+G,+M,-n,-L,-C,+j,-K,+h,+b,+S,-p,+e,+R,-F,+d,+D,+f,-r,-E,+P,+s,-B,-H,-k,+J,+c,+l,-N,+m,+g,+A) ;
	if(i==119) MV = GA5_4_1(+a,+G,+n,+M,-L,-C,+K,+j,+h,+b,+S,-R,+F,-p,+e,+d,+D,+E,-P,+f,-r,+s,-B,-H,-J,-k,+c,+l,-m,-N,+g,+A) ;
	if(i==120) MV = GA5_4_1(+a,-M,+c,+n,+L,-g,-E,+S,-p,+J,+f,-D,+b,+R,+k,-H,+h,+K,-r,+B,+d,+F,-j,-P,-s,+e,-G,-l,+N,+C,+m,+A) ;
	if(i==121) MV = GA5_4_1(+a,-n,+c,-L,-M,-g,+f,+R,+S,-k,+D,+E,+b,-p,+H,+J,+j,+h,-P,+B,+e,+d,-K,+s,+r,+F,-G,-m,-l,+C,-N,+A) ;
	if(i==122) MV = GA5_4_1(+a,-M,+c,-L,+n,-g,-E,+p,+S,+J,+D,+f,+b,+R,+H,+k,-K,+h,-r,+B,-F,+d,-j,+s,-P,+e,-G,-N,-l,+C,+m,+A) ;
	if(i==123) MV = GA5_4_1(+a,-M,+J,+P,+s,-E,-g,+l,-N,+c,+h,+K,+b,+R,+d,+F,+f,-D,-r,+B,+k,-H,+C,-n,-L,-G,+e,-S,+p,-j,+m,+A) ;
	if(i==124) MV = GA5_4_1(+a,-n,-k,-s,-r,+f,-g,+m,+l,+c,+j,+h,+b,-p,+e,+d,+D,+E,-P,+B,+H,+J,+C,+L,+M,-G,+F,-R,-S,-K,-N,+A) ;
	if(i==125) MV = GA5_4_1(+a,-M,+J,-s,+P,-E,-g,+N,+l,+c,-K,+h,+b,+R,-F,+d,+D,+f,-r,+B,+H,+k,+C,+L,-n,-G,+e,-p,-S,-j,+m,+A) ;
	if(i==126) MV = GA5_4_1(+a,-r,-k,-n,-s,+f,-h,-l,+p,+d,+g,+j,+b,+m,+c,+e,-E,+C,-M,+B,-J,+G,-D,+P,+L,+H,+F,+S,-N,-K,+R,+A) ;
	if(i==127) MV = GA5_4_1(+a,-P,-H,+L,+r,+D,-j,-m,-S,+e,+g,+K,+b,-N,+c,+F,+f,+C,-n,+B,+k,+G,+E,-s,-M,-J,+d,+R,+l,-h,-p,+A) ;
	if(i==128) MV = GA5_4_1(+a,-r,-H,+L,-P,+D,+K,-N,-S,-F,+g,+j,+b,+m,+c,+e,-E,+C,-M,+B,-J,+G,+f,-s,+n,-k,+d,+p,+l,-h,+R,+A) ;
	if(i==129) MV = GA5_4_1(+a,-s,-k,-r,-n,+f,-j,-p,-m,+e,+h,+g,+b,-l,+d,+c,-C,-D,-L,+B,-G,-H,+E,+M,+P,-J,+F,+N,+R,-K,-S,+A) ;
	if(i==130) MV = GA5_4_1(+a,-s,+J,+P,-M,-E,+K,+R,-N,-F,+h,+g,+b,-l,+d,+c,-C,-D,-L,+B,-G,-H,+f,-n,+r,-k,+e,-m,+p,-j,-S,+A) ;
	if(i==131) MV = GA5_4_1(+a,-P,-H,-r,+L,+D,-j,+S,-m,+e,-K,+g,+b,-N,-F,+c,-C,+f,-n,+B,-G,+k,+E,+M,-s,-J,+d,-l,+R,-h,-p,+A) ;
	if(i==132) MV = GA5_4_1(+a,-r,+d,-P,-L,-h,+f,-S,+N,-k,-E,+C,+b,+m,-J,+G,+g,+j,-M,+B,+c,+e,-K,+n,+s,+F,+H,+l,-p,-D,+R,+A) ;
	if(i==133) MV = GA5_4_1(+a,-P,+e,+s,+M,-j,+D,-R,-l,-H,+f,+C,+b,-N,+k,+G,+g,+K,-n,+B,+c,+F,-h,-L,-r,+d,-J,+m,+S,+E,-p,+A) ;
	if(i==134) MV = GA5_4_1(+a,-r,-F,+s,-n,+K,+D,-p,-l,-H,-E,+C,+b,+m,-J,+G,+g,+j,-M,+B,+c,+e,-h,-L,+P,+d,-k,+N,+S,+f,+R,+A) ;
	if(i==135) MV = GA5_4_1(+a,-s,+e,-M,-P,-j,+f,-N,-R,-k,-C,-D,+b,-l,-G,-H,+h,+g,-L,+B,+d,+c,-K,+r,+n,+F,-J,+p,+m,+E,-S,+A) ;
	if(i==136) MV = GA5_4_1(+a,-s,-F,+n,-r,+K,-E,+m,-p,+J,-C,-D,+b,-l,-G,-H,+h,+g,-L,+B,+d,+c,-j,-P,+M,+e,-k,-R,+N,+f,-S,+A) ;
	if(i==137) MV = GA5_4_1(+a,-P,+e,-M,+s,-j,+D,+l,-R,-H,-C,+f,+b,-N,-G,+k,-K,+g,-n,+B,-F,+c,-h,+r,-L,+d,-J,-S,+m,+E,-p,+A) ;
	if(i==138) MV = GA5_4_1(+a,-L,+d,+r,+P,-h,-C,+N,+m,+G,+f,+E,+b,-S,+k,+J,+j,+K,-s,+B,+e,+F,-g,-M,-n,+c,+H,-p,-R,-D,-l,+A) ;
	if(i==139) MV = GA5_4_1(+a,-L,+d,-P,+r,-h,-C,-m,+N,+G,-E,+f,+b,-S,-J,+k,-K,+j,-s,+B,-F,+e,-g,+n,-M,+c,+H,+R,-p,-D,-l,+A) ;
	if(i==140) MV = GA5_4_1(+a,-n,-F,+r,-s,+K,-C,-l,+m,+G,+D,+E,+b,-p,+H,+J,+j,+h,-P,+B,+e,+d,-g,-M,+L,+c,-k,+S,-R,+f,-N,+A) ;
	if(i==141) MV = GA5_4_1(+a,-L,+G,+M,+n,-C,-h,+p,+R,+d,+j,+K,+b,-S,+e,+F,+f,+E,-s,+B,+k,+J,-D,-r,-P,+H,+c,-N,-m,-g,-l,+A) ;
	if(i==142) MV = GA5_4_1(+a,-L,+G,-n,+M,-C,-h,-R,+p,+d,-K,+j,+b,-S,-F,+e,-E,+f,-s,+B,-J,+k,-D,+P,-r,+H,+c,+m,-N,-g,-l,+A) ;
	if(i==143) MV = GA5_4_1(+a,-n,+G,+M,-L,-C,+K,-S,+R,-F,+j,+h,+b,-p,+e,+d,+D,+E,-P,+B,+H,+J,+f,-r,+s,-k,+c,+l,-m,-g,-N,+A) ;
	if(i==144) MV = GA5_4_1(+a,-n,-M,+c,-L,-g,+S,+f,+R,-k,-E,+p,+J,+D,+b,+H,-h,-B,+d,+j,-P,+e,-K,-r,-F,+s,-G,+l,-C,-m,-N,+A) ;
	if(i==145) MV = GA5_4_1(+a,+M,-n,+c,-L,-g,+S,+E,-p,-J,+f,+R,-k,+D,+b,+H,-h,-B,+d,+K,+r,+F,+j,-P,+e,+s,-G,+l,-C,+N,-m,+A) ;
	if(i==146) MV = GA5_4_1(+a,-L,-M,+c,+n,-g,-p,-D,+R,+H,-E,+S,+J,+f,+b,+k,-K,-B,+F,+j,-s,+e,+h,-r,+d,-P,-G,-N,-C,-m,-l,+A) ;
	if(i==147) MV = GA5_4_1(+a,-n,-r,-k,-s,+f,+l,-g,+m,+c,-h,+p,+d,+j,+b,+e,-E,-B,+J,+D,-P,+H,+C,-M,+G,+L,+F,+S,+K,-R,-N,+A) ;
	if(i==148) MV = GA5_4_1(+a,-M,+s,+J,+P,-E,-N,-g,+l,+c,-K,-R,+F,+h,+b,+d,+D,-B,-H,+f,-r,+k,+C,+L,+G,-n,+e,-p,+j,-S,+m,+A) ;
	if(i==149) MV = GA5_4_1(+a,+L,-P,-H,+r,+D,+m,-g,-N,+c,-j,-S,+e,+K,+b,+F,+f,-B,-k,-E,+s,-J,+C,-n,+G,-M,+d,+R,+h,+p,+l,+A) ;
	if(i==150) MV = GA5_4_1(+a,+P,-M,+J,+s,-E,-l,-h,+R,+d,-g,-N,+c,+K,+b,+F,+f,-B,-k,-C,+n,-G,-D,-r,-H,-L,+e,-S,+j,-m,+p,+A) ;
	if(i==151) MV = GA5_4_1(+a,-s,-n,-k,-r,+f,-m,-j,-p,+e,-g,+l,+c,+h,+b,+d,+D,-B,-H,-C,-L,-G,+E,-P,+J,+M,+F,-R,+K,+N,-S,+A) ;
	if(i==152) MV = GA5_4_1(+a,-s,-M,+J,+P,-E,-N,+K,+R,-F,-g,+l,+c,+h,+b,+d,+D,-B,-H,-C,-L,-G,+f,-r,+k,-n,+e,-p,+j,-m,-S,+A) ;
	if(i==153) MV = GA5_4_1(+a,-r,-s,-k,-n,+f,+p,-h,-l,+d,-j,-m,+e,+g,+b,+c,-C,-B,+G,-E,-M,-J,-D,-L,-H,+P,+F,+N,+K,+S,+R,+A) ;
	if(i==154) MV = GA5_4_1(+a,-P,+r,-H,+L,+D,-S,-j,-m,+e,-K,+N,+F,+g,+b,+c,-C,-B,+G,+f,-n,+k,+E,+M,+J,-s,+d,-l,+h,+R,-p,+A) ;
	if(i==155) MV = GA5_4_1(+a,-r,-P,-H,+L,+D,-S,+K,-N,-F,-j,-m,+e,+g,+b,+c,-C,-B,+G,-E,-M,-J,+f,-n,+k,-s,+d,-l,+h,+p,+R,+A) ;
	if(i==156) MV = GA5_4_1(+a,-s,-P,+e,-M,-j,-R,+f,-N,-k,+D,+l,-H,-C,+b,-G,-g,-B,+c,+h,-L,+d,-K,-n,-F,+r,-J,-m,-E,+p,-S,+A) ;
	if(i==157) MV = GA5_4_1(+a,-s,+r,+F,-n,-K,+p,+E,-m,-J,+D,+l,-H,-C,+b,-G,-g,-B,+c,+h,-L,+d,+j,+M,+e,+P,+k,+N,+f,+R,-S,+A) ;
	if(i==158) MV = GA5_4_1(+a,+P,-s,+e,-M,-j,-R,-D,-l,+H,+f,-N,-k,-C,+b,-G,-g,-B,+c,+K,+n,+F,+h,-L,+d,+r,-J,-m,-E,+S,+p,+A) ;
	if(i==159) MV = GA5_4_1(+a,-r,-L,+d,-P,-h,+N,+f,-S,-k,-C,-m,+G,-E,+b,-J,-j,-B,+e,+g,-M,+c,-K,-s,-F,+n,+H,+p,+D,+l,+R,+A) ;
	if(i==160) MV = GA5_4_1(+a,-P,-L,+d,+r,-h,+m,+E,-S,-J,-C,+N,+G,+f,+b,+k,-K,-B,+F,+g,-n,+c,+j,-s,+e,-M,+H,+R,+D,+l,-p,+A) ;
	if(i==161) MV = GA5_4_1(+a,-r,+n,+F,-s,-K,+l,-D,+p,+H,-C,-m,+G,-E,+b,-J,-j,-B,+e,+g,-M,+c,+h,+P,+d,+L,+k,+S,+f,-N,+R,+A) ;
	if(i==162) MV = GA5_4_1(+a,+L,-r,+d,-P,-h,+N,+C,+m,-G,+f,-S,-k,-E,+b,-J,-j,-B,+e,+K,+s,+F,+g,-M,+c,+n,+H,+p,+D,-R,+l,+A) ;
	if(i==163) MV = GA5_4_1(+a,-n,+s,+F,-r,-K,-m,+C,+l,-G,-E,+p,+J,+D,+b,+H,-h,-B,+d,+j,-P,+e,+g,+L,+c,+M,+k,-R,+f,-S,-N,+A) ;
	if(i==164) MV = GA5_4_1(+a,-M,-P,+e,+s,-j,-l,+C,-N,-G,+D,-R,-H,+f,+b,+k,-K,-B,+F,+h,-r,+d,+g,-n,+c,-L,-J,-S,-E,+p,+m,+A) ;
	if(i==165) MV = GA5_4_1(+a,-L,+n,+G,+M,-C,+R,-h,+p,+d,-K,+S,+F,+j,+b,+e,-E,-B,+J,+f,-s,+k,-D,+P,-H,-r,+c,+m,+g,-N,-l,+A) ;
	if(i==166) MV = GA5_4_1(+a,+M,-L,+G,+n,-C,-p,-j,-S,+e,-h,+R,+d,+K,+b,+F,+f,-B,-k,+D,+r,+H,+E,-s,+J,-P,+c,-N,+g,+l,-m,+A) ;
	if(i==167) MV = GA5_4_1(+a,-n,-L,+G,+M,-C,+R,+K,-S,-F,-h,+p,+d,+j,+b,+e,-E,-B,+J,+D,-P,+H,+f,-s,+k,-r,+c,+m,+g,+l,-N,+A) ;
	if(i==168) MV = GA5_4_1(+a,-n,-L,-M,+c,-g,+R,+S,+f,-k,-p,-D,+H,-E,+J,+b,+B,-j,+e,-h,+d,-P,-K,+F,-s,-r,-G,+C,+m,+l,-N,+A) ;
	if(i==169) MV = GA5_4_1(+a,+L,-n,-M,+c,-g,+R,+p,+D,-H,+S,+f,-k,-E,+J,+b,+B,-j,+e,+K,-F,+s,-h,+d,-P,-r,-G,+C,+m,+N,+l,+A) ;
	if(i==170) MV = GA5_4_1(+a,+L,+M,-n,+c,-g,-p,+R,+D,-H,+S,+E,-J,+f,-k,+b,+B,-K,+F,-j,+e,+s,-h,+d,+r,-P,-G,+C,-N,+m,+l,+A) ;
	if(i==171) MV = GA5_4_1(+a,-n,-s,-r,-k,+f,+m,+l,-g,+c,-p,-j,+e,-h,+d,+b,+B,-D,+H,-E,+J,-P,+C,-G,-L,-M,+F,-K,+R,+S,-N,+A) ;
	if(i==172) MV = GA5_4_1(+a,+L,-P,-r,-H,+D,+m,+N,-g,+c,+S,-j,+e,+K,-F,+b,+B,+f,-k,-E,+J,+s,+C,-G,-n,-M,+d,-h,+R,+p,+l,+A) ;
	if(i==173) MV = GA5_4_1(+a,+L,+r,-P,-H,+D,-N,+m,-g,+c,+S,-K,+F,-j,+e,+b,+B,+E,-J,+f,-k,+s,+C,-G,+M,-n,+d,-h,-p,+R,+l,+A) ;
	if(i==174) MV = GA5_4_1(+a,-r,-n,-s,-k,+f,-l,+p,-h,+d,+m,-g,+c,-j,+e,+b,+B,+E,-J,-C,+G,-M,-D,+H,-P,-L,+F,-K,-S,+N,+R,+A) ;
	if(i==175) MV = GA5_4_1(+a,+P,-M,-s,+J,-E,-l,-R,-h,+d,+N,-g,+c,+K,-F,+b,+B,+f,-k,-C,+G,+n,-D,+H,-r,-L,+e,-j,-S,-m,+p,+A) ;
	if(i==176) MV = GA5_4_1(+a,-r,+L,-P,-H,+D,-N,-S,+K,-F,+m,-g,+c,-j,+e,+b,+B,+E,-J,-C,+G,-M,+f,-k,+s,-n,+d,-h,-p,-l,+R,+A) ;
	if(i==177) MV = GA5_4_1(+a,+P,+s,-M,+J,-E,+R,-l,-h,+d,+N,-K,+F,-g,+c,+b,+B,+C,-G,+f,-k,+n,-D,+H,+L,-r,+e,-j,+m,-S,+p,+A) ;
	if(i==178) MV = GA5_4_1(+a,-s,-r,-n,-k,+f,-p,-m,-j,+e,-l,-h,+d,-g,+c,+b,+B,+C,-G,+D,-H,-L,+E,-J,-M,-P,+F,-K,-N,-R,-S,+A) ;
	if(i==179) MV = GA5_4_1(+a,-s,+P,-M,+J,-E,+R,-N,+K,-F,-l,-h,+d,-g,+c,+b,+B,+C,-G,+D,-H,-L,+f,-k,+n,-r,+e,-j,+m,-p,-S,+A) ;
	if(i==180) MV = GA5_4_1(+a,-r,-P,-L,+d,-h,-S,+N,+f,-k,+m,+E,-J,-C,+G,+b,+B,-g,+c,-j,+e,-M,-K,+F,-n,-s,+H,-D,-l,+p,+R,+A) ;
	if(i==181) MV = GA5_4_1(+a,+P,-r,-L,+d,-h,-S,-m,-E,+J,+N,+f,-k,-C,+G,+b,+B,-g,+c,+K,-F,+n,-j,+e,-M,-s,+H,-D,-l,-R,+p,+A) ;
	if(i==182) MV = GA5_4_1(+a,-r,+s,-n,-F,+K,-p,-l,+D,-H,+m,+E,-J,-C,+G,+b,+B,-g,+c,-j,+e,-M,-h,+d,+L,-P,-k,+f,-N,-S,+R,+A) ;
	if(i==183) MV = GA5_4_1(+a,+P,+L,-r,+d,-h,+m,-S,-E,+J,+N,+C,-G,+f,-k,+b,+B,-K,+F,-g,+c,+n,-j,+e,+s,-M,+H,-D,+R,-l,+p,+A) ;
	if(i==184) MV = GA5_4_1(+a,-s,-M,-P,+e,-j,-N,-R,+f,-k,-l,+C,-G,+D,-H,+b,+B,-h,+d,-g,+c,-L,-K,+F,-r,-n,-J,+E,-p,-m,-S,+A) ;
	if(i==185) MV = GA5_4_1(+a,-s,+n,-r,-F,+K,+m,-p,-E,+J,-l,+C,-G,+D,-H,+b,+B,-h,+d,-g,+c,-L,-j,+e,+P,-M,-k,+f,+R,-N,-S,+A) ;
	if(i==186) MV = GA5_4_1(+a,+M,-s,-P,+e,-j,-N,+l,-C,+G,-R,+f,-k,+D,-H,+b,+B,-h,+d,+K,-F,+r,-g,+c,-L,-n,-J,+E,-p,+S,-m,+A) ;
	if(i==187) MV = GA5_4_1(+a,+M,+P,-s,+e,-j,-l,-N,-C,+G,-R,-D,+H,+f,-k,+b,+B,-K,+F,-h,+d,+r,-g,+c,+n,-L,-J,+E,-S,-p,-m,+A) ;
	if(i==188) MV = GA5_4_1(+a,-n,+r,-s,-F,+K,-l,+m,-C,+G,-p,-D,+H,-E,+J,+b,+B,-j,+e,-h,+d,-P,-g,+c,+M,-L,-k,+f,-S,+R,-N,+A) ;
	if(i==189) MV = GA5_4_1(+a,+M,-L,-n,+G,-C,-p,+S,-j,+e,-R,-h,+d,+K,-F,+b,+B,+f,-k,+D,-H,+r,+E,-J,-s,-P,+c,-g,-N,+l,-m,+A) ;
	if(i==190) MV = GA5_4_1(+a,+M,+n,-L,+G,-C,-S,-p,-j,+e,-R,-K,+F,-h,+d,+b,+B,-D,+H,+f,-k,+r,+E,-J,+P,-s,+c,-g,-l,-N,-m,+A) ;
	if(i==191) MV = GA5_4_1(+a,-n,+M,-L,+G,-C,-S,+R,+K,-F,-p,-j,+e,-h,+d,+b,+B,-D,+H,-E,+J,-P,+f,-k,+r,-s,+c,-g,-l,+m,-N,+A) ;
	if(i==192) MV = GA5_4_1(+a,-n,+L,+c,+G,+m,-R,+f,+K,-s,+D,+h,-P,+B,+e,+J,+j,+E,+b,-p,+H,+d,+S,-k,-F,+r,-M,-g,-C,-l,-N,+A) ;
	if(i==193) MV = GA5_4_1(+a,-L,-n,+c,+G,+m,-R,-D,-h,+P,+f,+K,-s,+B,+e,+J,+j,+E,+b,-S,+k,+F,-p,+H,+d,+r,-M,-g,-C,+N,-l,+A) ;
	if(i==194) MV = GA5_4_1(+a,-L,-M,+c,+G,+N,-p,-D,-h,+r,-E,+j,-s,+B,-F,+k,-K,+f,+b,-S,-J,+e,+R,+H,+d,-P,+n,-g,-C,-m,-l,+A) ;
	if(i==195) MV = GA5_4_1(+a,-n,+M,-G,+c,+l,-S,-K,+f,-r,+j,+E,-P,+B,+H,+d,+D,+h,+b,-p,+e,+J,-R,+F,-k,-s,+L,+C,-g,+m,-N,+A) ;
	if(i==196) MV = GA5_4_1(+a,-M,-n,-G,+c,+l,-S,-j,-E,+P,-K,+f,-r,+B,+H,+d,+D,+h,+b,+R,-F,+k,-p,+e,+J,-s,+L,+C,-g,+N,+m,+A) ;
	if(i==197) MV = GA5_4_1(+a,-M,-L,-G,+c,-N,+p,-j,-E,+s,+h,-D,-r,+B,+k,+F,+f,+K,+b,+R,+d,-H,-S,+e,+J,+P,-n,+C,-g,+l,+m,+A) ;
	if(i==198) MV = GA5_4_1(+a,-r,+P,+d,-H,-l,+S,+f,+K,-n,-E,+j,-M,+B,+c,+G,+g,+C,+b,+m,-J,+e,+N,-k,-F,+s,-L,-h,+D,-p,+R,+A) ;
	if(i==199) MV = GA5_4_1(+a,-P,-r,+d,-H,-l,+S,+E,-j,+M,+f,+K,-n,+B,+c,+G,+g,+C,+b,-N,+k,+F,+m,-J,+e,+s,-L,-h,+D,-R,-p,+A) ;
	if(i==200) MV = GA5_4_1(+a,-r,+s,+F,-k,+N,-p,-D,-h,-L,-E,+j,-M,+B,+c,+G,+g,+C,+b,+m,-J,+e,+l,+H,+d,-P,+n,-K,+f,-S,+R,+A) ;
	if(i==201) MV = GA5_4_1(+a,-s,+P,-J,+e,-m,+R,-K,+f,-n,+h,-D,-L,+B,-G,+c,-C,+g,+b,-l,+d,-H,+N,+F,-k,-r,+M,+E,-j,-p,-S,+A) ;
	if(i==202) MV = GA5_4_1(+a,-s,+r,-k,-F,-N,+p,-j,-E,-M,+h,-D,-L,+B,-G,+c,-C,+g,+b,-l,+d,-H,-m,+e,+J,+P,-n,+f,+K,+R,-S,+A) ;
	if(i==203) MV = GA5_4_1(+a,-P,-s,-J,+e,-m,+R,-h,+D,+L,-K,+f,-n,+B,-G,+c,-C,+g,+b,-N,-F,+k,-l,+d,-H,-r,+M,+E,-j,+S,-p,+A) ;
	if(i==204) MV = GA5_4_1(+a,-r,+L,+H,+d,+p,-N,-K,+f,-s,+g,+C,-M,+B,-J,+e,-E,+j,+b,+m,+c,+G,+S,+F,-k,-n,+P,-D,-h,-l,+R,+A) ;
	if(i==205) MV = GA5_4_1(+a,-P,-M,-J,+e,-S,+l,-h,+D,+r,+g,+C,-n,+B,+k,+F,+f,+K,+b,-N,+c,+G,+R,+d,-H,+L,-s,+E,-j,-m,-p,+A) ;
	if(i==206) MV = GA5_4_1(+a,-r,+n,-k,-F,-S,+l,-h,+D,-P,+g,+C,-M,+B,-J,+e,-E,+j,+b,+m,+c,+G,+p,+d,-H,+L,-s,+f,+K,-N,+R,+A) ;
	if(i==207) MV = GA5_4_1(+a,-P,-L,+d,-H,-R,+m,+E,-j,+s,-C,+g,-n,+B,-F,+k,-K,+f,+b,-N,-G,+c,-S,-J,+e,-M,+r,-h,+D,+l,-p,+A) ;
	if(i==208) MV = GA5_4_1(+a,-s,+M,+e,+J,-p,+N,+f,+K,-r,-C,+g,-L,+B,+d,-H,+h,-D,+b,-l,-G,+c,-R,-k,-F,+n,-P,-j,-E,+m,-S,+A) ;
	if(i==209) MV = GA5_4_1(+a,-s,+n,+F,-k,-R,+m,+E,-j,-P,-C,+g,-L,+B,+d,-H,+h,-D,+b,-l,-G,+c,+p,-J,+e,-M,+r,-K,+f,-N,-S,+A) ;
	if(i==210) MV = GA5_4_1(+a,-L,-r,+H,+d,+p,-N,-g,-C,+M,-K,+f,-s,+B,-J,+e,-E,+j,+b,-S,-F,+k,+m,+c,+G,-n,+P,-D,-h,-R,-l,+A) ;
	if(i==211) MV = GA5_4_1(+a,-L,-P,+H,+d,+R,-m,-g,-C,+n,+j,+E,-s,+B,+k,+F,+f,+K,+b,-S,+e,+J,-N,+c,+G,+M,-r,-D,-h,+p,-l,+A) ;
	if(i==212) MV = GA5_4_1(+a,-n,+s,-k,-F,+R,-m,-g,-C,-L,+j,+E,-P,+B,+H,+d,+D,+h,+b,-p,+e,+J,+l,+c,+G,+M,-r,+f,+K,-S,-N,+A) ;
	if(i==213) MV = GA5_4_1(+a,-M,-s,+e,+J,-p,+N,+C,-g,+L,+f,+K,-r,+B,+d,-H,+h,-D,+b,+R,+k,+F,-l,-G,+c,+n,-P,-j,-E,+S,+m,+A) ;
	if(i==214) MV = GA5_4_1(+a,-M,-P,+e,+J,+S,-l,+C,-g,+n,+D,+h,-r,+B,-F,+k,-K,+f,+b,+R,+H,+d,-N,-G,+c,-L,+s,-j,-E,+p,+m,+A) ;
	if(i==215) MV = GA5_4_1(+a,-n,+r,+F,-k,+S,-l,+C,-g,-M,+D,+h,-P,+B,+e,+J,+j,+E,+b,-p,+H,+d,-m,-G,+c,-L,+s,-K,+f,+R,-N,+A) ;
	if(i==216) MV = GA5_4_1(+a,-n,+c,+M,-G,+l,+f,-S,-K,-r,-E,-B,+d,+j,-P,+H,-h,+p,+J,+D,+b,+e,-R,+k,+s,+F,+L,+g,-m,+C,-N,+A) ;
	if(i==217) MV = GA5_4_1(+a,-M,+c,-L,-G,-N,-E,+p,-j,+s,+D,-B,+F,+h,-r,+k,-K,-R,-H,+f,+b,+d,-S,-J,-P,+e,-n,+g,-l,+C,+m,+A) ;
	if(i==218) MV = GA5_4_1(+a,+L,+c,+n,-G,+m,+D,-R,-h,-P,+f,-B,+e,+K,+s,-J,-j,-S,-k,-E,+b,+F,+p,+H,-r,+d,-M,+g,+N,+C,+l,+A) ;
	if(i==219) MV = GA5_4_1(+a,-M,+G,-n,+c,+l,+j,-S,-E,+P,-K,-B,-H,+f,-r,+d,+D,-R,+F,+h,+b,+k,-p,+e,+s,+J,+L,+C,-N,-g,+m,+A) ;
	if(i==220) MV = GA5_4_1(+a,-n,+G,+L,+c,+m,+K,-R,+f,-s,-h,-B,+J,+D,-P,+e,-E,+p,+d,+j,+b,+H,+S,+F,-r,-k,-M,+C,+l,-g,-N,+A) ;
	if(i==221) MV = GA5_4_1(+a,+L,+G,-M,+c,-N,+h,+p,+D,+r,-j,-B,-k,-E,+s,+F,+f,-S,+e,+K,+b,-J,+R,+d,-P,-H,-n,+C,+m,-g,+l,+A) ;
	if(i==222) MV = GA5_4_1(+a,+P,+d,+r,+H,-l,-E,+S,-j,-M,+f,-B,+c,+K,+n,-G,-g,-N,-k,-C,+b,+F,-m,-J,-s,+e,-L,+h,-R,-D,+p,+A) ;
	if(i==223) MV = GA5_4_1(+a,-s,+e,+P,-J,-m,+f,+R,-K,-n,+D,-B,+c,+h,-L,-G,-g,+l,-H,-C,+b,+d,+N,+k,+r,+F,+M,+j,+p,+E,-S,+A) ;
	if(i==224) MV = GA5_4_1(+a,-s,-F,+r,-k,-N,-E,+p,-j,-M,+D,-B,+c,+h,-L,-G,-g,+l,-H,-C,+b,+d,-m,-J,-P,+e,-n,-K,-R,+f,-S,+A) ;
	if(i==225) MV = GA5_4_1(+a,-r,-H,+P,+d,-l,+K,+S,+f,-n,-j,-B,+G,-E,-M,+c,-C,-m,+e,+g,+b,-J,+N,+F,-s,-k,-L,-D,+p,-h,+R,+A) ;
	if(i==226) MV = GA5_4_1(+a,-P,+J,-s,+e,-m,+h,+R,+D,+L,-K,-B,+G,+f,-n,+c,-C,+N,+F,+g,+b,+k,-l,+d,+r,-H,+M,+E,-S,-j,-p,+A) ;
	if(i==227) MV = GA5_4_1(+a,-r,+k,-s,-F,-N,+h,+p,+D,+L,-j,-B,+G,-E,-M,+c,-C,-m,+e,+g,+b,-J,-l,+d,-P,-H,-n,+f,-S,+K,+R,+A) ;
	if(i==228) MV = GA5_4_1(+a,+P,-H,-L,+d,+R,+j,-m,-E,+s,-g,-B,-k,-C,+n,+F,+f,-N,+c,+K,+b,-G,-S,+e,-M,+J,-r,-D,-l,-h,+p,+A) ;
	if(i==229) MV = GA5_4_1(+a,-s,+J,+M,+e,-p,+K,+N,+f,-r,-g,-B,-H,-C,-L,+d,+D,+l,+c,+h,+b,-G,-R,+F,-n,-k,-P,+E,-m,-j,-S,+A) ;
	if(i==230) MV = GA5_4_1(+a,-s,+k,-n,-F,+R,+j,-m,-E,+P,-g,-B,-H,-C,-L,+d,+D,+l,+c,+h,+b,-G,-p,+e,-M,+J,-r,+f,-N,+K,-S,+A) ;
	if(i==231) MV = GA5_4_1(+a,-r,+d,+L,+H,+p,+f,-N,-K,-s,-C,-B,+e,+g,-M,-J,-j,-m,+G,-E,+b,+c,+S,+k,+n,+F,+P,+h,+l,-D,+R,+A) ;
	if(i==232) MV = GA5_4_1(+a,-P,+e,-M,-J,-S,+D,+l,-h,+r,-C,-B,+F,+g,-n,+k,-K,+N,+G,+f,+b,+c,+R,+H,-L,+d,-s,+j,+m,+E,-p,+A) ;
	if(i==233) MV = GA5_4_1(+a,-r,-F,+n,-k,-S,+D,+l,-h,-P,-C,-B,+e,+g,-M,-J,-j,-m,+G,-E,+b,+c,+p,+H,-L,+d,-s,-K,+N,+f,+R,+A) ;
	if(i==234) MV = GA5_4_1(+a,-L,-H,-r,+d,+p,+g,-N,-C,+M,-K,-B,+J,+f,-s,+e,-E,+S,+F,+j,+b,+k,+m,+c,+n,+G,+P,-D,+R,-h,-l,+A) ;
	if(i==235) MV = GA5_4_1(+a,+M,+J,-P,+e,-S,+g,+l,-C,+n,-h,-B,-k,+D,+r,+F,+f,+R,+d,+K,+b,+H,-N,+c,-L,+G,-s,+E,-p,-j,-m,+A) ;
	if(i==236) MV = GA5_4_1(+a,-n,+k,-r,-F,-S,+g,+l,-C,+M,-h,-B,+J,+D,-P,+e,-E,+p,+d,+j,+b,+H,+m,+c,-L,+G,-s,+f,+R,+K,-N,+A) ;
	if(i==237) MV = GA5_4_1(+a,-L,+d,-P,+H,+R,-C,-m,-g,+n,-E,-B,+F,+j,-s,+k,-K,+S,+J,+f,+b,+e,-N,-G,-M,+c,-r,+h,-p,-D,-l,+A) ;
	if(i==238) MV = GA5_4_1(+a,+M,+e,+s,-J,-p,-C,+N,-g,-L,+f,-B,+d,+K,+r,+H,-h,+R,-k,+D,+b,+F,+l,-G,-n,+c,-P,+j,+S,+E,-m,+A) ;
	if(i==239) MV = GA5_4_1(+a,-n,-F,+s,-k,+R,-C,-m,-g,-L,-E,-B,+d,+j,-P,+H,-h,+p,+J,+D,+b,+e,+l,-G,-M,+c,-r,-K,+S,+f,-N,+A) ;
	if(i==240) MV = GA5_4_1(+a,-n,+c,+G,+L,+m,+f,+K,-R,-s,+B,-D,+e,-h,+J,-P,-p,-j,+H,-E,+d,+b,+S,+r,+k,+F,-M,-l,+g,+C,-N,+A) ;
	if(i==241) MV = GA5_4_1(+a,+L,+c,+G,+n,+m,+D,+h,-R,-P,+B,+f,+e,+K,+J,+s,+S,-j,-k,-E,-F,+b,+p,+r,+H,+d,-M,-N,+g,+C,+l,+A) ;
	if(i==242) MV = GA5_4_1(+a,+L,+c,+G,-M,-N,+D,+h,+p,+r,+B,+E,+F,-j,-k,+s,+S,-K,-J,+f,+e,+b,+R,+P,+H,+d,-n,-m,+g,+C,+l,+A) ;
	if(i==243) MV = GA5_4_1(+a,-n,-G,+c,+M,+l,-K,+f,-S,-r,+B,-j,+H,-E,+d,-P,-p,-D,+e,-h,+J,+b,-R,-s,-F,+k,+L,+m,-C,+g,-N,+A) ;
	if(i==244) MV = GA5_4_1(+a,+M,-G,+c,+n,+l,+j,+E,-S,-P,+B,-K,+H,+f,+d,+r,-R,-D,+F,-h,-k,+b,+p,-s,+e,+J,+L,-N,-C,+g,-m,+A) ;
	if(i==245) MV = GA5_4_1(+a,+M,-G,+c,-L,+N,+j,+E,-p,+s,+B,-h,-k,+D,-F,+r,-R,+f,+d,+K,-H,+b,-S,-P,+e,+J,+n,+l,-C,+g,-m,+A) ;
	if(i==246) MV = GA5_4_1(+a,-r,+d,-H,+P,-l,+f,+K,+S,-n,+B,+E,+c,-j,+G,-M,+m,-g,-J,-C,+e,+b,+N,+s,+k,+F,-L,-p,+h,-D,+R,+A) ;
	if(i==247) MV = GA5_4_1(+a,+P,+d,-H,+r,-l,-E,+j,+S,-M,+B,+f,+c,+K,+G,+n,+N,-g,-k,-C,-F,+b,-m,+s,-J,+e,-L,+R,+h,-D,+p,+A) ;
	if(i==248) MV = GA5_4_1(+a,-r,-F,+k,-s,-N,+D,+h,+p,+L,+B,+E,+c,-j,+G,-M,+m,-g,-J,-C,+e,+b,-l,+P,+H,+d,-n,+S,-K,+f,+R,+A) ;
	if(i==249) MV = GA5_4_1(+a,-s,-J,+e,+P,-m,-K,+f,+R,-n,+B,-h,-G,+D,+c,-L,-l,+C,+d,-g,-H,+b,+N,-r,-F,+k,+M,-p,-E,+j,-S,+A) ;
	if(i==250) MV = GA5_4_1(+a,-s,+k,+F,-r,+N,+j,+E,-p,+M,+B,-h,-G,+D,+c,-L,-l,+C,+d,-g,-H,+b,+m,-P,+e,+J,+n,-R,+f,+K,-S,+A) ;
	if(i==251) MV = GA5_4_1(+a,+P,-J,+e,+s,-m,+h,-D,+R,-L,+B,-K,-G,+f,+c,+n,+N,+C,+F,-g,-k,+b,+l,-r,+d,-H,+M,-S,-E,+j,+p,+A) ;
	if(i==252) MV = GA5_4_1(+a,-r,+H,+d,+L,+p,-K,+f,-N,-s,+B,-g,-J,-C,+e,-M,+m,+E,+c,-j,+G,+b,+S,-n,-F,+k,+P,-l,+D,+h,+R,+A) ;
	if(i==253) MV = GA5_4_1(+a,+P,-J,+e,-M,+S,+h,-D,-l,+r,+B,-g,-k,-C,-F,+n,+N,+f,+c,+K,+G,+b,+R,-L,+d,-H,+s,-m,-E,+j,+p,+A) ;
	if(i==254) MV = GA5_4_1(+a,-r,+k,+F,-n,+S,+h,-D,-l,+P,+B,-g,-J,-C,+e,-M,+m,+E,+c,-j,+G,+b,-p,-L,+d,-H,+s,+N,+f,+K,+R,+A) ;
	if(i==255) MV = GA5_4_1(+a,+P,+d,-H,-L,+R,-E,+j,-m,+s,+B,+C,+F,-g,-k,+n,+N,-K,-G,+f,+c,+b,-S,+M,-J,+e,-r,+l,+h,-D,+p,+A) ;
	if(i==256) MV = GA5_4_1(+a,-s,+e,+J,+M,-p,+f,+K,+N,-r,+B,+C,+d,-g,-H,-L,-l,-h,-G,+D,+c,+b,-R,+n,+k,+F,-P,+m,+j,+E,-S,+A) ;
	if(i==257) MV = GA5_4_1(+a,-s,-F,+k,-n,+R,-E,+j,-m,+P,+B,+C,+d,-g,-H,-L,-l,-h,-G,+D,+c,+b,-p,+M,-J,+e,-r,+N,-K,+f,-S,+A) ;
	if(i==258) MV = GA5_4_1(+a,+L,+H,+d,+r,+p,+g,+C,-N,-M,+B,-K,-J,+f,+e,+s,+S,+E,+F,-j,-k,+b,-m,-n,+c,+G,+P,+R,+D,+h,+l,+A) ;
	if(i==259) MV = GA5_4_1(+a,+L,+H,+d,-P,-R,+g,+C,+m,+n,+B,-j,-k,-E,-F,+s,+S,+f,+e,+K,+J,+b,-N,-M,+c,+G,+r,+p,+D,+h,+l,+A) ;
	if(i==260) MV = GA5_4_1(+a,-n,+k,+F,-s,-R,+g,+C,+m,+L,+B,-j,+H,-E,+d,-P,-p,-D,+e,-h,+J,+b,-l,-M,+c,+G,+r,+S,+f,+K,-N,+A) ;
	if(i==261) MV = GA5_4_1(+a,+M,+e,+J,+s,-p,-C,+g,+N,-L,+B,+f,+d,+K,-H,+r,-R,-h,-k,+D,-F,+b,+l,+n,-G,+c,-P,-S,+j,+E,-m,+A) ;
	if(i==262) MV = GA5_4_1(+a,+M,+e,+J,-P,-S,-C,+g,+l,+n,+B,-D,+F,-h,-k,+r,-R,-K,+H,+f,+d,+b,-N,+L,-G,+c,-s,+p,+j,+E,-m,+A) ;
	if(i==263) MV = GA5_4_1(+a,-n,-F,+k,-r,-S,-C,+g,+l,+M,+B,-D,+e,-h,+J,-P,-p,-j,+H,-E,+d,+b,+m,+L,-G,+c,-s,-R,-K,+f,-N,+A) ;
	if(i==264) MV = GA5_4_1(+a,+c,+n,-M,+G,+l,+f,+E,+B,+d,-S,-K,+r,+j,+P,-H,+h,-p,+J,-R,-k,-s,-D,+b,+e,+F,+L,-g,+m,-N,-C,+A) ;
	if(i==265) MV = GA5_4_1(+a,+c,+M,+n,+G,+l,-E,+f,+B,+d,-S,-j,-P,-K,+r,-H,+h,+R,+k,-p,+J,-s,-D,+b,-F,+e,+L,-g,+N,+m,-C,+A) ;
	if(i==266) MV = GA5_4_1(+a,+c,-L,-M,+G,+N,+D,+E,+B,-F,-p,-h,+r,+j,-s,+k,-K,+S,+J,-R,-H,-P,+f,+b,+e,+d,+n,-g,+m,+l,-C,+A) ;
	if(i==267) MV = GA5_4_1(+a,-G,-n,+L,+c,+m,+K,+h,+B,-J,-R,+f,-s,+D,-P,+e,-E,+p,+d,-S,-F,+r,+j,+b,+H,-k,-M,+C,+l,+N,-g,+A) ;
	if(i==268) MV = GA5_4_1(+a,-G,+M,+L,+c,-N,-j,+h,+B,+k,+p,+E,-s,+D,+r,+F,+f,+R,+d,-S,+e,+P,+K,+b,+H,-J,-n,+C,+l,+m,-g,+A) ;
	if(i==269) MV = GA5_4_1(+a,-G,-L,-n,+c,+m,-h,+K,+B,-J,-R,-D,+P,+f,-s,+e,-E,+S,+F,+p,+d,+r,+j,+b,+k,+H,-M,+C,-N,+l,-g,+A) ;
	if(i==270) MV = GA5_4_1(+a,+e,+s,-P,+J,-m,+f,-D,+B,+c,+R,-K,+n,+h,+L,+G,+g,-l,-H,+N,-k,-r,+C,+b,+d,+F,+M,-j,-p,-S,-E,+A) ;
	if(i==271) MV = GA5_4_1(+a,+e,+P,+s,+J,-m,+D,+f,+B,+c,+R,-h,-L,-K,+n,+G,+g,-N,+k,-l,-H,-r,+C,+b,-F,+d,+M,-j,+S,-p,-E,+A) ;
	if(i==272) MV = GA5_4_1(+a,+F,-r,+s,-k,+N,+D,+E,+B,+c,-p,-h,-L,+j,-M,+G,+g,-m,+J,-l,-H,-P,+C,+b,+e,+d,+n,-K,+S,-R,+f,+A) ;
	if(i==273) MV = GA5_4_1(+a,+H,-r,+P,+d,-l,+K,+j,+B,-G,+S,+f,-n,-E,-M,+c,-C,-m,+e,-N,-F,+s,+g,+b,-J,-k,-L,-D,+p,-R,-h,+A) ;
	if(i==274) MV = GA5_4_1(+a,+H,-P,-r,+d,-l,-j,+K,+B,-G,+S,+E,+M,+f,-n,+c,-C,+N,+F,-m,+e,+s,+g,+b,+k,-J,-L,-D,+R,+p,-h,+A) ;
	if(i==275) MV = GA5_4_1(+a,-k,+s,-r,-F,-N,-j,+h,+B,-G,+p,+E,+M,+D,+L,+c,-C,-l,+d,-m,+e,+P,+g,+b,+H,-J,-n,+f,+R,-S,+K,+A) ;
	if(i==276) MV = GA5_4_1(+a,-J,-s,+M,+e,-p,+K,+g,+B,+H,+N,+f,-r,-C,-L,+d,+D,+l,+c,+R,-F,+n,+h,+b,-G,-k,-P,+E,-m,+S,-j,+A) ;
	if(i==277) MV = GA5_4_1(+a,-J,+P,+M,+e,-S,-h,+g,+B,+k,+l,-D,-r,-C,+n,+F,+f,-N,+c,+R,+d,+L,+K,+b,-G,+H,-s,+E,-m,-p,-j,+A) ;
	if(i==278) MV = GA5_4_1(+a,-k,+r,-n,-F,-S,-h,+g,+B,-J,+l,-D,+P,-C,+M,+e,-E,+m,+c,+p,+d,+L,+j,+b,-G,+H,-s,+f,-N,+R,+K,+A) ;
	if(i==279) MV = GA5_4_1(+a,+H,+L,+P,+d,+R,-g,+j,+B,+k,-m,+C,-n,-E,+s,+F,+f,-S,+e,-N,+c,+M,+K,+b,-J,-G,-r,-D,+p,-l,-h,+A) ;
	if(i==280) MV = GA5_4_1(+a,-J,-M,-s,+e,-p,-g,+K,+B,+H,+N,+C,+L,+f,-r,+d,+D,-R,+F,+l,+c,+n,+h,+b,+k,-G,-P,+E,-S,-m,-j,+A) ;
	if(i==281) MV = GA5_4_1(+a,-k,+n,-s,-F,+R,-g,+j,+B,+H,-m,+C,+L,-E,+P,+d,+D,-p,+e,+l,+c,+M,+h,+b,-J,-G,-r,+f,-S,-N,+K,+A) ;
	if(i==282) MV = GA5_4_1(+a,+d,+r,-L,-H,+p,+f,+C,+B,+e,-N,-K,+s,+g,+M,+J,+j,+m,+G,+S,-k,-n,+E,+b,+c,+F,+P,-h,-l,+R,+D,+A) ;
	if(i==283) MV = GA5_4_1(+a,+d,-P,-L,-H,-R,-E,+C,+B,-F,+m,-j,+s,+g,-n,+k,-K,+N,+G,+S,+J,-M,+f,+b,+c,+e,+r,-h,-l,+p,+D,+A) ;
	if(i==284) MV = GA5_4_1(+a,+F,-s,+n,-k,-R,-E,+C,+B,+d,+m,-j,-P,+g,-L,-H,+h,+l,+G,-p,+J,-M,-D,+b,+c,+e,+r,-K,+N,+S,+f,+A) ;
	if(i==285) MV = GA5_4_1(+a,+d,+L,+r,-H,+p,-C,+f,+B,+e,-N,-g,-M,-K,+s,+J,+j,-S,+k,+m,+G,-n,+E,+b,-F,+c,+P,-h,-R,-l,+D,+A) ;
	if(i==286) MV = GA5_4_1(+a,+e,-M,-P,+J,+S,-C,-D,+B,-F,-l,-g,+n,+h,-r,+k,-K,-R,-H,+N,+G,-L,+f,+b,+d,+c,+s,-j,-p,-m,-E,+A) ;
	if(i==287) MV = GA5_4_1(+a,+F,-n,+r,-k,+S,-C,-D,+B,+e,-l,-g,-M,+h,-P,+J,+j,+p,-H,+m,+G,-L,+E,+b,+d,+c,+s,-K,-R,+N,+f,+A) ;
	if(i==288) MV = GA5_4_1(+a,+c,+M,-G,+n,+l,-E,-B,+f,+d,+j,-S,-P,-K,+H,+r,-R,+h,+k,-p,+s,+J,-D,+F,+b,+e,+L,-N,-g,+m,-C,+A) ;
	if(i==289) MV = GA5_4_1(+a,+c,+n,-G,-L,+m,+f,-B,+D,+e,+K,-R,+s,-h,-J,+P,+p,+j,+H,+S,-r,-k,+E,+d,+b,+F,-M,+l,-g,-N,-C,+A) ;
	if(i==290) MV = GA5_4_1(+a,+c,+M,-G,-L,+N,-E,-B,+D,-F,+j,-p,+s,-h,-k,+r,-R,-K,+H,+S,+P,+J,+f,+d,+b,+e,+n,+l,-g,+m,-C,+A) ;
	if(i==291) MV = GA5_4_1(+a,+G,-n,+c,+M,+l,-K,-B,+j,-H,+f,-S,-r,-E,+d,-P,-p,-D,+e,+R,+s,+F,-h,+J,+b,+k,+L,+m,-C,+N,+g,+A) ;
	if(i==292) MV = GA5_4_1(+a,+G,+L,+c,+n,+m,-h,-B,+K,+J,+D,-R,-P,+f,+e,+s,+S,+E,+F,-p,-r,+d,-j,-k,+b,+H,-M,-N,-C,-l,+g,+A) ;
	if(i==293) MV = GA5_4_1(+a,+G,+L,+c,+M,+N,-h,-B,+j,+k,+D,-p,-r,-E,-F,+s,+S,+f,+e,+R,+P,+d,+K,+J,+b,+H,+n,+m,-C,-l,+g,+A) ;
	if(i==294) MV = GA5_4_1(+a,+d,+r,+H,-P,-l,+f,-B,-E,+c,+K,+S,+n,-j,-G,+M,-m,+g,-J,+N,-s,-k,+C,+e,+b,+F,-L,+p,-h,+R,+D,+A) ;
	if(i==295) MV = GA5_4_1(+a,+F,-s,+k,-r,+N,-E,-B,+D,+c,+j,-p,+M,-h,-G,-L,-l,+g,+H,-m,+P,+J,+C,+d,+b,+e,+n,-R,-K,+S,+f,+A) ;
	if(i==296) MV = GA5_4_1(+a,+e,+P,-J,+s,-m,+D,-B,+f,+c,+h,+R,-L,-K,-G,+n,+N,+g,+k,-l,+r,-H,+C,+F,+b,+d,+M,-S,-j,-p,-E,+A) ;
	if(i==297) MV = GA5_4_1(+a,-H,+P,+d,+r,-l,-j,-B,+K,+G,-E,+S,-M,+f,+c,+n,+N,+C,+F,+m,-s,+e,-g,-k,+b,-J,-L,+R,+D,-p,+h,+A) ;
	if(i==298) MV = GA5_4_1(+a,+J,-s,+e,+P,-m,-K,-B,+h,+G,+f,+R,-n,+D,+c,-L,-l,+C,+d,-N,+r,+F,-g,-H,+b,+k,+M,-p,-E,+S,+j,+A) ;
	if(i==299) MV = GA5_4_1(+a,-k,+r,+F,-s,+N,-h,-B,+j,+G,+D,-p,+L,-E,+c,+M,-m,+C,+e,+l,+P,+d,-g,+J,+b,+H,+n,+S,+f,+R,+K,+A) ;
	if(i==300) MV = GA5_4_1(+a,-H,-r,+d,+L,+p,-K,-B,+g,+J,+f,-N,-s,-C,+e,-M,+m,+E,+c,-S,+n,+F,-j,+G,+b,+k,+P,-l,+D,-R,+h,+A) ;
	if(i==301) MV = GA5_4_1(+a,-H,+P,+d,+L,-R,-j,-B,+g,+k,-E,+m,-s,-C,-F,+n,+N,+f,+c,-S,+M,+e,+K,+G,+b,-J,+r,-l,+D,-p,+h,+A) ;
	if(i==302) MV = GA5_4_1(+a,-k,+s,+F,-n,-R,-j,-B,+g,-H,-E,+m,+P,-C,+d,+L,+l,-D,+c,+p,+M,+e,-h,+G,+b,-J,+r,+N,+f,-S,+K,+A) ;
	if(i==303) MV = GA5_4_1(+a,+J,+M,+e,+s,-p,-g,-B,+K,-H,-C,+N,-L,+f,+d,+r,-R,-D,+F,-l,-n,+c,-h,-k,+b,-G,-P,-S,-E,+m,+j,+A) ;
	if(i==304) MV = GA5_4_1(+a,+J,+M,+e,+P,+S,-g,-B,+h,+k,-C,-l,-n,+D,-F,+r,-R,+f,+d,-N,+L,+c,+K,-H,+b,-G,+s,-p,-E,+m,+j,+A) ;
	if(i==305) MV = GA5_4_1(+a,-k,+n,+F,-r,+S,-g,-B,+h,+J,-C,-l,+M,+D,+e,+P,+p,+E,+d,-m,+L,+c,-j,-H,+b,-G,+s,-R,+f,-N,+K,+A) ;
	if(i==306) MV = GA5_4_1(+a,+e,+s,-J,-M,-p,+f,-B,-C,+d,+K,+N,+r,-g,+H,+L,+l,+h,-G,-R,-n,-k,-D,+c,+b,+F,-P,-m,-j,-S,-E,+A) ;
	if(i==307) MV = GA5_4_1(+a,+e,+P,-J,-M,+S,+D,-B,-C,-F,+h,-l,+r,-g,-k,+n,+N,-K,-G,-R,+L,-H,+f,+c,+b,+d,+s,-m,-j,-p,-E,+A) ;
	if(i==308) MV = GA5_4_1(+a,+F,-r,+k,-n,+S,+D,-B,-C,+e,+h,-l,+P,-g,-J,-M,+m,+j,-G,+p,+L,-H,+E,+c,+b,+d,+s,+N,-K,-R,+f,+A) ;
	if(i==309) MV = GA5_4_1(+a,+d,+L,+H,+r,+p,-C,-B,+f,+e,+g,-N,-M,-K,-J,+s,+S,+j,+k,+m,+n,+G,+E,+F,+b,+c,+P,+R,-h,-l,+D,+A) ;
	if(i==310) MV = GA5_4_1(+a,+d,+L,+H,-P,-R,-C,-B,-E,-F,+g,+m,+n,-j,-k,+s,+S,-K,-J,+N,+M,+G,+f,+e,+b,+c,+r,+p,-h,-l,+D,+A) ;
	if(i==311) MV = GA5_4_1(+a,+F,-n,+k,-s,-R,-C,-B,-E,+d,+g,+m,+L,-j,+H,-P,-p,+h,-J,+l,+M,+G,-D,+e,+b,+c,+r,+S,-K,+N,+f,+A) ;
	if(i==312) MV = GA5_4_1(+a,+c,+G,+n,-L,+m,+B,+f,+D,+e,+K,+h,+J,-R,+s,+P,+p,-S,+r,+j,+H,-k,+E,+d,-F,+b,-M,+l,+N,-g,-C,+A) ;
	if(i==313) MV = GA5_4_1(+a,+c,+G,+L,+n,+m,+B,-D,+f,+e,-h,+K,+J,-R,-P,+s,+S,+p,+r,+j,+k,+H,+E,+F,+d,+b,-M,-N,+l,-g,-C,+A) ;
	if(i==314) MV = GA5_4_1(+a,+c,+G,+L,+M,+N,+B,-D,-E,-F,-h,+j,+k,-p,-r,+s,+S,-R,-P,-K,-J,+H,+f,+e,+d,+b,+n,+m,+l,-g,-C,+A) ;
	if(i==315) MV = GA5_4_1(+a,-G,+c,+n,-M,+l,+B,-K,+j,+H,+f,+E,+d,-S,+r,+P,+p,+R,-s,+D,+e,+F,+h,+J,-k,+b,+L,-m,+N,+C,-g,+A) ;
	if(i==316) MV = GA5_4_1(+a,-G,+c,+M,+n,+l,+B,-j,-K,+H,-E,+f,+d,-S,-P,+r,-R,+p,-s,+D,-F,+e,+h,+k,+J,+b,+L,-N,-m,+C,-g,+A) ;
	if(i==317) MV = GA5_4_1(+a,-G,+c,+M,+L,-N,+B,-j,+h,+k,-E,-D,+F,+p,-s,+r,-R,+S,+P,+f,+d,+e,+K,-H,+J,+b,-n,-l,-m,+C,-g,+A) ;
	if(i==318) MV = GA5_4_1(+a,+d,-H,+r,-P,-l,+B,+f,-E,+c,+K,+j,+G,+S,+n,+M,-m,-N,+s,+g,-J,-k,+C,+e,-F,+b,-L,+p,-R,-h,+D,+A) ;
	if(i==319) MV = GA5_4_1(+a,+d,-H,+P,+r,-l,+B,+E,+f,+c,-j,+K,+G,+S,-M,+n,+N,-m,+s,+g,+k,-J,+C,+F,+e,+b,-L,+R,+p,-h,+D,+A) ;
	if(i==320) MV = GA5_4_1(+a,+F,-k,+r,-s,+N,+B,-D,-E,+c,-h,+j,+G,-p,+L,+M,-m,-l,-P,+g,-J,+H,+C,+e,+d,+b,+n,+S,-R,-K,+f,+A) ;
	if(i==321) MV = GA5_4_1(+a,-J,+e,+s,-P,-m,+B,-K,+h,-G,+f,-D,+c,+R,+n,+L,+l,-N,-r,-C,+d,+F,+g,-H,-k,+b,+M,+p,+S,+E,-j,+A) ;
	if(i==322) MV = GA5_4_1(+a,-k,-F,+s,-r,-N,+B,-j,+h,-G,-E,-D,+c,+p,+M,+L,+l,+m,+P,-C,+d,+e,+g,-H,+J,+b,-n,-R,+S,+f,+K,+A) ;
	if(i==323) MV = GA5_4_1(+a,-J,+e,+P,+s,-m,+B,-h,-K,-G,+D,+f,+c,+R,-L,+n,+N,+l,-r,-C,-F,+d,+g,+k,-H,+b,+M,-S,+p,+E,-j,+A) ;
	if(i==324) MV = GA5_4_1(+a,+H,+d,+r,-L,+p,+B,-K,+g,-J,+f,+C,+e,-N,+s,+M,-m,-S,-n,-E,+c,+F,+j,+G,-k,+b,+P,+l,-R,-D,-h,+A) ;
	if(i==325) MV = GA5_4_1(+a,-J,+e,+P,+M,-S,+B,-h,+g,+k,+D,+C,+F,+l,-r,+n,+N,-R,+L,+f,+c,+d,+K,+G,-H,+b,-s,+m,+p,+E,-j,+A) ;
	if(i==326) MV = GA5_4_1(+a,-k,-F,+r,-n,-S,+B,-h,+g,-J,+D,+C,+e,+l,+P,+M,-m,-p,+L,-E,+c,+d,+j,+G,-H,+b,-s,+N,-R,+f,+K,+A) ;
	if(i==327) MV = GA5_4_1(+a,+H,+d,+L,+r,+p,+B,-g,-K,-J,-C,+f,+e,-N,-M,+s,+S,-m,-n,-E,-F,+c,+j,+k,+G,+b,+P,+R,+l,-D,-h,+A) ;
	if(i==328) MV = GA5_4_1(+a,+H,+d,+L,+P,+R,+B,-g,+j,+k,-C,+E,+F,-m,-n,+s,+S,+N,+M,+f,+e,+c,+K,+J,+G,+b,-r,-p,+l,-D,-h,+A) ;
	if(i==329) MV = GA5_4_1(+a,-k,-F,+n,-s,+R,+B,-g,+j,+H,-C,+E,+d,-m,+L,+P,+p,-l,+M,+D,+e,+c,+h,+J,+G,+b,-r,+S,+N,+f,+K,+A) ;
	if(i==330) MV = GA5_4_1(+a,+d,-H,+P,+L,-R,+B,+E,-C,-F,-j,+g,+k,+m,-s,+n,+N,+S,-M,-K,-G,-J,+f,+c,+e,+b,+r,-l,+p,-h,+D,+A) ;
	if(i==331) MV = GA5_4_1(+a,+e,+J,+s,-M,-p,+B,+f,-C,+d,+K,+g,-H,+N,+r,+L,+l,+R,+n,+h,-G,-k,-D,+c,-F,+b,-P,-m,+S,-j,-E,+A) ;
	if(i==332) MV = GA5_4_1(+a,+F,-k,+s,-n,-R,+B,+E,-C,+d,-j,+g,-H,+m,+P,+L,+l,-p,-M,+h,-G,-J,-D,+c,+e,+b,+r,+N,+S,-K,+f,+A) ;
	if(i==333) MV = GA5_4_1(+a,+e,+J,+M,+s,-p,+B,+C,+f,+d,-g,+K,-H,+N,-L,+r,-R,+l,+n,+h,+k,-G,-D,+F,+c,+b,-P,-S,-m,-j,-E,+A) ;
	if(i==334) MV = GA5_4_1(+a,+e,+J,+M,+P,+S,+B,+C,+D,-F,-g,+h,+k,-l,-n,+r,-R,+N,-L,-K,+H,-G,+f,+d,+c,+b,+s,-p,-m,-j,-E,+A) ;
	if(i==335) MV = GA5_4_1(+a,+F,-k,+n,-r,+S,+B,+C,+D,+e,-g,+h,+J,-l,+M,+P,+p,+m,-L,+j,+H,-G,+E,+d,+c,+b,+s,-R,+N,-K,+f,+A) ;
	if(i==336) MV = GA5_4_1(+a,+c,+e,+d,-F,+B,+m,+l,+N,+G,-p,+S,+J,-R,-H,+k,-K,-h,+r,+j,-s,-P,-g,+n,-L,-M,+b,+f,+D,+E,-C,+A) ;
	if(i==337) MV = GA5_4_1(+a,+c,+e,+F,+d,+B,+m,-N,+l,+G,-S,-p,+J,-R,-k,-H,+h,-K,+r,+j,+P,-s,-g,+L,+n,-M,+b,-D,+f,+E,-C,+A) ;
	if(i==338) MV = GA5_4_1(+a,+c,-F,+e,+d,+B,+N,+m,+l,+G,-S,+R,+k,-p,+J,-H,+h,-j,-P,-K,+r,-s,-g,+L,+M,+n,+b,-D,-E,+f,-C,+A) ;
	if(i==339) MV = GA5_4_1(+a,+d,+c,+e,-F,+B,-l,+p,-R,-H,+m,+N,+G,+S,+J,+k,-K,-j,+s,+g,-n,-M,-h,+r,-P,-L,+b,+f,-E,+C,+D,+A) ;
	if(i==340) MV = GA5_4_1(+a,+d,+c,+F,+e,+B,-l,+R,+p,-H,-N,+m,+G,+S,-k,+J,+j,-K,+s,+g,+M,-n,-h,+P,+r,-L,+b,+E,+f,+C,+D,+A) ;
	if(i==341) MV = GA5_4_1(+a,+F,+c,+e,+d,+B,+N,+S,-R,-k,+m,+l,+G,-p,+J,-H,+h,-j,-P,+g,-L,-M,-K,+r,-s,+n,+b,-D,-E,+C,+f,+A) ;
	if(i==342) MV = GA5_4_1(+a,+d,-F,+c,+e,+B,-R,-l,+p,-H,-N,-S,+k,+m,+G,+J,+j,-g,-M,-K,+s,-n,-h,+P,+L,+r,+b,+E,-C,+f,+D,+A) ;
	if(i==343) MV = GA5_4_1(+a,+e,+d,+c,-F,+B,-p,-m,+S,+J,-l,-R,-H,+N,+G,+k,-K,-g,+n,+h,-r,-L,-j,+s,-M,-P,+b,+f,-C,-D,-E,+A) ;
	if(i==344) MV = GA5_4_1(+a,+F,+d,+c,+e,+B,-R,+N,+S,-k,-l,+p,-H,+m,+G,+J,+j,-g,-M,+h,-P,-L,-K,+s,-n,+r,+b,+E,-C,-D,+f,+A) ;
	if(i==345) MV = GA5_4_1(+a,+e,+d,+F,+c,+B,-p,-S,-m,+J,+R,-l,-H,+N,-k,+G,+g,-K,+n,+h,+L,-r,-j,+M,+s,-P,+b,+C,+f,-D,-E,+A) ;
	if(i==346) MV = GA5_4_1(+a,+e,-F,+d,+c,+B,+S,-p,-m,+J,+R,-N,+k,-l,-H,+G,+g,-h,-L,-K,+n,-r,-j,+M,+P,+s,+b,+C,+D,+f,-E,+A) ;
	if(i==347) MV = GA5_4_1(+a,+F,+e,+d,+c,+B,+S,-R,+N,-k,-p,-m,+J,-l,-H,+G,+g,-h,-L,+j,-M,-P,-K,+n,-r,+s,+b,+C,+D,+E,+f,+A) ;
	if(i==348) MV = GA5_4_1(+a,-G,+H,-J,+k,+B,+l,+m,-N,+c,+p,+R,+d,-S,+e,+F,+f,+E,-s,+D,+r,+P,+C,-n,+M,+L,+b,+K,-j,+h,-g,+A) ;
	if(i==349) MV = GA5_4_1(+a,-G,+H,-k,-J,+B,+l,+N,+m,+c,-R,+p,+d,-S,-F,+e,-E,+f,-s,+D,-P,+r,+C,-M,-n,+L,+b,+j,+K,+h,-g,+A) ;
	if(i==350) MV = GA5_4_1(+a,-G,+k,+H,-J,+B,-N,+l,+m,+c,-R,+S,+F,+p,+d,+e,-E,-D,+P,+f,-s,+r,+C,-M,-L,-n,+b,+j,-h,+K,-g,+A) ;
	if(i==351) MV = GA5_4_1(+a,-J,-G,+H,+k,+B,-m,-p,-S,+e,+l,-N,+c,+R,+d,+F,+f,-D,-r,-C,+n,+L,+E,-s,+P,+M,+b,+K,-h,+g,-j,+A) ;
	if(i==352) MV = GA5_4_1(+a,-k,-G,+H,-J,+B,-N,+R,-S,-F,+l,+m,+c,+p,+d,+e,-E,-D,+P,-C,+M,+L,+f,-s,+r,-n,+b,+j,-h,+g,+K,+A) ;
	if(i==353) MV = GA5_4_1(+a,-J,-G,-k,+H,+B,-m,+S,-p,+e,+N,+l,+c,+R,-F,+d,+D,+f,-r,-C,-L,+n,+E,-P,-s,+M,+b,+h,+K,+g,-j,+A) ;
	if(i==354) MV = GA5_4_1(+a,+H,-J,-G,+k,+B,+p,-l,+R,+d,-m,-S,+e,-N,+c,+F,+f,+C,-n,-E,+s,+M,-D,-r,+L,+P,+b,+K,-g,+j,-h,+A) ;
	if(i==355) MV = GA5_4_1(+a,-J,+k,-G,+H,+B,-S,-m,-p,+e,+N,-R,+F,+l,+c,+d,+D,+C,+L,+f,-r,+n,+E,-P,-M,-s,+b,+h,-g,+K,-j,+A) ;
	if(i==356) MV = GA5_4_1(+a,-k,-J,-G,+H,+B,-S,-N,+R,-F,-m,-p,+e,+l,+c,+d,+D,+C,+L,-E,+P,+M,+f,-r,+n,-s,+b,+h,-g,+j,+K,+A) ;
	if(i==357) MV = GA5_4_1(+a,+H,-J,-k,-G,+B,+p,-R,-l,+d,+S,-m,+e,-N,-F,+c,-C,+f,-n,-E,-M,+s,-D,-L,-r,+P,+b,+g,+K,+j,-h,+A) ;
	if(i==358) MV = GA5_4_1(+a,+H,+k,-J,-G,+B,+R,+p,-l,+d,+S,+N,+F,-m,+e,+c,-C,+E,+M,+f,-n,+s,-D,-L,-P,-r,+b,+g,-j,+K,-h,+A) ;
	if(i==359) MV = GA5_4_1(+a,-k,+H,-J,-G,+B,+R,-S,-N,-F,+p,-l,+d,-m,+e,+c,-C,+E,+M,+D,+L,+P,+f,-n,+s,-r,+b,+g,-j,+h,+K,+A) ;

	return MV;
}

//////////////////////////////

GA5_4_1 Comp(GA5_4_1 r, int i)
{
	GA5_4_1 s;
	long Mask = 1;
	long index;

	unsigned long Matches[64] = {
		0x00000000, 0x03FFFFC0, 0x0CDE7B30, 0x0F2184F0, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 
		0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x30E1870C, 0x331E78CC, 0x3C3FFC3C, 0x3FC003FC, 
		0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x5119988A, 0x52E6674A, 0x5DC7E3BA, 0x5E381C7A, 
		0x61F81F86, 0x6207E046, 0x6D2664B6, 0x6ED99B76, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 
		0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x91266489, 0x92D99B49, 0x9DF81FB9, 0x9E07E079, 
		0xA1C7E385, 0xA2381C45, 0xAD1998B5, 0xAEE66775, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 
		0xC03FFC03, 0xC3C003C3, 0xCCE18733, 0xCF1E78F3, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 
		0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xF0DE7B0F, 0xF32184CF, 0xFC00003F, 0xFFFFFFFF
	} ;

	index = Matches[i&63];	// the &63 is to keep i inside range 0..63

	if((index & Mask) == 0) s.q     = r.q    ; else s.q     = -r.q    ;	Mask <<= 1;

	if((index & Mask) == 0) s.w     = r.w    ; else s.w     = -r.w    ;	Mask <<= 1;
	if((index & Mask) == 0) s.x     = r.x    ; else s.x     = -r.x    ;	Mask <<= 1;
	if((index & Mask) == 0) s.y     = r.y    ; else s.y     = -r.y    ;	Mask <<= 1;
	if((index & Mask) == 0) s.z     = r.z    ; else s.z     = -r.z    ;	Mask <<= 1;
	if((index & Mask) == 0) s.t     = r.t    ; else s.t     = -r.t    ;	Mask <<= 1;

	if((index & Mask) == 0) s.wx    = r.wx   ; else s.wx    = -r.wx   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wy    = r.wy   ; else s.wy    = -r.wy   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wz    = r.wz   ; else s.wz    = -r.wz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wt    = r.wt   ; else s.wt    = -r.wt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xy    = r.xy   ; else s.xy    = -r.xy   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xz    = r.xz   ; else s.xz    = -r.xz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xt    = r.xt   ; else s.xt    = -r.xt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.yz    = r.yz   ; else s.yz    = -r.yz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.yt    = r.yt   ; else s.yt    = -r.yt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.zt    = r.zt   ; else s.zt    = -r.zt   ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxy   = r.wxy  ; else s.wxy   = -r.wxy  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxz   = r.wxz  ; else s.wxz   = -r.wxz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxt   = r.wxt  ; else s.wxt   = -r.wxt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyz   = r.wyz  ; else s.wyz   = -r.wyz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyt   = r.wyt  ; else s.wyt   = -r.wyt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wzt   = r.wzt  ; else s.wzt   = -r.wzt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyz   = r.xyz  ; else s.xyz   = -r.xyz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyt   = r.xyt  ; else s.xyt   = -r.xyt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xzt   = r.xzt  ; else s.xzt   = -r.xzt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.yzt   = r.yzt  ; else s.yzt   = -r.yzt  ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxyz  = r.wxyz ; else s.wxyz  = -r.wxyz ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxyt  = r.wxyt ; else s.wxyt  = -r.wxyt ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxzt  = r.wxzt ; else s.wxzt  = -r.wxzt ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyzt  = r.wyzt ; else s.wyzt  = -r.wyzt ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyzt  = r.xyzt ; else s.xyzt  = -r.xyzt ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxyzt = r.wxyzt; else s.wxyzt = -r.wxyzt;	Mask <<= 1;

	return s;
}

//////////////////////////////////////////////////////

matrix GA5_4_1_To_Matrix(GA5_4_1 W)
{
	matrix q(4,4), x(4,4),y(4,4),z(4,4), xy(4,4),yz(4,4),xz(4,4), xyz(4,4);
	matrix t(4,4), xt(4,4),yt(4,4),zt(4,4), xyt(4,4),yzt(4,4),xzt(4,4), xyzt(4,4);
	matrix w(4,4), wx(4,4),wy(4,4),wz(4,4), wxy(4,4),wyz(4,4),wxz(4,4), wxyz(4,4);
	matrix wt(4,4), wxt(4,4),wyt(4,4),wzt(4,4), wxyt(4,4),wyzt(4,4),wxzt(4,4), wxyzt(4,4);

	matrix MV(4,4);
	matrix Zero(4,4);
/*

      q             wxyzt  
[ 1  0  0  0]  [ I  0  0  0]  
[ 0  1  0  0]  [ 0  I  0  0]
[ 0  0  1  0]  [ 0  0  I  0] 
[ 0  0  0  1]  [ 0  0  0  I] 

      w              x              y              z              t  
[ 0  0  0 -I]  [ 0  0  1  0]  [ 1  0  0  0]  [ 0  0  0  1]  [ 0  0 -1  0]
[ 0  0 -I  0]  [ 0  0  0  1]  [ 0  1  0  0]  [ 0  0 -1  0]  [ 0  0  0  1]
[ 0  I  0  0]  [ 1  0  0  0]  [ 0  0 -1  0]  [ 0 -1  0  0]  [ 1  0  0  0]
[ I  0  0  0]  [ 0  1  0  0]  [ 0  0  0 -1]  [ 1  0  0  0]  [ 0 -1  0  0]

      wx             wy             wz             wt             xy  
[ 0 -I  0  0]  [ 0  0  0  I]  [-I  0  0  0]  [ 0  I  0  0]  [ 0  0 -1  0]
[-I  0  0  0]  [ 0  0  I  0]  [ 0  I  0  0]  [-I  0  0  0]  [ 0  0  0 -1]
[ 0  0  0  I]  [ 0  I  0  0]  [ 0  0 -I  0]  [ 0  0  0  I]  [ 1  0  0  0]
[ 0  0  I  0]  [ I  0  0  0]  [ 0  0  0  I]  [ 0  0 -I  0]  [ 0  1  0  0]

      xz             xt             yz             yt             zt  
[ 0 -1  0  0]  [ 1  0  0  0]  [ 0  0  0  1]  [ 0  0 -1  0]  [ 0 -1  0  0]
[ 1  0  0  0]  [ 0 -1  0  0]  [ 0  0 -1  0]  [ 0  0  0  1]  [-1  0  0  0]
[ 0  0  0  1]  [ 0  0 -1  0]  [ 0  1  0  0]  [-1  0  0  0]  [ 0  0  0 -1]
[ 0  0 -1  0]  [ 0  0  0  1]  [-1  0  0  0]  [ 0  1  0  0]  [ 0  0 -1  0]

     wxy            wxz            wxt            wyz            wyt  
[ 0 -I  0  0]  [ 0  0  I  0]  [ 0  0  0 -I]  [ I  0  0  0]  [ 0 -I  0  0]
[-I  0  0  0]  [ 0  0  0 -I]  [ 0  0  I  0]  [ 0 -I  0  0]  [ I  0  0  0]
[ 0  0  0 -I]  [ I  0  0  0]  [ 0 -I  0  0]  [ 0  0 -I  0]  [ 0  0  0  I]
[ 0  0 -I  0]  [ 0 -I  0  0]  [ I  0  0  0]  [ 0  0  0  I]  [ 0  0 -I  0]

     wzt            xyz            xyt            xzt            yzt  
[ 0  0  I  0]  [ 0  1  0  0]  [-1  0  0  0]  [ 0  0  0 -1]  [ 0 -1  0  0]
[ 0  0  0  I]  [-1  0  0  0]  [ 0  1  0  0]  [ 0  0 -1  0]  [-1  0  0  0]
[-I  0  0  0]  [ 0  0  0  1]  [ 0  0 -1  0]  [ 0 -1  0  0]  [ 0  0  0  1]
[ 0 -I  0  0]  [ 0  0 -1  0]  [ 0  0  0  1]  [-1  0  0  0]  [ 0  0  1  0]

     wxyz           wxyt           wxzt           wyzt           xyzt  
[ 0  0  I  0]  [ 0  0  0 -I]  [ I  0  0  0]  [ 0  0 -I  0]  [ 0  0  0  1]
[ 0  0  0 -I]  [ 0  0  I  0]  [ 0  I  0  0]  [ 0  0  0 -I]  [ 0  0  1  0]
[-I  0  0  0]  [ 0  I  0  0]  [ 0  0 -I  0]  [-I  0  0  0]  [ 0 -1  0  0]
[ 0  I  0  0]  [-I  0  0  0]  [ 0  0  0 -I]  [ 0 -I  0  0]  [-1  0  0  0]

*/
	Zero = 	 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0  ;

	q = 	 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1  ;

	x = 	 0, 0, 1, 0,
		 0, 0, 0, 1,
		 1, 0, 0, 0,
		 0, 1, 0, 0  ;

	y = 	 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0,-1, 0,
		 0, 0, 0,-1  ;

	z = 	 0, 0, 0, 1,
		 0, 0,-1, 0,
		 0,-1, 0, 0,
		 1, 0, 0, 0  ;

	t = 	 0, 0,-1, 0,
		 0, 0, 0, 1,
		 1, 0, 0, 0,
		 0,-1, 0, 0  ;

	wxyzt =  I, 0, 0, 0,
		 0, I, 0, 0,
		 0, 0, I, 0,
		 0, 0, 0, I  ;

/*  // new ginac calling convention

	Zero = matrix({	{ 0, 0, 0, 0},
			{ 0, 0, 0, 0},
			{ 0, 0, 0, 0},
			{ 0, 0, 0, 0}  });

	q = matrix({	{ 1, 0, 0, 0},
			{ 0, 1, 0, 0},
			{ 0, 0, 1, 0},
			{ 0, 0, 0, 1}  });

	x = matrix({	{ 0, 0, 1, 0},
			{ 0, 0, 0, 1},
			{ 1, 0, 0, 0},
			{ 0, 1, 0, 0}  });

	y = matrix({	{ 1, 0, 0, 0},
			{ 0, 1, 0, 0},
			{ 0, 0,-1, 0},
			{ 0, 0, 0,-1}  });

	z = matrix({	{ 0, 0, 0, 1},
			{ 0, 0,-1, 0},
			{ 0,-1, 0, 0},
			{ 1, 0, 0, 0}  });

	t = matrix({	{ 0, 0,-1, 0},
			{ 0, 0, 0, 1},
			{ 1, 0, 0, 0},
			{ 0,-1, 0, 0}  });

	wxyzt = matrix({	{ I, 0, 0, 0},
				{ 0, I, 0, 0},
				{ 0, 0, I, 0},
				{ 0, 0, 0, I}  });

*/
//   calculate remaining terms. Synthesize up for xyzt

	xy = x.mul(y);
	xz = x.mul(z);
	yz = y.mul(z);
	xyz = xy.mul(z);

	xt = x.mul(t);
	yt = y.mul(t);
	zt = z.mul(t);

	xyt = xy.mul(t);
	xzt = xz.mul(t);
	yzt = yz.mul(t);
	xyzt = xyz.mul(t);

//  I assigned wxyzt to be I. Synthesize down for the remaining w terms

	wxyz = Zero.sub(wxyzt.mul(t));	// t^2 = -1
	wxyt = Zero.sub(wxyzt.mul(z));
	wxzt = wxyzt.mul(y);
	wyzt = Zero.sub(wxyzt.mul(x));

	wxy = wxyz.mul(z);
	wxz = Zero.sub(wxyz.mul(y));
	wxt = Zero.sub(wxyt.mul(y));
	wyz = Zero.sub(wyzt.mul(t));
	wyt = wxyt.mul(x);
	wzt = wyzt.mul(y);

	wx = wxy.mul(y);
	wy = Zero.sub(wxy.mul(x));
	wz = Zero.sub(wxz.mul(x));
	wt = Zero.sub(wxt.mul(x));

	w = wx.mul(x);

//	Build our matrix

	MV = Zero;
	MV = MV.add(q.mul_scalar(W.q));

	MV = MV.add(w.mul_scalar(W.w));
	MV = MV.add(x.mul_scalar(W.x));
	MV = MV.add(y.mul_scalar(W.y));
	MV = MV.add(z.mul_scalar(W.z));
	MV = MV.add(t.mul_scalar(W.t));

	MV = MV.add(wx.mul_scalar(W.wx));
	MV = MV.add(wy.mul_scalar(W.wy));
	MV = MV.add(wz.mul_scalar(W.wz));
	MV = MV.add(wt.mul_scalar(W.wt));
	MV = MV.add(xy.mul_scalar(W.xy));
	MV = MV.add(xz.mul_scalar(W.xz));
	MV = MV.add(xt.mul_scalar(W.xt));
	MV = MV.add(yz.mul_scalar(W.yz));
	MV = MV.add(yt.mul_scalar(W.yt));
	MV = MV.add(zt.mul_scalar(W.zt));

	MV = MV.add(wxy.mul_scalar(W.wxy));
	MV = MV.add(wxz.mul_scalar(W.wxz));
	MV = MV.add(wxt.mul_scalar(W.wxt));
	MV = MV.add(wyz.mul_scalar(W.wyz));
	MV = MV.add(wyt.mul_scalar(W.wyt));
	MV = MV.add(wzt.mul_scalar(W.wzt));
	MV = MV.add(xyz.mul_scalar(W.xyz));
	MV = MV.add(xyt.mul_scalar(W.xyt));
	MV = MV.add(xzt.mul_scalar(W.xzt));
	MV = MV.add(yzt.mul_scalar(W.yzt));

	MV = MV.add(wxyz.mul_scalar(W.wxyz));
	MV = MV.add(wxyt.mul_scalar(W.wxyt));
	MV = MV.add(wxzt.mul_scalar(W.wxzt));
	MV = MV.add(wyzt.mul_scalar(W.wyzt));
	MV = MV.add(xyzt.mul_scalar(W.xyzt));

	MV = MV.add(wxyzt.mul_scalar(W.wxyzt));

	return MV;

}

//////////////////////////////////////////////////////

GA5_4_1  Matrix_To_GA5_4_1(matrix W)
{
	GA5_4_1 MV;

	ex a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A;

	a = real_part( + W(0,0) + W(1,1) + W(2,2) + W(3,3) )/4;
	
	b = imag_part( - W(0,3) - W(1,2) + W(2,1) + W(3,0) )/4;
	c = real_part( + W(0,2) + W(1,3) + W(2,0) + W(3,1) )/4;
	d = real_part( + W(0,0) + W(1,1) - W(2,2) - W(3,3) )/4;
	e = real_part( + W(0,3) - W(1,2) - W(2,1) + W(3,0) )/4;
	f = real_part( - W(0,2) + W(1,3) + W(2,0) - W(3,1) )/4;
	
	g = imag_part( - W(0,1) - W(1,0) + W(2,3) + W(3,2) )/4;
	h = imag_part( + W(0,3) + W(1,2) + W(2,1) + W(3,0) )/4;
	j = imag_part( - W(0,0) + W(1,1) - W(2,2) + W(3,3) )/4;
	k = imag_part( + W(0,1) - W(1,0) + W(2,3) - W(3,2) )/4;
	l = real_part( - W(0,2) - W(1,3) + W(2,0) + W(3,1) )/4;
	m = real_part( - W(0,1) + W(1,0) + W(2,3) - W(3,2) )/4;
	n = real_part( + W(0,0) - W(1,1) - W(2,2) + W(3,3) )/4;
	p = real_part( + W(0,3) - W(1,2) + W(2,1) - W(3,0) )/4;
	r = real_part( - W(0,2) + W(1,3) - W(2,0) + W(3,1) )/4;
	s = real_part( - W(0,1) - W(1,0) - W(2,3) - W(3,2) )/4;
	
	S = imag_part( - W(0,1) - W(1,0) - W(2,3) - W(3,2) )/4;
	R = imag_part( + W(0,2) - W(1,3) + W(2,0) - W(3,1) )/4;
	P = imag_part( - W(0,3) + W(1,2) - W(2,1) + W(3,0) )/4;
	N = imag_part( + W(0,0) - W(1,1) - W(2,2) + W(3,3) )/4;
	M = imag_part( - W(0,1) + W(1,0) + W(2,3) - W(3,2) )/4;
	L = imag_part( + W(0,2) + W(1,3) - W(2,0) - W(3,1) )/4;
	K = real_part( + W(0,1) - W(1,0) + W(2,3) - W(3,2) )/4;
	J = real_part( - W(0,0) + W(1,1) - W(2,2) + W(3,3) )/4;
	H = real_part( - W(0,3) - W(1,2) - W(2,1) - W(3,0) )/4;
	G = real_part( - W(0,1) - W(1,0) + W(2,3) + W(3,2) )/4;
	
	F = imag_part( + W(0,2) - W(1,3) - W(2,0) + W(3,1) )/4;
	E = imag_part( - W(0,3) + W(1,2) + W(2,1) - W(3,0) )/4;
	D = imag_part( + W(0,0) + W(1,1) - W(2,2) - W(3,3) )/4;
	C = imag_part( - W(0,2) - W(1,3) - W(2,0) - W(3,1) )/4;
	B = real_part( + W(0,3) + W(1,2) - W(2,1) - W(3,0) )/4;
	
	A = imag_part( + W(0,0) + W(1,1) + W(2,2) + W(3,3) )/4;

	MV = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	return MV;

}


// 	int	main(void) { return 0; }

