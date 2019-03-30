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

GA5_4_1 Involution(GA5_4_1 a)	
// Corresponds to WPT parity transform: w -> -w,  x -> -x, y -> -y, z ->-z, t -> -t
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =   a.wxyz;	b.wxyt =   a.wxyt;	b.wxzt =   a.wxzt;	b.wyzt =   a.wyzt;	b.xyzt =   a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Transpose(const GA5_4_1 &a)	// verified
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

GA5_4_1 Conjugation(const GA5_4_1 &a)	
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

GA5_4_1 Hermitian(const GA5_4_1 &a)	
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

GA5_4_1 Conjugate(const GA5_4_1 &a)	
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

	det = g + G*I;
	
	return det;
}

/////////////////////////////////////////////////////////////////////

//  Adjugate left out at this time  Adjugate = Reciprocol*Determinant

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
	
	Result = MV1*MV2*MV3;	// need to divide by det
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

/////////////////////////// routines for printing equations //////////////////////////////////////////


///////////////////////////////////

void PrintSimpleBlade(GA5_4_1 a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q     ");	if(a.q == -1) printf("-q     ");	if(a.q ==  0) Count++;

	if(a.w ==  1) printf(" w     ");	if(a.w == -1) printf("-w     ");	if(a.w ==  0) Count++;
	if(a.x ==  1) printf(" x     ");	if(a.x == -1) printf("-x     ");	if(a.x ==  0) Count++;
	if(a.y ==  1) printf(" y     ");	if(a.y == -1) printf("-y     ");	if(a.y ==  0) Count++;
	if(a.z ==  1) printf(" z     ");	if(a.z == -1) printf("-z     ");	if(a.z ==  0) Count++;
	if(a.t ==  1) printf(" t     ");	if(a.t == -1) printf("-t     ");	if(a.t ==  0) Count++;

	if(a.wx ==  1) printf(" wx    ");	if(a.wx == -1) printf("-wx    ");	if(a.wx ==  0) Count++;
	if(a.wy ==  1) printf(" wy    ");	if(a.wy == -1) printf("-wy    ");	if(a.wy ==  0) Count++;
	if(a.wz ==  1) printf(" wz    ");	if(a.wz == -1) printf("-wz    ");	if(a.wz ==  0) Count++;
	if(a.wt ==  1) printf(" wt    ");	if(a.wt == -1) printf("-wt    ");	if(a.wt ==  0) Count++;
	if(a.xy ==  1) printf(" xy    ");	if(a.xy == -1) printf("-xy    ");	if(a.xy ==  0) Count++;
	if(a.xz ==  1) printf(" xz    ");	if(a.xz == -1) printf("-xz    ");	if(a.xz ==  0) Count++;
	if(a.xt ==  1) printf(" xt    ");	if(a.xt == -1) printf("-xt    ");	if(a.xt ==  0) Count++;
	if(a.yz ==  1) printf(" yz    ");	if(a.yz == -1) printf("-yz    ");	if(a.yz ==  0) Count++;
	if(a.yt ==  1) printf(" yt    ");	if(a.yt == -1) printf("-yt    ");	if(a.yt ==  0) Count++;
	if(a.zt ==  1) printf(" zt    ");	if(a.zt == -1) printf("-zt    ");	if(a.zt ==  0) Count++;

	if(a.wxy ==  1) printf(" wxy   ");	if(a.wxy == -1) printf("-wxy   ");	if(a.wxy ==  0) Count++;
	if(a.wxz ==  1) printf(" wxz   ");	if(a.wxz == -1) printf("-wxz   ");	if(a.wxz ==  0) Count++;
	if(a.wxt ==  1) printf(" wxt   ");	if(a.wxt == -1) printf("-wxt   ");	if(a.wxt ==  0) Count++;
	if(a.wyz ==  1) printf(" wyz   ");	if(a.wyz == -1) printf("-wyz   ");	if(a.wyz ==  0) Count++;
	if(a.wyt ==  1) printf(" wyt   ");	if(a.wyt == -1) printf("-wyt   ");	if(a.wyt ==  0) Count++;
	if(a.wzt ==  1) printf(" wzt   ");	if(a.wzt == -1) printf("-wzt   ");	if(a.wzt ==  0) Count++;
	if(a.xyz ==  1) printf(" xyz   ");	if(a.xyz == -1) printf("-xyz   ");	if(a.xyz ==  0) Count++;
	if(a.xyt ==  1) printf(" xyt   ");	if(a.xyt == -1) printf("-xyt   ");	if(a.xyt ==  0) Count++;
	if(a.xzt ==  1) printf(" xzt   ");	if(a.xzt == -1) printf("-xzt   ");	if(a.xzt ==  0) Count++;
	if(a.yzt ==  1) printf(" yzt   ");	if(a.yzt == -1) printf("-yzt   ");	if(a.yzt ==  0) Count++;


	if(a.wxyz ==  1) printf(" wxyz  ");	if(a.wxyz == -1) printf("-xwyz  ");	if(a.wxyz ==  0) Count++;
	if(a.wxyt ==  1) printf(" wxyt  ");	if(a.wxyt == -1) printf("-xwyt  ");	if(a.wxyt ==  0) Count++;
	if(a.wxzt ==  1) printf(" wxzt  ");	if(a.wxzt == -1) printf("-wxzt  ");	if(a.wxzt ==  0) Count++;
	if(a.wyzt ==  1) printf(" wyzt  ");	if(a.wyzt == -1) printf("-wyzt  ");	if(a.wyzt ==  0) Count++;
	if(a.xyzt ==  1) printf(" xyzt  ");	if(a.xyzt == -1) printf("-xyzt  ");	if(a.xyzt ==  0) Count++;

	if(a.wxyzt ==  1) printf(" wxyzt ");	if(a.wxyzt == -1) printf("-wxyzt ");	if(a.wxyzt ==  0) Count++;

	if(Count == 32) printf(" 0     ");
	if(Count < 31)  printf("error  ");

}


///////////////////////////////////

int	GetWeight(GA5_4_1 r) // assume only one non-zero component
{
	int weight = 0;		// default

	if (r.q ==  1) weight =  1;	if (r.q == -1) weight = -1;

	if (r.w ==  1) weight =  1;	if (r.w == -1) weight = -1;
	if (r.x ==  1) weight =  1;	if (r.x == -1) weight = -1;
	if (r.y ==  1) weight =  1;	if (r.y == -1) weight = -1;
	if (r.z ==  1) weight =  1;	if (r.z == -1) weight = -1;
	if (r.t ==  1) weight =  1;	if (r.t == -1) weight = -1;

	if (r.wx ==  1) weight =  1;	if (r.wx == -1) weight = -1;
	if (r.wy ==  1) weight =  1;	if (r.wy == -1) weight = -1;
	if (r.wz ==  1) weight =  1;	if (r.wz == -1) weight = -1;
	if (r.wt ==  1) weight =  1;	if (r.wt == -1) weight = -1;
	if (r.xy ==  1) weight =  1;	if (r.xy == -1) weight = -1;
	if (r.xz ==  1) weight =  1;	if (r.xz == -1) weight = -1;
	if (r.xt ==  1) weight =  1;	if (r.xt == -1) weight = -1;
	if (r.yz ==  1) weight =  1;	if (r.yz == -1) weight = -1;
	if (r.yt ==  1) weight =  1;	if (r.yt == -1) weight = -1;
	if (r.zt ==  1) weight =  1;	if (r.zt == -1) weight = -1;

	if (r.wxy ==  1) weight =  1;	if (r.wxy == -1) weight = -1;
	if (r.wxz ==  1) weight =  1;	if (r.wxz == -1) weight = -1;
	if (r.wxt ==  1) weight =  1;	if (r.wxt == -1) weight = -1;
	if (r.wyz ==  1) weight =  1;	if (r.wyz == -1) weight = -1;
	if (r.wyt ==  1) weight =  1;	if (r.wyt == -1) weight = -1;
	if (r.wzt ==  1) weight =  1;	if (r.wzt == -1) weight = -1;
	if (r.xyz ==  1) weight =  1;	if (r.xyz == -1) weight = -1;
	if (r.xyt ==  1) weight =  1;	if (r.xyt == -1) weight = -1;
	if (r.xzt ==  1) weight =  1;	if (r.xzt == -1) weight = -1;
	if (r.yzt ==  1) weight =  1;	if (r.yzt == -1) weight = -1;

	if (r.wxyz ==  1) weight =  1;	if (r.wxyz == -1) weight = -1;
	if (r.wxyt ==  1) weight =  1;	if (r.wxyt == -1) weight = -1;
	if (r.wxzt ==  1) weight =  1;	if (r.wxzt == -1) weight = -1;
	if (r.wyzt ==  1) weight =  1;	if (r.wyzt == -1) weight = -1;
	if (r.xyzt ==  1) weight =  1;	if (r.xyzt == -1) weight = -1;

	if (r.wxyzt ==  1) weight =  1;	if (r.wxyzt == -1) weight = -1;

	return weight;
}

/////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////

int	Higher(int a, int b)	// return higher of a or b
{
	int c;
	
	c = a;
	if (b>a) c = b;

	return c;

}

/////////////////////////////////////////////////////////////////////


GA5_4_1 AuntieWedge(const GA5_4_1 &a, const GA5_4_1 &b) {

	GA5_4_1 c;


	return c;

}


/////////////////////////////////////////////////////////////////////

int main(void)
{


	symbol A("A");
	symbol B("B"), C("C"), D("D"), E("E");
	symbol F("F"), G("G"), H("H"), J("J"), K("K"), L("L");
	symbol M("M"), N("N"), P("P"), R("R");
	symbol S("S");


	symbol a("a");
	symbol b("b"), c("c"), d("d"), e("e");
	symbol f("f"), g("g"), h("h"), j("j"), k("k"), l("l");
	symbol m("m"), n("n"), p("p"), r("r");
	symbol s("s");

	symbol a_q("a.q");
	symbol a_w("a.w"), a_x("a.x"), a_y("a.y"), a_z("a.z"), a_t("a.t");
	symbol a_wx("a.wx"), a_wy("a.wy"), a_wz("a.wz"), a_wt("a.wt"), a_xy("a.xy");
	symbol a_xz("a.xz"), a_xt("a.xt"), a_yz("a.yz"), a_yt("a.yt"), a_zt("a.zt");
	symbol a_wxy("a.wxy"), a_wxz("a.wxz"), a_wxt("a.wxt"), a_wyz("a.wyz"), a_wyt("a.wyt");
	symbol a_wzt("a.wzt"), a_xyz("a.xyz"), a_xyt("a.xyt"), a_xzt("a.xzt"), a_yzt("a.yzt");
	symbol a_wxyz("a.wxyz"), a_wxyt("a.wxyt"), a_wxzt("a.wxzt"), a_wyzt("a.wyzt"), a_xyzt("a.xyzt");
	symbol a_wxyzt("a.wxyzt");

	symbol b_q("b.q");
	symbol b_w("b.w"), b_x("b.x"), b_y("b.y"), b_z("b.z"), b_t("b.t");
	symbol b_wx("b.wx"), b_wy("b.wy"), b_wz("b.wz"), b_wt("b.wt"), b_xy("b.xy");
	symbol b_xz("b.xz"), b_xt("b.xt"), b_yz("b.yz"), b_yt("b.yt"), b_zt("b.zt");
	symbol b_wxy("b.wxy"), b_wxz("b.wxz"), b_wxt("b.wxt"), b_wyz("b.wyz"), b_wyt("b.wyt");
	symbol b_wzt("b.wzt"), b_xyz("b.xyz"), b_xyt("b.xyt"), b_xzt("b.xzt"), b_yzt("b.yzt");
	symbol b_wxyz("b.wxyz"), b_wxyt("b.wxyt"), b_wxzt("b.wxzt"), b_wyzt("b.wyzt"), b_xyzt("b.xyzt");
	symbol b_wxyzt("b.wxyzt");


	int i;	// loop counters

	GA5_4_1 q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA5_4_1 wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA5_4_1 wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;

	// initialize the convenience variables above.
	
	q.q = 1; 
	w.w = 1; 
	x.x = 1;
	y.y = 1;
	z.z = 1;
	t.t = 1;
	wx.wx = 1;
	wy.wy = 1;
	wz.wz = 1;
	wt.wt = 1;
	xy.xy = 1;
	xz.xz = 1;
	xt.xt = 1;
	yz.yz = 1;
	yt.yt = 1;
	zt.zt = 1;
	wxy.wxy = 1;
	wxz.wxz = 1;
	wxt.wxt = 1;
	wyz.wyz = 1;
	wyt.wyt = 1;
	wzt.wzt = 1;
	xyz.xyz = 1;
	xyt.xyt = 1;
	xzt.xzt = 1;
	yzt.yzt = 1;
	wxyz.wxyz = 1;
	wxyt.wxyt = 1;
	wxzt.wxzt = 1;
	wyzt.wyzt = 1;
	xyzt.xyzt = 1;
	wxyzt.wxyzt = 1;

	GA5_4_1 Blade[32];

	i = 0;
	Blade[i++] = q; 
	Blade[i++] = w; 
	Blade[i++] = x;
	Blade[i++] = y;
	Blade[i++] = z;
	Blade[i++] = t;
	Blade[i++] = wx;
	Blade[i++] = wy;
	Blade[i++] = wz;
	Blade[i++] = wt;
	Blade[i++] = xy;
	Blade[i++] = xz;
	Blade[i++] = xt;
	Blade[i++] = yz;
	Blade[i++] = yt;
	Blade[i++] = zt;
	Blade[i++] = wxy;
	Blade[i++] = wxz;
	Blade[i++] = wxt;
	Blade[i++] = wyz;
	Blade[i++] = wyt;
	Blade[i++] = wzt;
	Blade[i++] = xyz;
	Blade[i++] = xyt;
	Blade[i++] = xzt;
	Blade[i++] = yzt;
	Blade[i++] = wxyz;
	Blade[i++] = wxyt;
	Blade[i++] = wxzt;
	Blade[i++] = wyzt;
	Blade[i++] = xyzt;
	Blade[i++] = wxyzt;

/*	cout << q << "\n";
	cout << w << "\n";
	cout << x << "\n";
	cout << y << "\n";
	cout << z << "\n";
	cout << t << "\n";
	cout << wx << "\n";
	cout << wy << "\n";
	cout << wz << "\n";
	cout << wt << "\n";
	cout << xy << "\n";
	cout << xz << "\n";
	cout << xt << "\n";
	cout << yz << "\n";
	cout << yt << "\n";
	cout << zt << "\n";
	cout << wxy << "\n";
	cout << wxz << "\n";
	cout << wxt << "\n";
	cout << wyz << "\n";
	cout << wyt << "\n";
	cout << wzt << "\n";
	cout << xyz << "\n";
	cout << xyt << "\n";
	cout << xzt << "\n";
	cout << yzt << "\n";
	cout << wxyz << "\n";
	cout << wxyt << "\n";
	cout << wxzt << "\n";
	cout << wyzt << "\n";
	cout << xyzt << "\n";
	cout << wxyzt << "\n";
*/

	GA5_4_1 MV1, MV2, MV3, MV4, MV5;

//	for (i=0;i<32;i++) cout << Blade[i] << "\n";
//	cout << "\n\n";

	// verify OverBar and UnderBar definitions

	for (i=0;i<32;i++) {
		MV1 = Blade[i];
		MV2 = OverBar(MV1); //  MV1^OverBar(MV1) = pseudovector blade wxyzt
		MV3 = MV1^MV2;	// 
		if (MV3 == wxyzt); else {
			cout << MV1 << "\n";  cout << MV2 << "\n";  cout << MV3 << "\n\n";
		}
	}

	for (i=0;i<32;i++) {
		MV1 = Blade[i];
		MV2 = UnderBar(MV1); //  MV1^OverBar(MV1) = pseudovector blade wxyzt
		MV3 = MV2^MV1;	// 
		if (MV3 == wxyzt); else {
			cout << MV1 << "\n";  cout << MV2 << "\n";  cout << MV3 << "\n\n";
		}
	}

//	MV1 = wxy;   MV2 = zt;  MV3 = MV1^MV2; cout << MV1 << "\n";  cout << MV2 << "\n";  cout << MV3 << "\n\n";

//	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
//	cout << "MV1       = " << MV1 << "\n";

//	printout Dual
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	MV2 = Dual(MV1);
	cout << "Dual(MV1)  = " << MV2 << "\n";

/*	MV1 = GA5_4_1(a_q, a_w, a_x, a_y, a_z, a_t, 
			a_wx,  a_wy,  a_wz,  a_wt,  a_xy,  a_xz,  a_xt,  a_yz,  a_yt,  a_zt,
			a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
			a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt,  a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w, b_x, b_y, b_z, b_t, 
			b_wx,  b_wy,  b_wz,  b_wt,  b_xy,  b_xz,  b_xt,  b_yz,  b_yt,  b_zt,
			b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
			b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt,  b_wxyzt);

	MV3 = Product(MV1,MV2) - MV1*MV2;

	cout << "delta products = " << MV3 << " good!\n";
*/


//	printout DorstDual
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	MV2 = DorstDual(MV1);
	cout << "DorstDual(MV1)  = " << MV2 << "\n";

//	printout DorstUnDual
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	MV2 = DorstUnDual(MV1);
	cout << "DorstUnDual(MV1)  = " << MV2 << "\n";

//	printout DorstDual(DorstUnDual(MV1))
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	MV2 = DorstDual(DorstUnDual(MV1));
	cout << "DorstDual(DorstUnDual(MV1))  = " << MV2 << "\n";

	MV1 = GA5_4_1(a_q, a_w, a_x, a_y, a_z, a_t, 
			a_wx,  a_wy,  a_wz,  a_wt,  a_xy,  a_xz,  a_xt,  a_yz,  a_yt,  a_zt,
			a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
			a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt,  a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w, b_x, b_y, b_z, b_t, 
			b_wx,  b_wy,  b_wz,  b_wt,  b_xy,  b_xz,  b_xt,  b_yz,  b_yt,  b_zt,
			b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
			b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt,  b_wxyzt);

//	MV3 = AntiWedge(MV1,MV2);
//	cout << "AntiWedge = " << MV3 << "\n\n";

	MV3 = Regressive(MV1,MV2) - RegressiveViaFormula(MV1,MV2);
	cout << "Regressive  - RegressiveViaFormula(MV1,MV2) = " << MV3 << "\n\n";

////////////////////////////////////

/*	
	char Field[33][10];
	int NumTerms = 32;
	int EqnBasis, ii, jj, kk, weight;

	for (ii=0;ii<33;ii++) for (jj=0;jj<10;jj++) Field[ii][jj] = 0;	// organized by grade

	strcpy(Field[ 0],"q");

	strcpy(Field[ 1],"w");
	strcpy(Field[ 2],"x");
	strcpy(Field[ 3],"y");
	strcpy(Field[ 4],"z");
	strcpy(Field[ 5],"t");

	strcpy(Field[ 6],"wx");
	strcpy(Field[ 7],"wy");
	strcpy(Field[ 8],"wz");
	strcpy(Field[ 9],"wt");
	strcpy(Field[10],"xy");
	strcpy(Field[11],"xz");
	strcpy(Field[12],"xt");
	strcpy(Field[13],"yz");
	strcpy(Field[14],"yt");
	strcpy(Field[15],"zt");

	strcpy(Field[16],"wxy");
	strcpy(Field[17],"wxz");
	strcpy(Field[18],"wxt");
	strcpy(Field[19],"wyz");
	strcpy(Field[20],"wyt");
	strcpy(Field[21],"wzt");
	strcpy(Field[22],"xyz");
	strcpy(Field[23],"xyt");
	strcpy(Field[24],"xzt");
	strcpy(Field[25],"yzt");

	strcpy(Field[26],"wxyz");
	strcpy(Field[27],"wxyt");
	strcpy(Field[28],"wxzt");
	strcpy(Field[29],"wyzt");
	strcpy(Field[30],"xyzt");

	strcpy(Field[31],"wxyzt");

	strcpy(Field[32],"0    ");

/////////////////////////////////


	int Rank[32]  = {0, 1,1,1,1,1, 2,2,2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3,3,3,  4,4,4,4,4, 5} ;

	long Index[32] = {0,  1,2,4,8,16,  3,5,9,17,6,10,18,12,20,24, 7,11,19,13,21,25,14,22,26,28,  15,23,27,29,30,  31};
	// these are the blade bitmap, w=1, x=2,y=4,z=8,t=16


	int CountSetBits[32] = {0, 1,1,2,1,2, 2,3,1,2,2,3,2,3,3,4, 1,2,2,3,2,3,3,4,2,3, 3,4,3,4,4, 5};
	// number of set bits per the binary index  0 1 2 3 4 etc

/////////////////////////////////


// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format

	printf("Product\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					MV1 = Product(Blade[ii],Blade[jj]);
					weight = GetWeight(MV1);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


	printf("AntiWedge\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (31^( 31^Index[ii]) ^ (31^Index[jj]) )   ) {	// regressive and dual
					MV1 = AntiWedge(Blade[ii],Blade[jj]);
					weight = GetWeight(MV1);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");

	printf("Regressive\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (31^( 31^Index[ii]) ^ (31^Index[jj]) )   ) {	// regressive and dual
					MV1 = Regressive(Blade[ii],Blade[jj]);
					weight = GetWeight(MV1);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");

	printf("Expander equation set\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[ii]^Index[jj]] > Higher(Rank[ii],Rank[jj])  ) {
						MV1 = Blade[ii]*Blade[jj];
						weight = GetWeight(MV1);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");

	printf("Conserver equation set\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[ii]^Index[jj]] == Higher(Rank[ii],Rank[jj])  ) {
						MV1 = Blade[ii]*Blade[jj];
						weight = GetWeight(MV1);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");

	printf("Shrinker equation set\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[ii]^Index[jj]] < Higher(Rank[ii],Rank[jj])  ) {
						MV1 = Blade[ii]*Blade[jj];
						weight = GetWeight(MV1);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


	printf("Symmetric Product\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					MV1 = SymmetricProductViaFormula(Blade[ii],Blade[jj]);
					weight = GetWeight(MV1);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


	printf("AntiSymmetric Product\n");
	for (kk=0;kk<NumTerms;kk++){	// equation number
		EqnBasis =  kk;
		printf("c.%s = ", Field[EqnBasis]);
		for (ii=0;ii<NumTerms;ii++) {
			for (jj=0;jj<NumTerms;jj++) {  
				if(Index[kk] == (Index[ii]^Index[jj])) {	// potential term for this equation - Clifford product as done
					MV1 = AntiSymmetricProductViaFormula(Blade[ii],Blade[jj]);
					weight = GetWeight(MV1);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[ii],Field[jj]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");

*/
	MV1 = GA5_4_1(a_q, a_w, a_x, a_y, a_z, a_t, 
			a_wx,  a_wy,  a_wz,  a_wt,  a_xy,  a_xz,  a_xt,  a_yz,  a_yt,  a_zt,
			a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
			a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt,  a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w, b_x, b_y, b_z, b_t, 
			b_wx,  b_wy,  b_wz,  b_wt,  b_xy,  b_xz,  b_xt,  b_yz,  b_yt,  b_zt,
			b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
			b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt,  b_wxyzt);

//	MV3 = AntiWedge(MV1,MV2);
//	cout << "AntiWedge = " << MV3 << "\n\n";

//	MV3 = AntiWedge(MV1,MV2) - AuntieWedge(MV1,MV2);
//	cout << "AntiWedge(MV1,MV2)  - AuntieWedge(MV1,MV2) = " << MV3 << "\n\n";

//	MV3 = Regressive(MV1,MV2) - AuntieWedge(MV1,MV2);
//	cout << "Regressive(MV1,MV2)  - AuntieWedge(MV1,MV2) = " << MV3 << "\n\n";

	MV3 = Symmetric(MV1,MV2) - SymmetricProductViaFormula(MV1,MV2);
	cout << "Symmetric(MV1,MV2) - SymmetricProductViaFormula(MV1,MV2) = " << MV3 << "\n\n";

	MV3 = AntiSymmetric(MV1,MV2) - AntiSymmetricProductViaFormula(MV1,MV2);
	cout << "AntiSymmetric(MV1,MV2) - AntiSymmetricProductViaFormula(MV1,MV2) = " << MV3 << "\n\n";





////////////////////////////////////

	return 0;
}

/*
GA5_4_1 Parity(const GA5_4_1 &a)	// Corresponds to parity transform: w -> -w, x -> -x, y -> -y, z ->-z, t -> -t
GA5_4_1 Dual(const GA5_4_1 &a)	// multiply by pseudoscalar. LeftDual = RightDual = Dual
GA5_4_1 ComplexConjugate(const GA5_4_1 &a)	// negate every term including w
GA5_4_1 Transpose(const GA5_4_1 &a)	
GA5_4_1 Hermitian(const GA5_4_1 &a)	
GA5_4_1 Reverse(const GA5_4_1 &a)

*/





