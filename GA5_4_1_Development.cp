#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

///////////////////////////////////

	matrix zero(4,4);
	matrix Q(4,4);	// scalar
	
	matrix X(4,4);
	matrix Y(4,4);
	matrix Z(4,4);
	matrix W(4,4);
	matrix T(4,4);

	matrix  WX(4,4);
	matrix  WY(4,4);
	matrix  WZ(4,4);
	matrix  WT(4,4);
	matrix  XY(4,4);
	matrix  XZ(4,4);
	matrix  XT(4,4);
	matrix  YZ(4,4);
	matrix  YT(4,4);
	matrix  ZT(4,4);

	matrix  WXY(4,4);
	matrix  WXZ(4,4);
	matrix  WXT(4,4);
	matrix  WYZ(4,4);
	matrix  WYT(4,4);
	matrix  WZT(4,4);
	matrix  XYZ(4,4);
	matrix  XYT(4,4);
	matrix  XZT(4,4);
	matrix  YZT(4,4);

	matrix  WXYZ(4,4);
	matrix  WXYT(4,4);
	matrix  WXZT(4,4);
	matrix  WYZT(4,4);
	matrix  XYZT(4,4);

	matrix  WXYZT(4,4);

///////////////////////////////////


//********************************* GA5_4_1 ***************************************

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
		
			w = w1;
			x = x1;
			y = y1;
			z = z1;
			t = t1;
		
			wx = wx1;
			wy = wy1;
			wz = wz1;
			wt = wt1;
			xy = xy1;
			xz = xz1;
			xt = xt1;
			yz = yz1;
			yt = yt1;
			zt = zt1;
		
			wxy = wxy1;
			wxz = wxz1;
			wxt = wxt1;
			wyz = wyz1;
			wyt = wyt1;
			wzt = wzt1;
			xyz = xyz1;
			xyt = xyt1;
			xzt = xzt1;
			yzt = yzt1;
		
			wxyz = wxyz1;
			wxyt = wxyt1;
			wxzt = wxzt1;
			wyzt = wyzt1;
			xyzt = xyzt1;
		
			wxyzt = wxyzt1;
		}
	};


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




GA5_4_1 operator*(const GA5_4_1 &a, const GA5_4_1 &b) {
	GA5_4_1 c;


c.q    =  + a.q   *b.q    + a.w   *b.w    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    - a.t   *b.t    - a.wx  *b.wx   - a.wy  *b.wy   - a.wz  *b.wz   + a.wt  *b.wt   - a.xy  *b.xy   - a.xz  *b.xz   + a.xt  *b.xt   - a.yz  *b.yz   + a.yt  *b.yt   + a.zt  *b.zt   - a.wxy *b.wxy  - a.wxz *b.wxz  + a.wxt *b.wxt  - a.wyz *b.wyz  + a.wyt *b.wyt  + a.wzt *b.wzt  - a.xyz *b.xyz  + a.xyt *b.xyt  + a.xzt *b.xzt  + a.yzt *b.yzt  + a.wxyz*b.wxyz - a.wxyt*b.wxyt - a.wxzt*b.wxzt - a.wyzt*b.wyzt - a.xyzt*b.xyzt - a.wxyzt*b.wxyzt ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    - a.x   *b.wx   - a.y   *b.wy   - a.z   *b.wz   + a.t   *b.wt   + a.wx  *b.x    + a.wy  *b.y    + a.wz  *b.z    - a.wt  *b.t    - a.xy  *b.wxy  - a.xz  *b.wxz  + a.xt  *b.wxt  - a.yz  *b.wyz  + a.yt  *b.wyt  + a.zt  *b.wzt  - a.wxy *b.xy   - a.wxz *b.xz   + a.wxt *b.xt   - a.wyz *b.yz   + a.wyt *b.yt   + a.wzt *b.zt   + a.xyz *b.wxyz - a.xyt *b.wxyt - a.xzt *b.wxzt - a.yzt *b.wyzt - a.wxyz*b.xyz  + a.wxyt*b.xyt  + a.wxzt*b.xzt  + a.wyzt*b.yzt  - a.xyzt*b.wxyzt - a.wxyzt*b.xyzt ; 
c.x    =  + a.q   *b.x    + a.w   *b.wx   + a.x   *b.q    - a.y   *b.xy   - a.z   *b.xz   + a.t   *b.xt   - a.wx  *b.w    + a.wy  *b.wxy  + a.wz  *b.wxz  - a.wt  *b.wxt  + a.xy  *b.y    + a.xz  *b.z    - a.xt  *b.t    - a.yz  *b.xyz  + a.yt  *b.xyt  + a.zt  *b.xzt  + a.wxy *b.wy   + a.wxz *b.wz   - a.wxt *b.wt   - a.wyz *b.wxyz + a.wyt *b.wxyt + a.wzt *b.wxzt - a.xyz *b.yz   + a.xyt *b.yt   + a.xzt *b.zt   - a.yzt *b.xyzt + a.wxyz*b.wyz  - a.wxyt*b.wyt  - a.wxzt*b.wzt  + a.wyzt*b.wxyzt + a.xyzt*b.yzt  + a.wxyzt*b.wyzt ; 
c.y    =  + a.q   *b.y    + a.w   *b.wy   + a.x   *b.xy   + a.y   *b.q    - a.z   *b.yz   + a.t   *b.yt   - a.wx  *b.wxy  - a.wy  *b.w    + a.wz  *b.wyz  - a.wt  *b.wyt  - a.xy  *b.x    + a.xz  *b.xyz  - a.xt  *b.xyt  + a.yz  *b.z    - a.yt  *b.t    + a.zt  *b.yzt  - a.wxy *b.wx   + a.wxz *b.wxyz - a.wxt *b.wxyt + a.wyz *b.wz   - a.wyt *b.wt   + a.wzt *b.wyzt + a.xyz *b.xz   - a.xyt *b.xt   + a.xzt *b.xyzt + a.yzt *b.zt   - a.wxyz*b.wxz  + a.wxyt*b.wxt  - a.wxzt*b.wxyzt - a.wyzt*b.wzt  - a.xyzt*b.xzt  - a.wxyzt*b.wxzt ; 
c.z    =  + a.q   *b.z    + a.w   *b.wz   + a.x   *b.xz   + a.y   *b.yz   + a.z   *b.q    + a.t   *b.zt   - a.wx  *b.wxz  - a.wy  *b.wyz  - a.wz  *b.w    - a.wt  *b.wzt  - a.xy  *b.xyz  - a.xz  *b.x    - a.xt  *b.xzt  - a.yz  *b.y    - a.yt  *b.yzt  - a.zt  *b.t    - a.wxy *b.wxyz - a.wxz *b.wx   - a.wxt *b.wxzt - a.wyz *b.wy   - a.wyt *b.wyzt - a.wzt *b.wt   - a.xyz *b.xy   - a.xyt *b.xyzt - a.xzt *b.xt   - a.yzt *b.yt   + a.wxyz*b.wxy  + a.wxyt*b.wxyzt + a.wxzt*b.wxt  + a.wyzt*b.wyt  + a.xyzt*b.xyt  + a.wxyzt*b.wxyt ; 
c.t    =  + a.q   *b.t    + a.w   *b.wt   + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   + a.t   *b.q    - a.wx  *b.wxt  - a.wy  *b.wyt  - a.wz  *b.wzt  - a.wt  *b.w    - a.xy  *b.xyt  - a.xz  *b.xzt  - a.xt  *b.x    - a.yz  *b.yzt  - a.yt  *b.y    - a.zt  *b.z    - a.wxy *b.wxyt - a.wxz *b.wxzt - a.wxt *b.wx   - a.wyz *b.wyzt - a.wyt *b.wy   - a.wzt *b.wz   - a.xyz *b.xyzt - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   + a.wxyz*b.wxyzt + a.wxyt*b.wxy  + a.wxzt*b.wxz  + a.wyzt*b.wyz  + a.xyzt*b.xyz  + a.wxyzt*b.wxyz ; 
c.wx   =  + a.q   *b.wx   + a.w   *b.x    - a.x   *b.w    + a.y   *b.wxy  + a.z   *b.wxz  - a.t   *b.wxt  + a.wx  *b.q    - a.wy  *b.xy   - a.wz  *b.xz   + a.wt  *b.xt   + a.xy  *b.wy   + a.xz  *b.wz   - a.xt  *b.wt   - a.yz  *b.wxyz + a.yt  *b.wxyt + a.zt  *b.wxzt + a.wxy *b.y    + a.wxz *b.z    - a.wxt *b.t    - a.wyz *b.xyz  + a.wyt *b.xyt  + a.wzt *b.xzt  + a.xyz *b.wyz  - a.xyt *b.wyt  - a.xzt *b.wzt  + a.yzt *b.wxyzt - a.wxyz*b.yz   + a.wxyt*b.yt   + a.wxzt*b.zt   - a.wyzt*b.xyzt + a.xyzt*b.wyzt + a.wxyzt*b.yzt  ; 
c.wy   =  + a.q   *b.wy   + a.w   *b.y    - a.x   *b.wxy  - a.y   *b.w    + a.z   *b.wyz  - a.t   *b.wyt  + a.wx  *b.xy   + a.wy  *b.q    - a.wz  *b.yz   + a.wt  *b.yt   - a.xy  *b.wx   + a.xz  *b.wxyz - a.xt  *b.wxyt + a.yz  *b.wz   - a.yt  *b.wt   + a.zt  *b.wyzt - a.wxy *b.x    + a.wxz *b.xyz  - a.wxt *b.xyt  + a.wyz *b.z    - a.wyt *b.t    + a.wzt *b.yzt  - a.xyz *b.wxz  + a.xyt *b.wxt  - a.xzt *b.wxyzt - a.yzt *b.wzt  + a.wxyz*b.xz   - a.wxyt*b.xt   + a.wxzt*b.xyzt + a.wyzt*b.zt   - a.xyzt*b.wxzt - a.wxyzt*b.xzt  ; 
c.wz   =  + a.q   *b.wz   + a.w   *b.z    - a.x   *b.wxz  - a.y   *b.wyz  - a.z   *b.w    - a.t   *b.wzt  + a.wx  *b.xz   + a.wy  *b.yz   + a.wz  *b.q    + a.wt  *b.zt   - a.xy  *b.wxyz - a.xz  *b.wx   - a.xt  *b.wxzt - a.yz  *b.wy   - a.yt  *b.wyzt - a.zt  *b.wt   - a.wxy *b.xyz  - a.wxz *b.x    - a.wxt *b.xzt  - a.wyz *b.y    - a.wyt *b.yzt  - a.wzt *b.t    + a.xyz *b.wxy  + a.xyt *b.wxyzt + a.xzt *b.wxt  + a.yzt *b.wyt  - a.wxyz*b.xy   - a.wxyt*b.xyzt - a.wxzt*b.xt   - a.wyzt*b.yt   + a.xyzt*b.wxyt + a.wxyzt*b.xyt  ; 
c.wt   =  + a.q   *b.wt   + a.w   *b.t    - a.x   *b.wxt  - a.y   *b.wyt  - a.z   *b.wzt  - a.t   *b.w    + a.wx  *b.xt   + a.wy  *b.yt   + a.wz  *b.zt   + a.wt  *b.q    - a.xy  *b.wxyt - a.xz  *b.wxzt - a.xt  *b.wx   - a.yz  *b.wyzt - a.yt  *b.wy   - a.zt  *b.wz   - a.wxy *b.xyt  - a.wxz *b.xzt  - a.wxt *b.x    - a.wyz *b.yzt  - a.wyt *b.y    - a.wzt *b.z    + a.xyz *b.wxyzt + a.xyt *b.wxy  + a.xzt *b.wxz  + a.yzt *b.wyz  - a.wxyz*b.xyzt - a.wxyt*b.xy   - a.wxzt*b.xz   - a.wyzt*b.yz   + a.xyzt*b.wxyz + a.wxyzt*b.xyz  ; 
c.xy   =  + a.q   *b.xy   + a.w   *b.wxy  + a.x   *b.y    - a.y   *b.x    + a.z   *b.xyz  - a.t   *b.xyt  - a.wx  *b.wy   + a.wy  *b.wx   - a.wz  *b.wxyz + a.wt  *b.wxyt + a.xy  *b.q    - a.xz  *b.yz   + a.xt  *b.yt   + a.yz  *b.xz   - a.yt  *b.xt   + a.zt  *b.xyzt + a.wxy *b.w    - a.wxz *b.wyz  + a.wxt *b.wyt  + a.wyz *b.wxz  - a.wyt *b.wxt  + a.wzt *b.wxyzt + a.xyz *b.z    - a.xyt *b.t    + a.xzt *b.yzt  - a.yzt *b.xzt  - a.wxyz*b.wz   + a.wxyt*b.wt   - a.wxzt*b.wyzt + a.wyzt*b.wxzt + a.xyzt*b.zt   + a.wxyzt*b.wzt  ; 
c.xz   =  + a.q   *b.xz   + a.w   *b.wxz  + a.x   *b.z    - a.y   *b.xyz  - a.z   *b.x    - a.t   *b.xzt  - a.wx  *b.wz   + a.wy  *b.wxyz + a.wz  *b.wx   + a.wt  *b.wxzt + a.xy  *b.yz   + a.xz  *b.q    + a.xt  *b.zt   - a.yz  *b.xy   - a.yt  *b.xyzt - a.zt  *b.xt   + a.wxy *b.wyz  + a.wxz *b.w    + a.wxt *b.wzt  - a.wyz *b.wxy  - a.wyt *b.wxyzt - a.wzt *b.wxt  - a.xyz *b.y    - a.xyt *b.yzt  - a.xzt *b.t    + a.yzt *b.xyt  + a.wxyz*b.wy   + a.wxyt*b.wyzt + a.wxzt*b.wt   - a.wyzt*b.wxyt - a.xyzt*b.yt   - a.wxyzt*b.wyt  ; 
c.xt   =  + a.q   *b.xt   + a.w   *b.wxt  + a.x   *b.t    - a.y   *b.xyt  - a.z   *b.xzt  - a.t   *b.x    - a.wx  *b.wt   + a.wy  *b.wxyt + a.wz  *b.wxzt + a.wt  *b.wx   + a.xy  *b.yt   + a.xz  *b.zt   + a.xt  *b.q    - a.yz  *b.xyzt - a.yt  *b.xy   - a.zt  *b.xz   + a.wxy *b.wyt  + a.wxz *b.wzt  + a.wxt *b.w    - a.wyz *b.wxyzt - a.wyt *b.wxy  - a.wzt *b.wxz  - a.xyz *b.yzt  - a.xyt *b.y    - a.xzt *b.z    + a.yzt *b.xyz  + a.wxyz*b.wyzt + a.wxyt*b.wy   + a.wxzt*b.wz   - a.wyzt*b.wxyz - a.xyzt*b.yz   - a.wxyzt*b.wyz  ; 
c.yz   =  + a.q   *b.yz   + a.w   *b.wyz  + a.x   *b.xyz  + a.y   *b.z    - a.z   *b.y    - a.t   *b.yzt  - a.wx  *b.wxyz - a.wy  *b.wz   + a.wz  *b.wy   + a.wt  *b.wyzt - a.xy  *b.xz   + a.xz  *b.xy   + a.xt  *b.xyzt + a.yz  *b.q    + a.yt  *b.zt   - a.zt  *b.yt   - a.wxy *b.wxz  + a.wxz *b.wxy  + a.wxt *b.wxyzt + a.wyz *b.w    + a.wyt *b.wzt  - a.wzt *b.wyt  + a.xyz *b.x    + a.xyt *b.xzt  - a.xzt *b.xyt  - a.yzt *b.t    - a.wxyz*b.wx   - a.wxyt*b.wxzt + a.wxzt*b.wxyt + a.wyzt*b.wt   + a.xyzt*b.xt   + a.wxyzt*b.wxt  ; 
c.yt   =  + a.q   *b.yt   + a.w   *b.wyt  + a.x   *b.xyt  + a.y   *b.t    - a.z   *b.yzt  - a.t   *b.y    - a.wx  *b.wxyt - a.wy  *b.wt   + a.wz  *b.wyzt + a.wt  *b.wy   - a.xy  *b.xt   + a.xz  *b.xyzt + a.xt  *b.xy   + a.yz  *b.zt   + a.yt  *b.q    - a.zt  *b.yz   - a.wxy *b.wxt  + a.wxz *b.wxyzt + a.wxt *b.wxy  + a.wyz *b.wzt  + a.wyt *b.w    - a.wzt *b.wyz  + a.xyz *b.xzt  + a.xyt *b.x    - a.xzt *b.xyz  - a.yzt *b.z    - a.wxyz*b.wxzt - a.wxyt*b.wx   + a.wxzt*b.wxyz + a.wyzt*b.wz   + a.xyzt*b.xz   + a.wxyzt*b.wxz  ; 
c.zt   =  + a.q   *b.zt   + a.w   *b.wzt  + a.x   *b.xzt  + a.y   *b.yzt  + a.z   *b.t    - a.t   *b.z    - a.wx  *b.wxzt - a.wy  *b.wyzt - a.wz  *b.wt   + a.wt  *b.wz   - a.xy  *b.xyzt - a.xz  *b.xt   + a.xt  *b.xz   - a.yz  *b.yt   + a.yt  *b.yz   + a.zt  *b.q    - a.wxy *b.wxyzt - a.wxz *b.wxt  + a.wxt *b.wxz  - a.wyz *b.wyt  + a.wyt *b.wyz  + a.wzt *b.w    - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xzt *b.x    + a.yzt *b.y    + a.wxyz*b.wxyt - a.wxyt*b.wxyz - a.wxzt*b.wx   - a.wyzt*b.wy   - a.xyzt*b.xy   - a.wxyzt*b.wxy  ; 
c.wxy  =  + a.q   *b.wxy  + a.w   *b.xy   - a.x   *b.wy   + a.y   *b.wx   - a.z   *b.wxyz + a.t   *b.wxyt + a.wx  *b.y    - a.wy  *b.x    + a.wz  *b.xyz  - a.wt  *b.xyt  + a.xy  *b.w    - a.xz  *b.wyz  + a.xt  *b.wyt  + a.yz  *b.wxz  - a.yt  *b.wxt  + a.zt  *b.wxyzt + a.wxy *b.q    - a.wxz *b.yz   + a.wxt *b.yt   + a.wyz *b.xz   - a.wyt *b.xt   + a.wzt *b.xyzt - a.xyz *b.wz   + a.xyt *b.wt   - a.xzt *b.wyzt + a.yzt *b.wxzt + a.wxyz*b.z    - a.wxyt*b.t    + a.wxzt*b.yzt  - a.wyzt*b.xzt  + a.xyzt*b.wzt  + a.wxyzt*b.zt   ; 
c.wxz  =  + a.q   *b.wxz  + a.w   *b.xz   - a.x   *b.wz   + a.y   *b.wxyz + a.z   *b.wx   + a.t   *b.wxzt + a.wx  *b.z    - a.wy  *b.xyz  - a.wz  *b.x    - a.wt  *b.xzt  + a.xy  *b.wyz  + a.xz  *b.w    + a.xt  *b.wzt  - a.yz  *b.wxy  - a.yt  *b.wxyzt - a.zt  *b.wxt  + a.wxy *b.yz   + a.wxz *b.q    + a.wxt *b.zt   - a.wyz *b.xy   - a.wyt *b.xyzt - a.wzt *b.xt   + a.xyz *b.wy   + a.xyt *b.wyzt + a.xzt *b.wt   - a.yzt *b.wxyt - a.wxyz*b.y    - a.wxyt*b.yzt  - a.wxzt*b.t    + a.wyzt*b.xyt  - a.xyzt*b.wyt  - a.wxyzt*b.yt   ; 
c.wxt  =  + a.q   *b.wxt  + a.w   *b.xt   - a.x   *b.wt   + a.y   *b.wxyt + a.z   *b.wxzt + a.t   *b.wx   + a.wx  *b.t    - a.wy  *b.xyt  - a.wz  *b.xzt  - a.wt  *b.x    + a.xy  *b.wyt  + a.xz  *b.wzt  + a.xt  *b.w    - a.yz  *b.wxyzt - a.yt  *b.wxy  - a.zt  *b.wxz  + a.wxy *b.yt   + a.wxz *b.zt   + a.wxt *b.q    - a.wyz *b.xyzt - a.wyt *b.xy   - a.wzt *b.xz   + a.xyz *b.wyzt + a.xyt *b.wy   + a.xzt *b.wz   - a.yzt *b.wxyz - a.wxyz*b.yzt  - a.wxyt*b.y    - a.wxzt*b.z    + a.wyzt*b.xyz  - a.xyzt*b.wyz  - a.wxyzt*b.yz   ; 
c.wyz  =  + a.q   *b.wyz  + a.w   *b.yz   - a.x   *b.wxyz - a.y   *b.wz   + a.z   *b.wy   + a.t   *b.wyzt + a.wx  *b.xyz  + a.wy  *b.z    - a.wz  *b.y    - a.wt  *b.yzt  - a.xy  *b.wxz  + a.xz  *b.wxy  + a.xt  *b.wxyzt + a.yz  *b.w    + a.yt  *b.wzt  - a.zt  *b.wyt  - a.wxy *b.xz   + a.wxz *b.xy   + a.wxt *b.xyzt + a.wyz *b.q    + a.wyt *b.zt   - a.wzt *b.yt   - a.xyz *b.wx   - a.xyt *b.wxzt + a.xzt *b.wxyt + a.yzt *b.wt   + a.wxyz*b.x    + a.wxyt*b.xzt  - a.wxzt*b.xyt  - a.wyzt*b.t    + a.xyzt*b.wxt  + a.wxyzt*b.xt   ; 
c.wyt  =  + a.q   *b.wyt  + a.w   *b.yt   - a.x   *b.wxyt - a.y   *b.wt   + a.z   *b.wyzt + a.t   *b.wy   + a.wx  *b.xyt  + a.wy  *b.t    - a.wz  *b.yzt  - a.wt  *b.y    - a.xy  *b.wxt  + a.xz  *b.wxyzt + a.xt  *b.wxy  + a.yz  *b.wzt  + a.yt  *b.w    - a.zt  *b.wyz  - a.wxy *b.xt   + a.wxz *b.xyzt + a.wxt *b.xy   + a.wyz *b.zt   + a.wyt *b.q    - a.wzt *b.yz   - a.xyz *b.wxzt - a.xyt *b.wx   + a.xzt *b.wxyz + a.yzt *b.wz   + a.wxyz*b.xzt  + a.wxyt*b.x    - a.wxzt*b.xyz  - a.wyzt*b.z    + a.xyzt*b.wxz  + a.wxyzt*b.xz   ; 
c.wzt  =  + a.q   *b.wzt  + a.w   *b.zt   - a.x   *b.wxzt - a.y   *b.wyzt - a.z   *b.wt   + a.t   *b.wz   + a.wx  *b.xzt  + a.wy  *b.yzt  + a.wz  *b.t    - a.wt  *b.z    - a.xy  *b.wxyzt - a.xz  *b.wxt  + a.xt  *b.wxz  - a.yz  *b.wyt  + a.yt  *b.wyz  + a.zt  *b.w    - a.wxy *b.xyzt - a.wxz *b.xt   + a.wxt *b.xz   - a.wyz *b.yt   + a.wyt *b.yz   + a.wzt *b.q    + a.xyz *b.wxyt - a.xyt *b.wxyz - a.xzt *b.wx   - a.yzt *b.wy   - a.wxyz*b.xyt  + a.wxyt*b.xyz  + a.wxzt*b.x    + a.wyzt*b.y    - a.xyzt*b.wxy  - a.wxyzt*b.xy   ; 
c.xyz  =  + a.q   *b.xyz  + a.w   *b.wxyz + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.t   *b.xyzt - a.wx  *b.wyz  + a.wy  *b.wxz  - a.wz  *b.wxy  - a.wt  *b.wxyzt + a.xy  *b.z    - a.xz  *b.y    - a.xt  *b.yzt  + a.yz  *b.x    + a.yt  *b.xzt  - a.zt  *b.xyt  + a.wxy *b.wz   - a.wxz *b.wy   - a.wxt *b.wyzt + a.wyz *b.wx   + a.wyt *b.wxzt - a.wzt *b.wxyt + a.xyz *b.q    + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   - a.wxyz*b.w    - a.wxyt*b.wzt  + a.wxzt*b.wyt  - a.wyzt*b.wxt  - a.xyzt*b.t    - a.wxyzt*b.wt   ; 
c.xyt  =  + a.q   *b.xyt  + a.w   *b.wxyt + a.x   *b.yt   - a.y   *b.xt   + a.z   *b.xyzt + a.t   *b.xy   - a.wx  *b.wyt  + a.wy  *b.wxt  - a.wz  *b.wxyzt - a.wt  *b.wxy  + a.xy  *b.t    - a.xz  *b.yzt  - a.xt  *b.y    + a.yz  *b.xzt  + a.yt  *b.x    - a.zt  *b.xyz  + a.wxy *b.wt   - a.wxz *b.wyzt - a.wxt *b.wy   + a.wyz *b.wxzt + a.wyt *b.wx   - a.wzt *b.wxyz + a.xyz *b.zt   + a.xyt *b.q    - a.xzt *b.yz   + a.yzt *b.xz   - a.wxyz*b.wzt  - a.wxyt*b.w    + a.wxzt*b.wyz  - a.wyzt*b.wxz  - a.xyzt*b.z    - a.wxyzt*b.wz   ; 
c.xzt  =  + a.q   *b.xzt  + a.w   *b.wxzt + a.x   *b.zt   - a.y   *b.xyzt - a.z   *b.xt   + a.t   *b.xz   - a.wx  *b.wzt  + a.wy  *b.wxyzt + a.wz  *b.wxt  - a.wt  *b.wxz  + a.xy  *b.yzt  + a.xz  *b.t    - a.xt  *b.z    - a.yz  *b.xyt  + a.yt  *b.xyz  + a.zt  *b.x    + a.wxy *b.wyzt + a.wxz *b.wt   - a.wxt *b.wz   - a.wyz *b.wxyt + a.wyt *b.wxyz + a.wzt *b.wx   - a.xyz *b.yt   + a.xyt *b.yz   + a.xzt *b.q    - a.yzt *b.xy   + a.wxyz*b.wyt  - a.wxyt*b.wyz  - a.wxzt*b.w    + a.wyzt*b.wxy  + a.xyzt*b.y    + a.wxyzt*b.wy   ; 
c.yzt  =  + a.q   *b.yzt  + a.w   *b.wyzt + a.x   *b.xyzt + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   - a.wx  *b.wxyzt - a.wy  *b.wzt  + a.wz  *b.wyt  - a.wt  *b.wyz  - a.xy  *b.xzt  + a.xz  *b.xyt  - a.xt  *b.xyz  + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    - a.wxy *b.wxzt + a.wxz *b.wxyt - a.wxt *b.wxyz + a.wyz *b.wt   - a.wyt *b.wz   + a.wzt *b.wy   + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.yzt *b.q    - a.wxyz*b.wxt  + a.wxyt*b.wxz  - a.wxzt*b.wxy  - a.wyzt*b.w    - a.xyzt*b.x    - a.wxyzt*b.wx   ; 
c.wxyz =  + a.q   *b.wxyz + a.w   *b.xyz  - a.x   *b.wyz  + a.y   *b.wxz  - a.z   *b.wxy  - a.t   *b.wxyzt + a.wx  *b.yz   - a.wy  *b.xz   + a.wz  *b.xy   + a.wt  *b.xyzt + a.xy  *b.wz   - a.xz  *b.wy   - a.xt  *b.wyzt + a.yz  *b.wx   + a.yt  *b.wxzt - a.zt  *b.wxyt + a.wxy *b.z    - a.wxz *b.y    - a.wxt *b.yzt  + a.wyz *b.x    + a.wyt *b.xzt  - a.wzt *b.xyt  - a.xyz *b.w    - a.xyt *b.wzt  + a.xzt *b.wyt  - a.yzt *b.wxt  + a.wxyz*b.q    + a.wxyt*b.zt   - a.wxzt*b.yt   + a.wyzt*b.xt   - a.xyzt*b.wt   - a.wxyzt*b.t    ; 
c.wxyt =  + a.q   *b.wxyt + a.w   *b.xyt  - a.x   *b.wyt  + a.y   *b.wxt  - a.z   *b.wxyzt - a.t   *b.wxy  + a.wx  *b.yt   - a.wy  *b.xt   + a.wz  *b.xyzt + a.wt  *b.xy   + a.xy  *b.wt   - a.xz  *b.wyzt - a.xt  *b.wy   + a.yz  *b.wxzt + a.yt  *b.wx   - a.zt  *b.wxyz + a.wxy *b.t    - a.wxz *b.yzt  - a.wxt *b.y    + a.wyz *b.xzt  + a.wyt *b.x    - a.wzt *b.xyz  - a.xyz *b.wzt  - a.xyt *b.w    + a.xzt *b.wyz  - a.yzt *b.wxz  + a.wxyz*b.zt   + a.wxyt*b.q    - a.wxzt*b.yz   + a.wyzt*b.xz   - a.xyzt*b.wz   - a.wxyzt*b.z    ; 
c.wxzt =  + a.q   *b.wxzt + a.w   *b.xzt  - a.x   *b.wzt  + a.y   *b.wxyzt + a.z   *b.wxt  - a.t   *b.wxz  + a.wx  *b.zt   - a.wy  *b.xyzt - a.wz  *b.xt   + a.wt  *b.xz   + a.xy  *b.wyzt + a.xz  *b.wt   - a.xt  *b.wz   - a.yz  *b.wxyt + a.yt  *b.wxyz + a.zt  *b.wx   + a.wxy *b.yzt  + a.wxz *b.t    - a.wxt *b.z    - a.wyz *b.xyt  + a.wyt *b.xyz  + a.wzt *b.x    + a.xyz *b.wyt  - a.xyt *b.wyz  - a.xzt *b.w    + a.yzt *b.wxy  - a.wxyz*b.yt   + a.wxyt*b.yz   + a.wxzt*b.q    - a.wyzt*b.xy   + a.xyzt*b.wy   + a.wxyzt*b.y    ; 
c.wyzt =  + a.q   *b.wyzt + a.w   *b.yzt  - a.x   *b.wxyzt - a.y   *b.wzt  + a.z   *b.wyt  - a.t   *b.wyz  + a.wx  *b.xyzt + a.wy  *b.zt   - a.wz  *b.yt   + a.wt  *b.yz   - a.xy  *b.wxzt + a.xz  *b.wxyt - a.xt  *b.wxyz + a.yz  *b.wt   - a.yt  *b.wz   + a.zt  *b.wy   - a.wxy *b.xzt  + a.wxz *b.xyt  - a.wxt *b.xyz  + a.wyz *b.t    - a.wyt *b.z    + a.wzt *b.y    - a.xyz *b.wxt  + a.xyt *b.wxz  - a.xzt *b.wxy  - a.yzt *b.w    + a.wxyz*b.xt   - a.wxyt*b.xz   + a.wxzt*b.xy   + a.wyzt*b.q    - a.xyzt*b.wx   - a.wxyzt*b.x    ; 
c.xyzt =  + a.q   *b.xyzt + a.w   *b.wxyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  - a.wx  *b.wyzt + a.wy  *b.wxzt - a.wz  *b.wxyt + a.wt  *b.wxyz + a.xy  *b.zt   - a.xz  *b.yt   + a.xt  *b.yz   + a.yz  *b.xt   - a.yt  *b.xz   + a.zt  *b.xy   + a.wxy *b.wzt  - a.wxz *b.wyt  + a.wxt *b.wyz  + a.wyz *b.wxt  - a.wyt *b.wxz  + a.wzt *b.wxy  + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    - a.wxyz*b.wt   + a.wxyt*b.wz   - a.wxzt*b.wy   + a.wyzt*b.wx   + a.xyzt*b.q    + a.wxyzt*b.w    ; 
c.wxyzt =  + a.q   *b.wxyzt + a.w   *b.xyzt - a.x   *b.wyzt + a.y   *b.wxzt - a.z   *b.wxyt + a.t   *b.wxyz + a.wx  *b.yzt  - a.wy  *b.xzt  + a.wz  *b.xyt  - a.wt  *b.xyz  + a.xy  *b.wzt  - a.xz  *b.wyt  + a.xt  *b.wyz  + a.yz  *b.wxt  - a.yt  *b.wxz  + a.zt  *b.wxy  + a.wxy *b.zt   - a.wxz *b.yt   + a.wxt *b.yz   + a.wyz *b.xt   - a.wyt *b.xz   + a.wzt *b.xy   - a.xyz *b.wt   + a.xyt *b.wz   - a.xzt *b.wy   + a.yzt *b.wx   + a.wxyz*b.t    - a.wxyt*b.z    + a.wxzt*b.y    - a.wyzt*b.x    + a.xyzt*b.w    + a.wxyzt*b.q    ; 


	c.q = expand(c.q);

	c.w = expand(c.w);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);

	c.wx = expand(c.wx);
	c.wy = expand(c.wy);
	c.wz = expand(c.wz);
	c.wt = expand(c.wt);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);

	c.wxy = expand(c.wxy);
	c.wxz = expand(c.wxz);
	c.wxt = expand(c.wxt);
	c.wyz = expand(c.wyz);
	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);

	c.wxyz = expand(c.wxyz);
	c.wxyt = expand(c.wxyt);
	c.wxzt = expand(c.wxzt);
	c.wyzt = expand(c.wyzt);
	c.xyzt = expand(c.xyzt);

	c.wxyzt = expand(c.wxyzt);

	return c;
}

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

	c.w = expand(c.w);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);
	c.t = expand(c.t);

	c.wx = expand(c.wx);
	c.wy = expand(c.wy);
	c.wz = expand(c.wz);
	c.wt = expand(c.wt);
	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xt = expand(c.xt);
	c.yt = expand(c.yt);
	c.zt = expand(c.zt);

	c.wxy = expand(c.wxy);
	c.wxz = expand(c.wxz);
	c.wxt = expand(c.wxt);
	c.wyz = expand(c.wyz);
	c.wyt = expand(c.wyt);
	c.wzt = expand(c.wzt);
	c.xyz = expand(c.xyz);
	c.xyt = expand(c.xyt);
	c.xzt = expand(c.xzt);
	c.yzt = expand(c.yzt);

	c.wxyz = expand(c.wxyz);
	c.wxyt = expand(c.wxyt);
	c.wxzt = expand(c.wxzt);
	c.wyzt = expand(c.wyzt);
	c.xyzt = expand(c.xyzt);

	c.wxyzt = expand(c.wxyzt);

	return c;
}

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


GA5_4_1 S(const GA5_4_1 &u, const GA5_4_1 &v)	// symmetric product
{
	GA5_4_1 a;
	a = u*v/2 + v*u/2;
	return a;
}

GA5_4_1 A(const GA5_4_1 &u, const GA5_4_1 &v)	// antisymmetric product
{
	GA5_4_1 a;
	a = u*v/2 - v*u/2;
	return a;
}

GA5_4_1 Dual(const GA5_4_1 &a)	// dual =  a*I_inv 
{
	GA5_4_1 I_inv,c;	// I_inv created initially zero
	I_inv.wxyzt = -1;	// pure pseudoscalar
	c = a*I_inv;	// form dual
	return c;
}

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


	int ii,jj;	// loop counters

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

	GA5_4_1 FiveVect;
	GA5_4_1 MV5, MV6, MV7, MV8, MV9, MV10;
//	FiveVect.q = a;     FiveVect.w = b;     FiveVect.x = c;     FiveVect.y = d;     FiveVect.z = e;     FiveVect.t = f;
//	FiveVect.wx = g;    FiveVect.wy = h;    FiveVect.wz = j;    FiveVect.wt = k;    FiveVect.xy = l;
//	FiveVect.xz = m;    FiveVect.xt = n;    FiveVect.yz = p;    FiveVect.yt = r;    FiveVect.zt = s;
//
//	FiveVect.wxy = S;   FiveVect.wxz = R;   FiveVect.wxt = P;   FiveVect.wyz = N;   FiveVect.wyt = M;
//	FiveVect.wzt = L;   FiveVect.xyz = K;   FiveVect.xyt = J;   FiveVect.xzt = H;   FiveVect.yzt = G;
//
//	FiveVect.wxyz = F;  FiveVect.wxyt = E;  FiveVect.wxzt = D;  FiveVect.wyzt = C;  FiveVect.xyzt = B;  FiveVect.wxyzt = A;

	FiveVect = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "FiveVect       = " << FiveVect << "\n";

	MV5 = Dual(FiveVect);
	cout << "Dual(FiveVect) = " << MV5 << "\n\n";

	MV6 = Conjugate(FiveVect);
	cout << "Conjugate(FiveVect) = " << MV6 << "\n\n";

	MV7 = Parity(FiveVect);
	cout << "Parity(FiveVect) = " << MV7 << "\n\n";

	MV8 = Transpose(FiveVect);
	cout << "Transpose(FiveVect) = " << MV8 << "\n\n";

	MV9 = Hermitian(FiveVect);
	cout << "Hermitian(FiveVect) = " << MV9 << "\n\n";
	MV9 = FiveVect*Hermitian(FiveVect);
	cout << "FiveVect*Hermitian(FiveVect) = " << MV9 << "\n\n";

	MV10 = Reverse(FiveVect);
	cout << "Reverse(FiveVect) = " << MV10 << "\n\n";

	MV5 = FiveVect*Reverse(FiveVect);
	cout << "FiveVect*Reverse(FiveVect) = " << MV5 << "\n\n";

//	MV6 = MV5*Conjugate(MV5);
//	cout << "MV5*Conjugate(MV5) = " << MV6 << "\n\n";

/*	ex sam;
	sam = a + I*b;
	cout << "sam = " << sam << "\n";
	cout << "sam.real() = " << sam.real_part() << "\n";
	cout << "sam.imag() = " << sam.imag_part() << "\n";
		sam = a+I*b
		sam.real() = -imag_part(b)+real_part(a)
		sam.imag() = real_part(b)+imag_part(a)
*/
////////////////////////////////////  begin matrix land


	zero  =  0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0;

	
	Q =  	 1,  0,  0,  0,
		 0,  1,  0,  0,
		 0,  0,  1,  0,
		 0,  0,  0,  1;

	W = 	 0,  0,  0, -I,
		 0,  0,  I,  0,
		 0, -I,  0,  0,
		 I,  0,  0,  0;

	X = 	 1,  0,  0,  0,
		 0, -1,  0,  0,
		 0,  0,  1,  0,
		 0,  0,  0, -1;

	Y = 	 0,  1,  0,  0, 
		 1,  0,  0,  0,  
		 0,  0,  0,  1,
		 0,  0,  1,  0;

	Z = 	 0,  0,  0, -1,
		 0,  0,  1,  0,
		 0,  1,  0,  0, 
		-1,  0,  0,  0;

	T = 	 0, -1,  0,  0,  
		 1,  0,  0,  0,
		 0,  0,  0,  1,
		 0,  0, -1,  0;
/*
	GA5_4_1 q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA5_4_1 wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA5_4_1 wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;
*/



//	create our remaining terms


	WX = W.mul(X);	// bivectors
	WY = W.mul(Y);
	WZ = W.mul(Z);
	WT = W.mul(T);
	XY = X.mul(Y);
	XZ = X.mul(Z);
	XT = X.mul(T);
	YZ = Y.mul(Z);
	YT = Y.mul(T);
	ZT = Z.mul(T);


	WXY = WX.mul(Y); // Trivectors
	WXZ = WX.mul(Z);
	WXT = WX.mul(T);
	WYZ = WY.mul(Z);
	WYT = WY.mul(T);
	WZT = WZ.mul(T);
	XYZ = XY.mul(Z);
	XYT = XY.mul(T);
	XZT = XZ.mul(T);
	YZT = YZ.mul(T);


	WXYZ = WXY.mul(Z);  // Tetra vectors
	WXYT = WXY.mul(T);
	WXZT = WXZ.mul(T);
	WYZT = WYZ.mul(T);
	XYZT = XYZ.mul(T);

	WXYZT = WXYZ.mul(T);	// penta vector

//	build generic MV in matrix format

	matrix spacetime(4,4);

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ g*WX(ii,jj) + h*WY(ii,jj) + j*WZ(ii,jj) + k*WT(ii,jj) + l*XY(ii,jj)
					+ m*XZ(ii,jj) + n*XT(ii,jj) + p*YZ(ii,jj) + r*YT(ii,jj) + s*ZT(ii,jj)

					+ S*WXY(ii,jj) + R*WXZ(ii,jj) + P*WXT(ii,jj) + N*WYZ(ii,jj) + M*WYT(ii,jj)
					+ L*WZT(ii,jj) + K*XYZ(ii,jj) + J*XYT(ii,jj) + H*XZT(ii,jj) + G*YZT(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

	cout << "spacetime = " << spacetime << "\n";
	cout << "spacetime.transpose() = " << spacetime.transpose() << "\n";
	matrix hmm(4,4);
	matrix arg(4,4);

	arg = spacetime.transpose();
	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			hmm(ii,jj) = expand(spacetime(ii,jj) - arg(ii,jj));
		}
	}

	cout << "spacetime - spacetime.transpose() = " << hmm << "\n";

/*
spacetime = 
[[I*j+J+c-I*C+a+I*A+r-I*R,-f+I*F-n-I*N+d+I*D+l-I*L,s+I*S+I*k+K+p-I*P-I*h+H,-I*M+I*g+G-e+I*E-I*b+B-m],
 [ f-I*F-n-I*N+d+I*D-l+I*L,I*j+J-c+I*C+a+I*A-r+I*R,-I*M+I*g+G+e-I*E+I*b-B-m,s+I*S+I*k+K-p+I*P+I*h-H],
 [ s+I*S-I*k-K-p+I*P-I*h+H,I*M+I*g+G+e-I*E-I*b+B+m,-I*j-J+c-I*C+a+I*A-r+I*R,f-I*F+n+I*N+d+I*D+l-I*L],
 [ I*M+I*g+G-e+I*E+I*b-B+m,s+I*S-I*k-K+p-I*P+I*h-H,-f+I*F+n+I*N+d+I*D-l+I*L,-I*j-J-c+I*C+a+I*A+r-I*R]]
spacetime.transpose() = 
[[I*j+J+c-I*C+a+I*A+r-I*R,f-I*F-n-I*N+d+I*D-l+I*L,s+I*S-I*k-K-p+I*P-I*h+H,I*M+I*g+G-e+I*E+I*b-B+m],
 [-f+I*F-n-I*N+d+I*D+l-I*L,I*j+J-c+I*C+a+I*A-r+I*R,I*M+I*g+G+e-I*E-I*b+B+m,s+I*S-I*k-K+p-I*P+I*h-H],
 [ s+I*S+I*k+K+p-I*P-I*h+H,-I*M+I*g+G+e-I*E+I*b-B-m,-I*j-J+c-I*C+a+I*A-r+I*R,-f+I*F+n+I*N+d+I*D-l+I*L],
 [-I*M+I*g+G-e+I*E-I*b+B-m,s+I*S+I*k+K-p+I*P+I*h-H,f-I*F+n+I*N+d+I*D+l-I*L,-I*j-J-c+I*C+a+I*A+r-I*R]]
spacetime - spacetime.transpose() = 
[[0,2*l-(2*I)*L-2*f+(2*I)*F,2*p-(2*I)*P+(2*I)*k+2*K,-(2*I)*b+2*B-2*m-(2*I)*M],   b f k l m p P M L K F B
 [-2*l+(2*I)*L+2*f-(2*I)*F,0,(2*I)*b-2*B-2*m-(2*I)*M,-2*p+(2*I)*P+(2*I)*k+2*K],  b f k l m p
 [-2*p+(2*I)*P-(2*I)*k-2*K,-(2*I)*b+2*B+2*m+(2*I)*M,0,2*l-(2*I)*L+2*f-(2*I)*F],  b f k l m p
 [(2*I)*b-2*B+2*m+(2*I)*M,2*p-(2*I)*P-(2*I)*k-2*K,-2*l+(2*I)*L-2*f+(2*I)*F,0]]

	FiveVect = GA5_4_1(+a, +b,+c,+d,+e,+f, +g,+h,+j,+k,+l, +m,+n,+p,+r,+s, +S,+R,+P,+N,+M, +L,+K,+J,+H,+G, +F,+E,+D,+C,+B, +A);

	FiveVect = GA5_4_1(+a, -b,+c,+d,+e,-f, +g,+h,+j,-k,-l, -m,+n,-p,+r,+s, +S,+R,-P,+N,-M, -L,-K,+J,+H,+G, -F,+E,+D,+C,-B, +A);
Transpose

*/

	MV5 = wxyzt*wxyzt;
	cout << "wxyzt*wxyzt = " << MV5 << "\n";

	cout << "check determinates of spacetime and it's transpose\n\n";
	cout << "spacetime.determinant() = " << expand(spacetime.determinant()) << "\n\n";
	arg = spacetime.transpose();
	cout << "arg.determinant() = " << expand(arg.determinant()) << "\n\n";
	cout << "delta determinants = " << expand(spacetime.determinant() - arg.determinant()) << " expect 0\n\n";


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ g*WX(ii,jj) + h*WY(ii,jj) + j*WZ(ii,jj) + k*WT(ii,jj) + l*XY(ii,jj)
					+ m*XZ(ii,jj) + n*XT(ii,jj) + p*YZ(ii,jj) + r*YT(ii,jj) + s*ZT(ii,jj)

					+ S*WXY(ii,jj) + R*WXZ(ii,jj) + P*WXT(ii,jj) + N*WYZ(ii,jj) + M*WYT(ii,jj)
					+ L*WZT(ii,jj) + K*XYZ(ii,jj) + J*XYT(ii,jj) + H*XZT(ii,jj) + G*YZT(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

	matrix reversed(4,4);
	matrix product(4,4);


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			reversed(ii,jj) = expand(  // Reverse(spacetime) negate bivectors, trivectors
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)

					- g*WX(ii,jj) - h*WY(ii,jj) - j*WZ(ii,jj) - k*WT(ii,jj) - l*XY(ii,jj)
					- m*XZ(ii,jj) - n*XT(ii,jj) - p*YZ(ii,jj) - r*YT(ii,jj) - s*ZT(ii,jj)

					- S*WXY(ii,jj) - R*WXZ(ii,jj) - P*WXT(ii,jj) - N*WYZ(ii,jj) - M*WYT(ii,jj)
					- L*WZT(ii,jj) - K*XYZ(ii,jj) - J*XYT(ii,jj) - H*XZT(ii,jj) - G*YZT(ii,jj)

					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

	product = spacetime.mul(reversed);

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			product(ii,jj) = expand(  product(ii,jj) );
		}
	}

//	cout << "spacetime = " << spacetime << "\n";
//	cout << "reversed = " << reversed << "\n";
//	cout << "product = " << product << "\n";
//
//	ex det1, det2, delta;
//	det1 = spacetime.determinant();
//	det2 = reversed.determinant();
//	delta = expand(expand (det1) - expand(det2));

//	cout << "det1 - det2 = " << delta << "\n";
//	The MV and its reverse have the same determinant.

////////////////////////////////////

	FiveVect = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "FiveVect       = " << FiveVect << "\n";
	MV10 = Reverse(FiveVect);
	cout << "Reverse(FiveVect) = " << MV10 << "\n\n";

	MV5 = FiveVect*Reverse(FiveVect);
	cout << "FiveVect*Reverse(FiveVect) = " << MV5 << "\n\n";
/*
index = 2046  count = 30 
s = 
(         3,
         -5,        -7,       -11,       -13,       -17, 
          0,         0,         0,         0,         0,         0,         0,         0,         0,         0, 
          0,         0,         0,         0,         0,         0,         0,         0,         0,         0, 
        -19,       -23,       -29,       -31,       -37, 
         41) 

t = 
(      1592,
          0,         0,         0,         0,         0, 
          0,         0,         0,         0,         0,         0,         0,         0,         0,         0, 
          0,         0,         0,         0,         0,         0,         0,         0,         0,         0, 
          0,         0,         0,         0,         0, 
       -376) 
*/
////////////////////////////////////

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ 3*Q(ii,jj)
					+ 5*W(ii,jj) + 7*X(ii,jj) + 11*Y(ii,jj) + 13*Z(ii,jj) + 17*T(ii,jj)
					+ 19*WXYZ(ii,jj) + 23*WXYT(ii,jj) + 29*WXZT(ii,jj) + 31*WYZT(ii,jj) + 37*XYZT(ii,jj)
					+ 41*WXYZT(ii,jj) 
                        );
		}
	}

	cout << "Test vector determinant = " << spacetime.determinant() << "\n";
// Test vector determinant = 2393088-1197184*I


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}
	ex	det1, det2, delta;
	det1 = spacetime.determinant();

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			hmm(ii,jj) = expand(
					+ a*Q(ii,jj)
					- b*W(ii,jj) - c*X(ii,jj) - d*Y(ii,jj) - e*Z(ii,jj) - f*T(ii,jj)
					- F*WXYZ(ii,jj) - E*WXYT(ii,jj) - D*WXZT(ii,jj) - C*WYZT(ii,jj) - B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}
	det2 = hmm.determinant();

	cout << "det1 = " << det1 << "\n\n";
	cout << "det2 = " << det2 << "\n\n";
	delta = expand(expand (det1) - expand(det2));

	cout << "det1 - det2 = " << delta << " (match! good!)\n";

	arg = spacetime.mul(hmm);
	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			arg(ii,jj) = expand(arg(ii,jj));
		}
	}

	cout << "arg = " << arg << "\n";

	GA5_4_1 MV1, MV2, MV3;
	MV1 = GA5_4_1(a, b,c,d,e,f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, F,E,D,C,B, A);
	MV2 = GA5_4_1(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			hmm(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}
	det1 = hmm.determinant();
	cout << "Matrix Determinant(a + A*wxyzt) = " << det1 << "\n";
	delta = expand( det1 - ( (a + I*A)*(a + I*A)*(a + I*A)*(a + I*A) ) ) ;
	cout << " det1 - ( (a + I*A)^4 = " << delta << "\n";

	MV1 = GA5_4_1(a, 0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,  A);
	MV2 = GA5_4_1(a, 0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0, -A);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			hmm(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}
	det1 = hmm.determinant();

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			arg(ii,jj) = expand(
					+ a*Q(ii,jj)
					- A*WXYZT(ii,jj) 
                        );
		}
	}
	det2 = arg.determinant();
	cout << "Matrix Determinant(a + A*wxyzt) = " << det1 << "\n";
	cout << "Matrix Determinant(a - A*wxyzt) = " << det2 << "\n";
	delta = expand( det1 - det2 ) ;
	cout << " det1 - det2 = " << delta << "\n";

// now do the two full determinants

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ g*WX(ii,jj) + h*WY(ii,jj) + j*WZ(ii,jj) + k*WT(ii,jj) + l*XY(ii,jj)
					+ m*XZ(ii,jj) + n*XT(ii,jj) + p*YZ(ii,jj) + r*YT(ii,jj) + s*ZT(ii,jj)

					+ S*WXY(ii,jj) + R*WXZ(ii,jj) + P*WXT(ii,jj) + N*WYZ(ii,jj) + M*WYT(ii,jj)
					+ L*WZT(ii,jj) + K*XYZ(ii,jj) + J*XYT(ii,jj) + H*XZT(ii,jj) + G*YZT(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

/*	cout << "\nNow we compare matrix versus multivector formula\n";
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det1 = spacetime.determinant();
	det2 = Determinant(MV1);
	delta = expand( det1 - det2 ) ;
	cout << " det1 - det2 = " << delta << " Agreement is so nice!\n";
*/
//  Now we build faster, specialized formulas for determinant
	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
	a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
	a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
	a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
	-a_wx, -a_wy, -a_wz, -a_wt, -a_xy, -a_xz, -a_xt, -a_yz, -a_yt, -a_zt, 
	-a_wxy, -a_wxz, -a_wxt, -a_wyz, -a_wyt, -a_wzt, -a_xyz, -a_xyt, -a_xzt, -a_yzt,
	a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV3 = MV1*MV2;
/*	cout << "MV3 = " << MV3 << "\n\n";

	cout << "a = " << MV3.q << " ; \n";
	cout << "b = " << MV3.w << " ; \n";
	cout << "c = " << MV3.x << " ; \n";
	cout << "d = " << MV3.y << " ; \n";
	cout << "e = " << MV3.z << " ; \n";
	cout << "f = " << MV3.t << " ; \n";

	cout << "F = " << MV3.wxyz  << " ; \n";
	cout << "E = " << MV3.wxyt  << " ; \n";
	cout << "D = " << MV3.wxzt  << " ; \n";
	cout << "C = " << MV3.wyzt  << " ; \n";
	cout << "B = " << MV3.xyzt  << " ; \n";
	cout << "A = " << MV3.wxyzt << " ; \n";
*/
/*	MV1 = GA5_4_1(a, b,c,d,e,f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, F,E,D,C,B, A);
	MV2 = GA5_4_1(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "g = " << MV3.q << " ; \n";

	cout << "G = " << MV3.wxyzt  << " ; \n";
*/
/*	cout << "\nNow we compare matrix versus multivector formula\n";
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det1 = spacetime.determinant();
	det2 = Determinant(MV1);
	delta = expand( det1 - det2 ) ;
	cout << " det1 - det2 = " << delta << " Agreement is so nice!\n";
*/
////////////////////////////////////

//	do some time calculations

/*	time_t start_time, end_time;
	double deltat;
	
	time(&start_time);
	det1 = spacetime.determinant();	//10x
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	det1 = spacetime.determinant();
	time(&end_time);
	deltat = difftime(end_time, start_time);
	cout << "deltat for 10X matrix determinant = " << deltat << "\n";

	int i;
	time(&start_time);
	for(i=0;i<10000;i++ ) det2 = Determinant(MV1);
	time(&end_time);
	deltat = difftime(end_time, start_time);
	cout << "deltat for 10000X MV determinant = " << deltat << "\n";
*/

/*
Now we compare matrix versus multivector formula
 det1 - det2 = 0 Agreement is so nice!
deltat for 10X matrix determinant = 45    4.5 seconds per determinant
deltat for 10000X MV determinant = 38     3.8 msec per determinant!

*/

	// do a prime integer MV determinant, checking for pythagorean  set

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	cout << "MV1 = " << MV1 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
//Determinant(MV1) = -98748240-54398496*I
// cc = 112740459
// cc = sqrt(12710411270159616)


/*	ex aa,bb,cc,dd,ff,gg;
	aa = -98748240;	bb = -54398496;
	cc = sqrt(aa*aa + bb*bb);
	cout << "cc = " << cc << "\n";
	dd = cc*cc - aa*aa - bb*bb;
	cout << "dd = " << dd << "\n";

	ff = 112740459;
	dd = ff*ff - (aa*aa + bb*bb);
	cout << "dd = 112740459*112740459 - 12710411270159616 = " << dd << "\n";
*/
	MV2 = Reciprocal(MV1);
	cout << "MV2 = " << MV2 << "\n";
	MV3 = MV1*MV2;
	cout << "MV3 = " << MV3 << "\n";



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
	r = Set(  3,    5,  7, 11, 13, 17,
		 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 
		 19, 23, 29, 31, 37,    41) ; 

*/





