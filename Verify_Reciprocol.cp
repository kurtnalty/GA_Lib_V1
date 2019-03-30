// compare matrix product and geometric product

// compile by
// g++ Demo_Mink.cp -l ginac -l cln


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;


//********************************* Mink ***************************************

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

ostream &operator<<(ostream &ff, Mink &v) {
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



void PrintMink(Mink &v) {
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

Mink Product_Clifford(const Mink &a, const Mink &b)
{

	Mink c;
c.q    =  + a.q   *b.q    + a.x   *b.x    + a.y   *b.y    + a.z   *b.z    - a.t   *b.t    - a.xy  *b.xy   - a.xz  *b.xz   - a.yz  *b.yz   + a.xt  *b.xt   + a.yt  *b.yt   + a.zt  *b.zt   - a.xyz *b.xyz  + a.xyt *b.xyt  + a.xzt *b.xzt  + a.yzt *b.yzt  - a.xyzt*b.xyzt ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    - a.y   *b.xy   - a.z   *b.xz   + a.t   *b.xt   + a.xy  *b.y    + a.xz  *b.z    - a.yz  *b.xyz  - a.xt  *b.t    + a.yt  *b.xyt  + a.zt  *b.xzt  - a.xyz *b.yz   + a.xyt *b.yt   + a.xzt *b.zt   - a.yzt *b.xyzt + a.xyzt*b.yzt  ; 
c.y    =  + a.q   *b.y    + a.x   *b.xy   + a.y   *b.q    - a.z   *b.yz   + a.t   *b.yt   - a.xy  *b.x    + a.xz  *b.xyz  + a.yz  *b.z    - a.xt  *b.xyt  - a.yt  *b.t    + a.zt  *b.yzt  + a.xyz *b.xz   - a.xyt *b.xt   + a.xzt *b.xyzt + a.yzt *b.zt   - a.xyzt*b.xzt  ; 
c.z    =  + a.q   *b.z    + a.x   *b.xz   + a.y   *b.yz   + a.z   *b.q    + a.t   *b.zt   - a.xy  *b.xyz  - a.xz  *b.x    - a.yz  *b.y    - a.xt  *b.xzt  - a.yt  *b.yzt  - a.zt  *b.t    - a.xyz *b.xy   - a.xyt *b.xyzt - a.xzt *b.xt   - a.yzt *b.yt   + a.xyzt*b.xyt  ; 
c.t    =  + a.q   *b.t    + a.x   *b.xt   + a.y   *b.yt   + a.z   *b.zt   + a.t   *b.q    - a.xy  *b.xyt  - a.xz  *b.xzt  - a.yz  *b.yzt  - a.xt  *b.x    - a.yt  *b.y    - a.zt  *b.z    - a.xyz *b.xyzt - a.xyt *b.xy   - a.xzt *b.xz   - a.yzt *b.yz   + a.xyzt*b.xyz  ; 
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.z   *b.xyz  - a.t   *b.xyt  + a.xy  *b.q    - a.xz  *b.yz   + a.yz  *b.xz   + a.xt  *b.yt   - a.yt  *b.xt   + a.zt  *b.xyzt + a.xyz *b.z    - a.xyt *b.t    + a.xzt *b.yzt  - a.yzt *b.xzt  + a.xyzt*b.zt   ; 
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.y   *b.xyz  - a.z   *b.x    - a.t   *b.xzt  + a.xy  *b.yz   + a.xz  *b.q    - a.yz  *b.xy   + a.xt  *b.zt   - a.yt  *b.xyzt - a.zt  *b.xt   - a.xyz *b.y    - a.xyt *b.yzt  - a.xzt *b.t    + a.yzt *b.xyt  - a.xyzt*b.yt   ; 
c.yz   =  + a.q   *b.yz   + a.x   *b.xyz  + a.y   *b.z    - a.z   *b.y    - a.t   *b.yzt  - a.xy  *b.xz   + a.xz  *b.xy   + a.yz  *b.q    + a.xt  *b.xyzt + a.yt  *b.zt   - a.zt  *b.yt   + a.xyz *b.x    + a.xyt *b.xzt  - a.xzt *b.xyt  - a.yzt *b.t    + a.xyzt*b.xt   ; 
c.xt   =  + a.q   *b.xt   + a.x   *b.t    - a.y   *b.xyt  - a.z   *b.xzt  - a.t   *b.x    + a.xy  *b.yt   + a.xz  *b.zt   - a.yz  *b.xyzt + a.xt  *b.q    - a.yt  *b.xy   - a.zt  *b.xz   - a.xyz *b.yzt  - a.xyt *b.y    - a.xzt *b.z    + a.yzt *b.xyz  - a.xyzt*b.yz   ; 
c.yt   =  + a.q   *b.yt   + a.x   *b.xyt  + a.y   *b.t    - a.z   *b.yzt  - a.t   *b.y    - a.xy  *b.xt   + a.xz  *b.xyzt + a.yz  *b.zt   + a.xt  *b.xy   + a.yt  *b.q    - a.zt  *b.yz   + a.xyz *b.xzt  + a.xyt *b.x    - a.xzt *b.xyz  - a.yzt *b.z    + a.xyzt*b.xz   ; 
c.zt   =  + a.q   *b.zt   + a.x   *b.xzt  + a.y   *b.yzt  + a.z   *b.t    - a.t   *b.z    - a.xy  *b.xyzt - a.xz  *b.xt   - a.yz  *b.yt   + a.xt  *b.xz   + a.yt  *b.yz   + a.zt  *b.q    - a.xyz *b.xyt  + a.xyt *b.xyz  + a.xzt *b.x    + a.yzt *b.y    - a.xyzt*b.xy   ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.t   *b.xyzt + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    - a.xt  *b.yzt  + a.yt  *b.xzt  - a.zt  *b.xyt  + a.xyz *b.q    + a.xyt *b.zt   - a.xzt *b.yt   + a.yzt *b.xt   - a.xyzt*b.t    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.z   *b.xyzt + a.t   *b.xy   + a.xy  *b.t    - a.xz  *b.yzt  + a.yz  *b.xzt  - a.xt  *b.y    + a.yt  *b.x    - a.zt  *b.xyz  + a.xyz *b.zt   + a.xyt *b.q    - a.xzt *b.yz   + a.yzt *b.xz   - a.xyzt*b.z    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.y   *b.xyzt - a.z   *b.xt   + a.t   *b.xz   + a.xy  *b.yzt  + a.xz  *b.t    - a.yz  *b.xyt  - a.xt  *b.z    + a.yt  *b.xyz  + a.zt  *b.x    - a.xyz *b.yt   + a.xyt *b.yz   + a.xzt *b.q    - a.yzt *b.xy   + a.xyzt*b.y    ; 
c.yzt  =  + a.q   *b.yzt  + a.x   *b.xyzt + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   - a.xy  *b.xzt  + a.xz  *b.xyt  + a.yz  *b.t    - a.xt  *b.xyz  - a.yt  *b.z    + a.zt  *b.y    + a.xyz *b.xt   - a.xyt *b.xz   + a.xzt *b.xy   + a.yzt *b.q    - a.xyzt*b.x    ; 
c.xyzt =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.yz  *b.xt   + a.xt  *b.yz   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q    ; 

	return c;
}

///////////////////////////////////


//********************************* GA4E ***************************************

struct GA4E{	// 0 . 1 2 4 8 . 3 5 6 9 10 12 . 7 11 13 14 . 15
	ex q,  x,y,z,t, xy,xz,yz,xt,yt,zt,  xyz,xyt,xzt,yzt,  xyzt;
	GA4E() {q = 0; x = 0; y = 0; z = 0; t = 0; 
		xy = 0; xz = 0; yz = 0;  xt = 0; yt = 0; zt = 0;
		xyz = 0; xyt = 0; xzt = 0; yzt = 0; xyzt=0;}
	GA4E(ex qq, ex xx, ex yy, ex zz, ex tt, 
		ex xxyy, ex xxzz, ex yyzz, ex xxtt, ex yytt, ex zztt,
		ex xxyyzz, ex xxyytt, ex xxzztt, ex yyzztt,   ex xxyyzztt) 

		{q = qq; x = xx; y = yy; z = zz; t = tt;
		 xy = xxyy; xz = xxzz; yz = yyzz; xt = xxtt; yt = yytt; zt = zztt; 
		 xyz = xxyyzz;  xyt = xxyytt; xzt = xxzztt; yzt = yyzztt; xyzt=xxyyzztt;}
};


ostream &operator<<(ostream &ff, GA4E &v) {
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


void PrintGA4E(GA4E &v) {
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

GA4E Product_GA4E(const GA4E &a, const GA4E &b)
{

	GA4E c;

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


/////////////////////////////////// 

Mink Dual_Mink(Mink a, int i)
{
	Mink b;
	int SignArray[16];
	int j,k;

	k = 1;
	for (j=0;j<16;j++) {
		if ( (k & i) != 0) SignArray[j] = -1; else SignArray[j] = +1;
		k <<= 1;
	}

	b.q = SignArray[ 0]*a.xyzt;

	b.x = SignArray[ 1]*a.yzt;
	b.y = SignArray[ 2]*a.xzt;
	b.z = SignArray[ 3]*a.xyt;
	b.t = SignArray[ 4]*a.xyz;

	b.xy = SignArray[ 5]*a.zt;
	b.xz = SignArray[ 6]*a.yt;
	b.yz = SignArray[ 7]*a.xt;
	b.xt = SignArray[ 8]*a.yz;
	b.yt = SignArray[ 9]*a.xz;
	b.zt = SignArray[10]*a.xy;

	b.xyz = SignArray[11]*a.t;
	b.xyt = SignArray[12]*a.z;
	b.xzt = SignArray[13]*a.y;
	b.yzt = SignArray[14]*a.x;

	b.xyzt = SignArray[15]*a.q;


	return b;

}


/////////////////////////////////// 

Mink Complement_Mink(Mink a, int i)
{
	Mink b;
	int SignArray[16];
	int j,k;

	k = 1;
	for (j=0;j<16;j++) {
		if ( (k & i) != 0) SignArray[j] = -1; else SignArray[j] = +1;
		k <<= 1;
	}

	b.q = SignArray[ 0]*a.q;

	b.x = SignArray[ 1]*a.x;
	b.y = SignArray[ 2]*a.y;
	b.z = SignArray[ 3]*a.z;
	b.t = SignArray[ 4]*a.t;

	b.xy = SignArray[ 5]*a.xy;
	b.xz = SignArray[ 6]*a.xz;
	b.yz = SignArray[ 7]*a.yz;
	b.xt = SignArray[ 8]*a.xt;
	b.yt = SignArray[ 9]*a.yt;
	b.zt = SignArray[10]*a.zt;

	b.xyz = SignArray[11]*a.xyz;
	b.xyt = SignArray[12]*a.xyt;
	b.xzt = SignArray[13]*a.xzt;
	b.yzt = SignArray[14]*a.yzt;

	b.xyzt = SignArray[15]*a.xyzt;


	return b;

}


/////////////////////////////////// 

GA4E Complement_GA4E(GA4E a, int i)
{
	GA4E b;
	int SignArray[16];
	int j,k;

	k = 1;
	for (j=0;j<16;j++) {
		if ( (k & i) != 0) SignArray[j] = -1; else SignArray[j] = +1;
		k <<= 1;
	}

	b.q = SignArray[ 0]*a.xyzt;

	b.x = SignArray[ 1]*a.yzt;
	b.y = SignArray[ 2]*a.xzt;
	b.z = SignArray[ 3]*a.xyt;
	b.t = SignArray[ 4]*a.xyz;

	b.xy = SignArray[ 5]*a.zt;
	b.xz = SignArray[ 6]*a.yt;
	b.yz = SignArray[ 7]*a.xt;
	b.xt = SignArray[ 8]*a.yz;
	b.yt = SignArray[ 9]*a.xz;
	b.zt = SignArray[10]*a.xy;

	b.xyz = SignArray[11]*a.t;
	b.xyt = SignArray[12]*a.z;
	b.xzt = SignArray[13]*a.y;
	b.yzt = SignArray[14]*a.x;

	b.xyzt = SignArray[15]*a.q;


	return b;

}



/////////////////////////////////// 

GA4E GA4E_Unary(GA4E a, int i)
{
	GA4E b;
	int SignArray[16];
	int j,k;

	k = 1;
	for (j=0;j<16;j++) {
		if ( (k & i) != 0) SignArray[j] = -1; else SignArray[j] = +1;
		k <<= 1;
	}

	b.q = SignArray[ 0]*a.q;

	b.x = SignArray[ 1]*a.x;
	b.y = SignArray[ 2]*a.y;
	b.z = SignArray[ 3]*a.z;
	b.t = SignArray[ 4]*a.t;

	b.xy = SignArray[ 5]*a.xy;
	b.xz = SignArray[ 6]*a.xz;
	b.yz = SignArray[ 7]*a.yz;
	b.xt = SignArray[ 8]*a.xt;
	b.yt = SignArray[ 9]*a.yt;
	b.zt = SignArray[10]*a.zt;

	b.xyz = SignArray[11]*a.xyz;
	b.xyt = SignArray[12]*a.xyt;
	b.xzt = SignArray[13]*a.xzt;
	b.yzt = SignArray[14]*a.yzt;

	b.xyzt = SignArray[15]*a.xyzt;


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

ex Determinant_Mink(Mink A) {
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
	C = Product_Clifford(B,A);

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

Mink Conjugate(Mink a) {

	Mink u;

	u.q =  a.q;
	u.x = -a.x;
	u.y = -a.y;
	u.z = -a.z;
	u.t = -a.t;
	u.xy = -a.xy;
	u.xz = -a.xz;
	u.yz = -a.yz;
	u.xt = -a.xt;
	u.yt = -a.yt;
	u.zt = -a.zt;
	u.xyz = -a.xyz;
	u.xyt = -a.xyt;
	u.xzt = -a.xzt;
	u.yzt = -a.yzt;
	u.xyzt = -a.xyzt;

	return u;
}


Mink Adjugate(Mink V)
{
	Mink W;

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

Mink Reciprocal(Mink a)
{
	ex b;
	Mink c,d;

	b = Determinant_Mink(a);

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

ex Determinant_GA4E(GA4E A) {
/*


*/	
	GA4E B, C;
	ex a, b, c, d, e, s, det;

	B = Reverse(A);
	C = Product_GA4E(B,A);

	a = C.q;
	b = C.x;
	c = C.y;
	d = C.z;
	e = C.t;

	s = C.xyzt;
	det = expand(a*a - b*b - c*c - d*d - e*e - s*s);

	return(det);
}

//////////////////////////////////////////////////////

GA4E Conjugate(GA4E a) {

	GA4E u;

	u.q =  a.q;
	u.x = -a.x;
	u.y = -a.y;
	u.z = -a.z;
	u.t = -a.t;
	u.xy = -a.xy;
	u.xz = -a.xz;
	u.yz = -a.yz;
	u.xt = -a.xt;
	u.yt = -a.yt;
	u.zt = -a.zt;
	u.xyz = -a.xyz;
	u.xyt = -a.xyt;
	u.xzt = -a.xzt;
	u.yzt = -a.yzt;
	u.xyzt = -a.xyzt;

	return u;
}


//////////////////////////////////////////////////////

GA4E Adjugate(GA4E a) {

	GA4E u;

//Product_Clifford(Reverse(r),Conjugate(Product_Clifford(r,Reverse(r)))) = Adjugate(r);

	u = Product_GA4E(Reverse(a),Conjugate(Product_GA4E(a,Reverse(a))));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.t = expand(u.t);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xt = expand(u.xt);
	u.yt = expand(u.yt);
	u.zt = expand(u.zt);
	u.xyz = expand(u.xyz);
	u.xyt = expand(u.xyt);
	u.xzt = expand(u.xzt);
	u.yzt = expand(u.yzt);
	u.xyzt = expand(u.xyzt);

	return u;
}


GA4E Reciprocal(GA4E a)
{
	ex b;
	GA4E c,d;

	b = Determinant_GA4E(a);

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


///////////////////////////////// 


int main(void)
{
	int i, j, k, ZeroCount, SignArray[16];

	Mink r,s,t,u,w;
	GA4E R,S,T,U,W;

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d"),   a_t("e");
	symbol a_xy("f"),  a_xz("g"),  a_yz("h"),  a_xt("i"), a_yt("j"), a_zt("k");
	symbol a_xyz("l"), a_xyt("m"), a_xzt("n"), a_yzt("o");
	symbol a_xyzt("p");

	symbol b_q("A");
	symbol b_x("B") ,  b_y("C") ,  b_z("D"),   b_t("E");
	symbol b_xy("F"),  b_xz("G"),  b_yz("H"),  b_xt("I"), b_yt("J"), b_zt("K");
	symbol b_xyz("L"), b_xyt("M"), b_xzt("N"), b_yzt("O");
	symbol b_xyzt("P");

	ex det_ref, det_candidate, check;
	int det_matches;

///////////////////////////////// 

	r = Mink(a_q, a_x,a_y,a_z,a_t, a_xy,a_xz,a_yz,a_xt,a_yt,a_zt, a_xyz,a_xyt,a_xzt,a_yzt, a_xyzt);  
	s = Reverse(r);
	
	printf("r = ");	PrintMink(r); printf("\n");
	printf("s = ");	PrintMink(s); printf("\n");

	t = Product_Clifford(r,s);
	w = Product_Clifford(s,r);

	cout << "r*s = " << t << "\n";
	cout << "s*r = " << w << "\n";

	u = t - w;

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.t = expand(u.t);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xt = expand(u.xt);
	u.yt = expand(u.yt);
	u.zt = expand(u.zt);
	u.xyz = expand(u.xyz);
	u.xyt = expand(u.xyt);
	u.xzt = expand(u.xzt);
	u.yzt = expand(u.yzt);
	u.xyzt = expand(u.xyzt);

	if (u==Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0)) 
		cout << "r and s=Reverse(r) commute"  << endl;	
	else cout << "r and s do not commute.\nu = r*s - s*r = " << u << "\n";

	t = Product_Clifford(r,s);

//	w = Mink(A,  -B,-C,-D,-E, 0,0,0,0,0,0,  0,0,0,0, -P)

	w.q = t.q;

	w.x = -t.x;
	w.y = -t.y;
	w.z = -t.z;
	w.t = -t.t;

	w.xy = 0;
	w.xz = 0;
	w.yz = 0;
	w.xt = 0;
	w.yt = 0;
	w.zt = 0;

	w.xyz = 0;
	w.xyt = 0;
	w.xzt = 0;
	w.yzt = 0;

	w.xyzt = -t.xyzt;


	u = Product_Clifford(t,w) - Product_Clifford(w,t);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.t = expand(u.t);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xt = expand(u.xt);
	u.yt = expand(u.yt);
	u.zt = expand(u.zt);
	u.xyz = expand(u.xyz);
	u.xyt = expand(u.xyt);
	u.xzt = expand(u.xzt);
	u.yzt = expand(u.yzt);
	u.xyzt = expand(u.xyzt);

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
		cout << "w and t commute"  << endl;	
	else cout << "u = w*t - t*w = " << u << "\n";



///////////////////////////////// 

/*



r = Mink(a, b,c,d,e, f,g,h,i,j,k, l,m,n,o, p)

Form the Reverse Conjugate, (index 32736 in my codes below.)

s = Mink(a, b,c,d,e, -f,-g,-h,-i,-j,-k, -l,-m,-n,-o, p));  

Form Clifford product of r and s. Ten terms disappear in the product.

t = r*s = Mink(A,  B,C,D,E, 0,0,0,0,0,0,  0,0,0,0, P)

A = g^2+l^2+b^2-j^2-o^2-e^2+h^2-m^2+c^2-p^2+f^2-k^2-i^2-n^2+a^2+d^2
B = -2*e*i-2*n*k+2*g*d-2*o*p+2*f*c+2*h*l-2*m*j+2*a*b
C = -2*k*o+2*m*i+2*a*c+2*n*p-2*g*l-2*j*e-2*b*f+2*d*h
D =  2*j*o+2*n*i-2*m*p-2*h*c-2*k*e-2*g*b+2*d*a+2*l*f
E = -2*d*k+2*m*f+2*h*o+2*a*e+2*g*n-2*j*c-2*l*p-2*b*i
P = -2*b*o+2*j*g+2*n*c+2*a*p-2*m*d-2*h*i-2*k*f+2*e*l

We now form a complex conjugate-like term from t

w = Mink(A,  -B,-C,-D,-E, 0,0,0,0,0,0,  0,0,0,0, -P)

Multiply t and w

t*w = Mink(A^2-B^2-C^2-D^2+E^2+P^2,   0,0,0,0,   0,0,0,0,0,0,   0,0,0,0, 0)



Start with   t*w = w*t = det(r) = |r|^4

w  = det(r) *(1/t) = det(r)*(1/(r*s))

w*s = det(r)*(1/r)

(1/r) = (w*s)/det(r)


*/
	r = Mink(a_q, a_x,a_y,a_z,a_t, a_xy,a_xz,a_yz,a_xt,a_yt,a_zt, a_xyz,a_xyt,a_xzt,a_yzt, a_xyzt);  
	s = Reverse(r);
	
	printf("r = ");	PrintMink(r); printf("\n");
	printf("s = ");	PrintMink(s); printf("\n");

	t = Product_Clifford(r,s);

//	w = Mink(A,  -B,-C,-D,-E, 0,0,0,0,0,0,  0,0,0,0, -P)

	w.q = t.q;

	w.x = -t.x;
	w.y = -t.y;
	w.z = -t.z;
	w.t = -t.t;

	w.xy = 0;
	w.xz = 0;
	w.yz = 0;
	w.xt = 0;
	w.yt = 0;
	w.zt = 0;

	w.xyz = 0;
	w.xyt = 0;
	w.xzt = 0;
	w.yzt = 0;

	w.xyzt = -t.xyzt;

//	w = det(r)*(1/s)*(1/r)   => s*w = det(r)*(1/r)
	// (w*s)*r = det(r)

	u = Product_Clifford(r,Product_Clifford(s,w));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.t = expand(u.t);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xt = expand(u.xt);
	u.yt = expand(u.yt);
	u.zt = expand(u.zt);
	u.xyz = expand(u.xyz);
	u.xyt = expand(u.xyt);
	u.xzt = expand(u.xzt);
	u.yzt = expand(u.yzt);
	u.xyzt = expand(u.xyzt);

	cout << "u = " << u << "\nExpect scalar = det(r)\n\n";

	check = expand(Determinant_Mink(r) - u.q);
	cout << "check = " << check << " expect zero.\n";

///////////////////////////////// Repeat for GA4E

	R = GA4E(a_q, a_x,a_y,a_z,a_t, a_xy,a_xz,a_yz,a_xt,a_yt,a_zt, a_xyz,a_xyt,a_xzt,a_yzt, a_xyzt);  
	S.q =  R.q;

	S.x = R.x;
	S.y = R.y;
	S.z = R.z;
	S.t = R.t;

	S.xy = -R.xy;
	S.xz = -R.xz;
	S.yz = -R.yz;
	S.xt = -R.xt;
	S.yt = -R.yt;
	S.zt = -R.zt;

	S.xyz = -R.xyz;
	S.xyt = -R.xyt;
	S.xzt = -R.xzt;
	S.yzt = -R.yzt;

	S.xyzt =  R.xyzt;
	
	printf("R = ");	PrintGA4E(R); printf("\n");
	printf("S = ");	PrintGA4E(S); printf("\n");

	T = Product_GA4E(R,S);

//	W = GA4E(A,  -B,-C,-D,-E, 0,0,0,0,0,0,  0,0,0,0, -P)

	W.q = T.q;

	W.x = -T.x;
	W.y = -T.y;
	W.z = -T.z;
	W.t = -T.t;

	W.xy = 0;
	W.xz = 0;
	W.yz = 0;
	W.xt = 0;
	W.yt = 0;
	W.zt = 0;

	W.xyz = 0;
	W.xyt = 0;
	W.xzt = 0;
	W.yzt = 0;

	W.xyzt = -T.xyzt;

//	w = det(r)*(1/s)*(1/r)   => s*w = det(r)*(1/r)
	// (w*s)*r = det(r)

	U = Product_GA4E(R,Product_GA4E(S,W));

	U.q = expand(U.q);
	U.x = expand(U.x);
	U.y = expand(U.y);
	U.z = expand(U.z);
	U.t = expand(U.t);
	U.xy = expand(U.xy);
	U.xz = expand(U.xz);
	U.yz = expand(U.yz);
	U.xt = expand(U.xt);
	U.yt = expand(U.yt);
	U.zt = expand(U.zt);
	U.xyz = expand(U.xyz);
	U.xyt = expand(U.xyt);
	U.xzt = expand(U.xzt);
	U.yzt = expand(U.yzt);
	U.xyzt = expand(U.xyzt);

	cout << "U = " << U << "\nExpect scalar = det(r)\n\n";

	check = expand(Determinant_GA4E(R) - U.q);
	cout << "check = " << check << "\n";

/////////////////////////////////

	R = GA4E(b_q, b_x,b_y,b_z,b_t, 0,0,0,0,0,0, 0,0,0,0, b_xyzt); 
	S = GA4E(b_q, -b_x,-b_y,-b_z,-b_t, 0,0,0,0,0,0, 0,0,0,0, -b_xyzt); 

	U = Product_GA4E(R,S);

	U.q = expand(U.q);
	U.x = expand(U.x);
	U.y = expand(U.y);
	U.z = expand(U.z);
	U.t = expand(U.t);
	U.xy = expand(U.xy);
	U.xz = expand(U.xz);
	U.yz = expand(U.yz);
	U.xt = expand(U.xt);
	U.yt = expand(U.yt);
	U.zt = expand(U.zt);
	U.xyz = expand(U.xyz);
	U.xyt = expand(U.xyt);
	U.xzt = expand(U.xzt);
	U.yzt = expand(U.yzt);
	U.xyzt = expand(U.xyzt);

	cout << "U = " << U << "\nExpect scalar = det(R)\n\n";


/////////////////////////////////


/*	Mink X, Y, Z;

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	Y = Mink( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	Z = Mink(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 
*/

	r = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Reciprocal(r);
	t = Product_Clifford(r,s);

	cout << "r = " << r << "\n";
	cout << "det(r) = << " << Determinant_Mink(r) << "\n";
	cout << "s = (1/r) = " << s << "\n";
	cout << "t = r*(1/r) = " << t << "\n";

//////////////////////////////////

//test mink adjugate

	r = Mink(a_q, a_x,a_y,a_z,a_t, a_xy,a_xz,a_yz,a_xt,a_yt,a_zt, a_xyz,a_xyt,a_xzt,a_yzt, a_xyzt);  

	u = Product_Clifford(Reverse(r),Conjugate(Product_Clifford(r,Reverse(r)))) - Adjugate(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.t = expand(u.t);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xt = expand(u.xt);
	u.yt = expand(u.yt);
	u.zt = expand(u.zt);
	u.xyz = expand(u.xyz);
	u.xyt = expand(u.xyt);
	u.xzt = expand(u.xzt);
	u.yzt = expand(u.yzt);
	u.xyzt = expand(u.xyzt);

	cout << "u = " << u << "\nExpect 0\n\n";


/////////////////////////////////

// test GA4E reciprocol, determinant, adjugate

	R = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	S = Reciprocal(R);
	T = Product_GA4E(R,S);

	cout << "R = " << R << "\n";
	cout << "det(R) = << " << Determinant_GA4E(R) << "\n";
	cout << "S = (1/R) = " << S << "\n";
	cout << "T = R*(1/R) = " << T << "\n";

	U = Adjugate(R);
	cout << "U = Adjugate(r) = " << U << "\n\n";

/////////////////////////////////

	return 0;
}

