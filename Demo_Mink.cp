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

int SignArray[16];


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
	return ff << "(" 
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




Mink operator*(const Mink &u, const Mink &v) {
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


Mink Clifford(const Mink &a, const Mink &b)
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

Mink CliffordViaPython(const Mink &a, const Mink &b)
{

	Mink c;

c.q    = + (a.q*b.q - a.t*b.t + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt);

c.x    = + (a.q*b.x + a.t*b.xt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt);
c.y    = + (a.q*b.y + a.t*b.yt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt);
c.z    = + (a.q*b.z + a.t*b.zt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t);
c.t    = + (a.q*b.t + a.t*b.q + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z);

c.xy   = + (a.q*b.xy - a.t*b.xyt + a.x*b.y + a.xt*b.yt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt - a.xz*b.yz + a.xzt*b.yzt - a.y*b.x - a.yt*b.xt + a.yz*b.xz - a.yzt*b.xzt + a.z*b.xyz + a.zt*b.xyzt);
c.xz   = + (a.q*b.xz - a.t*b.xzt + a.x*b.z + a.xt*b.zt + a.xy*b.yz - a.xyt*b.yzt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt - a.yz*b.xy + a.yzt*b.xyt - a.z*b.x - a.zt*b.xt);
c.xt   = + (a.q*b.xt - a.t*b.x + a.x*b.t + a.xt*b.q + a.xy*b.yt - a.xyt*b.y - a.xyz*b.yzt - a.xyzt*b.yz + a.xz*b.zt - a.xzt*b.z - a.y*b.xyt - a.yt*b.xy - a.yz*b.xyzt + a.yzt*b.xyz - a.z*b.xzt - a.zt*b.xz);
c.yz   = + (a.q*b.yz - a.t*b.yzt + a.x*b.xyz + a.xt*b.xyzt - a.xy*b.xz + a.xyt*b.xzt + a.xyz*b.x + a.xyzt*b.xt + a.xz*b.xy - a.xzt*b.xyt + a.y*b.z + a.yt*b.zt + a.yz*b.q - a.yzt*b.t - a.z*b.y - a.zt*b.yt);
c.yt   = + (a.q*b.yt - a.t*b.y + a.x*b.xyt + a.xt*b.xy - a.xy*b.xt + a.xyt*b.x + a.xyz*b.xzt + a.xyzt*b.xz + a.xz*b.xyzt - a.xzt*b.xyz + a.y*b.t + a.yt*b.q + a.yz*b.zt - a.yzt*b.z - a.z*b.yzt - a.zt*b.yz);
c.zt   = + (a.q*b.zt - a.t*b.z + a.x*b.xzt + a.xt*b.xz - a.xy*b.xyzt + a.xyt*b.xyz - a.xyz*b.xyt - a.xyzt*b.xy - a.xz*b.xt + a.xzt*b.x + a.y*b.yzt + a.yt*b.yz - a.yz*b.yt + a.yzt*b.y + a.z*b.t + a.zt*b.q);

c.xyz  = + (a.q*b.xyz + a.t*b.xyzt + a.x*b.yz - a.xt*b.yzt + a.xy*b.z + a.xyt*b.zt + a.xyz*b.q - a.xyzt*b.t - a.xz*b.y - a.xzt*b.yt - a.y*b.xz + a.yt*b.xzt + a.yz*b.x + a.yzt*b.xt + a.z*b.xy - a.zt*b.xyt);
c.xyt  = + (a.q*b.xyt + a.t*b.xy + a.x*b.yt - a.xt*b.y + a.xy*b.t + a.xyt*b.q + a.xyz*b.zt - a.xyzt*b.z - a.xz*b.yzt - a.xzt*b.yz - a.y*b.xt + a.yt*b.x + a.yz*b.xzt + a.yzt*b.xz + a.z*b.xyzt - a.zt*b.xyz);
c.xzt  = + (a.q*b.xzt + a.t*b.xz + a.x*b.zt - a.xt*b.z + a.xy*b.yzt + a.xyt*b.yz - a.xyz*b.yt + a.xyzt*b.y + a.xz*b.t + a.xzt*b.q - a.y*b.xyzt + a.yt*b.xyz - a.yz*b.xyt - a.yzt*b.xy - a.z*b.xt + a.zt*b.x);
c.yzt  = + (a.q*b.yzt + a.t*b.yz + a.x*b.xyzt - a.xt*b.xyz - a.xy*b.xzt - a.xyt*b.xz + a.xyz*b.xt - a.xyzt*b.x + a.xz*b.xyt + a.xzt*b.xy + a.y*b.zt - a.yt*b.z + a.yz*b.t + a.yzt*b.q - a.z*b.yt + a.zt*b.y);

c.xyzt = + (a.q*b.xyzt - a.t*b.xyz + a.x*b.yzt + a.xt*b.yz + a.xy*b.zt - a.xyt*b.z + a.xyz*b.t + a.xyzt*b.q - a.xz*b.yt + a.xzt*b.y - a.y*b.xzt - a.yt*b.xz + a.yz*b.xt - a.yzt*b.x + a.z*b.xyt + a.zt*b.xy);


	return c;
}


Mink WedgeViaPython(const Mink &a, const Mink &b)
{

	Mink c;

c.q    = + (a.q*b.q);

c.x    = + (a.q*b.x + a.x*b.q);
c.y    = + (a.q*b.y + a.y*b.q);
c.z    = + (a.q*b.z + a.z*b.q);
c.t    = + (a.q*b.t + a.t*b.q);

c.xy   = + (a.q*b.xy + a.x*b.y + a.xy*b.q - a.y*b.x);
c.xz   = + (a.q*b.xz + a.x*b.z + a.xz*b.q - a.z*b.x);
c.xt   = + (a.q*b.xt - a.t*b.x + a.x*b.t + a.xt*b.q);
c.yz   = + (a.q*b.yz + a.y*b.z + a.yz*b.q - a.z*b.y);
c.yt   = + (a.q*b.yt - a.t*b.y + a.y*b.t + a.yt*b.q);
c.zt   = + (a.q*b.zt - a.t*b.z + a.z*b.t + a.zt*b.q);

c.xyz  = + (a.q*b.xyz + a.x*b.yz + a.xy*b.z + a.xyz*b.q - a.xz*b.y - a.y*b.xz + a.yz*b.x + a.z*b.xy);
c.xyt  = + (a.q*b.xyt + a.t*b.xy + a.x*b.yt - a.xt*b.y + a.xy*b.t + a.xyt*b.q - a.y*b.xt + a.yt*b.x);
c.xzt  = + (a.q*b.xzt + a.t*b.xz + a.x*b.zt - a.xt*b.z + a.xz*b.t + a.xzt*b.q - a.z*b.xt + a.zt*b.x);
c.yzt  = + (a.q*b.yzt + a.t*b.yz + a.y*b.zt - a.yt*b.z + a.yz*b.t + a.yzt*b.q - a.z*b.yt + a.zt*b.y);

c.xyzt = + (a.q*b.xyzt - a.t*b.xyz + a.x*b.yzt + a.xt*b.yz + a.xy*b.zt - a.xyt*b.z + a.xyz*b.t + a.xyzt*b.q - a.xz*b.yt + a.xzt*b.y - a.y*b.xzt - a.yt*b.xz + a.yz*b.xt - a.yzt*b.x + a.z*b.xyt + a.zt*b.xy);

	return c;
}


Mink InnerViaPython(const Mink &a, const Mink &b)
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

	return c;
}



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


	return c;
}


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

	return c;
}



Mink operator^(const Mink &u, const Mink &v) {
	Mink w;

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


Mink Wedge(const Mink &u, const Mink &v) {
	Mink w;

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


Mink AntiWedge(const Mink a, const Mink b)
{
	Mink c;

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

Mink OverBar(Mink a)   // Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	Mink b;		// a blade wedge b blade = pseudovector blade

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


Mink UnderBar(Mink a)   // Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	Mink b;		// b blade wedge a blade = pseudovector blade

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


Mink Involution(Mink w)	// Corresponds to PT parity transform: x -> -x, y -> -y, z ->-z, t -> -t
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


Mink Regressive(Mink a, Mink b) 
{
	Mink c;

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



Mink RegressiveViaFormula(Mink a, Mink b)
{
	Mink c;

	Mink I,I_inv,d,e,f;
	
	I.xyzt = 1;	I_inv.xyzt = -1;

	d = a*I_inv^b*I_inv;
	c = d*I;

	return c;
}


Mink LowerRightViaFormula(Mink a, Mink b)
{

// -Wedge(Blade[i]*xyzt,xyzt*Blade[j])

	Mink c;

	Mink I,I_inv,d,e,f;
	
	I.xyzt = 1;	I_inv.xyzt = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}



Mink Expander(const Mink &a, const Mink &b)
{

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

	return c;
}


Mink Conserver(const Mink &a, const Mink &b)
{

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


	return c;
}



Mink Shrinker(const Mink &a, const Mink &b)
{

	Mink c;

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

	return c;
}


///////////////////////////////////////////////////////


Mink SymmetricViaFormula(Mink a, Mink b)
{
	Mink  c,d,e,f;
	
	c = (a*b + b*a)/2;

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


Mink AntiSymmetricViaFormula(Mink a, Mink b)
{
	Mink  c,d,e,f;
	
	c = (a*b - b*a)/2;

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


Mink Symmetric(const Mink &a, const Mink &b)
{

	Mink c;

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

	return c;
}


Mink AntiSymmetric(const Mink &a, const Mink &b)
{

	Mink c;

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

	return c;
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

///////////////////////////////////

void PrintSimpleBlade(Mink a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q    ");
	if(a.q == -1) printf("-q    ");
	if(a.q ==  0) Count++;


	if(a.x ==  1) printf(" x    ");
	if(a.x == -1) printf("-x    ");
	if(a.x ==  0) Count++;

	if(a.y ==  1) printf(" y    ");
	if(a.y == -1) printf("-y    ");
	if(a.y ==  0) Count++;

	if(a.z ==  1) printf(" z    ");
	if(a.z == -1) printf("-z    ");
	if(a.z ==  0) Count++;

	if(a.t ==  1) printf(" t    ");
	if(a.t == -1) printf("-t    ");
	if(a.t ==  0) Count++;


	if(a.xy ==  1) printf(" xy   ");
	if(a.xy == -1) printf("-xy   ");
	if(a.xy ==  0) Count++;

	if(a.xz ==  1) printf(" xz   ");
	if(a.xz == -1) printf("-xz   ");
	if(a.xz ==  0) Count++;

	if(a.yz ==  1) printf(" yz   ");
	if(a.yz == -1) printf("-yz   ");
	if(a.yz ==  0) Count++;

	if(a.xt ==  1) printf(" xt   ");
	if(a.xt == -1) printf("-xt   ");
	if(a.xt ==  0) Count++;

	if(a.yt ==  1) printf(" yt   ");
	if(a.yt == -1) printf("-yt   ");
	if(a.yt ==  0) Count++;

	if(a.zt ==  1) printf(" zt   ");
	if(a.zt == -1) printf("-zt   ");
	if(a.zt ==  0) Count++;


	if(a.xyz ==  1) printf(" xyz  ");
	if(a.xyz == -1) printf("-xyz  ");
	if(a.xyz ==  0) Count++;

	if(a.xyt ==  1) printf(" xyt  ");
	if(a.xyt == -1) printf("-xyt  ");
	if(a.xyt ==  0) Count++;

	if(a.xzt ==  1) printf(" xzt  ");
	if(a.xzt == -1) printf("-xzt  ");
	if(a.xzt ==  0) Count++;

	if(a.yzt ==  1) printf(" yzt  ");
	if(a.yzt == -1) printf("-yzt  ");
	if(a.yzt ==  0) Count++;


	if(a.xyzt ==  1) printf(" xyzt ");
	if(a.xyzt == -1) printf("-xyzt ");
	if(a.xyzt ==  0) Count++;

	if(Count == 16) printf(" 0    ");
	if(Count < 15)  printf("error ");


}

///////////////////////////////////

int	GetWeight(Mink r) // assume only one non-zero component
{
	int weight = 0;		// default

	if (r.q ==  1) weight =  1;
	if (r.q == -1) weight = -1;

	if (r.x ==  1) weight =  1;
	if (r.x == -1) weight = -1;
	if (r.y ==  1) weight =  1;
	if (r.y == -1) weight = -1;
	if (r.z ==  1) weight =  1;
	if (r.z == -1) weight = -1;
	if (r.t ==  1) weight =  1;
	if (r.t == -1) weight = -1;

	if (r.xy ==  1) weight =  1;
	if (r.xy == -1) weight = -1;
	if (r.xz ==  1) weight =  1;
	if (r.xz == -1) weight = -1;
	if (r.yz ==  1) weight =  1;
	if (r.yz == -1) weight = -1;
	if (r.xt ==  1) weight =  1;
	if (r.xt == -1) weight = -1;
	if (r.yt ==  1) weight =  1;
	if (r.yt == -1) weight = -1;
	if (r.zt ==  1) weight =  1;
	if (r.zt == -1) weight = -1;

	if (r.xyz ==  1) weight =  1;
	if (r.xyz == -1) weight = -1;
	if (r.xyt ==  1) weight =  1;
	if (r.xyt == -1) weight = -1;
	if (r.xzt ==  1) weight =  1;
	if (r.xzt == -1) weight = -1;
	if (r.yzt ==  1) weight =  1;
	if (r.yzt == -1) weight = -1;

	if (r.xyzt ==  1) weight =  1;
	if (r.xyzt == -1) weight = -1;

	return weight;
}

///////////////////////////////////

int	Rank[16]  = {0, 1,1,1,1, 2,2,2,2,2,2, 3,3,3,3,  4} ;
long	Index[16] = {0,1,2,4,8,3,5,6,9,10,12,7,11,13,14,15};	// organized by grade
int	CountSetBits[16] = {0, 1,1,2,1, 2,2,3,1,2,2, 3,2,3,3, 4};

int	Higher(int a, int b)	// return higher of a or b
{
	int c;
	
	c = a;
	if (b>a) c = b;

	return c;

}

///////////////////////////////////

Mink DorstDual(Mink a)
{
	Mink b, I_inv;

	I_inv.xyzt = -1;

	b = LeftContraction(a,I_inv);
	return b;

}


///////////////////////////////////

Mink DorstUnDual(Mink a)
{
	Mink b, I;

	I.xyzt = 1;

	b = LeftContraction(a,I);
	return b;

}

///////////////////////////////////


int main(void)
{

//	long SwapArray[16]={0,1,3,2,6,7,5,4,12,13,15,14,10,11,9,8};	// gray code
//	long SwapArray[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	// straight count
//	long SwapArray[16]={0,3,5,6,9,10,12,15,1,2,4,7,8,11,13,14};	// even subset first
	long SwapArray[16]={0,1,2,4,8,3,5,6,9,10,12,7,11,13,14,15};	// organized by grade
	
	int i,j,k;

/*	char Field[17][10];

	for (i=0;i<17;i++) for (j=0;j<10;j++) Field[i][j] = 0;
	strcpy(Field[0], "q   ");
	strcpy(Field[1], "x   ");
	strcpy(Field[2], "y   ");
	strcpy(Field[3], "xy  ");
	strcpy(Field[4], "z   ");
	strcpy(Field[5], "xz  ");
	strcpy(Field[6], "yz  ");
	strcpy(Field[7], "xyz ");
	strcpy(Field[8], "t   ");
	strcpy(Field[9], "xt  ");
	strcpy(Field[10],"yt  ");
	strcpy(Field[11],"xyt ");
	strcpy(Field[12],"zt  ");
	strcpy(Field[13],"xzt ");
	strcpy(Field[14],"yzt ");
	strcpy(Field[15],"xyzt");
	strcpy(Field[16],"0   ");
*/
	int mm;


	Mink first;
	Mink second;
	Mink third;
	Mink check;
	Mink Zero;   Zero = Mink();

/*	symbol a("a");
	symbol b("b"), c("c"), d("d"), e("e");
	symbol f("f"), g("g"), h("h"), j("j"), k("k"), l("l");
	symbol m("m"), n("n"), p("p"), r("r");
	symbol s("s");

	symbol A("A");
	symbol B("B"), C("C"), D("D"), E("E");
	symbol F("F"), G("G"), H("H"), J("J"), K("K"), L("L");
	symbol M("M"), N("N"), P("P"), R("R");
	symbol S("S");
*/

	Mink q, x,y,z,t, xy,xz,yz,xt,yt,zt, xyz,xyt,xzt,yzt, xyzt;
	Mink r,s,u;

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

	symbol c_q("c.a");
	symbol c_x("c.b") ,  c_y("c.c") ,  c_z("c.d"),   c_t("c.e");
	symbol c_xy("c.f"),  c_xz("c.g"),  c_yz("c.h"),  c_xt("c.i"), c_yt("c.j"), c_zt("c.k");
	symbol c_xyz("c.l"), c_xyt("c.m"), c_xzt("c.n"), c_yzt("c.o");
	symbol c_xyzt("c.p");

	symbol d_q("d.a");
	symbol d_x("d.b") ,  d_y("d.c") ,  d_z("d.d"),   d_t("d.e");
	symbol d_xy("d.f"),  d_xz("d.g"),  d_yz("d.h"),  d_xt("d.i"), d_yt("d.j"), d_zt("d.k");
	symbol d_xyz("d.l"), d_xyt("d.m"), d_xzt("d.n"), d_yzt("d.o");
	symbol d_xyzt("d.p");


/////////////////////////////////

//	Test OverBar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = OverBar(r);
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

	cout << "r = " << r << "\n";
	cout << "OverBar(r) = " << u << "\n";
	cout << "\n";



/////////////////////////////////

//	Test Double OverBar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = OverBar(r);
	u = OverBar(s);

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

	cout << "r = " << r << "\n";
	cout << "OverBar(OverBar(r))) = " << u << "\n";
	cout << "\n";


/////////////////////////////////

//	Test UnderBar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = UnderBar(r);
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

	cout << "r = " << r << "\n";
	cout << "UnderBar(r) = " << u << "\n";
	cout << "\n";



/////////////////////////////////

//	Test Double UnderBar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = UnderBar(r);
	u = UnderBar(s);

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

	cout << "r = " << r << "\n";
	cout << "UnderBar(UnderBar(r)) = " << u << "\n";
	cout << "\n";


/////////////////////////////////

//	Test Over Under Bar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = OverBar(r);
	u = UnderBar(s);

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

	cout << "r = " << r << "\n";
	cout << "OverBar(UnderBar(r)) = " << u << "\n";
	cout << "\n";


/////////////////////////////////

//	Test Under Over Bar

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = UnderBar(r);
	u = OverBar(s);

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

	cout << "r = " << r << "\n";
	cout << "UnderBar(OverBar(r)) = " << u << "\n";
	cout << "\n";


////////////////////////////////////

//	Test AntiWedge Equations

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	t = OverBar((UnderBar(r)^UnderBar(s)));
	u = AntiWedge(r,s) - t;

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Formula is Consistent with OverBar((UnderBar(r)^UnderBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test AntiWedge Equations

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	t = UnderBar((OverBar(r)^OverBar(s)));
	u = AntiWedge(r,s) - t;

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Formula is Consistent with UnderBar((OverBar(r)^OverBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test Regressive Equations

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Regressive(r,s) - RegressiveViaFormula(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Regressive is Consistent with RegressiveViaFormula"  << endl;	
	}else cout << "Check RegressiveProductViaEquations \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Clifford

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = (r*s)*t - r*(s*t);
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Clifford Product is Associative "  << endl;	
	}else cout << "Clifford Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Wedge

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = ((r^s)^t) - (r^(s^t));
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Wedge Product is Associative "  << endl;	
	}else cout << "Wedge Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of AntiWedge

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = AntiWedge(AntiWedge(r,s),t) - AntiWedge(r,AntiWedge(s,t));
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Product is Associative "  << endl;	
	}else cout << "AntiWedge Product is Non-associative \n";
	cout << "\n";


///////////////////////////////////////

// print multiplication tables
//	ex q,  x,y,z,t, xy,xz,yz,xt,yt,zt,  xyz,xyt,xzt,yzt,  xyzt;

	q    = Mink(1, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);

	x    = Mink(0, 1,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	y    = Mink(0, 0,1,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	z    = Mink(0, 0,0,1,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	t    = Mink(0, 0,0,0,1, 0,0,0,0,0,0, 0,0,0,0, 0);

	xy   = Mink(0, 0,0,0,0, 1,0,0,0,0,0, 0,0,0,0, 0);
	xz   = Mink(0, 0,0,0,0, 0,1,0,0,0,0, 0,0,0,0, 0);
	yz   = Mink(0, 0,0,0,0, 0,0,1,0,0,0, 0,0,0,0, 0);
	xt   = Mink(0, 0,0,0,0, 0,0,0,1,0,0, 0,0,0,0, 0);
	yt   = Mink(0, 0,0,0,0, 0,0,0,0,1,0, 0,0,0,0, 0);
	zt   = Mink(0, 0,0,0,0, 0,0,0,0,0,1, 0,0,0,0, 0);

	xyz  = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 1,0,0,0, 0);
	xyt  = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,1,0,0, 0);
	xzt  = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,0, 0);
	yzt  = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,1, 0);

	xyzt = Mink(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 1);

	Mink Blade[16];

	Blade[ 0] = q;

	Blade[ 1] = x;
	Blade[ 2] = y;
	Blade[ 3] = z;
	Blade[ 4] = t;

	Blade[ 5] = xy;
	Blade[ 6] = xz;
	Blade[ 7] = yz;
	Blade[ 8] = xt;
	Blade[ 9] = yt;
	Blade[10] = zt;

	Blade[11] = xyz;
	Blade[12] = xyt;
	Blade[13] = xzt;
	Blade[14] = yzt;

	Blade[15] = xyzt;

//	print geometric product multiplication table

	printf(" *    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade(Blade[i]*Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Wedge product multiplication table

	printf(" ^    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade(Blade[i]^Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print AntiWedge product multiplication table

	printf("Lengyel\n");
	printf(" V    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade(AntiWedge(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Regressive product multiplication table

	printf("Hestenes\n");
	printf(" V    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade(Regressive(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print Clifford - Wedge

	printf("Blade[i]*Blade[j] - (Blade[i]^Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Blade[i]*Blade[j] - (Blade[i]^Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print ?? Contraction product multiplication table

	printf("AntiWedge(OverBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = AntiWedge(OverBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(UnderBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = AntiWedge(UnderBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = AntiWedge(Blade[i],OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = AntiWedge(Blade[i],UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(OverBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(OverBar(Blade[i]),UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(UnderBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(UnderBar(Blade[i]),UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////   Now we try various Hestenes forms


///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(Blade[i]*xyzt,Blade[j]*xyzt)\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(Blade[i]*xyzt,Blade[j]*xyzt);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(Blade[i]*xyzt,xyzt*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(Blade[i]*xyzt,xyzt*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyzt*Blade[i],Blade[j]*xyzt)\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(xyzt*Blade[i],Blade[j]*xyzt);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyzt*Blade[i],xyzt*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Wedge(xyzt*Blade[i],xyzt*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////////// -Wedge(Blade[i]*xyzt,xyzt*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula

	printf("Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


/////////////////////////////////////// -Wedge(Blade[i]*xyzt,xyzt*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula - Wedge

	printf("Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			u = Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////

//	Test Associativity of LowerRightViaFormula

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = LowerRightViaFormula(LowerRightViaFormula(r,s),t) - LowerRightViaFormula(r,LowerRightViaFormula(s,t));
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "LowerRightViaFormula Product is Associative "  << endl;	
	}else cout << "LowerRightViaFormula Product is Non-associative \n";
	cout << "\n";

////////////////////////////////

//	print geometric product increased rank terms

	printf("Terms with increased rank\n");
	printf(" <    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			k = CountSetBits[Index[i]^Index[j]];	// clifford product index
			if (k > Higher(Rank[i],Rank[j])) PrintSimpleBlade(Blade[i]*Blade[j]); 
					else PrintSimpleBlade(Zero);
		}
		printf("\n");
	}
	printf("\n\n\n");



////////////////////////////////

//	print geometric product equal rank terms

	printf("Terms with preserved rank\n");
	printf(" <    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			k = CountSetBits[Index[i]^Index[j]];	// clifford product index
			if (k == Higher(Rank[i],Rank[j])) PrintSimpleBlade(Blade[i]*Blade[j]); 
					else PrintSimpleBlade(Zero);
		}
		printf("\n");
	}
	printf("\n\n\n");



////////////////////////////////

//	print geometric product reduced rank terms

	printf("Terms with reduced rank\n");
	printf(" <    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			k = CountSetBits[Index[i]^Index[j]];	// clifford product index
			if (k < Higher(Rank[i],Rank[j])) PrintSimpleBlade(Blade[i]*Blade[j]); 
					else PrintSimpleBlade(Zero);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////////////////


//	long SwapArray[16]={0, 1,2,4,8, 3,5,6,9,10,12, 7,11,13,14, 15};	// organized by grade
	int NumTerms = 16;
	int EqnBasis, ii, jj, kk, bitmap, weight;
	
	char Field[17][10];

//	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

	for (i=0;i<17;i++) for (j=0;j<10;j++) Field[i][j] = 0;	// organized by grade
	strcpy(Field[0], "q   ");

	strcpy(Field[1], "x   ");
	strcpy(Field[2], "y   ");
	strcpy(Field[3], "z   ");
	strcpy(Field[4], "t   ");

	strcpy(Field[5], "xy  ");
	strcpy(Field[6], "xz  ");
	strcpy(Field[7], "yz  ");
	strcpy(Field[8], "xt  ");
	strcpy(Field[9], "yt  ");
	strcpy(Field[10],"zt  ");

	strcpy(Field[11],"xyz ");
	strcpy(Field[12],"xyt ");
	strcpy(Field[13],"xzt ");
	strcpy(Field[14],"yzt ");

	strcpy(Field[15],"xyzt");

	strcpy(Field[16],"0   ");



/////////////////////////////////


// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format

	printf("Clifford equation set for test purposes\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					r = Blade[i]*Blade[j];
					weight = GetWeight(r);
		//			if(weight == 0) printf("   0      ");
					if(weight  < 0) printf(" - "); 
					if(weight >  0) printf(" + ");
					if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////

//	Test Clifford equation set

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = r*s - Clifford(r,s);;
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Clifford products agree "  << endl;	
	}else cout << "Clifford equations are not correct \n";
	cout << "\n";

	cout << "u = " << u << "\n";


/////////////////////////////////

// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format

	printf("Expander equation set for test purposes\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[i]^Index[j]] > Higher(Rank[i],Rank[j])  ) {
						r = Blade[i]*Blade[j];
						weight = GetWeight(r);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////

// wedge - expander

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Wedge(r,s) - Expander(r,s);;
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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Wedge and expander products agree "  << endl;	
	}else cout << "Wedge - Expander differ (as expected)\n";
	cout << "\n";

	cout << "u = " << u << "\n\n\n";



/////////////////////////////////

//	Test Associativity of Expander

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Expander(Expander(r,s),t) - Expander(r,Expander(s,t));

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Expander Product is Associative "  << endl;	
	}else cout << "Expander Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////


/////////////////////////////////

// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format

	printf("Conserver equation set for test purposes\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[i]^Index[j]] == Higher(Rank[i],Rank[j])  ) {
						r = Blade[i]*Blade[j];
						weight = GetWeight(r);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////


/////////////////////////////////

//	Test Associativity of Conserver

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Conserver(Conserver(r,s),t) - Conserver(r,Conserver(s,t));

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Conserver Product is Associative "  << endl;	
	}else cout << "Conserver Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format

	printf("Shrinker equation set for test purposes\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					if(CountSetBits[Index[i]^Index[j]] < Higher(Rank[i],Rank[j])  ) {
						r = Blade[i]*Blade[j];
						weight = GetWeight(r);
		//				if(weight == 0) printf("   0      ");
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////

// Shrinker


/////////////////////////////////

//	Test Associativity of Shrinker

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Shrinker(Shrinker(r,s),t) - Shrinker(r,Shrinker(s,t));

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Shrinker Product is Associative "  << endl;	
	}else cout << "Shrinker Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

// test composition: Clifford = Expander + Conserver + Shrinker

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Clifford(r,s) - Expander(r,s) - Conserver(r,s) - Shrinker(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Clifford = Expander + Conserver + Shrinker"  << endl;	
	}else cout << "Clifford != Expander + Conserver + Shrinker \n";
	cout << "\n";

////////////////////////////////

//	print symmetric product terms

	printf("Symmetric Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( SymmetricViaFormula(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print antisymmetric product terms

	printf("AntiSymmetric Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( AntiSymmetricViaFormula(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




/////////////////////////////////

// we now print out the constituent equations for Symmetric product

	printf("Symmetric Product Equations\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					r = SymmetricViaFormula(Blade[i],Blade[j]);
					weight = GetWeight(r);
					if(weight != 0) {
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////

// test Symmetric

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = SymmetricViaFormula(r,s) - Symmetric(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Symmetric and SymmetricViaFormula Agree"  << endl;	
	}else cout << "Symmetric and SymmetricViaFormula Disagree \n";
	cout << "\n";




/////////////////////////////////

// we now print out the constituent equations for AntiSymmetric product

	printf("AntiSymmetric Product Equations\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					r = AntiSymmetricViaFormula(Blade[i],Blade[j]);
					weight = GetWeight(r);
					if(weight != 0) {
						if(weight  < 0) printf(" - "); 
						if(weight >  0) printf(" + ");
						if(weight != 0) printf("a.%s*b.%s",Field[i],Field[j]);
					}
				}
			}
		}
		printf(" ; \n");
	}
	printf("\n\n");


/////////////////////////////////

// test AntiSymmetric

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = AntiSymmetricViaFormula(r,s) - AntiSymmetric(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiSymmetric and AntiSymmetricViaFormula Agree"  << endl;	
	}else cout << "AntiSymmetric and AntiSymmetricViaFormula Disagree \n";
	cout << "\n";




/////////////////////////////////

// test CliffordViaPython

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Clifford(r,s) - CliffordViaPython(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Clifford and CliffordViaPython Agree"  << endl;	
	}else cout << "Clifford and CliffordViaPython Disagree \n";
	cout << "\n";


/////////////////////////////////

// test WedgeViaPython

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Wedge(r,s) - WedgeViaPython(r,s);

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

	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Wedge and WedgeViaPython Agree"  << endl;	
	}else cout << "Wedge and WedgeViaPython Disagree \n";
	cout << "\n";


////////////////////////////////

//	print InnerViaPython product terms

	printf("InnerViaPython Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( InnerViaPython(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print LeftContraction product terms

	printf("LeftContraction Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( LeftContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print RightContraction product terms

	printf("RightContraction Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( RightContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




/////////////////////////////////

// test Adjugate

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = r*Adjugate(r);

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

	cout << "Test Adjugate = Inverse*determinant\n";
	cout << "u = r*Adjugate(r) = " << u << "\nExpect scalar only\n\n";


/////////////////////////////////

// test Reciprocal

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	r = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Reciprocal(r);
	u = r*s;

	cout << "Test Reciprocal numerically\n";
	cout << "r = " << r << "\n";
	cout << "s = 1/r = " << s << "\n";
	cout << "u = r*s = " << u << "\n";

	u = s*r;

	cout << "Test Reciprocal numerically cummutative\n";
	cout << "r = " << r << "\n";
	cout << "s = 1/r = " << s << "\n";
	cout << "u = s*r = " << u << "\n";



/////////////////////////////////

// test Dual

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = Mink(0,  0,0,0,0,  0,0,0,0,0,0, 0,0,0,0, -1);   // s = I_inv
	u = r*s;
	t = Dual(r);

	cout << "u = r*s = " << u << "\n";

	s = t - u;
	cout << "u - Dual(r) = " << s << " expect 0\n";

	s = Dual(Dual(r));
	cout << "Dual(Dual(r)) = " << s << " expect -r\n";




///////////////////////////////// p. 78 Dorst

//	Test LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));
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

//	cout << "u = " << u << "\n";
	if (u==Mink(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "LeftContraction(Wedge(r,s), t) == LeftContraction(r,LeftContraction(s,t)) "  << endl;	
	}else cout << "LeftContraction(Wedge(r,s), t) != LeftContraction(r,LeftContraction(s,t)) \n";
	cout << "\n";


/////////////////////////////////  Dorst p. 80, 81

//	Test DorstDual;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Dual(r);
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

	cout << "     Dual(r) = u = " << u << "\n";

	u = DorstDual(r);
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

	cout << "DorstDual(r) = u = " << u << "\n";

	u = DorstUnDual(r);
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

	cout << "DorstUnDual(r) = u = " << u << "\n";

	u = DorstDual(DorstUnDual(r));
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

	cout << "DorstDual(DorsrUnDual(r)) = u = " << u << "\n";

	cout << "\n";

/////////////////////////////////  


	return 0;
}



