// Routines for CHO Algebra and GiNaC
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 6 September 2018
// License: Freeware
// Alternative License: BSD
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

struct C
{
	ex q,  xyz;
	C() {q = 0;  xyz = 0;} // initialize to zero
	C(ex qq) {q = qq;   xyz = 0;}  // scalar to real component
	C(ex a1, ex a2) { q = a1;    xyz = a2;}
};

struct Q
{
	ex q,xy,yz,xz;
	Q(ex qq, ex xyxy, ex yzyz, ex xzxz) {q = qq; xy = xyxy; yz = yzyz; xz = xzxz;}
	Q(ex qq) {q = qq; xy = 0; yz = 0; xz = 0;}
	Q() {q = 0; xy = 0; yz = 0; xz = 0;}
};

struct O{	
	ex q,  i,j,ij, E, iE,jE,ijE;
	O() {q = 0; i = 0; j = 0; ij = 0; E = 0; iE = 0; jE = 0; ijE = 0;}
	O(ex qq) {q = qq; i = 0; j = 0; ij = 0; E = 0; iE = 0; jE = 0; ijE = 0;}
	O(ex qq, ex ii, ex jj, ex ijij,  ex EE, ex iEiE, ex jEjE, ex ijEijE) 
	 {q = qq; i = ii; j = jj; ij = ijij; E = EE; iE = iEiE; jE = jEjE; ijE = ijEijE;}
};

struct CH{
	ex	q,	x,	y,	xy,	
		z,	xz,	yz,	xyz;
	CH() {	q = 0;	x = 0;	y = 0;	xy = 0;
		z = 0;	xz = 0;	yz = 0;	xyz = 0; }
	CH(ex qq) {	q = qq;	x = 0;	y = 0;	xy = 0;
			z = 0;	xz = 0;	yz = 0;	xyz = 0; }
	CH( ex a0,  ex a1,  ex a2,  ex a3,  ex b0,  ex b1,  ex b2,  ex b3 )
	{	q = a0 ; 	x = a1 ; 	y = a2 ; 	xy = a3 ;
	 	z = b0 ; 	xz = b1 ; 	yz = b2 ; 	xyz = b3 ; }
	
};

struct CO{
	ex	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
		xyz,	xyzi,	xyzj,	xyzij,	xyzE,	xyziE,	xyzjE,	xyzijE;	
	CO() {
	   q = 0;    i = 0;      j = 0;     ij = 0;     E = 0;     iE = 0;     jE = 0;     ijE = 0;
	   xyz = 0;  xyzi = 0;   xyzj = 0;  xyzij = 0;  xyzE = 0;  xyziE = 0;  xyzjE = 0;  xyzijE = 0;
	}
	CO(ex qq) {
	   q = qq;   i = 0;      j = 0;     ij = 0;     E = 0;     iE = 0;     jE = 0;     ijE = 0;
	   xyz = 0;  xyzi = 0;   xyzj = 0;  xyzij = 0;  xyzE = 0;  xyziE = 0;  xyzjE = 0;  xyzijE = 0;
	}
	CO (	ex a0,  ex a1,  ex a2,  ex a3,  ex a4,  ex a5,  ex a6,  ex a7, 
		ex b0,  ex b1,  ex b2,  ex b3,  ex b4,  ex b5,  ex b6,  ex b7 ) {
			q = a0 ; 	i = a1 ; 	j = a2 ; 	ij = a3 ;
		 	E = a4 ; 	iE = a5 ; 	jE = a6 ; 	ijE = a7 ; 
			xyz = b0 ; 	xyzi = b1 ; 	xyzj = b2 ; 	xyzij = b3 ;
		 	xyzE = b4 ; 	xyziE = b5 ; 	xyzjE = b6 ; 	xyzijE = b7 ; }
};

struct HO{ex
	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
	xy,	xyi,	xyj,	xyij,	xyE,	xyiE,	xyjE,	xyijE,
	yz,	yzi,	yzj,	yzij,	yzE,	yziE,	yzjE,	yzijE,
	xz,	xzi,	xzj,	xzij,	xzE,	xziE,	xzjE,	xzijE;
	HO() {  q = 0;    i = 0;      j = 0;      ij = 0;	
		E = 0;    iE = 0;     jE = 0;     ijE = 0;
		xy = 0;   xyi = 0;    xyj = 0;    xyij = 0;
		xyE = 0;  xyiE = 0;   xyjE = 0;   xyijE = 0;
		yz = 0;   yzi = 0;    yzj = 0;    yzij = 0;
		yzE = 0;  yziE = 0;   yzjE = 0;   yzijE = 0;
		xz = 0;   xzi = 0;    xzj = 0;    xzij = 0;
		xzE = 0;  xziE = 0;   xzjE = 0;   xzijE = 0; }
	HO(ex qq) {
		q = qq;   i = 0;      j = 0;      ij = 0;	
		E = 0;    iE = 0;     jE = 0;     ijE = 0;
		xy = 0;   xyi = 0;    xyj = 0;    xyij = 0;
		xyE = 0;  xyiE = 0;   xyjE = 0;   xyijE = 0;
		yz = 0;   yzi = 0;    yzj = 0;    yzij = 0;
		yzE = 0;  yziE = 0;   yzjE = 0;   yzijE = 0;
		xz = 0;   xzi = 0;    xzj = 0;    xzij = 0;
		xzE = 0;  xziE = 0;   xzjE = 0;   xzijE = 0; }

	HO( 	ex a0,  ex a1,  ex a2,  ex a3,  ex a4,  ex a5,  ex a6,  ex a7, 
		ex b0,  ex b1,  ex b2,  ex b3,  ex b4,  ex b5,  ex b6,  ex b7, 
		ex c0,  ex c1,  ex c2,  ex c3,  ex c4,  ex c5,  ex c6,  ex c7, 
		ex d0,  ex d1,  ex d2,  ex d3,  ex d4,  ex d5,  ex d6,  ex d7 )
	{
	q = a0 ; 	i = a1 ; 	j = a2 ; 	ij = a3 ; 	E = a4 ; 	iE = a5 ; 	jE = a6 ; 	ijE = a7 ; 
	xy = b0 ; 	xyi = b1 ; 	xyj = b2 ; 	xyij = b3 ; 	xyE = b4 ; 	xyiE = b5 ; 	xyjE = b6 ; 	xyijE = b7 ; 
	yz = c0 ; 	yzi = c1 ; 	yzj = c2 ; 	yzij = c3 ; 	yzE = c4 ; 	yziE = c5 ; 	yzjE = c6 ; 	yzijE = c7 ; 
	xz = d0 ; 	xzi = d1 ; 	xzj = d2 ; 	xzij = d3 ; 	xzE = d4 ; 	xziE = d5 ; 	xzjE = d6 ; 	xzijE = d7 ; }
};

struct CHO{ex
	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
	x,	xi,	xj,	xij,	xE,	xiE,	xjE,	xijE,
	y,	yi,	yj,	yij,	yE,	yiE,	yjE,	yijE,
	xy,	xyi,	xyj,	xyij,	xyE,	xyiE,	xyjE,	xyijE,
	z,	zi,	zj,	zij,	zE,	ziE,	zjE,	zijE,
	xz,	xzi,	xzj,	xzij,	xzE,	xziE,	xzjE,	xzijE,
	yz,	yzi,	yzj,	yzij,	yzE,	yziE,	yzjE,	yzijE,
	xyz,	xyzi,	xyzj,	xyzij,	xyzE,	xyziE,	xyzjE,	xyzijE;	
	CHO() {
	q = 0;		i = 0;		j = 0;		ij = 0;		E = 0;		iE = 0;		jE = 0;		ijE = 0;
	x = 0;		xi = 0;		xj = 0;		xij = 0;	xE = 0;		xiE = 0;	xjE = 0;	xijE = 0;
	y = 0;		yi = 0;		yj = 0;		yij = 0;	yE = 0;		yiE = 0;	yjE = 0;	yijE = 0;
	xy = 0;		xyi = 0;	xyj = 0;	xyij = 0;	xyE = 0;	xyiE = 0;	xyjE = 0;	xyijE = 0;
	z = 0;		zi = 0;		zj = 0;		zij = 0;	zE = 0;		ziE = 0;	zjE = 0;	zijE = 0;
	xz = 0;		xzi = 0;	xzj = 0;	xzij = 0;	xzE = 0;	xziE = 0;	xzjE = 0;	xzijE = 0;
	yz = 0;		yzi = 0;	yzj = 0;	yzij = 0;	yzE = 0;	yziE = 0;	yzjE = 0;	yzijE = 0;
	xyz = 0;	xyzi = 0;	xyzj = 0;	xyzij = 0;	xyzE = 0;	xyziE = 0;	xyzjE = 0;	xyzijE = 0; }
	CHO(ex qq) {
	q = qq;		i = 0;		j = 0;		ij = 0;		E = 0;		iE = 0;		jE = 0;		ijE = 0;
	x = 0;		xi = 0;		xj = 0;		xij = 0;	xE = 0;		xiE = 0;	xjE = 0;	xijE = 0;
	y = 0;		yi = 0;		yj = 0;		yij = 0;	yE = 0;		yiE = 0;	yjE = 0;	yijE = 0;
	xy = 0;		xyi = 0;	xyj = 0;	xyij = 0;	xyE = 0;	xyiE = 0;	xyjE = 0;	xyijE = 0;
	z = 0;		zi = 0;		zj = 0;		zij = 0;	zE = 0;		ziE = 0;	zjE = 0;	zijE = 0;
	xz = 0;		xzi = 0;	xzj = 0;	xzij = 0;	xzE = 0;	xziE = 0;	xzjE = 0;	xzijE = 0;
	yz = 0;		yzi = 0;	yzj = 0;	yzij = 0;	yzE = 0;	yziE = 0;	yzjE = 0;	yzijE = 0;
	xyz = 0;	xyzi = 0;	xyzj = 0;	xyzij = 0;	xyzE = 0;	xyziE = 0;	xyzjE = 0;	xyzijE = 0; }

	CHO( 	ex a0,  ex a1,  ex a2,  ex a3,  ex a4,  ex a5,  ex a6,  ex a7, 
		ex b0,  ex b1,  ex b2,  ex b3,  ex b4,  ex b5,  ex b6,  ex b7, 
		ex c0,  ex c1,  ex c2,  ex c3,  ex c4,  ex c5,  ex c6,  ex c7, 
		ex d0,  ex d1,  ex d2,  ex d3,  ex d4,  ex d5,  ex d6,  ex d7, 
		ex e0,  ex e1,  ex e2,  ex e3,  ex e4,  ex e5,  ex e6,  ex e7, 
		ex f0,  ex f1,  ex f2,  ex f3,  ex f4,  ex f5,  ex f6,  ex f7, 
		ex g0,  ex g1,  ex g2,  ex g3,  ex g4,  ex g5,  ex g6,  ex g7, 
		ex h0,  ex h1,  ex h2,  ex h3,  ex h4,  ex h5,  ex h6,  ex h7 )
{
	q = a0 ; 	i = a1 ; 	j = a2 ; 	ij = a3 ; 	E = a4 ; 	iE = a5 ; 	jE = a6 ; 	ijE = a7 ; 
	x = b0 ; 	xi = b1 ; 	xj = b2 ; 	xij = b3 ; 	xE = b4 ; 	xiE = b5 ; 	xjE = b6 ; 	xijE = b7 ; 
	y = c0 ; 	yi = c1 ; 	yj = c2 ; 	yij = c3 ; 	yE = c4 ; 	yiE = c5 ; 	yjE = c6 ; 	yijE = c7 ; 
	xy = d0 ; 	xyi = d1 ; 	xyj = d2 ; 	xyij = d3 ; 	xyE = d4 ; 	xyiE = d5 ; 	xyjE = d6 ; 	xyijE = d7 ; 
	z = e0 ; 	zi = e1 ; 	zj = e2 ; 	zij = e3 ; 	zE = e4 ; 	ziE = e5 ; 	zjE = e6 ; 	zijE = e7 ; 
	xz = f0 ; 	xzi = f1 ; 	xzj = f2 ; 	xzij = f3 ; 	xzE = f4 ; 	xziE = f5 ; 	xzjE = f6 ; 	xzijE = f7 ; 
	yz = g0 ; 	yzi = g1 ; 	yzj = g2 ; 	yzij = g3 ; 	yzE = g4 ; 	yziE = g5 ; 	yzjE = g6 ; 	yzijE = g7 ; 
	xyz = h0 ; 	xyzi = h1 ; 	xyzj = h2 ; 	xyzij = h3 ; 	xyzE = h4 ; 	xyziE = h5 ; 	xyzjE = h6 ; 	xyzijE = h7 ; }
};


void Program_Failed (const char* desc);

istream &operator>>(istream &ff, C &v);
istream &operator>>(istream &ff, Q &v) ;
istream &operator>>(istream &ff, O &v);
istream &operator>>(istream &ff, CH &v);
istream &operator>>(istream &ff, CO &v) ;
istream &operator>>(istream &ff, HO &v) ;
istream &operator>>(istream &ff, CHO &v) ;

ostream &operator<<(ostream &ff, C &v) ;
ostream &operator<<(ostream &ff, Q &v) ;
ostream &operator<<(ostream &ff, O &v) ;
ostream &operator<<(ostream &ff, CH &v) ;
ostream &operator<<(ostream &ff, CO &v) ;
ostream &operator<<(ostream &ff, HO &v) ;
ostream &operator<<(ostream &ff, CHO &v) ;

C operator+(const C &a, const C &b) ;
Q operator+(const Q &a, const Q &b) ;
O operator+(const O &a, const O &b) ;
CH operator+(const CH &a, const CH &b) ;
CO operator+(const CO &a, const CO &b) ;
HO operator+(const HO &a, const HO &b) ;
CHO operator+(const CHO &a, const CHO &b) ;

C operator-(const C &a, const C &b) ;
Q operator-(const Q &a, const Q &b);
O operator-(const O &a, const O &b) ;
CH operator-(const CH &a, const CH &b) ;
CO operator-(const CO &a, const CO &b) ;
HO operator-(const HO &a, const HO &b) ;
CHO operator-(const CHO &a, const CHO &b) ;

int operator==(const C &a, const C &b) ;
int operator==(const Q &a, const Q &b);
int operator==(const O &a, const O &b) ;
int operator==(const CH &a, const CH &b) ;
int operator==(const CO &a, const CO &b) ;
int operator==(const HO &a, const HO &b) ;
int operator==(const CHO &a, const CHO &b) ;

int operator!=(const C &a, const C &b) ;
int operator!=(const Q &a, const Q &b) ;
int operator!=(const O &a, const O &b);
int operator!=(const CH &a, const CH &b) ;
int operator!=(const CO &a, const CO &b) ;
int operator!=(const HO &a, const HO &b) ;
int operator!=(const CHO &a, const CHO &b);

void PrintSimpleC(const C &a);
void PrintSimpleQ(const Q &a);
void PrintSimpleO(const O &a);
void PrintSimpleCH(const CH &a);
void PrintSimpleCO(const CO &a);
void PrintSimpleHO(const HO &a);
void PrintSimpleCHO(const CHO &a);

void Fill_C_Basis(void) ;
void Fill_Q_Basis(void) ;
void Fill_O_Basis(void) ;
void Fill_CH_Basis(void) ;
void Fill_CO_Basis(void) ;
void Fill_HO_Basis(void);
void Fill_CHO_Basis(void) ;

C operator*(const C &a, const C &b);
Q operator*(const Q &a, const Q &b) ;
O operator*(const O &a, const O &b) ;
CH operator*(const CH &a, const CH &b) ;
CO operator*(const CO &a, const CO &b) ;
HO operator*(const HO &a, const HO &b) ;
CHO operator*(const CHO &a, const CHO &b) ;
CHO GA3EO_Product(const CHO &a, const CHO &b) ;

void CHO_Zero_Grid(const CHO &a) ;
void HO_Zero_Grid(const HO &a) ;
void CO_Zero_Grid(const CO &a) ;
void CH_Zero_Grid(const CH &a) ;
void O_Zero_Grid(const O &a) ;
