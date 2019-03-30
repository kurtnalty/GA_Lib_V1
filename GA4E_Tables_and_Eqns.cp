// compare matrix product and geometric product

// compile by
// g++ Demo_GA4E.cp -l ginac -l cln


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>

#include "GA4E_Routines.cp"

using namespace std;
using namespace GiNaC;

int SignArray[16];


///////////////////////////////////

void PrintSimpleBlade(GA4E a) 
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

int	GetWeight(GA4E r) // assume only one non-zero component
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

///////////////////////////////////////////////////////


GA4E SymmetricViaFormula(GA4E a, GA4E b)
{
	GA4E  c,d,e,f;
	
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


GA4E AntiSymmetricViaFormula(GA4E a, GA4E b)
{
	GA4E  c,d,e,f;
	
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


	GA4E first;
	GA4E second;
	GA4E third;
	GA4E check;
	GA4E Zero;   Zero = GA4E();

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

	GA4E q, x,y,z,t, xy,xz,yz,xt,yt,zt, xyz,xyt,xzt,yzt, xyzt;
	GA4E r,s,u;

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


////////////////////////////////////

//	Test AntiWedge Equations

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Formula is Consistent with OverBar((UnderBar(r)^UnderBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test AntiWedge Equations

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Formula is Consistent with UnderBar((OverBar(r)^OverBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test Regressive Equations

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Regressive is Consistent with RegressiveViaFormula"  << endl;	
	}else cout << "Check RegressiveProductViaEquations \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Clifford

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Clifford Product is Associative "  << endl;	
	}else cout << "Clifford Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Wedge

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Wedge Product is Associative "  << endl;	
	}else cout << "Wedge Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of AntiWedge

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiWedge Product is Associative "  << endl;	
	}else cout << "AntiWedge Product is Non-associative \n";
	cout << "\n";


///////////////////////////////////////

// print multiplication tables
//	ex q,  x,y,z,t, xy,xz,yz,xt,yt,zt,  xyz,xyt,xzt,yzt,  xyzt;

	q    = GA4E(1, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);

	x    = GA4E(0, 1,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	y    = GA4E(0, 0,1,0,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	z    = GA4E(0, 0,0,1,0, 0,0,0,0,0,0, 0,0,0,0, 0);
	t    = GA4E(0, 0,0,0,1, 0,0,0,0,0,0, 0,0,0,0, 0);

	xy   = GA4E(0, 0,0,0,0, 1,0,0,0,0,0, 0,0,0,0, 0);
	xz   = GA4E(0, 0,0,0,0, 0,1,0,0,0,0, 0,0,0,0, 0);
	yz   = GA4E(0, 0,0,0,0, 0,0,1,0,0,0, 0,0,0,0, 0);
	xt   = GA4E(0, 0,0,0,0, 0,0,0,1,0,0, 0,0,0,0, 0);
	yt   = GA4E(0, 0,0,0,0, 0,0,0,0,1,0, 0,0,0,0, 0);
	zt   = GA4E(0, 0,0,0,0, 0,0,0,0,0,1, 0,0,0,0, 0);

	xyz  = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 1,0,0,0, 0);
	xyt  = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 0,1,0,0, 0);
	xzt  = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,0, 0);
	yzt  = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,1, 0);

	xyzt = GA4E(0, 0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 1);

	GA4E Blade[16];

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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "LowerRightViaFormula Product is Associative "  << endl;	
	}else cout << "LowerRightViaFormula Product is Non-associative \n";
	cout << "\n";

////////////////////////////////

//	print geometric product increased rank terms

	printf("Terms with increased rank\n");
	printf(" >    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
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
	printf(" =    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
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

//	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = r*s - Product(r,s);;
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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Wedge and expander products agree "  << endl;	
	}else cout << "Wedge - Expander differ (as expected)\n";
	cout << "\n";

	cout << "u = " << u << "\n\n\n";



/////////////////////////////////

//	Test Associativity of Expander

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "Shrinker Product is Associative "  << endl;	
	}else cout << "Shrinker Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

// test composition: Clifford = Expander + Conserver + Shrinker

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	u = Product(r,s) - Expander(r,s) - Conserver(r,s) - Shrinker(r,s);

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "AntiSymmetric and AntiSymmetricViaFormula Agree"  << endl;	
	}else cout << "AntiSymmetric and AntiSymmetricViaFormula Disagree \n";
	cout << "\n";




////////////////////////////////

//	print Inner product terms

	printf("Inner Product\n");
	printf(" ?    | ");   for (i=0;i<16;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------------------------------------------------------\n");

	for (i=0;i<16;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<16;j++) {
			PrintSimpleBlade( Inner(Blade[i],Blade[j] ) ); 
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	r = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
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

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	s = GA4E(0,  0,0,0,0,  0,0,0,0,0,0, 0,0,0,0, 1);   // s = I_inv
	u = r*s;
	t = Dual(r);

	cout << "u = r*s = " << u << "\n";

	s = t - u;
	cout << "u - Dual(r) = " << s << " expect 0\n";

	s = Dual(Dual(r));
	cout << "Dual(Dual(r)) = " << s << " expect r\n";




///////////////////////////////// p. 78 Dorst

//	Test LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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
	if (u==GA4E(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) {
		cout << "LeftContraction(Wedge(r,s), t) == LeftContraction(r,LeftContraction(s,t)) "  << endl;	
	}else cout << "LeftContraction(Wedge(r,s), t) != LeftContraction(r,LeftContraction(s,t)) \n";
	cout << "\n";


/////////////////////////////////  Dorst p. 80, 81

//	Test DorstDual;

	r = GA4E(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = GA4E(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = GA4E(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = GA4E(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

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



