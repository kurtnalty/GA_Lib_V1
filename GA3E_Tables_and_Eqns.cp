// compare matrix product and geometric product

// compile by
// g++ Demo_GA3E.cp -l ginac -l cln


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>

#include "GA3E_Routines.cp"

using namespace std;
using namespace GiNaC;

int SignArray[8];


///////////////////////////////////

void PrintSimpleBlade(GA3E a) 
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


	if(a.xy ==  1) printf(" xy   ");
	if(a.xy == -1) printf("-xy   ");
	if(a.xy ==  0) Count++;

	if(a.xz ==  1) printf(" xz   ");
	if(a.xz == -1) printf("-xz   ");
	if(a.xz ==  0) Count++;

	if(a.yz ==  1) printf(" yz   ");
	if(a.yz == -1) printf("-yz   ");
	if(a.yz ==  0) Count++;


	if(a.xyz ==  1) printf(" xyz  ");
	if(a.xyz == -1) printf("-xyz  ");
	if(a.xyz ==  0) Count++;

	if(Count == 8) printf(" 0    ");
	if(Count < 7)  printf("error ");


}

///////////////////////////////////

int	GetWeight(GA3E r) // assume only one non-zero component
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

	if (r.xy ==  1) weight =  1;
	if (r.xy == -1) weight = -1;
	if (r.xz ==  1) weight =  1;
	if (r.xz == -1) weight = -1;
	if (r.yz ==  1) weight =  1;
	if (r.yz == -1) weight = -1;

	if (r.xyz ==  1) weight =  1;
	if (r.xyz == -1) weight = -1;

	return weight;
}

///////////////////////////////////////////////////////


GA3E SymmetricViaFormula(GA3E a, GA3E b)
{
	GA3E  c,d,e,f;
	
	c = (a*b + b*a)/2;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}



///////////////////////////////////////////////////////


GA3E AntiSymmetricViaFormula(GA3E a, GA3E b)
{
	GA3E  c,d,e,f;
	
	c = (a*b - b*a)/2;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.z = expand(c.z);

	c.xy = expand(c.xy);
	c.xz = expand(c.xz);
	c.yz = expand(c.yz);
	c.xyz = expand(c.xyz);

	return c;
}


///////////////////////////////////

int	Rank[8]  = {0, 1,1,1, 2,2,2, 3} ;
long	Index[8] = {0, 1,2,4, 3,5,6, 7};	// organized by grade
int	CountSetBits[8] = {0, 1,1,2,1, 2,2,3};

int	Higher(int a, int b)	// return higher of a or b
{
	int c;
	
	c = a;
	if (b>a) c = b;

	return c;

}

///////////////////////////////////

int main(void)
{

	
	int i,j,k;

	GA3E first;
	GA3E second;
	GA3E third;
	GA3E check;
	GA3E Zero;   Zero = GA3E();


	GA3E q, x,y,z, xy,xz,yz, xyz;
	GA3E r,s,t,u;

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d");
	symbol a_xy("e"),  a_xz("f"),  a_yz("g");
	symbol a_xyz("h");

	symbol b_q("A");
	symbol b_x("B") ,  b_y("C") ,  b_z("D");
	symbol b_xy("E"),  b_xz("F"),  b_yz("G");
	symbol b_xyz("H");

	symbol c_q("c.a");
	symbol c_x("c.b") ,  c_y("c.c") ,  c_z("c.d");
	symbol c_xy("c.e"),  c_xz("c.f"),  c_yz("c.g");
	symbol c_xyz("c.h");

	symbol d_q("d.a");
	symbol d_x("d.b") ,  d_y("d.c") ,  d_z("d.d");
	symbol d_xy("d.e"),  d_xz("d.f"),  d_yz("d.g");
	symbol d_xyz("d.h");


////////////////////////////////////

//	Test AntiWedge Equations

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	t = OverBar((UnderBar(r)^UnderBar(s)));
	u = AntiWedge(r,s) - t;

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "AntiWedge Formula is Consistent with OverBar((UnderBar(r)^UnderBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test AntiWedge Equations

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	t = UnderBar((OverBar(r)^OverBar(s)));
	u = AntiWedge(r,s) - t;

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "AntiWedge Formula is Consistent with UnderBar((OverBar(r)^OverBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test Regressive Equations

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Regressive(r,s) - RegressiveViaFormula(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Regressive is Consistent with RegressiveViaFormula"  << endl;	
	}else cout << "Check RegressiveProductViaEquations \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Clifford

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = (r*s)*t - r*(s*t);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Clifford Product is Associative "  << endl;	
	}else cout << "Clifford Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Wedge

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = ((r^s)^t) - (r^(s^t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Wedge Product is Associative "  << endl;	
	}else cout << "Wedge Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of AntiWedge

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = AntiWedge(AntiWedge(r,s),t) - AntiWedge(r,AntiWedge(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "AntiWedge Product is Associative "  << endl;	
	}else cout << "AntiWedge Product is Non-associative \n";
	cout << "\n";


///////////////////////////////////////

// print multiplication tables
//	ex q,  x,y,z, xy,xz,yz,  xyz;

	q    = GA3E(1, 0,0,0, 0,0,0, 0);

	x    = GA3E(0, 1,0,0, 0,0,0, 0);
	y    = GA3E(0, 0,1,0, 0,0,0, 0);
	z    = GA3E(0, 0,0,1, 0,0,0, 0);

	xy   = GA3E(0, 0,0,0, 1,0,0, 0);
	xz   = GA3E(0, 0,0,0, 0,1,0, 0);
	yz   = GA3E(0, 0,0,0, 0,0,1, 0);

	xyz  = GA3E(0, 0,0,0, 0,0,0, 1);

	GA3E Blade[8];

	Blade[ 0] = q;

	Blade[ 1] = x;
	Blade[ 2] = y;
	Blade[ 3] = z;

	Blade[ 4] = xy;
	Blade[ 5] = xz;
	Blade[ 6] = yz;

	Blade[ 7] = xyz;

//	print geometric product multiplication table

	printf(" *    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n"); 
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade(Blade[i]*Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Wedge product multiplication table

	printf(" ^    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade(Blade[i]^Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print AntiWedge product multiplication table

	printf("Lengyel\n");
	printf(" V    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade(AntiWedge(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Regressive product multiplication table

	printf("Hestenes\n");
	printf(" V    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade(Regressive(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print Clifford - Wedge

	printf("Blade[i]*Blade[j] - (Blade[i]^Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Blade[i]*Blade[j] - (Blade[i]^Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print ?? Contraction product multiplication table

	printf("AntiWedge(OverBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = AntiWedge(OverBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(UnderBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = AntiWedge(UnderBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = AntiWedge(Blade[i],OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = AntiWedge(Blade[i],UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(OverBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(OverBar(Blade[i]),UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(UnderBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(UnderBar(Blade[i]),UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////   Now we try various Hestenes forms


///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(Blade[i]*xyz,Blade[j]*xyz)\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(Blade[i]*xyz,Blade[j]*xyz);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(Blade[i]*xyz,xyz*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(Blade[i]*xyz,xyz*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyz*Blade[i],Blade[j]*xyz)\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(xyz*Blade[i],Blade[j]*xyz);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyz*Blade[i],xyz*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Wedge(xyz*Blade[i],xyz*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////////// -Wedge(Blade[i]*xyz,xyz*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula

	printf("Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


/////////////////////////////////////// -Wedge(Blade[i]*xyz,xyz*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula - Wedge

	printf("Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			u = Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////

//	Test Associativity of LowerRightViaFormula

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = LowerRightViaFormula(LowerRightViaFormula(r,s),t) - LowerRightViaFormula(r,LowerRightViaFormula(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "LowerRightViaFormula Product is Associative "  << endl;	
	}else cout << "LowerRightViaFormula Product is Non-associative \n";
	cout << "\n";

////////////////////////////////

//	print geometric product increased rank terms

	printf("Terms with increased rank\n");
	printf(" >    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
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
	printf(" =    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
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
	printf(" <    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			k = CountSetBits[Index[i]^Index[j]];	// clifford product index
			if (k < Higher(Rank[i],Rank[j])) PrintSimpleBlade(Blade[i]*Blade[j]); 
					else PrintSimpleBlade(Zero);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////////////////


	int NumTerms = 8;
	int EqnBasis, weight;
	
	char Field[9][10];

//	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

	for (i=0;i<9;i++) for (j=0;j<10;j++) Field[i][j] = 0;	// organized by grade
	strcpy(Field[0], "q  ");

	strcpy(Field[1], "x  ");
	strcpy(Field[2], "y  ");
	strcpy(Field[3], "z  ");

	strcpy(Field[4], "xy ");
	strcpy(Field[5], "xz ");
	strcpy(Field[6], "yz ");

	strcpy(Field[7],"xyz");

	strcpy(Field[8],"0  ");



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


// we now print out the constituent equations for C = A^B using s x y z xy yz xz xyz format

	printf("Wedge\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
					r = Wedge(Blade[i],Blade[j]);
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


// we now print out the constituent equations for C = A B using s x y z xy yz xz xyz format
// KN Fix This
	printf("AntiWedge\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
			//	if(Index[k] == (Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
				if(Index[k] == (7^Index[i]^Index[j])) {	// potential term for this equation - AntiWedge?
					r = AntiWedge(Blade[i],Blade[j]);
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = r*s - Product(r,s);;
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Wedge(r,s) - Expander(r,s);;
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Wedge and expander products agree "  << endl;	
	}else cout << "Wedge - Expander differ (as expected)\n";
	cout << "\n";

	cout << "u = " << u << "\n\n\n";



/////////////////////////////////

//	Test Associativity of Expander

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Expander(Expander(r,s),t) - Expander(r,Expander(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Conserver(Conserver(r,s),t) - Conserver(r,Conserver(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Shrinker(Shrinker(r,s),t) - Shrinker(r,Shrinker(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Shrinker Product is Associative "  << endl;	
	}else cout << "Shrinker Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

// test composition: Clifford = Expander + Conserver + Shrinker

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Product(r,s) - Expander(r,s) - Conserver(r,s) - Shrinker(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "Clifford = Expander + Conserver + Shrinker"  << endl;	
	}else cout << "Clifford != Expander + Conserver + Shrinker \n";
	cout << "\n";

////////////////////////////////

//	print symmetric product terms

	printf("Symmetric Product\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade( SymmetricViaFormula(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print antisymmetric product terms

	printf("AntiSymmetric Product\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = SymmetricViaFormula(r,s) - Symmetric(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = AntiSymmetricViaFormula(r,s) - AntiSymmetric(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	if (u==Zero) {
		cout << "AntiSymmetric and AntiSymmetricViaFormula Agree"  << endl;	
	}else cout << "AntiSymmetric and AntiSymmetricViaFormula Disagree \n";
	cout << "\n";




////////////////////////////////

//	print Inner product terms

	printf("Inner Product\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade( Inner(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print LeftContraction product terms

	printf("LeftContraction Product\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade( LeftContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print RightContraction product terms

	printf("RightContraction Product\n");
	printf(" ?    | ");   for (i=0;i<8;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("-------------------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleBlade( RightContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




/////////////////////////////////

// test Adjugate

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = r*Adjugate(r);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	cout << "Test Adjugate = Inverse*determinant\n";
	cout << "u = r*Adjugate(r) = " << u << "\nExpect scalar only\n\n";


/////////////////////////////////

// test Reciprocal

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	r = GA3E(  3,    5,  7, 11,     13, 17, 19,   23) ; 
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

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	s = GA3E(0,  0,0,0,  0,0,0, -1);   // s = I_inv
	u = r*s;
	t = Dual(r);

	cout << "u = r*s = " << u << "\n";

	s = t - u;
	cout << "u - Dual(r) = " << s << " expect 0\n";

	s = Dual(Dual(r));
	cout << "Dual(Dual(r)) = " << s << " expect -r\n";




///////////////////////////////// p. 78 Dorst

//	Test LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

//	cout << "u = " << u << "\n";
	if (u==Zero) {
		cout << "LeftContraction(Wedge(r,s), t) == LeftContraction(r,LeftContraction(s,t)) "  << endl;	
	}else cout << "LeftContraction(Wedge(r,s), t) != LeftContraction(r,LeftContraction(s,t)) \n";
	cout << "\n";


/////////////////////////////////  Dorst p. 80, 81

//	Test DorstDual;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	u = Dual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	cout << "     Dual(r) = u = " << u << "\n";

	u = DorstDual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	cout << "DorstDual(r) = u = " << u << "\n";

	u = DorstUnDual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	cout << "DorstUnDual(r) = u = " << u << "\n";

	u = DorstDual(DorstUnDual(r));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.z = expand(u.z);
	u.xy = expand(u.xy);
	u.xz = expand(u.xz);
	u.yz = expand(u.yz);
	u.xyz = expand(u.xyz);

	cout << "DorstDual(DorsrUnDual(r)) = u = " << u << "\n";

	cout << "\n";

/////////////////////////////////  


	return 0;
}



