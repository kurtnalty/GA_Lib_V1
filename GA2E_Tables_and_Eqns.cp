// compare matrix product and geometric product

// compile by
// g++ Demo_GA2E.cp -l ginac -l cln


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>

#include "GA2E_Routines.cp"

using namespace std;
using namespace GiNaC;

int SignArray[4];


///////////////////////////////////

void PrintSimpleBlade(GA2E a) 
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

	if(a.xy ==  1) printf(" xy   ");
	if(a.xy == -1) printf("-xy   ");
	if(a.xy ==  0) Count++;

	if(Count == 4) printf(" 0    ");
	if(Count < 3)  printf("error ");


}

///////////////////////////////////

int	GetWeight(GA2E r) // assume only one non-zero component
{
	int weight = 0;		// default

	if (r.q ==  1) weight =  1;
	if (r.q == -1) weight = -1;

	if (r.x ==  1) weight =  1;
	if (r.x == -1) weight = -1;
	if (r.y ==  1) weight =  1;
	if (r.y == -1) weight = -1;

	if (r.xy ==  1) weight =  1;
	if (r.xy == -1) weight = -1;

	return weight;
}

///////////////////////////////////////////////////////


GA2E SymmetricViaFormula(GA2E a, GA2E b)
{
	GA2E  c,d,e,f;
	
	c = (a*b + b*a)/2;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}



///////////////////////////////////////////////////////


GA2E AntiSymmetricViaFormula(GA2E a, GA2E b)
{
	GA2E  c,d,e,f;
	
	c = (a*b - b*a)/2;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}


///////////////////////////////////

int	Rank[4]  = {0, 1,1, 2} ;
long	Index[4] = {0, 1,2, 3};	// organized by grade
int	CountSetBits[4] = {0, 1,1, 2};

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

	GA2E first;
	GA2E second;
	GA2E third;
	GA2E check;
	GA2E Zero;   Zero = GA2E();

	GA2E q, x,y, xy;
	GA2E r,s,t,u;

	symbol a_q("a"), a_x("b"), a_y("c"), a_xy("d") ;
	symbol b_q("A"), b_x("B") ,  b_y("C"), b_xy("D") ;
	symbol c_q("c.a"), c_x("c.b"), c_y("c.c"), c_xy("c.d") ;
	symbol d_q("d.a"), d_x("d.b"), d_y("d.c"), d_xy("d.d");

////////////////////////////////////

//	Test AntiWedge Equations

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	t = OverBar((UnderBar(r)^UnderBar(s)));
	u = AntiWedge(r,s) - t;

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "AntiWedge Formula is Consistent with OverBar((UnderBar(r)^UnderBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test AntiWedge Equations

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	t = UnderBar((OverBar(r)^OverBar(s)));
	u = AntiWedge(r,s) - t;

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "AntiWedge Formula is Consistent with UnderBar((OverBar(r)^OverBar(s)))"  << endl;	
	}else cout << "Check AntiWedge Equations \n";
	cout << "\n";


////////////////////////////////////

//	Test Regressive Equations

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Regressive(r,s) - RegressiveViaFormula(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Regressive is Consistent with RegressiveViaFormula"  << endl;	
	}else cout << "Check RegressiveProductViaEquations \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Clifford

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = (r*s)*t - r*(s*t);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Clifford Product is Associative "  << endl;	
	}else cout << "Clifford Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of Wedge

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = ((r^s)^t) - (r^(s^t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Wedge Product is Associative "  << endl;	
	}else cout << "Wedge Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

//	Test Associativity of AntiWedge

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = AntiWedge(AntiWedge(r,s),t) - AntiWedge(r,AntiWedge(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "AntiWedge Product is Associative "  << endl;	
	}else cout << "AntiWedge Product is Non-associative \n";
	cout << "\n";


///////////////////////////////////////

// print multiplication tables
//	ex q,  x,y,z, xy,xz,yz,  xyz;

	q    = GA2E(1, 0,0, 0);

	x    = GA2E(0, 1,0, 0);
	y    = GA2E(0, 0,1, 0);

	xy   = GA2E(0, 0,0, 1);

	GA2E Blade[4];

	Blade[ 0] = q;
	Blade[ 1] = x;
	Blade[ 2] = y;
	Blade[ 3] = xy;

//	print geometric product multiplication table

	printf(" *    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n"); 
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade(Blade[i]*Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Wedge product multiplication table

	printf(" ^    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade(Blade[i]^Blade[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print AntiWedge product multiplication table

	printf("Lengyel\n");
	printf(" V    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade(AntiWedge(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print Regressive product multiplication table

	printf("Hestenes\n");
	printf(" V    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade(Regressive(Blade[i],Blade[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");

//	print Clifford - Wedge

	printf("Blade[i]*Blade[j] - (Blade[i]^Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Blade[i]*Blade[j] - (Blade[i]^Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print ?? Contraction product multiplication table

	printf("AntiWedge(OverBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = AntiWedge(OverBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(UnderBar(Blade[i]),Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = AntiWedge(UnderBar(Blade[i]),Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = AntiWedge(Blade[i],OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



//	print ?? Contraction product multiplication table

	printf("AntiWedge(Blade[i],UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = AntiWedge(Blade[i],UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(OverBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(OverBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(OverBar(Blade[i]),UnderBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),OverBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(UnderBar(Blade[i]),OverBar(Blade[j]));
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////

//	print ?? Contraction product multiplication table

	printf("Wedge(UnderBar(Blade[i]),UnderBar(Blade[j]))\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
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
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(Blade[i]*xy,Blade[j]*xy);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(Blade[i]*xyz,xyz*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(Blade[i]*xy,xy*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyz*Blade[i],Blade[j]*xyz)\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(xy*Blade[i],Blade[j]*xy);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("Wedge(xyz*Blade[i],xyz*Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Wedge(xy*Blade[i],xy*Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////  

//	print ?? Contraction product multiplication table

	printf("LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////////// -Wedge(Blade[i]*xyz,xyz*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula

	printf("Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Blade[i]*Blade[j] - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");


/////////////////////////////////////// -Wedge(Blade[i]*xyz,xyz*Blade[j]) LowerRightViaFormula(a,b)


//	print Clifford - LowerRightViaFormula - Wedge

	printf("Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j])\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			u = Blade[i]*Blade[j] - Wedge(Blade[i],Blade[j]) - LowerRightViaFormula(Blade[i],Blade[j]);
			PrintSimpleBlade(u);
		}
		printf("\n");
	}
	printf("\n\n\n");



/////////////////////////////////

//	Test Associativity of LowerRightViaFormula

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = LowerRightViaFormula(LowerRightViaFormula(r,s),t) - LowerRightViaFormula(r,LowerRightViaFormula(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "LowerRightViaFormula Product is Associative "  << endl;	
	}else cout << "LowerRightViaFormula Product is Non-associative \n";
	cout << "\n";

////////////////////////////////

//	print geometric product increased rank terms

	printf("Terms with increased rank\n");
	printf(" >    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
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
	printf(" =    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
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
	printf(" <    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			k = CountSetBits[Index[i]^Index[j]];	// clifford product index
			if (k < Higher(Rank[i],Rank[j])) PrintSimpleBlade(Blade[i]*Blade[j]); 
					else PrintSimpleBlade(Zero);
		}
		printf("\n");
	}
	printf("\n\n\n");


///////////////////////////////////////////////////


	int NumTerms = 4;
	int EqnBasis, weight;
	
	char Field[5][10];

//	r = GA2E(a_q, a_x,a_y, a_xy);  

	for (i=0;i<5;i++) for (j=0;j<10;j++) Field[i][j] = 0;	// organized by grade
	strcpy(Field[0], "q  ");

	strcpy(Field[1], "x  ");
	strcpy(Field[2], "y  ");

	strcpy(Field[3], "xy ");

	strcpy(Field[4],"0  ");



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

	printf("AntiWedge\n");
	for (k=0;k<NumTerms;k++){	// equation number
		EqnBasis =  k;
		printf("c.%s = ", Field[EqnBasis]);
		for (i=0;i<NumTerms;i++) {
			for (j=0;j<NumTerms;j++) {  
				if(Index[k] == (3^Index[i]^Index[j])) {	// potential term for this equation - Clifford product as done
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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = r*s - Product(r,s);;
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Wedge(r,s) - Expander(r,s);;
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Wedge and expander products agree "  << endl;	
	}else cout << "Wedge - Expander differ (as expected)\n";
	cout << "\n";

	cout << "u = " << u << "\n\n\n";



/////////////////////////////////

//	Test Associativity of Expander

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Expander(Expander(r,s),t) - Expander(r,Expander(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Conserver(Conserver(r,s),t) - Conserver(r,Conserver(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Shrinker(Shrinker(r,s),t) - Shrinker(r,Shrinker(s,t));

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Shrinker Product is Associative "  << endl;	
	}else cout << "Shrinker Product is Non-associative \n";
	cout << "\n";


/////////////////////////////////

// test composition: Clifford = Expander + Conserver + Shrinker

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Product(r,s) - Expander(r,s) - Conserver(r,s) - Shrinker(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "Clifford = Expander + Conserver + Shrinker"  << endl;	
	}else cout << "Clifford != Expander + Conserver + Shrinker \n";
	cout << "\n";

////////////////////////////////

//	print symmetric product terms

	printf("Symmetric Product\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade( SymmetricViaFormula(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print antisymmetric product terms

	printf("AntiSymmetric Product\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = SymmetricViaFormula(r,s) - Symmetric(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = AntiSymmetricViaFormula(r,s) - AntiSymmetric(r,s);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	if (u==Zero) {
		cout << "AntiSymmetric and AntiSymmetricViaFormula Agree"  << endl;	
	}else cout << "AntiSymmetric and AntiSymmetricViaFormula Disagree \n";
	cout << "\n";




////////////////////////////////

//	print Inner product terms

	printf("Inner Product\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade( Inner(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print LeftContraction product terms

	printf("LeftContraction Product\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade( LeftContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




////////////////////////////////

//	print RightContraction product terms

	printf("RightContraction Product\n");
	printf(" ?    | ");   for (i=0;i<4;i++) PrintSimpleBlade(Blade[i]);    printf("\n");
	printf("------------------------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleBlade(Blade[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleBlade( RightContraction(Blade[i],Blade[j] ) ); 
		}
		printf("\n");
	}
	printf("\n\n\n");




/////////////////////////////////

// test Adjugate

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = r*Adjugate(r);

	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	cout << "Test Adjugate = Inverse*determinant\n";
	cout << "u = r*Adjugate(r) = " << u << "\nExpect scalar only\n\n";


/////////////////////////////////

// test Reciprocal

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	r = GA2E(  3,    5,  7,  11) ; 
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

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	s = GA2E(0,  0,0, -1);   // s = I_inv
	u = r*s;
	t = Dual(r);

	cout << "u = r*s = " << u << "\n";

	s = t - u;
	cout << "u - Dual(r) = " << s << " expect 0\n";

	s = Dual(Dual(r));
	cout << "Dual(Dual(r)) = " << s << " expect -r\n";




///////////////////////////////// p. 78 Dorst

//	Test LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = LeftContraction(Wedge(r,s), t) - LeftContraction(r,LeftContraction(s,t));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

//	cout << "u = " << u << "\n";
	if (u==Zero) {
		cout << "LeftContraction(Wedge(r,s), t) == LeftContraction(r,LeftContraction(s,t)) "  << endl;	
	}else cout << "LeftContraction(Wedge(r,s), t) != LeftContraction(r,LeftContraction(s,t)) \n";
	cout << "\n";


/////////////////////////////////  Dorst p. 80, 81

//	Test DorstDual;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	u = Dual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	cout << "     Dual(r) = u = " << u << "\n";

	u = DorstDual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	cout << "DorstDual(r) = u = " << u << "\n";

	u = DorstUnDual(r);
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	cout << "DorstUnDual(r) = u = " << u << "\n";

	u = DorstDual(DorstUnDual(r));
	u.q = expand(u.q);
	u.x = expand(u.x);
	u.y = expand(u.y);
	u.xy = expand(u.xy);

	cout << "DorstDual(DorsrUnDual(r)) = u = " << u << "\n";

	cout << "\n";

/////////////////////////////////  


	return 0;
}



