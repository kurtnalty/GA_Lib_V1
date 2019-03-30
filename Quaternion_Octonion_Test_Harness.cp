// Test Routines for Vector, Quaternion and Octonion Algebra
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 9 August 2018
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

#include "Quaternion_Octonion_Routines.cp"


///////////////////////////////////

void PrintSimpleQuaternion(Quaternion a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q ");	if(a.q == -1) printf("-q ");	if(a.q ==  0) Count++;
	if(a.i ==  1) printf(" i ");	if(a.i == -1) printf("-i ");	if(a.i ==  0) Count++;
	if(a.j ==  1) printf(" j ");	if(a.j == -1) printf("-j ");	if(a.j ==  0) Count++;
	if(a.k ==  1) printf(" k ");	if(a.k == -1) printf("-k ");	if(a.k ==  0) Count++;

	if(Count == 4) printf(" 0 "); if(Count < 3)  printf("error ");
}


///////////////////////////////////

void PrintSimpleOctonion(Octonion a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q ");	if(a.q == -1) printf("-q ");	if(a.q ==  0) Count++;
	if(a.i ==  1) printf(" i ");	if(a.i == -1) printf("-i ");	if(a.i ==  0) Count++;
	if(a.j ==  1) printf(" j ");	if(a.j == -1) printf("-j ");	if(a.j ==  0) Count++;
	if(a.k ==  1) printf(" k ");	if(a.k == -1) printf("-k ");	if(a.k ==  0) Count++;
	if(a.E ==  1) printf(" E ");	if(a.E == -1) printf("-E ");	if(a.E ==  0) Count++;
	if(a.I ==  1) printf(" I ");	if(a.I == -1) printf("-I ");	if(a.I ==  0) Count++;
	if(a.J ==  1) printf(" J ");	if(a.J == -1) printf("-J ");	if(a.J ==  0) Count++;
	if(a.K ==  1) printf(" K ");	if(a.K == -1) printf("-K ");	if(a.K ==  0) Count++;

	if(Count == 8) printf(" 0 "); if(Count < 7)  printf("error ");
}

///////////////////////////////////

void PrintFancierOctonion(Octonion a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" e_0 ");	if(a.q == -1) printf("-e_0 ");	if(a.q ==  0) Count++;
	if(a.i ==  1) printf(" e_1 ");	if(a.i == -1) printf("-e_1 ");	if(a.i ==  0) Count++;
	if(a.j ==  1) printf(" e_2 ");	if(a.j == -1) printf("-e_2 ");	if(a.j ==  0) Count++;
	if(a.k ==  1) printf(" e_3 ");	if(a.k == -1) printf("-e_3 ");	if(a.k ==  0) Count++;
	if(a.E ==  1) printf(" e_4 ");	if(a.E == -1) printf("-e_4 ");	if(a.E ==  0) Count++;
	if(a.I ==  1) printf(" e_5 ");	if(a.I == -1) printf("-e_5 ");	if(a.I ==  0) Count++;
	if(a.J ==  1) printf(" e_6 ");	if(a.J == -1) printf("-e_6 ");	if(a.J ==  0) Count++;
	if(a.K ==  1) printf(" e_7 ");	if(a.K == -1) printf("-e_7 ");	if(a.K ==  0) Count++;

	if(Count == 8) printf(" 0 "); if(Count < 7)  printf("error ");
}

//////////////////////////////////////////////////////

int main(void)
{
	int i,j,k;

	Quaternion Quaternion_Basis[4];
	Octonion  Octonion_Basis[8];

	Quaternion_Basis[0] = Quaternion(1, 0,0,0);
	Quaternion_Basis[1] = Quaternion(0, 1,0,0);
	Quaternion_Basis[2] = Quaternion(0, 0,1,0);
	Quaternion_Basis[3] = Quaternion(0, 0,0,1);

	Octonion_Basis[0] = Octonion(1, 0,0,0, 0,0,0, 0);
	Octonion_Basis[1] = Octonion(0, 1,0,0, 0,0,0, 0);
	Octonion_Basis[2] = Octonion(0, 0,1,0, 0,0,0, 0);
	Octonion_Basis[3] = Octonion(0, 0,0,1, 0,0,0, 0);
	Octonion_Basis[4] = Octonion(0, 0,0,0, 1,0,0, 0);
	Octonion_Basis[5] = Octonion(0, 0,0,0, 0,1,0, 0);
	Octonion_Basis[6] = Octonion(0, 0,0,0, 0,0,1, 0);
	Octonion_Basis[7] = Octonion(0, 0,0,0, 0,0,0, 1);


//	print quaternion product multiplication table

	printf(" * | ");   for (i=0;i<4;i++) PrintSimpleQuaternion(Quaternion_Basis[i]);    printf("\n"); 
	printf("-----------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleQuaternion(Quaternion_Basis[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleQuaternion(Quaternion_Basis[i]*Quaternion_Basis[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");

/////////////////////////// Print Baez' Octonions //////////////////
	
	Baez = 1;

	printf("Baez's Octonion Table, Featuring Cyclic Indices\n");

//	print octonion product multiplication table

	printf(" * | ");   for (i=0;i<8;i++) PrintSimpleOctonion(Octonion_Basis[i]);    printf("\n"); 
	printf("-----------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleOctonion(Octonion_Basis[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleOctonion(Octonion_Basis[i]*Octonion_Basis[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print fancier octonion product multiplication table

	printf(" *   | ");   for (i=0;i<8;i++) PrintFancierOctonion(Octonion_Basis[i]);    printf("\n"); 
	printf("-----------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintFancierOctonion(Octonion_Basis[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintFancierOctonion(Octonion_Basis[i]*Octonion_Basis[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");

/////////////////////////// Print Cayley's Octonions //////////////////
	
	Baez = 0;

	printf("Cayley's Octonion Table, Showing Quaternion Inheritance\n");

//	print octonion product multiplication table

	printf(" * | ");   for (i=0;i<8;i++) PrintSimpleOctonion(Octonion_Basis[i]);    printf("\n"); 
	printf("-----------------------------\n");

	for (i=0;i<8;i++) {
		PrintSimpleOctonion(Octonion_Basis[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintSimpleOctonion(Octonion_Basis[i]*Octonion_Basis[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");


//	print fancier octonion product multiplication table

	printf(" *   | ");   for (i=0;i<8;i++) PrintFancierOctonion(Octonion_Basis[i]);    printf("\n"); 
	printf("-----------------------------------------------\n");

	for (i=0;i<8;i++) {
		PrintFancierOctonion(Octonion_Basis[i]);  printf("| ");
		for(j=0;j<8;j++) {
			PrintFancierOctonion(Octonion_Basis[i]*Octonion_Basis[j]);
		}
		printf("\n");
	}
	printf("\n\n\n");

///////////////////////////////// Verify eight square identity ///////////////////////

	ex A, B,C,D,E,  F,G,H;
	ex a, b,c,d,e,  f,g,h;
	ex t1, t2, t3, delta;

	A = symbol("A");	B = symbol("B");	C = symbol("C");	D = symbol("D");
	E = symbol("E");	F = symbol("F");	G = symbol("G");	H = symbol("H");

	a = symbol("a");	b = symbol("b");	c = symbol("c");	d = symbol("d");
	e = symbol("e");	f = symbol("f");	g = symbol("g");	h = symbol("h");

	Octonion Oct1, Oct2, Oct3, Oct4, Oct5, Oct6, Zero;

	Zero = Octonion();
	Oct1 = Octonion(a, b,c,d, e, f,g,h);
	Oct2 = Octonion(A, B,C,D, E, F,G,H);
	Oct3 = Oct1*Oct2;

	Oct4 = (Oct1*conj(Oct1));
	Oct5 = (Oct2*conj(Oct2));
	Oct6 = (Oct3*conj(Oct3));

	cout << "Oct1*conj(Oct1) = " << Oct4 << "\n\n";
	cout << "Oct2*conj(Oct2) = " << Oct5 << "\n\n";
	cout << "Oct3*conj(Oct3) = " << Oct6 << "\n\n";

	t1 = (Oct1*conj(Oct1)).q;
	t2 = (Oct2*conj(Oct2)).q;
	t3 = (Oct3*conj(Oct3)).q;
	delta = expand(t3 - t1*t2);
	cout << "delta = " << delta << " (expect 0) \n\n";

/////////////////////////////////////////////////////

//	Demonstrate Quaternion non-commutative, associative behaviors

	Quaternion Q1, Q2, Q3, Q4;

	Q1 = Quaternion(a, b,c,d);
	Q2 = Quaternion(A, B,C,D);
	Q3 = Quaternion(e, f,g,h);

	Q4 = Q1*Q2 - Q2*Q1;	// Quaternion Commutator
	cout << "Q1*Q2 - Q2*Q1 = " << Q4 << " (expect non-zero result) \n\n";
	Q4 = (Q1*Q2)*Q3 - Q1*(Q2*Q3);	// Quaternion Associator
	cout << "(Q1*Q2)*Q3 - Q1*(Q2*Q3) = " << Q4 << " (expect zero) \n\n";

/////////////////////////////////////////////////////

//	Demonstrate Octonion non-commutative, non-associative, triality behaviors

	ex M, N,P,R, S, T,U,V;
	M = symbol("M");	N = symbol("N");	P = symbol("P");	R = symbol("R");
	S = symbol("S");	T = symbol("T");	U = symbol("U");	V = symbol("V");


	Oct1 = Octonion(a, b,c,d, e, f,g,h);
	Oct2 = Octonion(A, B,C,D, E, F,G,H);
	Oct3 = Octonion(M, N,P,R, S, T,U,V);

	Baez = 0;

	cout << "Using Baez = 0 \n";
	Oct4 = Oct1*Oct2 - Oct2*Oct1;	// Octonion Commutator
	cout << "Oct1*Oct2 - Oct2*Oct1 = " << Oct4 << " (expect non-zero result) \n\n";
	Oct4 = (Oct1*Oct2)*Oct3 - Oct1*(Oct2*Oct3);	// Octonion Associator
	cout << "(Oct1*Oct2)*Oct3 - Oct1*(Oct2*Oct3) = " << Oct4 << " (expect non-zero result) \n\n";
	Oct4 = (Oct1*Oct2)*Oct3 + (Oct2*Oct3)*Oct1 + (Oct3*Oct1)*Oct2;
	cout << "(Oct1*Oct2)*Oct3 + (Oct2*Oct3)*Oct1 + (Oct3*Oct1)*Oct2 = " << Oct4 << "  \n\n";
	Oct4 = Oct1*(Oct2*Oct3) + Oct2*(Oct3*Oct1) + Oct3*(Oct1*Oct2);
	cout << "Oct1*(Oct2*Oct3) + Oct2*(Oct3*Oct1) + Oct3*(Oct1*Oct2) = " << Oct4 << "  \n\n";

	Baez = 1;

	cout << "Using Baez = 1 \n";
	Oct4 = Oct1*Oct2 - Oct2*Oct1;	// Octonion Commutator
	cout << "Oct1*Oct2 - Oct2*Oct1 = " << Oct4 << " (expect non-zero result) \n\n";
	Oct4 = (Oct1*Oct2)*Oct3 - Oct1*(Oct2*Oct3);	// Octonion Associator
	cout << "(Oct1*Oct2)*Oct3 - Oct1*(Oct2*Oct3) = " << Oct4 << " (expect non-zero result) \n\n";
	Oct4 = (Oct1*Oct2)*Oct3 + (Oct2*Oct3)*Oct1 + (Oct3*Oct1)*Oct2;
	cout << "(Oct1*Oct2)*Oct3 + (Oct2*Oct3)*Oct1 + (Oct3*Oct1)*Oct2 = " << Oct4 << "  \n\n";
	Oct4 = Oct1*(Oct2*Oct3) + Oct2*(Oct3*Oct1) + Oct3*(Oct1*Oct2);
	cout << "Oct1*(Oct2*Oct3) + Oct2*(Oct3*Oct1) + Oct3*(Oct1*Oct2) = " << Oct4 << "  \n\n";

	Oct4 = Jacobi(Oct1, Oct2, Oct3);
	cout << "Jacobi(Oct1, Oct2, Oct3) = " << Oct4 << "  \n\n";
	Oct4 = Jacobi(Oct1, Oct2, Oct2);
	cout << "Jacobi(Oct1, Oct2, Oct2) = " << Oct4 << "  \n\n";
	Oct4 = Jacobi(Oct1, Oct1, Oct2);
	cout << "Jacobi(Oct1, Oct1, Oct2) = " << Oct4 << "  \n\n";

	Baez = 0;
	cout << "Using Baez = 0 \n";
	Oct4 = Jacobi(Oct1, Oct2, Oct3);
	cout << "Jacobi(Oct1, Oct2, Oct3) = " << Oct4 << "  \n\n";
	Oct4 = Jacobi(Oct1, Oct2, Oct2);
	cout << "Jacobi(Oct1, Oct2, Oct2) = " << Oct4 << "  \n\n";
	Oct4 = Jacobi(Oct1, Oct1, Oct2);
	cout << "Jacobi(Oct1, Oct1, Oct2) = " << Oct4 << "  \n\n";

//	try alternator relations
	Oct4 = Associator(Oct1, Oct1, Oct2);
	cout << "Associator(Oct1, Oct1, Oct2) = " << Oct4 << " (expect 0) \n\n";
	Oct4 = Associator(Oct1, Oct2, Oct2);
	cout << "Associator(Oct1, Oct2, Oct2) = " << Oct4 << " (expect 0) \n\n";
	Oct4 = Associator(Oct1, Oct2, Oct1);
	cout << "Associator(Oct1, Oct2, Oct1) = " << Oct4 << " (expect 0) \n\n";

//	test flexible identity x*(y*x) = (x*y)*x 
	Oct4 = Oct1*(Oct2*Oct1) - (Oct1*Oct2)*Oct1;
	cout << "Flexible Identity => Oct1*(Oct2*Oct1) - (Oct1*Oct2)*Oct1 = " << Oct4 << " (expect 0) \n\n";

	Baez = 1;
	cout << "Using Baez = 1 \n";
	Oct4 = Jacobi(Oct1, Oct2, Oct3);
	cout << "Jacobi(Oct1, Oct2, Oct3) = " << Oct4 << "  \n\n";

//	try alternator relations
	Oct4 = Associator(Oct1, Oct1, Oct2);
	cout << "Associator(Oct1, Oct1, Oct2) = " << Oct4 << " (expect 0) \n\n";
	Oct4 = Associator(Oct1, Oct2, Oct2);
	cout << "Associator(Oct1, Oct2, Oct2) = " << Oct4 << " (expect 0) \n\n";
	Oct4 = Associator(Oct1, Oct2, Oct1);
	cout << "Associator(Oct1, Oct2, Oct1) = " << Oct4 << " (expect 0) \n\n";

//	test flexible identity x*(y*x) = (x*y)*x 
	Oct4 = Oct1*(Oct2*Oct1) - (Oct1*Oct2)*Oct1;
	cout << "Flexible Identity => Oct1*(Oct2*Oct1) - (Oct1*Oct2)*Oct1 = " << Oct4 << " (expect 0) \n\n";

////////////////////////// Test == and != for Octonions ///////////////////////////

	Oct4 = Jacobi(Oct1, Oct2, Oct3);	// not zero
	if (Oct4 == Zero) cout << "Jacobi = 0 \n\n";
	if (Oct4 != Zero) cout << "Jacobi != 0 (as expected) \n\n";
	Oct4 = Associator(Oct1, Oct2, Oct1);
	if (Oct4 == Zero) cout << "Associator(Oct1, Oct2, Oct1) = 0 (as expected) \n\n";
	if (Oct4 != Zero) cout << "Associator(Oct1, Oct2, Oct1) != 0 \n\n";

////////////////////////// Test == and != for Quaternions ///////////////////////////

	Q1 = Quaternion(1, 2,3,4);	// not zero
	Q2 = Quaternion(0, 0,0,0);	// zero
	Q3 = Quaternion(0, 0,0,0);	// zero
	if (Q1 == Q3) cout << "Q1 = 0 \n\n";
	if (Q1 != Q3) cout << "Q1 != 0 (as expected) \n\n";

	if (Q2 == Q3) cout << "Q2 = 0 (as expected) \n\n";
	if (Q2 != Q3) cout << "Q2 != 0 \n\n";

////////////////////////// Test == and != for Vectors ///////////////////////////

	Vector V1, V2, V3;

	V1 = Vector(1,2,3);	// not zero
	V2 = Vector(0,0,0);	// zero
	V3 = Vector(0,0,0);	// zero
	if (V1 == V3) cout << "V1 = 0 \n\n";
	if (V1 != V3) cout << "V1 != 0 (as expected) \n\n";

	if (V2 == V3) cout << "V2 = 0 (as expected) \n\n";
	if (V2 != V3) cout << "V2 != 0 \n\n";

/////////////////// Test Vector Cross Product //////////////////////////////////

	V1 = Vector(a, b, c);
	V2 = Vector(A, B, C);
	V3 = Cross(V1, V2);
	cout << "Cross(V1, V2) = " << V3 << "\n\n";

	V3 = Cross(V1,V1);	// expect zero
	cout << "Cross(V1, V1) = " << V3 << " (expect 0) \n\n";


/////////////////////////////////////////////////////

	return 0;

}


//////////////////////////////////////////////////////

/*

 * |  q  i  j  k 
-----------------
 q |  q  i  j  k 
 i |  i -q  k -j 
 j |  j -k -q  i 
 k |  k  j -i -q 



 * |  q  i  j  k  E  I  J  K 
-----------------------------
 q |  q  i  j  k  E  I  J  K 
 i |  i -q  k -j  I -E -K  J 
 j |  j -k -q  i  J  K -E -I 
 k |  k  j -i -q  K -J  I -E 
 E |  E -I -J -K -q  i  j  k 
 I |  I  E -K  J -i -q -k  j 
 J |  J  K  E -I -j  k -q -i 
 K |  K -J  I  E -k -j  i -q 



 *   |  e_0  e_1  e_2  e_3  e_4  e_5  e_6  e_7 
-----------------------------------------------
 e_0 |  e_0  e_1  e_2  e_3  e_4  e_5  e_6  e_7 
 e_1 |  e_1 -e_0  e_3 -e_2  e_5 -e_4 -e_7  e_6 
 e_2 |  e_2 -e_3 -e_0  e_1  e_6  e_7 -e_4 -e_5 
 e_3 |  e_3  e_2 -e_1 -e_0  e_7 -e_6  e_5 -e_4 
 e_4 |  e_4 -e_5 -e_6 -e_7 -e_0  e_1  e_2  e_3 
 e_5 |  e_5  e_4 -e_7  e_6 -e_1 -e_0 -e_3  e_2 
 e_6 |  e_6  e_7  e_4 -e_5 -e_2  e_3 -e_0 -e_1 
 e_7 |  e_7 -e_6  e_5  e_4 -e_3 -e_2  e_1 -e_0 

*/


