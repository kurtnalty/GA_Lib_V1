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

#include "Quaternion_Octonion_Routines.c"

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
	int i,j;

	Quaternion Quaternion_Basis[4];
	Octonion  Octonion_Basis[8];

	Quaternion_Basis[0] = Set_Quaternion(1, 0,0,0);
	Quaternion_Basis[1] = Set_Quaternion(0, 1,0,0);
	Quaternion_Basis[2] = Set_Quaternion(0, 0,1,0);
	Quaternion_Basis[3] = Set_Quaternion(0, 0,0,1);

	Octonion_Basis[0] = Set_Octonion(1, 0,0,0, 0,0,0, 0);
	Octonion_Basis[1] = Set_Octonion(0, 1,0,0, 0,0,0, 0);
	Octonion_Basis[2] = Set_Octonion(0, 0,1,0, 0,0,0, 0);
	Octonion_Basis[3] = Set_Octonion(0, 0,0,1, 0,0,0, 0);
	Octonion_Basis[4] = Set_Octonion(0, 0,0,0, 1,0,0, 0);
	Octonion_Basis[5] = Set_Octonion(0, 0,0,0, 0,1,0, 0);
	Octonion_Basis[6] = Set_Octonion(0, 0,0,0, 0,0,1, 0);
	Octonion_Basis[7] = Set_Octonion(0, 0,0,0, 0,0,0, 1);


//	print quaternion product multiplication table

	printf(" * | ");   for (i=0;i<4;i++) PrintSimpleQuaternion(Quaternion_Basis[i]);    printf("\n"); 
	printf("-----------------\n");

	for (i=0;i<4;i++) {
		PrintSimpleQuaternion(Quaternion_Basis[i]);  printf("| ");
		for(j=0;j<4;j++) {
			PrintSimpleQuaternion(Product_Quaternion(Quaternion_Basis[i],Quaternion_Basis[j]));
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
			PrintSimpleOctonion(Product_Octonion(Octonion_Basis[i],Octonion_Basis[j]));
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
			PrintFancierOctonion(Product_Octonion(Octonion_Basis[i],Octonion_Basis[j]));
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
			PrintSimpleOctonion(Product_Octonion(Octonion_Basis[i],Octonion_Basis[j]));
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
			PrintFancierOctonion(Product_Octonion(Octonion_Basis[i],Octonion_Basis[j]));
		}
		printf("\n");
	}
	printf("\n\n\n");

	printf("\n***********************************************\n\n");

/////////////// Test Vector Routines //////////////////////////

	Vector V1, V2, V3;
	V1 = Set_Vector(2,3, 6);  // mag 7
	V2 = Set_Vector(3,4,12);  // mag 13
	printf("V1 = "); Print_Vector(V1); printf("  Mag(V1) = %lf \n", Mag_Vector(V1));
	printf("V2 = "); Print_Vector(V2); printf("  Mag(V2) = %lf \n\n", Mag_Vector(V2));

	V3 = Add_Vector(V1, V2);
	Print_Vector(V1); printf(" + "); Print_Vector(V2); printf(" = "); Print_Vector(V3); printf("\n\n");

	V3 = Sub_Vector(V1, V2);
	Print_Vector(V1); printf(" - "); Print_Vector(V2); printf(" = "); Print_Vector(V3); printf("\n\n");

	V3 = Cross(V1, V2);
	Print_Vector(V1); printf(" x "); Print_Vector(V2); printf(" = "); Print_Vector(V3); printf("\n\n");

	V3 = Cross(V1, V1);
	Print_Vector(V1); printf(" x "); Print_Vector(V1); printf(" = "); Print_Vector(V3); printf(" (expect 0) \n\n");

	printf("V1.V2 = %lf \n\n", Dot_Vector(V1,V2) );

	printf("\n***********************************************\n\n");

/////////////// Test Quaternion Routines ///////////////

	Quaternion Q1, Q2, Q3, Q4;
	Q1 = Set_Quaternion(2,4,5, 6);  // mag 9
	Q2 = Set_Quaternion(5,6,8,10);  // mag 15
	printf("Q1 = "); Print_Quaternion(Q1); printf("  Mag(Q1) = %lf \n", Mag_Quaternion(Q1));
	printf("Q2 = "); Print_Quaternion(Q2); printf("  Mag(Q2) = %lf \n\n", Mag_Quaternion(Q2));

	Q3 = Add_Quaternion(Q1, Q2);
	Print_Quaternion(Q1); printf(" + "); Print_Quaternion(Q2); printf(" = "); Print_Quaternion(Q3); printf("\n\n");

	Q3 = Sub_Quaternion(Q1, Q2);
	Print_Quaternion(Q1); printf(" - "); Print_Quaternion(Q2); printf(" = "); Print_Quaternion(Q3); printf("\n\n");

	Q3 = Product_Quaternion(Q1, Q2);
	Print_Quaternion(Q1); printf(" x "); Print_Quaternion(Q2); printf(" = "); Print_Quaternion(Q3); printf("\n\n");

	printf("  Mag(Q3) = %lf (expect 135 = 9*15)\n\n", Mag_Quaternion(Q3));

	Q3 = Conj_Quaternion(Q1);
	printf("    (Q1) = "); Print_Quaternion(Q1); printf("  \n");
	printf("Conj(Q1) = "); Print_Quaternion(Q3); printf("  \n\n");

	Q4 = Product_Quaternion(Q1,Q3);
	Print_Quaternion(Q1); printf(" x "); Print_Quaternion(Q3); printf(" = "); Print_Quaternion(Q4); printf(" \n\n");

	printf("Q1.Q2 = %lf \n\n", Dot_Quaternion(Q1,Q2) );

	Q3 = Reciprocol_Quaternion(Q1);
	Q4 = Product_Quaternion(Q1,Q3);
	printf("     (Q1) = "); Print_Quaternion(Q1); printf("  \n");
	printf("   (1/Q1) = "); Print_Quaternion(Q3); printf("  \n");
	printf("Q1*(1/Q1) = "); Print_Quaternion(Q4); printf("  \n\n");

	printf("\n***********************************************\n\n");

///////////// Test Octonion Routines ////////////////

	Octonion O1, O2, O3, O4;
	O1 = Set_Octonion(1,2,3,4,5,7, 8,11);  // mag 17
	O2 = Set_Octonion(3,4,5,6,7,9,12,13);  // mag 23
	printf("O1 = "); Print_Octonion(O1); printf("  Mag(O1) = %lf (expect 17)\n", Mag_Octonion(O1));
	printf("O2 = "); Print_Octonion(O2); printf("  Mag(O2) = %lf (expect 23)\n\n", Mag_Octonion(O2));

	O3 = Add_Octonion(O1, O2);
	Print_Octonion(O1); printf(" + "); Print_Octonion(O2); printf(" = "); Print_Octonion(O3); printf("\n\n");

	O3 = Sub_Octonion(O1, O2);
	Print_Octonion(O1); printf(" - "); Print_Octonion(O2); printf(" = "); Print_Octonion(O3); printf("\n\n");

	O3 = Product_Octonion(O1, O2);
	Print_Octonion(O1); printf(" x "); Print_Octonion(O2); printf(" = "); Print_Octonion(O3); printf("\n\n");

	printf("  Mag(O3) = %lf (expect 391 = 17*23)\n\n", Mag_Octonion(O3));

	O3 = Conj_Octonion(O1);
	printf("    (O1) = "); Print_Octonion(O1); printf("  \n");
	printf("Conj(O1) = "); Print_Octonion(O3); printf("  \n\n");

	O4 = Product_Octonion(O1,O3);
	Print_Octonion(O1); printf(" x "); Print_Octonion(O3); printf(" = "); Print_Octonion(O4); printf(" \n\n");

	printf("O1.O2 = %lf \n\n", Dot_Octonion(O1,O2) );

	O3 = Reciprocol_Octonion(O1);
	O4 = Product_Octonion(O1,O3);
	printf("     (O1) = "); Print_Octonion(O1); printf("  \n");
	printf("   (1/O1) = "); Print_Octonion(O3); printf("  \n");
	printf("O1*(1/O1) = "); Print_Octonion(O4); printf("  \n\n");

	printf("\n***********************************************\n\n");

/////////////////////////////////////////////////////

	return 0;

}



