// Tests for CHO Algebra Software
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 4 September 2018
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

#include "CHO_Routines.cp"

//////////////////////////////////////////////////////

int main(void)
{
	int i;
	int CR_Count;

	C C1, C2, C3;
	Q H1, H2, H3;
	O O1, O2, O3;
	CH CH1, CH2, CH3;
	CO CO1, CO2, CO3;
	HO HO1, HO2, HO3;
	CHO CHO1, CHO2, CHO3;
	
	Fill_C_Basis() ;
	Fill_Q_Basis() ;
	Fill_O_Basis() ;
	Fill_CH_Basis() ;
	Fill_CO_Basis() ;
	Fill_HO_Basis() ;
	Fill_CHO_Basis() ;

//////////////////////////////////////////////////////

//
//	Test PrintSimple routines
//
	printf("\nComplex Basis\n");
	for (i=0; i< 2; i++) { PrintSimpleC  (  C_Basis[i]); printf("\t"); } 

	printf("\nQuaternion Basis\n");
	for (i=0; i< 4; i++) { PrintSimpleQ  (  Q_Basis[i]); printf("\t"); } 

	printf("\nOctonion Basis\n");
	for (i=0; i< 8; i++) { PrintSimpleO  (  O_Basis[i]); printf("\t"); } 

	printf("\nCxH Basis\n");
	for (i=0; i< 8; i++) { PrintSimpleCH ( CH_Basis[i]); printf("\t"); } 

	printf("\nCxO Basis\n");
	CR_Count = 0;
	for (i=0; i<16; i++) { PrintSimpleCO ( CO_Basis[i]); printf("\t"); if((++CR_Count &7) == 0) printf("\n"); } 

	printf("\nHxO Basis\n");
	CR_Count = 0;
	for (i=0; i<32; i++) { PrintSimpleHO ( HO_Basis[i]); printf("\t"); if((++CR_Count &7) == 0) printf("\n"); } 

	printf("\nCxHxO Basis\n");
	CR_Count = 0;
	for (i=0; i<64; i++) { PrintSimpleCHO(CHO_Basis[i]); printf("\t"); if((++CR_Count &7) == 0) printf("\n"); } 

//////////////////////////////////////////////////////

//
//	Test Set and print routines
//
	C1 = C(1, 2);

	H1 = Q(3,4,5,6);

	O1 = O(1,2,3,4,5,6,7,8);

	CH1 = CH(8,7,6,5,4,3,2,1);

	CO1 = CO(1, 2, 3, 4, 5, 6, 7, 8,
                     9,10,11,12,13,14,15,16);

	HO1 = HO(11,12,13,14,15,16,17,18,
                     19,20,21,22,23,24,25,26,
                     27,28,29,30,31,32,33,34,
                     35,36,37,38,39,40,41,42);

	CHO1 = CHO(11,12,13,14,15,16,17,18,
                       19,20,21,22,23,24,25,26,
                       27,28,29,30,31,32,33,34,
                       35,36,37,38,39,40,41,42,

		       43,44,45,46,47,48,49,50,
                       51,52,53,54,55,56,57,58,
                       59,60,61,62,63,64,65,66,
                       67,68,69,70,71,72,73,74);

	cout << "C1 = \n" << C1 << "\n\n";
	cout << "H1 = \n" << H1 << "\n\n";
	cout << "O1 = \n" << O1 << "\n\n";

	cout << "CH1 = \n" << CH1 << "\n\n";
	cout << "CO1 = \n" << CO1 << "\n\n";
	cout << "HO1 = \n" << HO1 << "\n\n";

	cout << "CHO1 = \n" << CHO1 << "\n\n";

//////////////////////////////////////////////////////

//
//	Test Zero, EQ and NE routines
//
	printf("\nTest Zero, EQ and NE routines\n\n");

	C2 = C();
	H2 = Q();
	O2 = O();

	CH2 = CH();
	CO2 = CO();
	HO2 = HO();

	CHO2 = CHO();

	if (C1 == C()) printf("C1 = 0\n");		if (C2 == C()) printf("C2 = 0\n");
	if (H1 == Q()) printf("H1 = 0\n");		if (H2 == Q()) printf("H2 = 0\n");
	if (O1 == O()) printf("O1 = 0\n");		if (O2 == O()) printf("O2 = 0\n");

	if (CH1 == CH()) printf("CH1 = 0\n");		if (CH2 == CH()) printf("CH2 = 0\n");
	if (CO1 == CO()) printf("CO1 = 0\n");		if (CO2 == CO()) printf("CO2 = 0\n");
	if (HO1 == HO()) printf("HO1 = 0\n");		if (HO2 == HO()) printf("HO2 = 0\n");

	if (CHO1 == CHO()) printf("CHO1 = 0\n");	if (CHO2 == CHO()) printf("CHO2 = 0\n");

	if (C1 != C()) printf("C1 != 0\n");		if (C2 != C()) printf("C2 != 0\n");
	if (H1 != Q()) printf("H1 != 0\n");		if (H2 != Q()) printf("H2 != 0\n");
	if (O1 != O()) printf("O1 != 0\n");		if (O2 != O()) printf("O2 != 0\n");

	if (CH1 != CH()) printf("CH1 != 0\n");		if (CH2 != CH()) printf("CH2 != 0\n");
	if (CO1 != CO()) printf("CO1 != 0\n");		if (CO2 != CO()) printf("CO2 != 0\n");
	if (HO1 != HO()) printf("HO1 != 0\n");		if (HO2 != HO()) printf("HO2 != 0\n");

	if (CHO1 != CHO()) printf("CHO1 != 0\n");	if (CHO2 != CHO()) printf("CHO2 != 0\n");

//////////////////////////////////////////////////////

//
//	Test Add, Subtract
//
	C2 = C1 + C1;	// C2 = 2*C1
	H2 = H1 + H1;
	O2 = O1 + O1;

	CH2 = CH1 + CH1;
	CO2 = CO1 + CO1;
	HO2 = HO1 + HO1;

	CHO2 = CHO1 + CHO1;

	cout << "C2 = \n" << C2 << "\n\n";
	cout << "H2 = \n" << H2 << "\n\n";
	cout << "O2 = \n" << O2 << "\n\n";

	cout << "CH2 = \n" << CH2 << "\n\n";
	cout << "CO2 = \n" << CO2 << "\n\n";
	cout << "HO2 = \n" << HO2 << "\n\n";

	cout << "CHO2 = \n" << CHO2 << "\n\n";

	C3 = C1 - C2;	// C3 = -C1
	H3 = H1 - H2;
	O3 = O1 - O2;

	CH3 = CH1 - CH2;
	CO3 = CO1 - CO2;
	HO3 = HO1 - HO2;

	CHO3 = CHO1 - CHO2;

	cout << "C3 = \n" << C3 << "\n\n";
	cout << "H3 = \n" << H3 << "\n\n";
	cout << "O3 = \n" << O3 << "\n\n";

	cout << "CH3 = \n" << CH3 << "\n\n";
	cout << "CO3 = \n" << CO3 << "\n\n";
	cout << "HO3 = \n" << HO3 << "\n\n";

	cout << "CHO3 = \n" << CHO3 << "\n\n";

//////////////////////////////////////////////////////

//
//	Test Assignment, EQ, NE for pairs
//
	C2 = C1;
	H2 = H1;
	O2 = O1;

	CH2 = CH1;
	CO2 = CO1;
	HO2 = HO1;

	CHO2 = CHO1;

	if(C1 == C2) printf("C1 = C2\n");		if(C1 == C3) printf("C1 = C3\n"); 
	if(H1 == H2) printf("H1 = H2\n");		if(H1 == H3) printf("H1 = H3\n"); 
	if(O1 == O2) printf("O1 = O2\n");		if(O1 == O3) printf("O1 = O3\n"); 

	if(CH1 == CH2) printf("CH1 = CH2\n");		if(CH1 == CH3) printf("CH1 = CH3\n"); 
	if(CO1 == CO2) printf("CO1 = CO2\n");		if(CO1 == CO3) printf("CO1 = CO3\n"); 
	if(HO1 == HO2) printf("HO1 = HO2\n");		if(HO1 == HO3) printf("HO1 = HO3\n"); 

	if(CHO1 == CHO2) printf("CHO1 = CHO2\n");	if(CHO1 == CHO3) printf("CHO1 = CHO3\n"); 

	printf("\n");

	if(C1 != C2) printf("C1 != C2\n");		if(C1 != C3) printf("C1 != C3\n"); 
	if(H1 != H2) printf("H1 != H2\n");		if(H1 != H3) printf("H1 != H3\n"); 
	if(O1 != O2) printf("O1 != O2\n");		if(O1 != O3) printf("O1 != O3\n"); 

	if(CH1 != CH2) printf("CH1 != CH2\n");		if(CH1 != CH3) printf("CH1 != CH3\n"); 
	if(CO1 != CO2) printf("CO1 != CO2\n");		if(CO1 != CO3) printf("CO1 != CO3\n"); 
	if(HO1 != HO2) printf("HO1 != HO2\n");		if(HO1 != HO3) printf("HO1 != HO3\n"); 

	if(CHO1 != CHO2) printf("CHO1 != CHO2\n");	if(CHO1 != CHO3) printf("CHO1 != CHO3\n"); 

//////////////////////////////////////////////////////

//
//	Test Product
//
	C2 = C1*C1;	// C2 = C1*C1
	H2 = H1*H1;
	O2 = O1*O1;

	CH2 = CH1*CH1;
	CO2 = CO1*CO1;
	HO2 = HO1*HO1;

	CHO2 = CHO2*CHO2;

	cout << "C2 = \n" << C2 << "\n\n";
	cout << "H2 = \n" << H2 << "\n\n";
	cout << "O2 = \n" << O2 << "\n\n";

	cout << "CH2 = \n" << CH2 << "\n\n";
	cout << "CO2 = \n" << CO2 << "\n\n";
	cout << "HO2 = \n" << HO2 << "\n\n";

	cout << "CHO2 = \n" << CHO2 << "\n\n";


//////////////////////////////////////////////////////

	return 0;
}

