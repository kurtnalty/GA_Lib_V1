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
#include <string.h>

//////////////////////////////////////////////////////

#include "CHO_Routines.c"

//////////////////////////////////////////////////////

int main(void)
{
	int i;
	int CR_Count;

	C C1, C2, C3;
	H H1, H2, H3;
	O O1, O2, O3;
	CH CH1, CH2, CH3;
	CO CO1, CO2, CO3;
	HO HO1, HO2, HO3;
	CHO CHO1, CHO2, CHO3;
	
	Fill_C_Basis() ;
	Fill_H_Basis() ;
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
	for (i=0; i< 4; i++) { PrintSimpleH  (  H_Basis[i]); printf("\t"); } 

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
	C1 = C_Set(1, 2);

	H1 = H_Set(3,4,5,6);

	O1 = O_Set(1,2,3,4,5,6,7,8);

	CH1 = CH_Set(8,7,6,5,4,3,2,1);

	CO1 = CO_Set(1, 2, 3, 4, 5, 6, 7, 8,
                     9,10,11,12,13,14,15,16);

	HO1 = HO_Set(11,12,13,14,15,16,17,18,
                     19,20,21,22,23,24,25,26,
                     27,28,29,30,31,32,33,34,
                     35,36,37,38,39,40,41,42);

	CHO1 = CHO_Set(11,12,13,14,15,16,17,18,
                       19,20,21,22,23,24,25,26,
                       27,28,29,30,31,32,33,34,
                       35,36,37,38,39,40,41,42,

		       43,44,45,46,47,48,49,50,
                       51,52,53,54,55,56,57,58,
                       59,60,61,62,63,64,65,66,
                       67,68,69,70,71,72,73,74);

	printf("C1 = \n"); Print_C(C1); printf("\n\n");
	printf("H1 = \n"); Print_H(H1); printf("\n\n");
	printf("O1 = \n"); Print_O(O1); printf("\n\n");

	printf("CH1 = \n"); Print_CH(CH1); printf("\n\n");
	printf("CO1 = \n"); Print_CO(CO1); printf("\n\n");
	printf("HO1 = \n"); Print_HO(HO1); printf("\n\n");

	printf("CHO1 = \n"); Print_CHO(CHO1); printf("\n\n");

//////////////////////////////////////////////////////

//
//	Test Zero, EQ and NE routines
//
	printf("\nTest Zero, EQ and NE routines\n\n");

	C2 = C_Zero();
	H2 = H_Zero();
	O2 = O_Zero();

	CH2 = CH_Zero();
	CO2 = CO_Zero();
	HO2 = HO_Zero();

	CHO2 = CHO_Zero();

	if (C_EQ_Zero(C1)) printf("C1 = 0\n");		if (C_EQ_Zero(C2)) printf("C2 = 0\n");
	if (H_EQ_Zero(H1)) printf("H1 = 0\n");		if (H_EQ_Zero(H2)) printf("H2 = 0\n");
	if (O_EQ_Zero(O1)) printf("O1 = 0\n");		if (O_EQ_Zero(O2)) printf("O2 = 0\n");

	if (CH_EQ_Zero(CH1)) printf("CH1 = 0\n");	if (CH_EQ_Zero(CH2)) printf("CH2 = 0\n");
	if (CO_EQ_Zero(CO1)) printf("CO1 = 0\n");	if (CO_EQ_Zero(CO2)) printf("CO2 = 0\n");
	if (HO_EQ_Zero(HO1)) printf("HO1 = 0\n");	if (HO_EQ_Zero(HO2)) printf("HO2 = 0\n");

	if (CHO_EQ_Zero(CHO1)) printf("CHO1 = 0\n");	if (CHO_EQ_Zero(CHO2)) printf("CHO2 = 0\n");

	if (C_NE_Zero(C1)) printf("C1 != 0\n");		if (C_NE_Zero(C2)) printf("C2 != 0\n");
	if (H_NE_Zero(H1)) printf("H1 != 0\n");		if (H_NE_Zero(H2)) printf("H2 != 0\n");
	if (O_NE_Zero(O1)) printf("O1 != 0\n");		if (O_NE_Zero(O2)) printf("O2 != 0\n");

	if (CH_NE_Zero(CH1)) printf("CH1 != 0\n");	if (CH_NE_Zero(CH2)) printf("CH2 != 0\n");
	if (CO_NE_Zero(CO1)) printf("CO1 != 0\n");	if (CO_NE_Zero(CO2)) printf("CO2 != 0\n");
	if (HO_NE_Zero(HO1)) printf("HO1 != 0\n");	if (HO_NE_Zero(HO2)) printf("HO2 != 0\n");

	if (CHO_NE_Zero(CHO1)) printf("CHO1 != 0\n");	if (CHO_NE_Zero(CHO2)) printf("CHO2 != 0\n");

//////////////////////////////////////////////////////

//
//	Test Add, Subtract
//
	C2 = C_Add(C1, C1);	// C2 = 2*C1
	H2 = H_Add(H1, H1);
	O2 = O_Add(O1, O1);

	CH2 = CH_Add(CH1, CH1);
	CO2 = CO_Add(CO1, CO1);
	HO2 = HO_Add(HO1, HO1);

	CHO2 = CHO_Add(CHO1, CHO1);

	printf("C2 = \n"); Print_C(C2); printf("\n\n");
	printf("H2 = \n"); Print_H(H2); printf("\n\n");
	printf("O2 = \n"); Print_O(O2); printf("\n\n");

	printf("CH2 = \n"); Print_CH(CH2); printf("\n\n");
	printf("CO2 = \n"); Print_CO(CO2); printf("\n\n");
	printf("HO2 = \n"); Print_HO(HO2); printf("\n\n");

	printf("CHO2 = \n"); Print_CHO(CHO2); printf("\n\n");

	C3 = C_Subtract(C1,C2);	// C3 = -C1
	H3 = H_Subtract(H1,H2);
	O3 = O_Subtract(O1,O2);

	CH3 = CH_Subtract(CH1,CH2);
	CO3 = CO_Subtract(CO1,CO2);
	HO3 = HO_Subtract(HO1,HO2);

	CHO3 = CHO_Subtract(CHO1,CHO2);

	printf("C3 = \n"); Print_C(C3); printf("\n\n");
	printf("H3 = \n"); Print_H(H3); printf("\n\n");
	printf("O3 = \n"); Print_O(O3); printf("\n\n");

	printf("CH3 = \n"); Print_CH(CH3); printf("\n\n");
	printf("CO3 = \n"); Print_CO(CO3); printf("\n\n");
	printf("HO3 = \n"); Print_HO(HO3); printf("\n\n");

	printf("CHO3 = \n"); Print_CHO(CHO3); printf("\n\n");

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

	if(C_EQ(C1,C2)) printf("C1 = C2\n");		if(C_EQ(C1,C3)) printf("C1 = C3\n"); 
	if(H_EQ(H1,H2)) printf("H1 = H2\n");		if(H_EQ(H1,H3)) printf("H1 = H3\n"); 
	if(O_EQ(O1,O2)) printf("O1 = O2\n");		if(O_EQ(O1,O3)) printf("O1 = O3\n"); 

	if(CH_EQ(CH1,CH2)) printf("CH1 = CH2\n");	if(CH_EQ(CH1,CH3)) printf("CH1 = CH3\n"); 
	if(CO_EQ(CO1,CO2)) printf("CO1 = CO2\n");	if(CO_EQ(CO1,CO3)) printf("CO1 = CO3\n"); 
	if(HO_EQ(HO1,HO2)) printf("HO1 = HO2\n");	if(HO_EQ(HO1,HO3)) printf("HO1 = HO3\n"); 

	if(CHO_EQ(CHO1,CHO2)) printf("CHO1 = CHO2\n");	if(CHO_EQ(CHO1,CHO3)) printf("CHO1 = CHO3\n"); 

	printf("\n");

	if(C_NE(C1,C2)) printf("C1 != C2\n");		if(C_NE(C1,C3)) printf("C1 != C3\n"); 
	if(H_NE(H1,H2)) printf("H1 != H2\n");		if(H_NE(H1,H3)) printf("H1 != H3\n"); 
	if(O_NE(O1,O2)) printf("O1 != O2\n");		if(O_NE(O1,O3)) printf("O1 != O3\n"); 

	if(CH_NE(CH1,CH2)) printf("CH1 != CH2\n");	if(CH_NE(CH1,CH3)) printf("CH1 != CH3\n"); 
	if(CO_NE(CO1,CO2)) printf("CO1 != CO2\n");	if(CO_NE(CO1,CO3)) printf("CO1 != CO3\n"); 
	if(HO_NE(HO1,HO2)) printf("HO1 != HO2\n");	if(HO_NE(HO1,HO3)) printf("HO1 != HO3\n"); 

	if(CHO_NE(CHO1,CHO2)) printf("CHO1 != CHO2\n");	if(CHO_NE(CHO1,CHO3)) printf("CHO1 != CHO3\n"); 

//////////////////////////////////////////////////////

//
//	Test Product
//
	C2 = C_Product(C1, C1);	// C2 = C1*C1
	H2 = H_Product(H1, H1);
	O2 = O_Product(O1, O1);

	CH2 = CH_Product(CH1, CH1);
	CO2 = CO_Product(CO1, CO1);
	HO2 = HO_Product(HO1, HO1);

	CHO2 = CHO_Product(CHO1, CHO1);

	printf("C2 = \n"); Print_C(C2); printf("\n\n");
	printf("H2 = \n"); Print_H(H2); printf("\n\n");
	printf("O2 = \n"); Print_O(O2); printf("\n\n");

	printf("CH2 = \n"); Print_CH(CH2); printf("\n\n");
	printf("CO2 = \n"); Print_CO(CO2); printf("\n\n");
	printf("HO2 = \n"); Print_HO(HO2); printf("\n\n");

	printf("CHO2 = \n"); Print_CHO(CHO2); printf("\n\n");

//////////////////////////////////////////////////////

	return 0;
}


//////////////////////////////////////////////////////


/*

C C_Product(C a, C b) 
H H_Product(H a, H b) 
O O_Product(O a, O b) 
CH CH_Product(CH a, CH b) 
CO CO_Product(CO a, CO b) 
HO HO_Product(HO a, HO b) 
CHO CHO_Product(CHO a, CHO b) 
CHO GA3EO_Product(CHO a, CHO b) 

*/



