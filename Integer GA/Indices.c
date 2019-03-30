#include <stdio.h>
#include <stdlib.h>

	long Matches[64] = {
		0x00000000, 0x03FFFFC0, 0x0CDE7B30, 0x0F2184F0, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 
		0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x30E1870C, 0x331E78CC, 0x3C3FFC3C, 0x3FC003FC, 
		0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x5119988A, 0x52E6674A, 0x5DC7E3BA, 0x5E381C7A, 
		0x61F81F86, 0x6207E046, 0x6D2664B6, 0x6ED99B76, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 
		0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x91266489, 0x92D99B49, 0x9DF81FB9, 0x9E07E079, 
		0xA1C7E385, 0xA2381C45, 0xAD1998B5, 0xAEE66775, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 
		0xC03FFC03, 0xC3C003C3, 0xCCE18733, 0xCF1E78F3, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 
		0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xF0DE7B0F, 0xF32184CF, 0xFC00003F, 0xFFFFFFFF
	} ;

	long MisMatches[64] = {
		0x044B2DDF, 0x07B4D21F, 0x089556EF, 0x0B6AA92F, 0x11269B77, 0x12D964B7, 0x1DF8E047, 0x1E071F87, 
		0x21C71C7B, 0x2238E3BB, 0x2D19674B, 0x2EE6988B, 0x34AAAAD3, 0x37555513, 0x3874D1E3, 0x3B8B2E23, 
		0x403F03FD, 0x43C0FC3D, 0x4CE178CD, 0x4F1E870D, 0x5552B555, 0x56AD4A95, 0x598CCE65, 0x5A7331A5, 
		0x65B33259, 0x664CCD99, 0x696D4969, 0x6A92B6A9, 0x70DE84F1, 0x73217B31, 0x7C00FFC1, 0x7FFF0001,
		0x8000FFFE, 0x83FF003E, 0x8CDE84CE, 0x8F217B0E, 0x956D4956, 0x9692B696, 0x99B33266, 0x9A4CCDA6, 
		0xA58CCE5A, 0xA673319A, 0xA952B56A, 0xAAAD4AAA, 0xB0E178F2, 0xB31E8732, 0xBC3F03C2, 0xBFC0FC02, 
		0xC474D1DC, 0xC78B2E1C, 0xC8AAAAEC, 0xCB55552C, 0xD1196774, 0xD2E698B4, 0xDDC71C44, 0xDE38E384, 
		0xE1F8E078, 0xE2071FB8, 0xED269B48, 0xEED96488, 0xF49556D0, 0xF76AA910, 0xF84B2DE0, 0xFBB4D220
	} ;

	long All_Indices[128] = {
	0x00000000, 0x03FFFFC0, 0x044B2DDF, 0x07B4D21F, 0x089556EF, 0x0B6AA92F, 0x0CDE7B30, 0x0F2184F0, 
	0x11269B77, 0x12D964B7, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 0x1DF8E047, 0x1E071F87, 
	0x21C71C7B, 0x2238E3BB, 0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x2D19674B, 0x2EE6988B, 
	0x30E1870C, 0x331E78CC, 0x34AAAAD3, 0x37555513, 0x3874D1E3, 0x3B8B2E23, 0x3C3FFC3C, 0x3FC003FC, 
	0x403F03FD, 0x43C0FC3D, 0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x4CE178CD, 0x4F1E870D, 
	0x5119988A, 0x52E6674A, 0x5552B555, 0x56AD4A95, 0x598CCE65, 0x5A7331A5, 0x5DC7E3BA, 0x5E381C7A, 
	0x61F81F86, 0x6207E046, 0x65B33259, 0x664CCD99, 0x696D4969, 0x6A92B6A9, 0x6D2664B6, 0x6ED99B76, 
	0x70DE84F1, 0x73217B31, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 0x7C00FFC1, 0x7FFF0001, 
	
	0x8000FFFE, 0x83FF003E, 0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x8CDE84CE, 0x8F217B0E, 
	0x91266489, 0x92D99B49, 0x956D4956, 0x9692B696, 0x99B33266, 0x9A4CCDA6, 0x9DF81FB9, 0x9E07E079, 
	0xA1C7E385, 0xA2381C45, 0xA58CCE5A, 0xA673319A, 0xA952B56A, 0xAAAD4AAA, 0xAD1998B5, 0xAEE66775, 
	0xB0E178F2, 0xB31E8732, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 0xBC3F03C2, 0xBFC0FC02, 
	0xC03FFC03, 0xC3C003C3, 0xC474D1DC, 0xC78B2E1C, 0xC8AAAAEC, 0xCB55552C, 0xCCE18733, 0xCF1E78F3, 
	0xD1196774, 0xD2E698B4, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 0xDDC71C44, 0xDE38E384, 
	0xE1F8E078, 0xE2071FB8, 0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xED269B48, 0xEED96488, 
	0xF0DE7B0F, 0xF32184CF, 0xF49556D0, 0xF76AA910, 0xF84B2DE0, 0xFBB4D220, 0xFC00003F, 0xFFFFFFFF
} ;

int	InGroup(long index)
{
	int i,Num_Match;

	Num_Match = 0;

	for (i=0; i<128; i++ ) if(All_Indices[i] == index) Num_Match++;

	if(Num_Match > 1) printf("Error in InGroup - Num_Match = %d \n", Num_Match);

	return (Num_Match == 1);
}


int	InMatches(long index)
{
	int i,Num_Match;

	Num_Match = 0;

	for (i=0; i<64; i++ ) if(Matches[i] == index) Num_Match++;

	if(Num_Match > 1) printf("Error in InGroup - Num_Match = %d \n", Num_Match);

	return (Num_Match == 1);
}


int	InMisMatches(long index)
{
	int i,Num_Match;

	Num_Match = 0;

	for (i=0; i<64; i++ ) if(MisMatches[i] == index) Num_Match++;

	if(Num_Match > 1) printf("Error in InGroup - Num_Match = %d \n", Num_Match);

	return (Num_Match == 1);
}

//////////////////////////////

int main(void)
{

	int i, j;

	long Index_20s[ 6] = {67108800, 1010826300, 1573381050, 1859754870, 2003457774, 2075405790 } ;
	long Index_16s[16] = { 129290783,  191539503,  316236983,  503783303,  566697083,  756639563,  883600083,  947180003, 
			      1077871613, 1289844941, 1431483733, 1502400101, 1716309401, 1787999913, 1931574065, 2147418113 } ;
	long Index_12s[10] = { 253854960, 378685800, 441201240, 629944740, 693258900, 820086540, 1148464674, 1219122450, 1360631946, 1644683334 } ;
	long Index_1s[16] = {  72035807,  144004847,  287742839,  502849607,  574153659,  786864267,  928339219,  998977059, 
			     1136720957, 1327400717, 1454197397, 1517498789, 1706242649, 1768769897, 1893631217, 2080440257 } ;
	long Index_0s[16] = {          0,  215907120,  359511720,  431213976,  645123684,  716027220,  857635020, 1069548540, 
			      1200345570, 1263905490, 1390831434, 1580735610, 1643650950, 1831232694, 1955965230, 2018234910 } ;

/*
 0s Conserved determinant (even indices)
 1s Conserved complement  ( odd indices)
12s Conserved determinant (even indices)
16s Conserved complement  ( odd indices)
20s Conserved determinant (even indices)
*/

// print the indices above in hex

	for (i=0; i< 6; i++) printf("%8lX ", Index_20s[i]);
	printf("\n\n");

	for (i=0; i<16; i++) printf("%8lX ", Index_16s[i]);
	printf("\n\n");

	for (i=0; i<10; i++) printf("%8lX ", Index_12s[i]);
	printf("\n\n");

	for (i=0; i<16; i++) printf("%8lX ", Index_1s[i]);
	printf("\n\n");

	for (i=0; i<16; i++) printf("%8lX ", Index_0s[i]);
	printf("\n\n");

	for (i=0; i < 128; i++) for (j=0; j<128; j++) 
		if(!InGroup(All_Indices[i]^All_Indices[j])) printf("Not in group at i= %d  j= %d \n",i,j);

	for (i=0; i < 64; i++) for (j=0; j<64; j++) 
		if(!InMatches(Matches[i]^Matches[j])) printf("Not in Matches group at i= %d  j= %d \n",i,j);

//	for (i=0; i < 64; i++) for (j=0; j<64; j++) 
//		if(!InMisMatches(MisMatches[i]^MisMatches[j])) printf("Not in MisMatchesgroup at i= %d  j= %d \n",i,j);

// as expect, a group for InGroup, InMatches

/*
	for (i=0;i<64;i++) {
		printf("0x%8lX, ",All_Indices[i]);
		if(((i+1)%8) == 0) printf("\n");
	}
	printf("\n");
	for (i=0;i<64;i++) {
		printf("0x%8lX, ",0xFFFFFFFF^All_Indices[63-i]);
		if(((i+1)%8) == 0) printf("\n");
	}
	printf("\n");
*/
	return 0;
}

/* some mirror symmetry. 
03FFFFC0 3C3FFC3C 5DC7E3BA 6ED99B76 776A56EE 7BB42DDE 

07B4D21F 0B6AA92F 12D964B7 1E071F87 21C71C7B 2D19674B 34AAAAD3 3874D1E3 403F03FD 4CE178CD 5552B555 598CCE65 664CCD99 6A92B6A9 73217B31 7FFF0001 

0F2184F0 16924968 1A4C3258 258C31A4 29524A94 30E1870C 44742E22 48AA5512 5119988A 6207E046 

044B2DDF 089556EF 11269B77 1DF8E047 2238E3BB 2EE6988B 37555513 3B8B2E23 43C0FC3D 4F1E870D 56AD4A95 5A7331A5 65B33259 696D4969 70DE84F1 7C00FFC1 

00000000 0CDE7B30 156DB6A8 19B3CD98 2673CE64 2AADB554 331E78CC 3FC003FC 478BD1E2 4B55AAD2 52E6674A 5E381C7A 61F81F86 6D2664B6 7495A92E 784BD21E 

*/

/*
	0x00000000, 0x03FFFFC0, 0x044B2DDF, 0x07B4D21F, 0x089556EF, 0x0B6AA92F, 0x0CDE7B30, 0x0F2184F0, 
	0x11269B77, 0x12D964B7, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 0x1DF8E047, 0x1E071F87, 
	0x21C71C7B, 0x2238E3BB, 0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x2D19674B, 0x2EE6988B, 
	0x30E1870C, 0x331E78CC, 0x34AAAAD3, 0x37555513, 0x3874D1E3, 0x3B8B2E23, 0x3C3FFC3C, 0x3FC003FC, 
	0x403F03FD, 0x43C0FC3D, 0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x4CE178CD, 0x4F1E870D, 
	0x5119988A, 0x52E6674A, 0x5552B555, 0x56AD4A95, 0x598CCE65, 0x5A7331A5, 0x5DC7E3BA, 0x5E381C7A, 
	0x61F81F86, 0x6207E046, 0x65B33259, 0x664CCD99, 0x696D4969, 0x6A92B6A9, 0x6D2664B6, 0x6ED99B76, 
	0x70DE84F1, 0x73217B31, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 0x7C00FFC1, 0x7FFF0001, 
	
	0x8000FFFE, 0x83FF003E, 0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x8CDE84CE, 0x8F217B0E, 
	0x91266489, 0x92D99B49, 0x956D4956, 0x9692B696, 0x99B33266, 0x9A4CCDA6, 0x9DF81FB9, 0x9E07E079, 
	0xA1C7E385, 0xA2381C45, 0xA58CCE5A, 0xA673319A, 0xA952B56A, 0xAAAD4AAA, 0xAD1998B5, 0xAEE66775, 
	0xB0E178F2, 0xB31E8732, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 0xBC3F03C2, 0xBFC0FC02, 
	0xC03FFC03, 0xC3C003C3, 0xC474D1DC, 0xC78B2E1C, 0xC8AAAAEC, 0xCB55552C, 0xCCE18733, 0xCF1E78F3, 
	0xD1196774, 0xD2E698B4, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 0xDDC71C44, 0xDE38E384, 
	0xE1F8E078, 0xE2071FB8, 0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xED269B48, 0xEED96488, 
	0xF0DE7B0F, 0xF32184CF, 0xF49556D0, 0xF76AA910, 0xF84B2DE0, 0xFBB4D220, 0xFC00003F, 0xFFFFFFFF


*/

