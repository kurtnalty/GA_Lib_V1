#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GA5_4_1_Routines.c"

/////////////////////////////////// Demonstrate GA5_4_1 Routines ////////////////////////////////


int main(void)
{
	
	GA5_4_1 r,s,t,u;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);

	Complex det_ref, det;
	det_ref = GA5_4_1_Determinant(r);
	printf("\nGA5_4_1 Routines Demo\n\n");

//////////////////////////////////////////////////////

// void GA5_4_1_PrintMV(GA5_4_1 v) 

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("GA5_4_1_PrintMV Demo\nr = ");  GA5_4_1_PrintMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void GA5_4_1_PrintlnMV(GA5_4_1 v) 

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("GA5_4_1_PrintlnMV Demo\nr = ");  GA5_4_1_PrintlnMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_OverBar(GA5_4_1 a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_OverBar(r);

	printf("GA5_4_1_OverBar Demo\n");
	printf("r =              ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_OverBar(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =       ");  GA5_4_1_PrintlnMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_UnderBar(GA5_4_1 a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_UnderBar(r);

	printf("GA5_4_1_UnderBar Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_UnderBar(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_OverBar(GA5_4_1_UnderBar(r));
	printf("t = GA5_4_1_OverBar(GA5_4_1_UnderBar(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Reverse(GA5_4_1 w) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Reverse(r);

	printf("GA5_4_1_Reverse Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Reverse(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Reverse(GA5_4_1_Reverse(r));
	printf("t = GA5_4_1_Reverse(GA5_4_1_Reverse(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Parity(GA5_4_1 w) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Parity(r);

	printf("GA5_4_1_Parity Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Parity(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Parity(GA5_4_1_Parity(r));
	printf("t = GA5_4_1_Parity(GA5_4_1_Parity(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Transpose(GA5_4_1 w) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Transpose(r);

	printf("GA5_4_1_Transpose Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Transpose(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Transpose(GA5_4_1_Transpose(r));
	printf("t = GA5_4_1_Transpose(GA5_4_1_Transpose(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Conjugation(GA5_4_1 w) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Conjugation(r);

	printf("GA5_4_1_Conjugation Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Conjugation(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Conjugation(GA5_4_1_Conjugation(r));
	printf("t = GA5_4_1_Conjugation(GA5_4_1_Conjugation(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Clifford_Conjugation(GA5_4_1 w) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Clifford_Conjugation(r);

	printf("GA5_4_1_Clifford_Conjugation Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Clifford_Conjugation(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Clifford_Conjugation(GA5_4_1_Clifford_Conjugation(r));
	printf("t = GA5_4_1_Clifford_Conjugation(GA5_4_1_Clifford_Conjugation(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Dual(GA5_4_1 w) ; 

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Dual(r);

	printf("GA5_4_1_Dual Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Dual(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	t = GA5_4_1_Dual(GA5_4_1_Dual(r));
	printf("t = GA5_4_1_Dual(GA5_4_1_Dual(r)) =        ");  GA5_4_1_PrintlnMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Dorst_Dual(GA5_4_1 a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Dorst_Dual(r);

	printf("GA5_4_1_Dorst_Dual Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Dorst_Dual(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Dorst_UnDual(GA5_4_1 a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	u = GA5_4_1_Dorst_UnDual(r);

	printf("GA5_4_1_Dorst_UnDual Demo\n");
	printf("r =               ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("u = GA5_4_1_Dorst_UnDual(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n");
	det = GA5_4_1_Determinant(u);
	if( (det.r == det_ref.r)&& (det.i == det_ref.i) ) printf("GA5_4_1_Determinant is conserved\n"); else printf("GA5_4_1_Determinant is not conserved\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*u =        ");  GA5_4_1_PrintlnMV(s);  printf("\n");

	s = GA5_4_1_Dorst_Dual(r);
	u = GA5_4_1_Dorst_UnDual(s);

	printf("GA5_4_1_Dorst_UnDual(GA5_4_1_Dorst_Dual(r)) = "); GA5_4_1_PrintlnMV(u); printf("\n");

	s = GA5_4_1_Dorst_UnDual(r);
	u = GA5_4_1_Dorst_Dual(s);

	printf("GA5_4_1_Dorst_Dual(GA5_4_1_Dorst_UnDual(r)) = "); GA5_4_1_PrintlnMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Add(GA5_4_1 u,GA5_4_1 v)

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	u = GA5_4_1_Add(r, s);

	printf("  r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("  s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	printf("r+s = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Subtract(GA5_4_1 u, GA5_4_1 v)

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	u = GA5_4_1_Subtract(r, s);

	printf("  r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("  s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	printf("r-s = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int GA5_4_1_Equal(GA5_4_1 u, GA5_4_1 v) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	t = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);

	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	printf("t = ");  GA5_4_1_PrintlnMV(t);  printf("\n");
	if (GA5_4_1_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (GA5_4_1_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int GA5_4_1_Not_Equal(GA5_4_1 u, GA5_4_1 v) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	t = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);

	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	printf("t = ");  GA5_4_1_PrintlnMV(t);  printf("\n");
	if (GA5_4_1_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (GA5_4_1_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Product(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Product(r,s);
	printf("u = r*s = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Divide_By_Constant(GA5_4_1 u, double a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	u = GA5_4_1_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Wedge(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Wedge(r,s);
	printf("u = GA5_4_1_Wedge(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_AntiWedge(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_AntiWedge(r,s);
	printf("u = GA5_4_1_AntiWedge(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Regressive(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Regressive(r,s);
	printf("u = GA5_4_1_Regressive(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Regressive_Via_Formula(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Regressive_Via_Formula(r,s);
	printf("u = GA5_4_1_Regressive_Via_Formula(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Lower_Right_Via_Formula(GA5_4_1 a, GA5_4_1 b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Lower_Right_Via_Formula(r,s);
	printf("u = GA5_4_1_Lower_Right_Via_Formula(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Expander(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Expander(r,s);
	printf("u = GA5_4_1_Expander(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Conserver(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Conserver(r,s);
	printf("u = GA5_4_1_Conserver(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Shrinker(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Shrinker(r,s);
	printf("u = GA5_4_1_Shrinker(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Symmetric(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Symmetric(r,s);
	printf("u = GA5_4_1_Symmetric(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_AntiSymmetric(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_AntiSymmetric(r,s);
	printf("u = GA5_4_1_AntiSymmetric(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Inner(const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Inner(r,s);
	printf("u = GA5_4_1_Inner(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Left_Contraction (const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Left_Contraction(r,s);
	printf("u = GA5_4_1_Left_Contraction(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Right_Contraction (const GA5_4_1 &a, const GA5_4_1 &b) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	printf("s = ");  GA5_4_1_PrintlnMV(s);  printf("\n");
	u = GA5_4_1_Right_Contraction(r,s);
	printf("u = GA5_4_1_Right_Contraction(r, s) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex GA5_4_1_Determinant(GA5_4_1 A) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	det = GA5_4_1_Determinant(r);
	printf("u = GA5_4_1_Determinant(r) = (%10.3e, %10.3e)\n\n",det.r, det.i); 


//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Adjugate(GA5_4_1 V) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	u = GA5_4_1_Adjugate(r);
	printf("u = GA5_4_1_Adjugate(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA5_4_1 GA5_4_1_Reciprocal(GA5_4_1 a) ;

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	printf("r = ");  GA5_4_1_PrintlnMV(r);  printf("\n");
	u = GA5_4_1_Reciprocal(r);
	printf("u = GA5_4_1_Reciprocal(r) = ");  GA5_4_1_PrintlnMV(u);  printf("\n\n");
	s = GA5_4_1_Product(r,u);
	printf("s = r*(1/r) = ");  GA5_4_1_PrintlnMV(s);  printf("\n\n");

////////////////////////////////////

//	verify the 64 determinant preserving GA5_4_1_Comp transforms

	GA5_4_1 MV1, MV2;
	int i;

	MV1 = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	det_ref = GA5_4_1_Determinant(MV1);

	printf("\nVerifying 64 GA5_4_1_Comp transforms\n");
	for (i=0;i<64;i++) 
{
		MV2 = GA5_4_1_Comp(MV1,i);
		det = GA5_4_1_Determinant(MV2);
		if((det.r != det_ref.r) || (det.i != det_ref.i)) printf("Failure at index %d \n", i);
		else printf("+ ");
}
	printf("\n");
	


////////////////////////////////////

//	verify the 360 determinant preserving GA5_4_1_Magic transforms

	MV1 = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	det_ref = GA5_4_1_Determinant(MV1);

	printf("\nVerifying 360 GA5_4_1_Magic transforms\n");
	for (i=0;i<360;i++) 
{
		MV2 = GA5_4_1_Magic(MV1,i);
		det = GA5_4_1_Determinant(MV2);
		if((det.r != det_ref.r) || (det.i != det_ref.i)) printf("Failure at index %d \n", i);
		else printf("+ ");
}
	printf("\n");
	

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = GA5_4_1_Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	s = GA5_4_1_Set(139,  149,151,157,163,167,
		173,179,181,191,193,   197,199,211,223,227,
		229,233,239,241,251,   257,263,269,271,277,
		281,283,293,307,311,   313);

) ; 

2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541

*/





