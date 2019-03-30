#include <stdio.h>
#include <stdlib.h>

#include "GA2E_Routines.c"

/////////////////////////////////// Demonstrate GA2E Routines ////////////////////////////////


int main(void)
{
	

	GA2E r,s,t,u;

	r = GA2E_Set( 3,  5,  7, 11) ; 
	s = GA2E_Set(13, 17, 19, 23) ;
	t = GA2E_Set(29, 31, 37, 41) ; 


	double det_ref, det;
	det_ref = GA2E_Determinant(r);
	printf("\nGA2E Routines Demo\n\n");

//////////////////////////////////////////////////////

// void GA2E_PrintlnMV(GA2E v) 

	r = GA2E_Set(  3,    5,  7,   11) ; 
	printf("GA2E_PrintlnMV Demo\nr = ");  GA2E_PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void GA2E_PrintMV(GA2E v) 

	r = GA2E_Set(  3,    5,  7,   11) ; 
	printf("GA2E_PrintMV Demo\nr = ");  GA2E_PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_OverBar(GA2E a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	u = GA2E_OverBar(r);

	printf("GA2E_OverBar Demo\n");
	printf("r =              ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_OverBar(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =       ");  GA2E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_UnderBar(GA2E a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	u = GA2E_UnderBar(r);

	printf("GA2E_UnderBar Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_UnderBar(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");

	r = GA2E_Set(  3,    5,  7,   11) ;  
	u = GA2E_UnderBar(r);

	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf(" expect zeroes in q, x, and y components\n");

	t = GA2E_OverBar(GA2E_UnderBar(r));
	printf("t = GA2E_OverBar(GA2E_UnderBar(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Reverse(GA2E w) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Reverse(r);

	printf("GA2E_Reverse Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Reverse(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Reverse(GA2E_Reverse(r));
	printf("t = GA2E_Reverse(GA2E_Reverse(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Parity(GA2E w) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Parity(r);

	printf("GA2E_Parity Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Parity(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Parity(GA2E_Parity(r));
	printf("t = GA2E_Parity(GA2E_Parity(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Transpose(GA2E w) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Transpose(r);

	printf("GA2E_Transpose Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Transpose(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Transpose(GA2E_Transpose(r));
	printf("t = GA2E_Transpose(GA2E_Transpose(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Conjugation(GA2E w) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Conjugation(r);

	printf("GA2E_Conjugation Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Conjugation(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Conjugation(GA2E_Conjugation(r));
	printf("t = GA2E_Conjugation(GA2E_Conjugation(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Clifford_Conjugation(GA2E w) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Clifford_Conjugation(r);

	printf("GA2E_Clifford_Conjugation Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Clifford_Conjugation(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Clifford_Conjugation(GA2E_Clifford_Conjugation(r));
	printf("t = GA2E_Clifford_Conjugation(GA2E_Clifford_Conjugation(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Dual(GA2E w) ; 

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Dual(r);

	printf("GA2E_Dual Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Dual(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	t = GA2E_Dual(GA2E_Dual(r));
	printf("t = GA2E_Dual(GA2E_Dual(r)) =        ");  GA2E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Dorst_Dual(GA2E a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Dorst_Dual(r);

	printf("GA2E_Dorst_Dual Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Reverse(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Dorst_UnDual(GA2E a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  

	u = GA2E_Dorst_UnDual(r);

	printf("GA2E_Dorst_UnDual Demo\n");
	printf("r =               ");  GA2E_PrintMV(r);  printf("\n");
	printf("u = GA2E_Dorst_UnDual(r) = ");  GA2E_PrintMV(u);  printf("\n");
	det = GA2E_Determinant(u);
	if(det == det_ref) printf("GA2E_Determinant is conserved\n"); else printf("GA2E_Determinant is not conserved\n");
	s = GA2E_Product(r,u);
	printf("s = r*u =        ");  GA2E_PrintMV(s);  printf("\n");

	s = GA2E_Dorst_Dual(r);
	u = GA2E_Dorst_UnDual(s);

	printf("GA2E_Dorst_UnDual(GA2E_Dorst_Dual(r)) = "); GA2E_PrintMV(u); printf("\n");

	s = GA2E_Dorst_UnDual(r);
	u = GA2E_Dorst_Dual(s);

	printf("GA2E_Dorst_Dual(GA2E_Dorst_UnDual(r)) = "); GA2E_PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA2E GA2E_Add(GA2E u,GA2E v)

	r = GA2E_Set(  3,    5,  7,   11) ; 
	s = GA2E_Set( 29,   31, 37,   41) ;
	u = GA2E_Add(r, s);

	printf("  r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA2E_PrintMV(s);  printf("\n");
	printf("r+s = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Subtract(GA2E u, GA2E v)

	r = GA2E_Set(  3,    5,  7,   11) ; 
	s = GA2E_Set( 29,   31, 37,   41) ;
	u = GA2E_Subtract(r, s);

	printf("  r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA2E_PrintMV(s);  printf("\n");
	printf("r-s = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int GA2E_Equal(GA2E u, GA2E v) ;

	r = GA2E_Set(  3,    5,  7,   11) ; 
	s = GA2E_Set( 29,   31, 37,   41) ;
	t = GA2E_Set( 29,   31, 37,   41) ; 

	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	printf("t = ");  GA2E_PrintMV(t);  printf("\n");
	if (GA2E_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (GA2E_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int GA2E_Not_Equal(GA2E u, GA2E v) ;

	r = GA2E_Set(  3,    5,  7,   11) ; 
	s = GA2E_Set( 29,   31, 37,   41) ;
	t = GA2E_Set( 29,   31, 37,   41) ; 

	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	printf("t = ");  GA2E_PrintMV(t);  printf("\n");
	if (GA2E_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (GA2E_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Product(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Product(r,s);
	printf("u = r*s = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Divide_By_Constant(GA2E u, double a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	u = GA2E_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Wedge(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Wedge(r,s);
	printf("u = GA2E_Wedge(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_AntiWedge(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_AntiWedge(r,s);
	printf("u = GA2E_AntiWedge(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Regressive(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Regressive(r,s);
	printf("u = GA2E_Regressive(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Regressive_Via_Formula(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Regressive_Via_Formula(r,s);
	printf("u = GA2E_Regressive_Via_Formula(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Lower_Right_Via_Formula(GA2E a, GA2E b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Lower_Right_Via_Formula(r,s);
	printf("u = GA2E_Lower_Right_Via_Formula(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Expander(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Expander(r,s);
	printf("u = GA2E_Expander(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Conserver(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Conserver(r,s);
	printf("u = GA2E_Conserver(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Shrinker(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Shrinker(r,s);
	printf("u = GA2E_Shrinker(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA2E GA2E_Symmetric(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Symmetric(r,s);
	printf("u = GA2E_Symmetric(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_AntiSymmetric(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_AntiSymmetric(r,s);
	printf("u = GA2E_AntiSymmetric(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Inner(const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Inner(r,s);
	printf("u = GA2E_Inner(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Left_Contraction (const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Left_Contraction(r,s);
	printf("u = GA2E_Left_Contraction(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Right_Contraction (const GA2E &a, const GA2E &b) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	s = GA2E_Set( 29,   31, 37,   41) ;
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	printf("s = ");  GA2E_PrintMV(s);  printf("\n");
	u = GA2E_Right_Contraction(r,s);
	printf("u = GA2E_Right_Contraction(r, s) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex GA2E_Determinant(GA2E A) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	det = GA2E_Determinant(r);
	printf("u = GA2E_Determinant(r) = %10.3e \n\n",det); 


//////////////////////////////////////////////////////

// GA2E GA2E_Adjugate(GA2E V) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	u = GA2E_Adjugate(r);
	printf("u = GA2E_Adjugate(r) = ");  GA2E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA2E GA2E_Reciprocal(GA2E a) ;

	r = GA2E_Set(  3,    5,  7,   11) ;  
	printf("r = ");  GA2E_PrintMV(r);  printf("\n");
	u = GA2E_Reciprocal(r);
	printf("u = GA2E_Reciprocal(r) = ");  GA2E_PrintMV(u);  printf("\n\n");
	s = GA2E_Product(r,u);
	printf("s = r*(1/r) = ");  GA2E_PrintMV(s);  printf("\n\n");

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = GA2E_Set(  3,    5,  7,   11) ; 
	s = GA2E_Set( 29,   31, 37,   41) ;

*/





