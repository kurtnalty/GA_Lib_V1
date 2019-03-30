#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GA3E_Int.c"

/////////////////////////////////// Demonstrate GA3E Routines ////////////////////////////////


int main(void)
{
	

	GA3E r,s,t,u;

//	GA3E X, Y, Z;

//	X = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
//	Y = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	Z = Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	t = Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	long det_ref, det;
	det_ref = Determinant(r);
	printf("\nGA3E Routines Demo\n\n");

//////////////////////////////////////////////////////

// void PrintlnMV(GA3E v) 

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("PrintlnMV Demo\nr = ");  PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void PrintMV(GA3E v) 

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("PrintMV Demo\nr = ");  PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E OverBar(GA3E a) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = OverBar(r);

	printf("OverBar Demo\n");
	printf("r =              ");  PrintMV(r);  printf("\n");
	printf("u = OverBar(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =       ");  PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E UnderBar(GA3E a) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = UnderBar(r);

	printf("UnderBar Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = UnderBar(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = OverBar(UnderBar(r));
	printf("t = OverBar(UnderBar(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Reverse(GA3E w) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = Reverse(r);

	printf("Reverse Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Reverse(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = Reverse(Reverse(r));
	printf("t = Reverse(Reverse(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Involution(GA3E w) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = Involution(r);

	printf("Involution Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Involution(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = Involution(Involution(r));
	printf("t = Involution(Involution(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Transpose(GA3E w) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = Transpose(r);

	printf("Transpose Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Transpose(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = Transpose(Transpose(r));
	printf("t = Transpose(Transpose(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Conjugation(GA3E w) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = Conjugation(r);

	printf("Conjugation Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Conjugation(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = Conjugation(Conjugation(r));
	printf("t = Conjugation(Conjugation(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E CliffordConjugation(GA3E w) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = CliffordConjugation(r);

	printf("CliffordConjugation Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = CliffordConjugation(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = CliffordConjugation(CliffordConjugation(r));
	printf("t = CliffordConjugation(CliffordConjugation(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Dual(GA3E w) ; 

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = Dual(r);

	printf("Dual Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Dual(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	t = Dual(Dual(r));
	printf("t = Dual(Dual(r)) =        ");  PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E DorstDual(GA3E a) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = DorstDual(r);

	printf("DorstDual Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Reverse(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E DorstUnDual(GA3E a) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = DorstUnDual(r);

	printf("DorstUnDual Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = DorstUnDual(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n");

	s = DorstDual(r);
	u = DorstUnDual(s);

	printf("DorstUnDual(DorstDual(r)) = "); PrintMV(u); printf("\n");

	s = DorstUnDual(r);
	u = DorstDual(s);

	printf("DorstDual(DorstUnDual(r)) = "); PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA3E Add(GA3E u,GA3E v)

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = Add(r, s);

	printf("  r = ");  PrintMV(r);  printf("\n");
	printf("  s = ");  PrintMV(s);  printf("\n");
	printf("r+s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Subtract(GA3E u, GA3E v)

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = Subtract(r, s);

	printf("  r = ");  PrintMV(r);  printf("\n");
	printf("  s = ");  PrintMV(s);  printf("\n");
	printf("r-s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int Equal(GA3E u, GA3E v) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	printf("t = ");  PrintMV(t);  printf("\n");
	if (Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int Not_Equal(GA3E u, GA3E v) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	printf("t = ");  PrintMV(t);  printf("\n");
	if (Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA3E Product(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Product(r,s);
	printf("u = r*s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Divide_By_Constant(GA3E u, long a) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  PrintMV(r);  printf("\n");
	u = Divide_By_Constant(r,2);
	printf("u = r/2 = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Wedge(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Wedge(r,s);
	printf("u = Wedge(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E AntiWedge(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = AntiWedge(r,s);
	printf("u = AntiWedge(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Regressive(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Regressive(r,s);
	printf("u = Regressive(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E RegressiveViaFormula(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = RegressiveViaFormula(r,s);
	printf("u = RegressiveViaFormula(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E LowerRightViaFormula(GA3E a, GA3E b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = LowerRightViaFormula(r,s);
	printf("u = LowerRightViaFormula(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Expander(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Expander(r,s);
	printf("u = Expander(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Conserver(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Conserver(r,s);
	printf("u = Conserver(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Shrinker(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Shrinker(r,s);
	printf("u = Shrinker(r, s) = ");  PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA3E Symmetric(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Symmetric(r,s);
	printf("u = Symmetric(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E AntiSymmetric(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = AntiSymmetric(r,s);
	printf("u = AntiSymmetric(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E Inner(const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Inner(r,s);
	printf("u = Inner(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E LeftContraction (const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = LeftContraction(r,s);
	printf("u = LeftContraction(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E RightContraction (const GA3E &a, const GA3E &b) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = RightContraction(r,s);
	printf("u = RightContraction(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex Determinant(GA3E A) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  PrintMV(r);  printf("\n");
	det = Determinant(r);
	printf("u = Determinant(r) = %10ld \n\n",det); 


//////////////////////////////////////////////////////

// GA3E Adjugate(GA3E V) ;

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  PrintMV(r);  printf("\n");
	u = Adjugate(r);
	printf("u = Adjugate(r) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

	return 0;
}

/* copy and paste bin

	r = Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

*/





