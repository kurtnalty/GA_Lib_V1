#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Mink_Int.c"

/////////////////////////////////// Demonstrate Mink Routines ////////////////////////////////


int main(void)
{
	
	Mink r,s,t,u;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	long det_ref, det;
	det_ref = Determinant(r);
	printf("\nMink Routines Demo\n\n");

//////////////////////////////////////////////////////

// void PrintlnMV(Mink v) 

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("PrintlnMV Demo\nr = ");  PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void PrintMV(Mink v) 

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("PrintMV Demo\nr = ");  PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink OverBar(Mink a) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = OverBar(r);

	printf("OverBar Demo\n");
	printf("r =              ");  PrintMV(r);  printf("\n");
	printf("u = OverBar(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =       ");  PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink UnderBar(Mink a) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Reverse(Mink w) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Involution(Mink w) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Transpose(Mink w) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Conjugation(Mink w) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink CliffordConjugation(Mink w) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Dual(Mink w) ; 

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink DorstDual(Mink a) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = DorstDual(r);

	printf("DorstDual Demo\n");
	printf("r =               ");  PrintMV(r);  printf("\n");
	printf("u = Reverse(r) = ");  PrintMV(u);  printf("\n");
	det = Determinant(u);
	if(det == det_ref) printf("Determinant is conserved\n"); else printf("Determinant is not conserved\n");
	s = Product(r,u);
	printf("s = r*u =        ");  PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink DorstUnDual(Mink a) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

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

// Mink Add(Mink u,Mink v)

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = Add(r, s);

	printf("  r = ");  PrintMV(r);  printf("\n");
	printf("  s = ");  PrintMV(s);  printf("\n");
	printf("r+s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Subtract(Mink u, Mink v)

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = Subtract(r, s);

	printf("  r = ");  PrintMV(r);  printf("\n");
	printf("  s = ");  PrintMV(s);  printf("\n");
	printf("r-s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int Equal(Mink u, Mink v) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	printf("t = ");  PrintMV(t);  printf("\n");
	if (Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int Not_Equal(Mink u, Mink v) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	printf("t = ");  PrintMV(t);  printf("\n");
	if (Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// Mink Product(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Product(r,s);
	printf("u = r*s = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Divide_By_Constant(Mink u, long a) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	u = Divide_By_Constant(r,2);
	printf("u = r/2 = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Wedge(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Wedge(r,s);
	printf("u = Wedge(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink AntiWedge(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = AntiWedge(r,s);
	printf("u = AntiWedge(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Regressive(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Regressive(r,s);
	printf("u = Regressive(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink RegressiveViaFormula(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = RegressiveViaFormula(r,s);
	printf("u = RegressiveViaFormula(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink LowerRightViaFormula(Mink a, Mink b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = LowerRightViaFormula(r,s);
	printf("u = LowerRightViaFormula(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Expander(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Expander(r,s);
	printf("u = Expander(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Conserver(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Conserver(r,s);
	printf("u = Conserver(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Shrinker(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Shrinker(r,s);
	printf("u = Shrinker(r, s) = ");  PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// Mink Symmetric(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Symmetric(r,s);
	printf("u = Symmetric(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink AntiSymmetric(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = AntiSymmetric(r,s);
	printf("u = AntiSymmetric(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Inner(const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = Inner(r,s);
	printf("u = Inner(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink LeftContraction (const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = LeftContraction(r,s);
	printf("u = LeftContraction(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink RightContraction (const Mink &a, const Mink &b) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	printf("s = ");  PrintMV(s);  printf("\n");
	u = RightContraction(r,s);
	printf("u = RightContraction(r, s) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex Determinant(Mink A) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	det = Determinant(r);
	printf("u = Determinant(r) = %10ld \n\n",det); 


//////////////////////////////////////////////////////

// Mink Adjugate(Mink V) ;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  PrintMV(r);  printf("\n");
	u = Adjugate(r);
	printf("u = Adjugate(r) = ");  PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

	return 0;
}

/* copy and paste bin

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

*/





