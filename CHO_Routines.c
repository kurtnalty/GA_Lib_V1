// Routines for CHO Algebra
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 4 September 2018
// License: Freeware
// Alternative License: BSD
//
//////////////////////////////////////////////////////

typedef struct {double
	q,  xyz;	
} C;

//////////////////////////////////////////////////////

typedef struct {double
	q,	xy,	yz,	xz;	
} H;

//////////////////////////////////////////////////////

typedef struct {double
	q,	i,	j,	ij,	E,	iE,	jE,	ijE;
} O;

//////////////////////////////////////////////////////

typedef struct {double
	q,	x,	y,	xy,	
	z,	xz,	yz,	xyz;	
} CH;

//////////////////////////////////////////////////////

typedef struct {double
	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
	xyz,	xyzi,	xyzj,	xyzij,	xyzE,	xyziE,	xyzjE,	xyzijE;	
} CO;

//////////////////////////////////////////////////////


typedef struct {double
	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
	xy,	xyi,	xyj,	xyij,	xyE,	xyiE,	xyjE,	xyijE,
	yz,	yzi,	yzj,	yzij,	yzE,	yziE,	yzjE,	yzijE,
	xz,	xzi,	xzj,	xzij,	xzE,	xziE,	xzjE,	xzijE;
} HO;

//////////////////////////////////////////////////////


typedef struct {double
	q,	i,	j,	ij,	E,	iE,	jE,	ijE,
	x,	xi,	xj,	xij,	xE,	xiE,	xjE,	xijE,
	y,	yi,	yj,	yij,	yE,	yiE,	yjE,	yijE,
	xy,	xyi,	xyj,	xyij,	xyE,	xyiE,	xyjE,	xyijE,
	z,	zi,	zj,	zij,	zE,	ziE,	zjE,	zijE,
	xz,	xzi,	xzj,	xzij,	xzE,	xziE,	xzjE,	xzijE,
	yz,	yzi,	yzj,	yzij,	yzE,	yziE,	yzjE,	yzijE,
	xyz,	xyzi,	xyzj,	xyzij,	xyzE,	xyziE,	xyzjE,	xyzijE;	
} CHO;

//////////////////////////////////////////////////////
/*   // These have been padded to 6 characters 
char C_Suffix[2][8] = {
		"q     ","xyz   "
	};

//////////////////////////////////////////////////////

char H_Suffix[4][8] = {
		"q     ","xy    ","yz    ","xz    "
	};

//////////////////////////////////////////////////////

char O_Suffix[8][8] = {
		"q     ","i     ","j     ","ij    ","E     ","iE    ","jE    ","ijE   "
	};

//////////////////////////////////////////////////////

char CH_Suffix[8][8] = {
		"q     ","x     ","y     ","xy    ",
		"z     ","xz    ","yz    ","xyz   "
	};

//////////////////////////////////////////////////////

char CO_Suffix[16][8] = {
		"q     ","i     ","j     ","ij    ","E     ","iE    ","jE    ","ijE   ",
		"xyz   ","xyzi  ","xyzj  ","xyzij ","xyzE  ","xyziE ","xyzjE ","xyzijE"
	};

//////////////////////////////////////////////////////

char HO_Suffix[32][8] = {
		"q     ","i     ","j     ","ij    ","E     ","iE    ","jE    ","ijE   ",
		"xy    ","xyi   ","xyj   ","xyij  ","xyE   ","xyiE  ","xyjE  ","xyijE ",
		"yz    ","yzi   ","yzj   ","yzij  ","yzE   ","yziE  ","yzjE  ","yzijE ",
		"xz    ","xzi   ","xzj   ","xzij  ","xzE   ","xziE  ","xzjE  ","xzijE "
	};

//////////////////////////////////////////////////////

char CHO_Suffix[64][8] = {
		"q     ","i     ","j     ","ij    ","E     ","iE    ","jE    ","ijE   ",
		"x     ","xi    ","xj    ","xij   ","xE    ","xiE   ","xjE   ","xijE  ",
		"y     ","yi    ","yj    ","yij   ","yE    ","yiE   ","yjE   ","yijE  ",
		"xy    ","xyi   ","xyj   ","xyij  ","xyE   ","xyiE  ","xyjE  ","xyijE ",
		"z     ","zi    ","zj    ","zij   ","zE    ","ziE   ","zjE   ","zijE  ",
		"xz    ","xzi   ","xzj   ","xzij  ","xzE   ","xziE  ","xzjE  ","xzijE ",
		"yz    ","yzi   ","yzj   ","yzij  ","yzE   ","yziE  ","yzjE  ","yzijE ",
		"xyz   ","xyzi  ","xyzj  ","xyzij ","xyzE  ","xyziE ","xyzjE ","xyzijE"
	};
*/
//////////////////////////////////////////////////////
//  no padding - tight format
char C_Suffix[2][8] = {
		"q","xyz"
	};

//////////////////////////////////////////////////////

char H_Suffix[4][8] = {
		"q","xy","yz","xz"
	};

//////////////////////////////////////////////////////

char O_Suffix[8][8] = {
		"q","i","j","ij","E","iE","jE","ijE"
	};

//////////////////////////////////////////////////////

char CH_Suffix[8][8] = {
		"q","x","y","xy",
		"z","xz","yz","xyz"
	};

//////////////////////////////////////////////////////

char CO_Suffix[16][8] = {
		"q","i","j","ij","E","iE","jE","ijE",
		"xyz","xyzi","xyzj","xyzij","xyzE","xyziE","xyzjE","xyzijE"
	};

//////////////////////////////////////////////////////

char HO_Suffix[32][8] = {
		"q","i","j","ij","E","iE","jE","ijE",
		"xy","xyi","xyj","xyij","xyE","xyiE","xyjE","xyijE",
		"yz","yzi","yzj","yzij","yzE","yziE","yzjE","yzijE",
		"xz","xzi","xzj","xzij","xzE","xziE","xzjE","xzijE"
	};

//////////////////////////////////////////////////////

char CHO_Suffix[64][8] = {
		"q","i","j","ij","E","iE","jE","ijE",
		"x","xi","xj","xij","xE","xiE","xjE","xijE",
		"y","yi","yj","yij","yE","yiE","yjE","yijE",
		"xy","xyi","xyj","xyij","xyE","xyiE","xyjE","xyijE",
		"z","zi","zj","zij","zE","ziE","zjE","zijE",
		"xz","xzi","xzj","xzij","xzE","xziE","xzjE","xzijE",
		"yz","yzi","yzj","yzij","yzE","yziE","yzjE","yzijE",
		"xyz","xyzi","xyzj","xyzij","xyzE","xyziE","xyzjE","xyzijE"
	};

//////////////////////////////////////////////////////

void Program_Failed (const char* desc)	//	Report error and quit
{
	printf("%s", desc);
	exit(0);
}

//////////////////////////////////////////////////////

C C_Set( double a0,  double a1 )
{
	C a;

	a.q = a0 ; 	a.xyz = a1 ; 

	return a;
}

//////////////////////////////////////////////////////

H H_Set( double a0,  double a1,  double a2,  double a3 )
{
	H a;

	a.q = a0 ; 	a.xy = a1 ; 	a.yz = a2 ; 	a.xz = a3 ; 

	return a;
}

//////////////////////////////////////////////////////

O O_Set( double a0,  double a1,  double a2,  double a3,  double a4,  double a5,  double a6,  double a7 )
{
	O a;

	a.q = a0 ; 	a.i = a1 ; 	a.j = a2 ; 	a.ij = a3 ;
 	a.E = a4 ; 	a.iE = a5 ; 	a.jE = a6 ; 	a.ijE = a7 ; 

	return a;
}

//////////////////////////////////////////////////////

CH CH_Set( double a0,  double a1,  double a2,  double a3,  double b0,  double b1,  double b2,  double b3 )
{

	CH a;

	a.q = a0 ; 	a.x = a1 ; 	a.y = a2 ; 	a.xy = a3 ;
 	a.z = b0 ; 	a.xz = b1 ; 	a.yz = b2 ; 	a.xyz = b3 ; 

	return a;
}

//////////////////////////////////////////////////////

CO CO_Set(	double a0,  double a1,  double a2,  double a3,  double a4,  double a5,  double a6,  double a7, 
		double b0,  double b1,  double b2,  double b3,  double b4,  double b5,  double b6,  double b7 )
{
	CO a;
	
	a.q = a0 ; 	a.i = a1 ; 	a.j = a2 ; 	a.ij = a3 ; 	a.E = a4 ; 	a.iE = a5 ; 	a.jE = a6 ; 	a.ijE = a7 ; 
	a.xyz = b0 ; 	a.xyzi = b1 ; 	a.xyzj = b2 ; 	a.xyzij = b3 ; 	a.xyzE = b4 ; 	a.xyziE = b5 ; 	a.xyzjE = b6 ; 	a.xyzijE = b7 ; 

	return a;
}

//////////////////////////////////////////////////////

HO HO_Set( 	double a0,  double a1,  double a2,  double a3,  double a4,  double a5,  double a6,  double a7, 
		double b0,  double b1,  double b2,  double b3,  double b4,  double b5,  double b6,  double b7, 
		double c0,  double c1,  double c2,  double c3,  double c4,  double c5,  double c6,  double c7, 
		double d0,  double d1,  double d2,  double d3,  double d4,  double d5,  double d6,  double d7 )
{
	HO a;

	a.q = a0 ; 	a.i = a1 ; 	a.j = a2 ; 	a.ij = a3 ; 	a.E = a4 ; 	a.iE = a5 ; 	a.jE = a6 ; 	a.ijE = a7 ; 
	a.xy = b0 ; 	a.xyi = b1 ; 	a.xyj = b2 ; 	a.xyij = b3 ; 	a.xyE = b4 ; 	a.xyiE = b5 ; 	a.xyjE = b6 ; 	a.xyijE = b7 ; 
	a.yz = c0 ; 	a.yzi = c1 ; 	a.yzj = c2 ; 	a.yzij = c3 ; 	a.yzE = c4 ; 	a.yziE = c5 ; 	a.yzjE = c6 ; 	a.yzijE = c7 ; 
	a.xz = d0 ; 	a.xzi = d1 ; 	a.xzj = d2 ; 	a.xzij = d3 ; 	a.xzE = d4 ; 	a.xziE = d5 ; 	a.xzjE = d6 ; 	a.xzijE = d7 ; 

	return a;
}

//////////////////////////////////////////////////////

CHO CHO_Set( 	double a0,  double a1,  double a2,  double a3,  double a4,  double a5,  double a6,  double a7, 
		double b0,  double b1,  double b2,  double b3,  double b4,  double b5,  double b6,  double b7, 
		double c0,  double c1,  double c2,  double c3,  double c4,  double c5,  double c6,  double c7, 
		double d0,  double d1,  double d2,  double d3,  double d4,  double d5,  double d6,  double d7, 
		double e0,  double e1,  double e2,  double e3,  double e4,  double e5,  double e6,  double e7, 
		double f0,  double f1,  double f2,  double f3,  double f4,  double f5,  double f6,  double f7, 
		double g0,  double g1,  double g2,  double g3,  double g4,  double g5,  double g6,  double g7, 
		double h0,  double h1,  double h2,  double h3,  double h4,  double h5,  double h6,  double h7 )
{
	CHO a;

	a.q = a0 ; 	a.i = a1 ; 	a.j = a2 ; 	a.ij = a3 ; 	a.E = a4 ; 	a.iE = a5 ; 	a.jE = a6 ; 	a.ijE = a7 ; 
	a.x = b0 ; 	a.xi = b1 ; 	a.xj = b2 ; 	a.xij = b3 ; 	a.xE = b4 ; 	a.xiE = b5 ; 	a.xjE = b6 ; 	a.xijE = b7 ; 
	a.y = c0 ; 	a.yi = c1 ; 	a.yj = c2 ; 	a.yij = c3 ; 	a.yE = c4 ; 	a.yiE = c5 ; 	a.yjE = c6 ; 	a.yijE = c7 ; 
	a.xy = d0 ; 	a.xyi = d1 ; 	a.xyj = d2 ; 	a.xyij = d3 ; 	a.xyE = d4 ; 	a.xyiE = d5 ; 	a.xyjE = d6 ; 	a.xyijE = d7 ; 
	a.z = e0 ; 	a.zi = e1 ; 	a.zj = e2 ; 	a.zij = e3 ; 	a.zE = e4 ; 	a.ziE = e5 ; 	a.zjE = e6 ; 	a.zijE = e7 ; 
	a.xz = f0 ; 	a.xzi = f1 ; 	a.xzj = f2 ; 	a.xzij = f3 ; 	a.xzE = f4 ; 	a.xziE = f5 ; 	a.xzjE = f6 ; 	a.xzijE = f7 ; 
	a.yz = g0 ; 	a.yzi = g1 ; 	a.yzj = g2 ; 	a.yzij = g3 ; 	a.yzE = g4 ; 	a.yziE = g5 ; 	a.yzjE = g6 ; 	a.yzijE = g7 ; 
	a.xyz = h0 ; 	a.xyzi = h1 ; 	a.xyzj = h2 ; 	a.xyzij = h3 ; 	a.xyzE = h4 ; 	a.xyziE = h5 ; 	a.xyzjE = h6 ; 	a.xyzijE = h7 ; 

	return a;

}

//////////////////////////////////////////////////////

C C_Zero(void)
{
	C a;

	a.q = 0.0;	a.xyz = 0.0;
	return a;
}

//////////////////////////////////////////////////////

H H_Zero(void)
{
	H a;

	a.q = 0.0;	a.xy = 0.0;	a.yz = 0.0;	a.xz = 0.0;
	return a;
}

//////////////////////////////////////////////////////

O O_Zero(void)
{
	O a;

	a.q = 0.0;	a.i = 0.0;	a.j = 0.0;	a.ij = 0.0;	a.E = 0.0;	a.iE = 0.0;	a.jE = 0.0;	a.ijE = 0.0;
	return a;
}

//////////////////////////////////////////////////////

CH CH_Zero(void)
{
	CH a;

	a.q = 0.0;	a.x = 0.0;	a.y = 0.0;	a.xy = 0.0;
	a.z = 0.0;	a.xz = 0.0;	a.yz = 0.0;	a.xyz = 0.0;
	return a;
}

//////////////////////////////////////////////////////

CO CO_Zero(void)
{
	CO a;

	a.q = 0.0;	a.i = 0.0;	a.j = 0.0;	a.ij = 0.0;	a.E = 0.0;	a.iE = 0.0;	a.jE = 0.0;	a.ijE = 0.0;
	a.xyz = 0.0;	a.xyzi = 0.0;	a.xyzj = 0.0;	a.xyzij = 0.0;	a.xyzE = 0.0;	a.xyziE = 0.0;	a.xyzjE = 0.0;	a.xyzijE = 0.0;

	return a;
}

//////////////////////////////////////////////////////

HO HO_Zero(void)
{
	HO a;

	a.q = 0.0;	a.i = 0.0;	a.j = 0.0;	a.ij = 0.0;	a.E = 0.0;	a.iE = 0.0;	a.jE = 0.0;	a.ijE = 0.0;
	a.xy = 0.0;	a.xyi = 0.0;	a.xyj = 0.0;	a.xyij = 0.0;	a.xyE = 0.0;	a.xyiE = 0.0;	a.xyjE = 0.0;	a.xyijE = 0.0;
	a.yz = 0.0;	a.yzi = 0.0;	a.yzj = 0.0;	a.yzij = 0.0;	a.yzE = 0.0;	a.yziE = 0.0;	a.yzjE = 0.0;	a.yzijE = 0.0;
	a.xz = 0.0;	a.xzi = 0.0;	a.xzj = 0.0;	a.xzij = 0.0;	a.xzE = 0.0;	a.xziE = 0.0;	a.xzjE = 0.0;	a.xzijE = 0.0;

	return a;
}

//////////////////////////////////////////////////////

CHO CHO_Zero(void)
{
	CHO a;

	a.q = 0.0;	a.i = 0.0;	a.j = 0.0;	a.ij = 0.0;	a.E = 0.0;	a.iE = 0.0;	a.jE = 0.0;	a.ijE = 0.0;
	a.x = 0.0;	a.xi = 0.0;	a.xj = 0.0;	a.xij = 0.0;	a.xE = 0.0;	a.xiE = 0.0;	a.xjE = 0.0;	a.xijE = 0.0;
	a.y = 0.0;	a.yi = 0.0;	a.yj = 0.0;	a.yij = 0.0;	a.yE = 0.0;	a.yiE = 0.0;	a.yjE = 0.0;	a.yijE = 0.0;
	a.xy = 0.0;	a.xyi = 0.0;	a.xyj = 0.0;	a.xyij = 0.0;	a.xyE = 0.0;	a.xyiE = 0.0;	a.xyjE = 0.0;	a.xyijE = 0.0;
	a.z = 0.0;	a.zi = 0.0;	a.zj = 0.0;	a.zij = 0.0;	a.zE = 0.0;	a.ziE = 0.0;	a.zjE = 0.0;	a.zijE = 0.0;
	a.xz = 0.0;	a.xzi = 0.0;	a.xzj = 0.0;	a.xzij = 0.0;	a.xzE = 0.0;	a.xziE = 0.0;	a.xzjE = 0.0;	a.xzijE = 0.0;
	a.yz = 0.0;	a.yzi = 0.0;	a.yzj = 0.0;	a.yzij = 0.0;	a.yzE = 0.0;	a.yziE = 0.0;	a.yzjE = 0.0;	a.yzijE = 0.0;
	a.xyz = 0.0;	a.xyzi = 0.0;	a.xyzj = 0.0;	a.xyzij = 0.0;	a.xyzE = 0.0;	a.xyziE = 0.0;	a.xyzjE = 0.0;	a.xyzijE = 0.0;

	return a;
}

//////////////////////////////////////////////////////

int C_EQ_Zero(C a) {
	int EQZ;
	EQZ = (
		(a.q    == 0.0) && (a.xyz  == 0.0)
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int H_EQ_Zero(H a) {
	int EQZ;
	EQZ = (
		(a.q == 0.0) && (a.xy == 0.0) && (a.yz == 0.0) && (a.xz == 0.0)
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int O_EQ_Zero(O a) {
	int EQZ;
	EQZ = (
		(a.q == 0.0) && (a.i  == 0.0) && (a.j  == 0.0) && (a.ij  == 0.0) &&
		(a.E == 0.0) && (a.iE == 0.0) && (a.jE == 0.0) && (a.ijE == 0.0)
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int CH_EQ_Zero(CH a) {
	int EQZ;
	EQZ = (
		(a.q == 0.0) && (a.x  == 0.0) && (a.y  == 0.0) && (a.xy  == 0.0) &&
		(a.z == 0.0) && (a.xz == 0.0) && (a.yz == 0.0) && (a.xyz == 0.0)
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int CO_EQ_Zero(CO a) {
	int EQZ;
	EQZ = (
		(a.q    == 0.0) && (a.i     == 0.0) && (a.j     == 0.0) && (a.ij     == 0.0) && 
		(a.E    == 0.0) && (a.iE    == 0.0) && (a.jE    == 0.0) && (a.ijE    == 0.0) && 
		(a.xyz  == 0.0) && (a.xyzi  == 0.0) && (a.xyzj  == 0.0) && (a.xyzij  == 0.0) && 
		(a.xyzE == 0.0) && (a.xyziE == 0.0) && (a.xyzjE == 0.0) && (a.xyzijE == 0.0) 
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int HO_EQ_Zero(HO a) {
	int EQZ;
	EQZ = (
		(a.q    == 0.0) && (a.i     == 0.0) && (a.j     == 0.0) && (a.ij     == 0.0) && 
		(a.E    == 0.0) && (a.iE    == 0.0) && (a.jE    == 0.0) && (a.ijE    == 0.0) && 
		(a.xy   == 0.0) && (a.xyi   == 0.0) && (a.xyj   == 0.0) && (a.xyij   == 0.0) && 
		(a.xyE  == 0.0) && (a.xyiE  == 0.0) && (a.xyjE  == 0.0) && (a.xyijE  == 0.0) && 
		(a.yz   == 0.0) && (a.yzi   == 0.0) && (a.yzj   == 0.0) && (a.yzij   == 0.0) && 
		(a.yzE  == 0.0) && (a.yziE  == 0.0) && (a.yzjE  == 0.0) && (a.yzijE  == 0.0) &&
		(a.xz   == 0.0) && (a.xzi   == 0.0) && (a.xzj   == 0.0) && (a.xzij   == 0.0) && 
		(a.xzE  == 0.0) && (a.xziE  == 0.0) && (a.xzjE  == 0.0) && (a.xzijE  == 0.0)  
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int CHO_EQ_Zero(CHO a) {
	int EQZ;
	EQZ = (
		(a.q    == 0.0) && (a.i     == 0.0) && (a.j     == 0.0) && (a.ij     == 0.0) && 
		(a.E    == 0.0) && (a.iE    == 0.0) && (a.jE    == 0.0) && (a.ijE    == 0.0) && 
		(a.x    == 0.0) && (a.xi    == 0.0) && (a.xj    == 0.0) && (a.xij    == 0.0) && 
		(a.xE   == 0.0) && (a.xiE   == 0.0) && (a.xjE   == 0.0) && (a.xijE   == 0.0) && 
		(a.y    == 0.0) && (a.yi    == 0.0) && (a.yj    == 0.0) && (a.yij    == 0.0) && 
		(a.yE   == 0.0) && (a.yiE   == 0.0) && (a.yjE   == 0.0) && (a.yijE   == 0.0) && 
		(a.xy   == 0.0) && (a.xyi   == 0.0) && (a.xyj   == 0.0) && (a.xyij   == 0.0) && 
		(a.xyE  == 0.0) && (a.xyiE  == 0.0) && (a.xyjE  == 0.0) && (a.xyijE  == 0.0) && 
		(a.z    == 0.0) && (a.zi    == 0.0) && (a.zj    == 0.0) && (a.zij    == 0.0) && 
		(a.zE   == 0.0) && (a.ziE   == 0.0) && (a.zjE   == 0.0) && (a.zijE   == 0.0) && 
		(a.xz   == 0.0) && (a.xzi   == 0.0) && (a.xzj   == 0.0) && (a.xzij   == 0.0) && 
		(a.xzE  == 0.0) && (a.xziE  == 0.0) && (a.xzjE  == 0.0) && (a.xzijE  == 0.0) && 
		(a.yz   == 0.0) && (a.yzi   == 0.0) && (a.yzj   == 0.0) && (a.yzij   == 0.0) && 
		(a.yzE  == 0.0) && (a.yziE  == 0.0) && (a.yzjE  == 0.0) && (a.yzijE  == 0.0) && 
		(a.xyz  == 0.0) && (a.xyzi  == 0.0) && (a.xyzj  == 0.0) && (a.xyzij  == 0.0) && 
		(a.xyzE == 0.0) && (a.xyziE == 0.0) && (a.xyzjE == 0.0) && (a.xyzijE == 0.0) 
	);
	return EQZ;
}

//////////////////////////////////////////////////////

int C_NE_Zero(C a) {
	int NEZ;
	NEZ = (
		(a.q    != 0.0) || (a.xyz     != 0.0)  
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int H_NE_Zero(H a) {
	int NEZ;
	NEZ = (
		(a.q != 0.0) || (a.xy != 0.0) || (a.yz != 0.0) || (a.xz != 0.0) 
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int O_NE_Zero(O a) {
	int NEZ;
	NEZ = (
		(a.q != 0.0) || (a.i  != 0.0) || (a.j  != 0.0) || (a.ij  != 0.0) ||
		(a.E != 0.0) || (a.iE != 0.0) || (a.jE != 0.0) || (a.ijE != 0.0) 
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int CH_NE_Zero(CH a) {
	int NEZ;
	NEZ = (
		(a.q != 0.0) || (a.x  != 0.0) || (a.y  != 0.0) || (a.xy  != 0.0) ||
		(a.z != 0.0) || (a.xz != 0.0) || (a.yz != 0.0) || (a.xyz != 0.0) 
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int CO_NE_Zero(CO a) {
	int NEZ;
	NEZ = (
		(a.q    != 0.0) || (a.i     != 0.0) || (a.j     != 0.0) || (a.ij     != 0.0) || 
		(a.E    != 0.0) || (a.iE    != 0.0) || (a.jE    != 0.0) || (a.ijE    != 0.0) || 
		(a.xyz  != 0.0) || (a.xyzi  != 0.0) || (a.xyzj  != 0.0) || (a.xyzij  != 0.0) || 
		(a.xyzE != 0.0) || (a.xyziE != 0.0) || (a.xyzjE != 0.0) || (a.xyzijE != 0.0) 
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int HO_NE_Zero(HO a) {
	int NEZ;
	NEZ = (
		(a.q    != 0.0) || (a.i     != 0.0) || (a.j     != 0.0) || (a.ij     != 0.0) || 
		(a.E    != 0.0) || (a.iE    != 0.0) || (a.jE    != 0.0) || (a.ijE    != 0.0) || 
		(a.xy   != 0.0) || (a.xyi   != 0.0) || (a.xyj   != 0.0) || (a.xyij   != 0.0) || 
		(a.xyE  != 0.0) || (a.xyiE  != 0.0) || (a.xyjE  != 0.0) || (a.xyijE  != 0.0) || 
		(a.yz   != 0.0) || (a.yzi   != 0.0) || (a.yzj   != 0.0) || (a.yzij   != 0.0) || 
		(a.yzE  != 0.0) || (a.yziE  != 0.0) || (a.yzjE  != 0.0) || (a.yzijE  != 0.0) ||
		(a.xz   != 0.0) || (a.xzi   != 0.0) || (a.xzj   != 0.0) || (a.xzij   != 0.0) || 
		(a.xzE  != 0.0) || (a.xziE  != 0.0) || (a.xzjE  != 0.0) || (a.xzijE  != 0.0)  
	);
	return NEZ;
}

//////////////////////////////////////////////////////

int CHO_NE_Zero(CHO a) {
	int NEZ;
	NEZ = (
		(a.q    != 0.0) || (a.i     != 0.0) || (a.j     != 0.0) || (a.ij     != 0.0) || 
		(a.E    != 0.0) || (a.iE    != 0.0) || (a.jE    != 0.0) || (a.ijE    != 0.0) || 
		(a.x    != 0.0) || (a.xi    != 0.0) || (a.xj    != 0.0) || (a.xij    != 0.0) || 
		(a.xE   != 0.0) || (a.xiE   != 0.0) || (a.xjE   != 0.0) || (a.xijE   != 0.0) || 
		(a.y    != 0.0) || (a.yi    != 0.0) || (a.yj    != 0.0) || (a.yij    != 0.0) || 
		(a.yE   != 0.0) || (a.yiE   != 0.0) || (a.yjE   != 0.0) || (a.yijE   != 0.0) || 
		(a.xy   != 0.0) || (a.xyi   != 0.0) || (a.xyj   != 0.0) || (a.xyij   != 0.0) || 
		(a.xyE  != 0.0) || (a.xyiE  != 0.0) || (a.xyjE  != 0.0) || (a.xyijE  != 0.0) || 
		(a.z    != 0.0) || (a.zi    != 0.0) || (a.zj    != 0.0) || (a.zij    != 0.0) || 
		(a.zE   != 0.0) || (a.ziE   != 0.0) || (a.zjE   != 0.0) || (a.zijE   != 0.0) || 
		(a.xz   != 0.0) || (a.xzi   != 0.0) || (a.xzj   != 0.0) || (a.xzij   != 0.0) || 
		(a.xzE  != 0.0) || (a.xziE  != 0.0) || (a.xzjE  != 0.0) || (a.xzijE  != 0.0) || 
		(a.yz   != 0.0) || (a.yzi   != 0.0) || (a.yzj   != 0.0) || (a.yzij   != 0.0) || 
		(a.yzE  != 0.0) || (a.yziE  != 0.0) || (a.yzjE  != 0.0) || (a.yzijE  != 0.0) || 
		(a.xyz  != 0.0) || (a.xyzi  != 0.0) || (a.xyzj  != 0.0) || (a.xyzij  != 0.0) || 
		(a.xyzE != 0.0) || (a.xyziE != 0.0) || (a.xyzjE != 0.0) || (a.xyzijE != 0.0) 
	);
	return NEZ;
}

//////////////////////////////////////////////////////

C C_Add(C a, C b) {

	C c;

	c.q   = a.q   + b.q;		
	c.xyz = a.xyz + b.xyz;
	return c;
}

//////////////////////////////////////////////////////

H H_Add(H a, H b) {

	H c;

	c.q = a.q + b.q;	c.xy = a.xy + b.xy;	c.yz = a.yz + b.yz;	c.xz = a.xz + b.xz;

	return c;
}

//////////////////////////////////////////////////////

O O_Add(O a, O b) {

	O c;

	c.q = a.q + b.q;	c.i  = a.i  + b.i;	c.j  = a.j  + b.j;	c.ij  = a.ij  + b.ij ;
	c.E = a.E + b.E;	c.iE = a.iE + b.iE;	c.jE = a.jE + b.jE;	c.ijE = a.ijE + b.ijE;

	return c;
}

//////////////////////////////////////////////////////

CH CH_Add(CH a, CH b) {

	CH c;

	c.q = a.q + b.q;	c.x  = a.x  + b.x;	c.y  = a.y  + b.y;	c.xy  = a.xy  + b.xy ;
	c.z = a.z + b.z;	c.xz = a.xz + b.xz;	c.yz = a.yz + b.yz;	c.xyz = a.xyz + b.xyz;

	return c;
}

//////////////////////////////////////////////////////

CO CO_Add(CO a, CO b) {

	CO c;

	c.q    = a.q    + b.q;		c.i     = a.i     + b.i;	c.j     = a.j     + b.j;	c.ij     = a.ij     + b.ij;
	c.E    = a.E    + b.E;		c.iE    = a.iE    + b.iE;	c.jE    = a.jE    + b.jE;	c.ijE    = a.ijE    + b.ijE;
	c.xyz  = a.xyz  + b.xyz;	c.xyzi  = a.xyzi  + b.xyzi;	c.xyzj  = a.xyzj  + b.xyzj;	c.xyzij  = a.xyzij  + b.xyzij;
	c.xyzE = a.xyzE + b.xyzE;	c.xyziE = a.xyziE + b.xyziE;	c.xyzjE = a.xyzjE + b.xyzjE;	c.xyzijE = a.xyzijE + b.xyzijE;

	return c;
}

//////////////////////////////////////////////////////

HO HO_Add(HO a, HO b) {

	HO c;

	c.q    = a.q    + b.q;		c.i     = a.i     + b.i;	c.j     = a.j     + b.j;	c.ij     = a.ij     + b.ij;
	c.E    = a.E    + b.E;		c.iE    = a.iE    + b.iE;	c.jE    = a.jE    + b.jE;	c.ijE    = a.ijE    + b.ijE;
	c.xy   = a.xy   + b.xy;		c.xyi   = a.xyi   + b.xyi;	c.xyj   = a.xyj   + b.xyj;	c.xyij   = a.xyij   + b.xyij;
	c.xyE  = a.xyE  + b.xyE;	c.xyiE  = a.xyiE  + b.xyiE;	c.xyjE  = a.xyjE  + b.xyjE;	c.xyijE  = a.xyijE  + b.xyijE;
	c.yz   = a.yz   + b.yz;		c.yzi   = a.yzi   + b.yzi;	c.yzj   = a.yzj   + b.yzj;	c.yzij   = a.yzij   + b.yzij;
	c.yzE  = a.yzE  + b.yzE;	c.yziE  = a.yziE  + b.yziE;	c.yzjE  = a.yzjE  + b.yzjE;	c.yzijE  = a.yzijE  + b.yzijE;
	c.xz   = a.xz   + b.xz;		c.xzi   = a.xzi   + b.xzi;	c.xzj   = a.xzj   + b.xzj;	c.xzij   = a.xzij   + b.xzij;
	c.xzE  = a.xzE  + b.xzE;	c.xziE  = a.xziE  + b.xziE;	c.xzjE  = a.xzjE  + b.xzjE;	c.xzijE  = a.xzijE  + b.xzijE;

	return c;
}

//////////////////////////////////////////////////////

CHO CHO_Add(CHO a, CHO b) {

	CHO c;

	c.q    = a.q    + b.q;		c.i     = a.i     + b.i;	c.j     = a.j     + b.j;	c.ij     = a.ij     + b.ij;
	c.E    = a.E    + b.E;		c.iE    = a.iE    + b.iE;	c.jE    = a.jE    + b.jE;	c.ijE    = a.ijE    + b.ijE;
	c.x    = a.x    + b.x;		c.xi    = a.xi    + b.xi;	c.xj    = a.xj    + b.xj;	c.xij    = a.xij    + b.xij;
	c.xE   = a.xE   + b.xE;		c.xiE   = a.xiE   + b.xiE;	c.xjE   = a.xjE   + b.xjE;	c.xijE   = a.xijE   + b.xijE;
	c.y    = a.y    + b.y;		c.yi    = a.yi    + b.yi;	c.yj    = a.yj    + b.yj;	c.yij    = a.yij    + b.yij;
	c.yE   = a.yE   + b.yE;		c.yiE   = a.yiE   + b.yiE;	c.yjE   = a.yjE   + b.yjE;	c.yijE   = a.yijE   + b.yijE;
	c.xy   = a.xy   + b.xy;		c.xyi   = a.xyi   + b.xyi;	c.xyj   = a.xyj   + b.xyj;	c.xyij   = a.xyij   + b.xyij;
	c.xyE  = a.xyE  + b.xyE;	c.xyiE  = a.xyiE  + b.xyiE;	c.xyjE  = a.xyjE  + b.xyjE;	c.xyijE  = a.xyijE  + b.xyijE;
	c.z    = a.z    + b.z;		c.zi    = a.zi    + b.zi;	c.zj    = a.zj    + b.zj;	c.zij    = a.zij    + b.zij;
	c.zE   = a.zE   + b.zE;		c.ziE   = a.ziE   + b.ziE;	c.zjE   = a.zjE   + b.zjE;	c.zijE   = a.zijE   + b.zijE;
	c.xz   = a.xz   + b.xz;		c.xzi   = a.xzi   + b.xzi;	c.xzj   = a.xzj   + b.xzj;	c.xzij   = a.xzij   + b.xzij;
	c.xzE  = a.xzE  + b.xzE;	c.xziE  = a.xziE  + b.xziE;	c.xzjE  = a.xzjE  + b.xzjE;	c.xzijE  = a.xzijE  + b.xzijE;
	c.yz   = a.yz   + b.yz;		c.yzi   = a.yzi   + b.yzi;	c.yzj   = a.yzj   + b.yzj;	c.yzij   = a.yzij   + b.yzij;
	c.yzE  = a.yzE  + b.yzE;	c.yziE  = a.yziE  + b.yziE;	c.yzjE  = a.yzjE  + b.yzjE;	c.yzijE  = a.yzijE  + b.yzijE;
	c.xyz  = a.xyz  + b.xyz;	c.xyzi  = a.xyzi  + b.xyzi;	c.xyzj  = a.xyzj  + b.xyzj;	c.xyzij  = a.xyzij  + b.xyzij;
	c.xyzE = a.xyzE + b.xyzE;	c.xyziE = a.xyziE + b.xyziE;	c.xyzjE = a.xyzjE + b.xyzjE;	c.xyzijE = a.xyzijE + b.xyzijE;

	return c;
}

//////////////////////////////////////////////////////

C C_Subtract(C a, C b) {

	C c;

	c.q   = a.q   - b.q;		
	c.xyz = a.xyz - b.xyz;
	return c;
}

//////////////////////////////////////////////////////

H H_Subtract(H a, H b) {

	H c;

	c.q = a.q - b.q;	c.xy = a.xy - b.xy;	c.yz = a.yz - b.yz;	c.xz = a.xz - b.xz;

	return c;
}

//////////////////////////////////////////////////////

O O_Subtract(O a, O b) {

	O c;

	c.q = a.q - b.q;	c.i  = a.i  - b.i;	c.j  = a.j  - b.j;	c.ij  = a.ij  - b.ij ;
	c.E = a.E - b.E;	c.iE = a.iE - b.iE;	c.jE = a.jE - b.jE;	c.ijE = a.ijE - b.ijE;

	return c;
}

//////////////////////////////////////////////////////

CH CH_Subtract(CH a, CH b) {

	CH c;

	c.q = a.q - b.q;	c.x  = a.x  - b.x;	c.y  = a.y  - b.y;	c.xy  = a.xy  - b.xy ;
	c.z = a.z - b.z;	c.xz = a.xz - b.xz;	c.yz = a.yz - b.yz;	c.xyz = a.xyz - b.xyz;

	return c;
}

//////////////////////////////////////////////////////

CO CO_Subtract(CO a, CO b) {

	CO c;

	c.q    = a.q    - b.q;		c.i     = a.i     - b.i;	c.j     = a.j     - b.j;	c.ij     = a.ij     - b.ij;
	c.E    = a.E    - b.E;		c.iE    = a.iE    - b.iE;	c.jE    = a.jE    - b.jE;	c.ijE    = a.ijE    - b.ijE;
	c.xyz  = a.xyz  - b.xyz;	c.xyzi  = a.xyzi  - b.xyzi;	c.xyzj  = a.xyzj  - b.xyzj;	c.xyzij  = a.xyzij  - b.xyzij;
	c.xyzE = a.xyzE - b.xyzE;	c.xyziE = a.xyziE - b.xyziE;	c.xyzjE = a.xyzjE - b.xyzjE;	c.xyzijE = a.xyzijE - b.xyzijE;

	return c;
}

//////////////////////////////////////////////////////

HO HO_Subtract(HO a, HO b) {

	HO c;

	c.q    = a.q    - b.q;		c.i     = a.i     - b.i;	c.j     = a.j     - b.j;	c.ij     = a.ij     - b.ij;
	c.E    = a.E    - b.E;		c.iE    = a.iE    - b.iE;	c.jE    = a.jE    - b.jE;	c.ijE    = a.ijE    - b.ijE;
	c.xy   = a.xy   - b.xy;		c.xyi   = a.xyi   - b.xyi;	c.xyj   = a.xyj   - b.xyj;	c.xyij   = a.xyij   - b.xyij;
	c.xyE  = a.xyE  - b.xyE;	c.xyiE  = a.xyiE  - b.xyiE;	c.xyjE  = a.xyjE  - b.xyjE;	c.xyijE  = a.xyijE  - b.xyijE;
	c.yz   = a.yz   - b.yz;		c.yzi   = a.yzi   - b.yzi;	c.yzj   = a.yzj   - b.yzj;	c.yzij   = a.yzij   - b.yzij;
	c.yzE  = a.yzE  - b.yzE;	c.yziE  = a.yziE  - b.yziE;	c.yzjE  = a.yzjE  - b.yzjE;	c.yzijE  = a.yzijE  - b.yzijE;
	c.xz   = a.xz   - b.xz;		c.xzi   = a.xzi   - b.xzi;	c.xzj   = a.xzj   - b.xzj;	c.xzij   = a.xzij   - b.xzij;
	c.xzE  = a.xzE  - b.xzE;	c.xziE  = a.xziE  - b.xziE;	c.xzjE  = a.xzjE  - b.xzjE;	c.xzijE  = a.xzijE  - b.xzijE;

	return c;
}

//////////////////////////////////////////////////////

CHO CHO_Subtract(CHO a, CHO b) {

	CHO c;

	c.q    = a.q    - b.q;		c.i     = a.i     - b.i;	c.j     = a.j     - b.j;	c.ij     = a.ij     - b.ij;
	c.E    = a.E    - b.E;		c.iE    = a.iE    - b.iE;	c.jE    = a.jE    - b.jE;	c.ijE    = a.ijE    - b.ijE;
	c.x    = a.x    - b.x;		c.xi    = a.xi    - b.xi;	c.xj    = a.xj    - b.xj;	c.xij    = a.xij    - b.xij;
	c.xE   = a.xE   - b.xE;		c.xiE   = a.xiE   - b.xiE;	c.xjE   = a.xjE   - b.xjE;	c.xijE   = a.xijE   - b.xijE;
	c.y    = a.y    - b.y;		c.yi    = a.yi    - b.yi;	c.yj    = a.yj    - b.yj;	c.yij    = a.yij    - b.yij;
	c.yE   = a.yE   - b.yE;		c.yiE   = a.yiE   - b.yiE;	c.yjE   = a.yjE   - b.yjE;	c.yijE   = a.yijE   - b.yijE;
	c.xy   = a.xy   - b.xy;		c.xyi   = a.xyi   - b.xyi;	c.xyj   = a.xyj   - b.xyj;	c.xyij   = a.xyij   - b.xyij;
	c.xyE  = a.xyE  - b.xyE;	c.xyiE  = a.xyiE  - b.xyiE;	c.xyjE  = a.xyjE  - b.xyjE;	c.xyijE  = a.xyijE  - b.xyijE;
	c.z    = a.z    - b.z;		c.zi    = a.zi    - b.zi;	c.zj    = a.zj    - b.zj;	c.zij    = a.zij    - b.zij;
	c.zE   = a.zE   - b.zE;		c.ziE   = a.ziE   - b.ziE;	c.zjE   = a.zjE   - b.zjE;	c.zijE   = a.zijE   - b.zijE;
	c.xz   = a.xz   - b.xz;		c.xzi   = a.xzi   - b.xzi;	c.xzj   = a.xzj   - b.xzj;	c.xzij   = a.xzij   - b.xzij;
	c.xzE  = a.xzE  - b.xzE;	c.xziE  = a.xziE  - b.xziE;	c.xzjE  = a.xzjE  - b.xzjE;	c.xzijE  = a.xzijE  - b.xzijE;
	c.yz   = a.yz   - b.yz;		c.yzi   = a.yzi   - b.yzi;	c.yzj   = a.yzj   - b.yzj;	c.yzij   = a.yzij   - b.yzij;
	c.yzE  = a.yzE  - b.yzE;	c.yziE  = a.yziE  - b.yziE;	c.yzjE  = a.yzjE  - b.yzjE;	c.yzijE  = a.yzijE  - b.yzijE;
	c.xyz  = a.xyz  - b.xyz;	c.xyzi  = a.xyzi  - b.xyzi;	c.xyzj  = a.xyzj  - b.xyzj;	c.xyzij  = a.xyzij  - b.xyzij;
	c.xyzE = a.xyzE - b.xyzE;	c.xyziE = a.xyziE - b.xyziE;	c.xyzjE = a.xyzjE - b.xyzjE;	c.xyzijE = a.xyzijE - b.xyzijE;

	return c;
}

//////////////////////////////////////////////////////

int C_EQ(C a, C b) {
	int CEQ = (
		(a.q    == b.q)    && (a.xyz  == b.xyz)
	);
	return CEQ;
}

//////////////////////////////////////////////////////

int H_EQ(H a, H b) {
	int CEQ = (
		(a.q == b.q) && (a.xy == b.xy) && (a.yz == b.yz) && (a.xz == b.xz) );
	return CEQ;
}

//////////////////////////////////////////////////////

int O_EQ(O a, O b) {
	int CEQ = (
		(a.q == b.q) && (a.i  == b.i ) && (a.j  == b.j ) && (a.ij  == b.ij ) &&
		(a.E == b.E) && (a.iE == b.iE) && (a.jE == b.jE) && (a.ijE == b.ijE) );
	return CEQ;
}

//////////////////////////////////////////////////////

int CH_EQ(CH a, CH b) {
	int CEQ = (
		(a.q == b.q) && (a.x  == b.x ) && (a.y  == b.y ) && (a.xy  == b.xy ) &&
		(a.z == b.z) && (a.xz == b.xz) && (a.yz == b.yz) && (a.xyz == b.xyz) );
	return CEQ;
}

//////////////////////////////////////////////////////

int CO_EQ(CO a, CO b) {
	int CEQ = (
		(a.q    == b.q)    && (a.i     == b.i)     && (a.j     == b.j)     && (a.ij     == b.ij)    && 
		(a.E    == b.E)    && (a.iE    == b.iE)    && (a.jE    == b.jE)    && (a.ijE    == b.ijE)   && 
		(a.xyz  == b.xyz)  && (a.xyzi  == b.xyzi)  && (a.xyzj  == b.xyzj)  && (a.xyzij  == b.xyzij) && 
		(a.xyzE == b.xyzE) && (a.xyziE == b.xyziE) && (a.xyzjE == b.xyzjE) && (a.xyzijE == b.xyzijE) );
	return CEQ;
}

//////////////////////////////////////////////////////

int HO_EQ(HO a, HO b) {
	int CEQ = (
		(a.q    == b.q)    && (a.i     == b.i)     && (a.j     == b.j)     && (a.ij     == b.ij)    && 
		(a.E    == b.E)    && (a.iE    == b.iE)    && (a.jE    == b.jE)    && (a.ijE    == b.ijE)   && 
		(a.xy   == b.xy)   && (a.xyi   == b.xyi)   && (a.xyj   == b.xyj)   && (a.xyij   == b.xyij)  && 
		(a.xyE  == b.xyE)  && (a.xyiE  == b.xyiE)  && (a.xyjE  == b.xyjE)  && (a.xyijE  == b.xyijE) && 
		(a.yz   == b.yz)   && (a.yzi   == b.yzi)   && (a.yzj   == b.yzj)   && (a.yzij   == b.yzij)  && 
		(a.yzE  == b.yzE)  && (a.yziE  == b.yziE)  && (a.yzjE  == b.yzjE)  && (a.yzijE  == b.yzijE) &&
		(a.xz   == b.xz)   && (a.xzi   == b.xzi)   && (a.xzj   == b.xzj)   && (a.xzij   == b.xzij)  && 
		(a.xzE  == b.xzE)  && (a.xziE  == b.xziE)  && (a.xzjE  == b.xzjE)  && (a.xzijE  == b.xzijE) ); 
	return CEQ;
}

//////////////////////////////////////////////////////

int CHO_EQ(CHO a, CHO b) {
	int CEQ = (
		(a.q    == b.q)    && (a.i     == b.i)     && (a.j     == b.j)     && (a.ij     == b.ij)    && 
		(a.E    == b.E)    && (a.iE    == b.iE)    && (a.jE    == b.jE)    && (a.ijE    == b.ijE)   && 
		(a.x    == b.x)    && (a.xi    == b.xi)    && (a.xj    == b.xj)    && (a.xij    == b.xij)   && 
		(a.xE   == b.xE)   && (a.xiE   == b.xiE)   && (a.xjE   == b.xjE)   && (a.xijE   == b.xijE)  && 
		(a.y    == b.y)    && (a.yi    == b.yi)    && (a.yj    == b.yj)    && (a.yij    == b.yij)   && 
		(a.yE   == b.yE)   && (a.yiE   == b.yiE)   && (a.yjE   == b.yjE)   && (a.yijE   == b.yijE)  && 
		(a.xy   == b.xy)   && (a.xyi   == b.xyi)   && (a.xyj   == b.xyj)   && (a.xyij   == b.xyij)  && 
		(a.xyE  == b.xyE)  && (a.xyiE  == b.xyiE)  && (a.xyjE  == b.xyjE)  && (a.xyijE  == b.xyijE) && 
		(a.z    == b.z)    && (a.zi    == b.zi)    && (a.zj    == b.zj)    && (a.zij    == b.zij)   && 
		(a.zE   == b.zE)   && (a.ziE   == b.ziE)   && (a.zjE   == b.zjE)   && (a.zijE   == b.zijE)  && 
		(a.xz   == b.xz)   && (a.xzi   == b.xzi)   && (a.xzj   == b.xzj)   && (a.xzij   == b.xzij)  && 
		(a.xzE  == b.xzE)  && (a.xziE  == b.xziE)  && (a.xzjE  == b.xzjE)  && (a.xzijE  == b.xzijE) && 
		(a.yz   == b.yz)   && (a.yzi   == b.yzi)   && (a.yzj   == b.yzj)   && (a.yzij   == b.yzij)  && 
		(a.yzE  == b.yzE)  && (a.yziE  == b.yziE)  && (a.yzjE  == b.yzjE)  && (a.yzijE  == b.yzijE) && 
		(a.xyz  == b.xyz)  && (a.xyzi  == b.xyzi)  && (a.xyzj  == b.xyzj)  && (a.xyzij  == b.xyzij) && 
		(a.xyzE == b.xyzE) && (a.xyziE == b.xyziE) && (a.xyzjE == b.xyzjE) && (a.xyzijE == b.xyzijE) );
	return CEQ;
}

//////////////////////////////////////////////////////

int C_NE(C a, C b) {
	int CNE = (
		(a.q != b.q) || (a.xyz != b.xyz) );
	return CNE;
}

//////////////////////////////////////////////////////

int H_NE(H a, H b) {
	int CNE = (
		(a.q != b.q) || (a.xy != b.xy) || (a.yz != b.yz) || (a.xz != b.xz) );
	return CNE;
}

//////////////////////////////////////////////////////

int O_NE(O a, O b) {
	int CNE = (
		(a.q != b.q) || (a.i  != b.i ) || (a.j  != b.j ) || (a.ij  != b.ij ) ||
		(a.E != b.E) || (a.iE != b.iE) || (a.jE != b.jE) || (a.ijE != b.ijE) );
	return CNE;
}

//////////////////////////////////////////////////////

int CH_NE(CH a, CH b) {
	int CNE = (
		(a.q != b.q) || (a.x  != b.x ) || (a.y  != b.y ) || (a.xy  != b.xy ) ||
		(a.z != b.z) || (a.xz != b.xz) || (a.yz != b.yz) || (a.xyz != b.xyz) );
	return CNE;
}

//////////////////////////////////////////////////////

int CO_NE(CO a, CO b) {
	int CNE = (
		(a.q    != b.q)    || (a.i     != b.i)     || (a.j     != b.j)     || (a.ij     != b.ij)    || 
		(a.E    != b.E)    || (a.iE    != b.iE)    || (a.jE    != b.jE)    || (a.ijE    != b.ijE)   || 
		(a.xyz  != b.xyz)  || (a.xyzi  != b.xyzi)  || (a.xyzj  != b.xyzj)  || (a.xyzij  != b.xyzij) || 
		(a.xyzE != b.xyzE) || (a.xyziE != b.xyziE) || (a.xyzjE != b.xyzjE) || (a.xyzijE != b.xyzijE) );
	return CNE;
}

//////////////////////////////////////////////////////

int HO_NE(HO a, HO b) {
	int CNE = (
		(a.q    != b.q)    || (a.i     != b.i)     || (a.j     != b.j)     || (a.ij     != b.ij)    || 
		(a.E    != b.E)    || (a.iE    != b.iE)    || (a.jE    != b.jE)    || (a.ijE    != b.ijE)   || 
		(a.xy   != b.xy)   || (a.xyi   != b.xyi)   || (a.xyj   != b.xyj)   || (a.xyij   != b.xyij)  || 
		(a.xyE  != b.xyE)  || (a.xyiE  != b.xyiE)  || (a.xyjE  != b.xyjE)  || (a.xyijE  != b.xyijE) || 
		(a.yz   != b.yz)   || (a.yzi   != b.yzi)   || (a.yzj   != b.yzj)   || (a.yzij   != b.yzij)  || 
		(a.yzE  != b.yzE)  || (a.yziE  != b.yziE)  || (a.yzjE  != b.yzjE)  || (a.yzijE  != b.yzijE) ||
		(a.xz   != b.xz)   || (a.xzi   != b.xzi)   || (a.xzj   != b.xzj)   || (a.xzij   != b.xzij)  || 
		(a.xzE  != b.xzE)  || (a.xziE  != b.xziE)  || (a.xzjE  != b.xzjE)  || (a.xzijE  != b.xzijE) ); 
	return CNE;
}

//////////////////////////////////////////////////////

int CHO_NE(CHO a, CHO b) {
	int CNE = (
		(a.q    != b.q)    || (a.i     != b.i)     || (a.j     != b.j)     || (a.ij     != b.ij)    || 
		(a.E    != b.E)    || (a.iE    != b.iE)    || (a.jE    != b.jE)    || (a.ijE    != b.ijE)   || 
		(a.x    != b.x)    || (a.xi    != b.xi)    || (a.xj    != b.xj)    || (a.xij    != b.xij)   || 
		(a.xE   != b.xE)   || (a.xiE   != b.xiE)   || (a.xjE   != b.xjE)   || (a.xijE   != b.xijE)  || 
		(a.y    != b.y)    || (a.yi    != b.yi)    || (a.yj    != b.yj)    || (a.yij    != b.yij)   || 
		(a.yE   != b.yE)   || (a.yiE   != b.yiE)   || (a.yjE   != b.yjE)   || (a.yijE   != b.yijE)  || 
		(a.xy   != b.xy)   || (a.xyi   != b.xyi)   || (a.xyj   != b.xyj)   || (a.xyij   != b.xyij)  || 
		(a.xyE  != b.xyE)  || (a.xyiE  != b.xyiE)  || (a.xyjE  != b.xyjE)  || (a.xyijE  != b.xyijE) || 
		(a.z    != b.z)    || (a.zi    != b.zi)    || (a.zj    != b.zj)    || (a.zij    != b.zij)   || 
		(a.zE   != b.zE)   || (a.ziE   != b.ziE)   || (a.zjE   != b.zjE)   || (a.zijE   != b.zijE)  || 
		(a.xz   != b.xz)   || (a.xzi   != b.xzi)   || (a.xzj   != b.xzj)   || (a.xzij   != b.xzij)  || 
		(a.xzE  != b.xzE)  || (a.xziE  != b.xziE)  || (a.xzjE  != b.xzjE)  || (a.xzijE  != b.xzijE) || 
		(a.yz   != b.yz)   || (a.yzi   != b.yzi)   || (a.yzj   != b.yzj)   || (a.yzij   != b.yzij)  || 
		(a.yzE  != b.yzE)  || (a.yziE  != b.yziE)  || (a.yzjE  != b.yzjE)  || (a.yzijE  != b.yzijE) || 
		(a.xyz  != b.xyz)  || (a.xyzi  != b.xyzi)  || (a.xyzj  != b.xyzj)  || (a.xyzij  != b.xyzij) || 
		(a.xyzE != b.xyzE) || (a.xyziE != b.xyziE) || (a.xyzjE != b.xyzjE) || (a.xyzijE != b.xyzijE) );
	return CNE;
}

//////////////////////////////////////////////////////

void PrintSimpleC(C a)
{
	int count = 0;

	if(a.q      == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.xyz    == +1) { printf("+xyz   "); count++;}	if(a.xyz    == -1) { printf("-xyz   "); count++;}
	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleC - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleH(H a)
{
	int count = 0;

	if(a.q      == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.xy     == +1) { printf("+xy    "); count++;}	if(a.xy     == -1) { printf("-xy    "); count++;}
	if(a.yz     == +1) { printf("+yz    "); count++;}	if(a.yz     == -1) { printf("-yz    "); count++;}
	if(a.xz     == +1) { printf("+xz    "); count++;}	if(a.xz     == -1) { printf("-xz    "); count++;}
	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleH - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleO(O a)
{
	int count = 0;

	if(a.q     == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.i     == +1) { printf("+i     "); count++;}	if(a.i      == -1) { printf("-i     "); count++;}
	if(a.j     == +1) { printf("+j     "); count++;}	if(a.j      == -1) { printf("-j     "); count++;}
	if(a.ij    == +1) { printf("+ij    "); count++;}	if(a.ij     == -1) { printf("-ij    "); count++;}
	if(a.E     == +1) { printf("+E     "); count++;}	if(a.E      == -1) { printf("-E     "); count++;}
	if(a.iE    == +1) { printf("+iE    "); count++;}	if(a.iE     == -1) { printf("-iE    "); count++;}
	if(a.jE    == +1) { printf("+jE    "); count++;}	if(a.jE     == -1) { printf("-jE    "); count++;}
	if(a.ijE   == +1) { printf("+ijE   "); count++;}	if(a.ijE    == -1) { printf("-ijE   "); count++;}
	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleO - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleCH(CH a)
{
	int count = 0;

	if(a.q     == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.x     == +1) { printf("+x     "); count++;}	if(a.x      == -1) { printf("-x     "); count++;}
	if(a.y     == +1) { printf("+y     "); count++;}	if(a.y      == -1) { printf("-y     "); count++;}
	if(a.xy    == +1) { printf("+xy    "); count++;}	if(a.xy     == -1) { printf("-xy    "); count++;}
	if(a.z     == +1) { printf("+z     "); count++;}	if(a.z      == -1) { printf("-z     "); count++;}
	if(a.xz    == +1) { printf("+xz    "); count++;}	if(a.xz     == -1) { printf("-xz    "); count++;}
	if(a.yz    == +1) { printf("+yz    "); count++;}	if(a.yz     == -1) { printf("-yz    "); count++;}
	if(a.xyz   == +1) { printf("+xyz   "); count++;}	if(a.xyz    == -1) { printf("-xyz   "); count++;}
	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleCH - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleCO(CO a)
{
	int count = 0;

	if(a.q      == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.i      == +1) { printf("+i     "); count++;}	if(a.i      == -1) { printf("-i     "); count++;}
	if(a.j      == +1) { printf("+j     "); count++;}	if(a.j      == -1) { printf("-j     "); count++;}
	if(a.ij     == +1) { printf("+ij    "); count++;}	if(a.ij     == -1) { printf("-ij    "); count++;}
	if(a.E      == +1) { printf("+E     "); count++;}	if(a.E      == -1) { printf("-E     "); count++;}
	if(a.iE     == +1) { printf("+iE    "); count++;}	if(a.iE     == -1) { printf("-iE    "); count++;}
	if(a.jE     == +1) { printf("+jE    "); count++;}	if(a.jE     == -1) { printf("-jE    "); count++;}
	if(a.ijE    == +1) { printf("+ijE   "); count++;}	if(a.ijE    == -1) { printf("-ijE   "); count++;}
	if(a.xyz    == +1) { printf("+xyz   "); count++;}	if(a.xyz    == -1) { printf("-xyz   "); count++;}
	if(a.xyzi   == +1) { printf("+xyzi  "); count++;}	if(a.xyzi   == -1) { printf("-xyzi  "); count++;}
	if(a.xyzj   == +1) { printf("+xyzj  "); count++;}	if(a.xyzj   == -1) { printf("-xyzj  "); count++;}
	if(a.xyzij  == +1) { printf("+xyzij "); count++;}	if(a.xyzij  == -1) { printf("-xyzij "); count++;}
	if(a.xyzE   == +1) { printf("+xyzE  "); count++;}	if(a.xyzE   == -1) { printf("-xyzE  "); count++;}
	if(a.xyziE  == +1) { printf("+xyziE "); count++;}	if(a.xyziE  == -1) { printf("-xyziE "); count++;}
	if(a.xyzjE  == +1) { printf("+xyzjE "); count++;}	if(a.xyzjE  == -1) { printf("-xyzjE "); count++;}
	if(a.xyzijE == +1) { printf("+xyzijE"); count++;}	if(a.xyzijE == -1) { printf("-xyzijE"); count++;}

	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleCO - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleHO(HO a)
{
	int count = 0;

	if(a.q      == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.i      == +1) { printf("+i     "); count++;}	if(a.i      == -1) { printf("-i     "); count++;}
	if(a.j      == +1) { printf("+j     "); count++;}	if(a.j      == -1) { printf("-j     "); count++;}
	if(a.ij     == +1) { printf("+ij    "); count++;}	if(a.ij     == -1) { printf("-ij    "); count++;}
	if(a.E      == +1) { printf("+E     "); count++;}	if(a.E      == -1) { printf("-E     "); count++;}
	if(a.iE     == +1) { printf("+iE    "); count++;}	if(a.iE     == -1) { printf("-iE    "); count++;}
	if(a.jE     == +1) { printf("+jE    "); count++;}	if(a.jE     == -1) { printf("-jE    "); count++;}
	if(a.ijE    == +1) { printf("+ijE   "); count++;}	if(a.ijE    == -1) { printf("-ijE   "); count++;}
	if(a.xy     == +1) { printf("+xy    "); count++;}	if(a.xy     == -1) { printf("-xy    "); count++;}
	if(a.xyi    == +1) { printf("+xyi   "); count++;}	if(a.xyi    == -1) { printf("-xyi   "); count++;}
	if(a.xyj    == +1) { printf("+xyj   "); count++;}	if(a.xyj    == -1) { printf("-xyj   "); count++;}
	if(a.xyij   == +1) { printf("+xyij  "); count++;}	if(a.xyij   == -1) { printf("-xyij  "); count++;}
	if(a.xyE    == +1) { printf("+xyE   "); count++;}	if(a.xyE    == -1) { printf("-xyE   "); count++;}
	if(a.xyiE   == +1) { printf("+xyiE  "); count++;}	if(a.xyiE   == -1) { printf("-xyiE  "); count++;}
	if(a.xyjE   == +1) { printf("+xyjE  "); count++;}	if(a.xyjE   == -1) { printf("-xyjE  "); count++;}
	if(a.xyijE  == +1) { printf("+xyijE "); count++;}	if(a.xyijE  == -1) { printf("-xyijE "); count++;}
	if(a.yz     == +1) { printf("+yz    "); count++;}	if(a.yz     == -1) { printf("-yz    "); count++;}
	if(a.yzi    == +1) { printf("+yzi   "); count++;}	if(a.yzi    == -1) { printf("-yzi   "); count++;}
	if(a.yzj    == +1) { printf("+yzj   "); count++;}	if(a.yzj    == -1) { printf("-yzj   "); count++;}
	if(a.yzij   == +1) { printf("+yzij  "); count++;}	if(a.yzij   == -1) { printf("-yzij  "); count++;}
	if(a.yzE    == +1) { printf("+yzE   "); count++;}	if(a.yzE    == -1) { printf("-yzE   "); count++;}
	if(a.yziE   == +1) { printf("+yziE  "); count++;}	if(a.yziE   == -1) { printf("-yziE  "); count++;}
	if(a.yzjE   == +1) { printf("+yzjE  "); count++;}	if(a.yzjE   == -1) { printf("-yzjE  "); count++;}
	if(a.yzijE  == +1) { printf("+yzijE "); count++;}	if(a.yzijE  == -1) { printf("-yzijE "); count++;}
	if(a.xz     == +1) { printf("+xz    "); count++;}	if(a.xz     == -1) { printf("-xz    "); count++;}
	if(a.xzi    == +1) { printf("+xzi   "); count++;}	if(a.xzi    == -1) { printf("-xzi   "); count++;}
	if(a.xzj    == +1) { printf("+xzj   "); count++;}	if(a.xzj    == -1) { printf("-xzj   "); count++;}
	if(a.xzij   == +1) { printf("+xzij  "); count++;}	if(a.xzij   == -1) { printf("-xzij  "); count++;}
	if(a.xzE    == +1) { printf("+xzE   "); count++;}	if(a.xzE    == -1) { printf("-xzE   "); count++;}
	if(a.xziE   == +1) { printf("+xziE  "); count++;}	if(a.xziE   == -1) { printf("-xziE  "); count++;}
	if(a.xzjE   == +1) { printf("+xzjE  "); count++;}	if(a.xzjE   == -1) { printf("-xzjE  "); count++;}
	if(a.xzijE  == +1) { printf("+xzijE "); count++;}	if(a.xzijE  == -1) { printf("-xzijE "); count++;}
	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleHO - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void PrintSimpleCHO(CHO a)
{
	int count = 0;

	if(a.q      == +1) { printf("+q     "); count++;}	if(a.q      == -1) { printf("-q     "); count++;}
	if(a.i      == +1) { printf("+i     "); count++;}	if(a.i      == -1) { printf("-i     "); count++;}
	if(a.j      == +1) { printf("+j     "); count++;}	if(a.j      == -1) { printf("-j     "); count++;}
	if(a.ij     == +1) { printf("+ij    "); count++;}	if(a.ij     == -1) { printf("-ij    "); count++;}
	if(a.E      == +1) { printf("+E     "); count++;}	if(a.E      == -1) { printf("-E     "); count++;}
	if(a.iE     == +1) { printf("+iE    "); count++;}	if(a.iE     == -1) { printf("-iE    "); count++;}
	if(a.jE     == +1) { printf("+jE    "); count++;}	if(a.jE     == -1) { printf("-jE    "); count++;}
	if(a.ijE    == +1) { printf("+ijE   "); count++;}	if(a.ijE    == -1) { printf("-ijE   "); count++;}
	if(a.x      == +1) { printf("+x     "); count++;}	if(a.x      == -1) { printf("-x     "); count++;}
	if(a.xi     == +1) { printf("+xi    "); count++;}	if(a.xi     == -1) { printf("-xi    "); count++;}
	if(a.xj     == +1) { printf("+xj    "); count++;}	if(a.xj     == -1) { printf("-xj    "); count++;}
	if(a.xij    == +1) { printf("+xij   "); count++;}	if(a.xij    == -1) { printf("-xij   "); count++;}
	if(a.xE     == +1) { printf("+xE    "); count++;}	if(a.xE     == -1) { printf("-xE    "); count++;}
	if(a.xiE    == +1) { printf("+xiE   "); count++;}	if(a.xiE    == -1) { printf("-xiE   "); count++;}
	if(a.xjE    == +1) { printf("+xjE   "); count++;}	if(a.xjE    == -1) { printf("-xjE   "); count++;}
	if(a.xijE   == +1) { printf("+xijE  "); count++;}	if(a.xijE   == -1) { printf("-xijE  "); count++;}
	if(a.y      == +1) { printf("+y     "); count++;}	if(a.y      == -1) { printf("-y     "); count++;}
	if(a.yi     == +1) { printf("+yi    "); count++;}	if(a.yi     == -1) { printf("-yi    "); count++;}
	if(a.yj     == +1) { printf("+yj    "); count++;}	if(a.yj     == -1) { printf("-yj    "); count++;}
	if(a.yij    == +1) { printf("+yij   "); count++;}	if(a.yij    == -1) { printf("-yij   "); count++;}
	if(a.yE     == +1) { printf("+yE    "); count++;}	if(a.yE     == -1) { printf("-yE    "); count++;}
	if(a.yiE    == +1) { printf("+yiE   "); count++;}	if(a.yiE    == -1) { printf("-yiE   "); count++;}
	if(a.yjE    == +1) { printf("+yjE   "); count++;}	if(a.yjE    == -1) { printf("-yjE   "); count++;}
	if(a.yijE   == +1) { printf("+yijE  "); count++;}	if(a.yijE   == -1) { printf("-yijE  "); count++;}
	if(a.xy     == +1) { printf("+xy    "); count++;}	if(a.xy     == -1) { printf("-xy    "); count++;}
	if(a.xyi    == +1) { printf("+xyi   "); count++;}	if(a.xyi    == -1) { printf("-xyi   "); count++;}
	if(a.xyj    == +1) { printf("+xyj   "); count++;}	if(a.xyj    == -1) { printf("-xyj   "); count++;}
	if(a.xyij   == +1) { printf("+xyij  "); count++;}	if(a.xyij   == -1) { printf("-xyij  "); count++;}
	if(a.xyE    == +1) { printf("+xyE   "); count++;}	if(a.xyE    == -1) { printf("-xyE   "); count++;}
	if(a.xyiE   == +1) { printf("+xyiE  "); count++;}	if(a.xyiE   == -1) { printf("-xyiE  "); count++;}
	if(a.xyjE   == +1) { printf("+xyjE  "); count++;}	if(a.xyjE   == -1) { printf("-xyjE  "); count++;}
	if(a.xyijE  == +1) { printf("+xyijE "); count++;}	if(a.xyijE  == -1) { printf("-xyijE "); count++;}
	if(a.z      == +1) { printf("+z     "); count++;}	if(a.z      == -1) { printf("-z     "); count++;}
	if(a.zi     == +1) { printf("+zi    "); count++;}	if(a.zi     == -1) { printf("-zi    "); count++;}
	if(a.zj     == +1) { printf("+zj    "); count++;}	if(a.zj     == -1) { printf("-zj    "); count++;}
	if(a.zij    == +1) { printf("+zij   "); count++;}	if(a.zij    == -1) { printf("-zij   "); count++;}
	if(a.zE     == +1) { printf("+zE    "); count++;}	if(a.zE     == -1) { printf("-zE    "); count++;}
	if(a.ziE    == +1) { printf("+ziE   "); count++;}	if(a.ziE    == -1) { printf("-ziE   "); count++;}
	if(a.zjE    == +1) { printf("+zjE   "); count++;}	if(a.zjE    == -1) { printf("-zjE   "); count++;}
	if(a.zijE   == +1) { printf("+zijE  "); count++;}	if(a.zijE   == -1) { printf("-zijE  "); count++;}
	if(a.xz     == +1) { printf("+xz    "); count++;}	if(a.xz     == -1) { printf("-xz    "); count++;}
	if(a.xzi    == +1) { printf("+xzi   "); count++;}	if(a.xzi    == -1) { printf("-xzi   "); count++;}
	if(a.xzj    == +1) { printf("+xzj   "); count++;}	if(a.xzj    == -1) { printf("-xzj   "); count++;}
	if(a.xzij   == +1) { printf("+xzij  "); count++;}	if(a.xzij   == -1) { printf("-xzij  "); count++;}
	if(a.xzE    == +1) { printf("+xzE   "); count++;}	if(a.xzE    == -1) { printf("-xzE   "); count++;}
	if(a.xziE   == +1) { printf("+xziE  "); count++;}	if(a.xziE   == -1) { printf("-xziE  "); count++;}
	if(a.xzjE   == +1) { printf("+xzjE  "); count++;}	if(a.xzjE   == -1) { printf("-xzjE  "); count++;}
	if(a.xzijE  == +1) { printf("+xzijE "); count++;}	if(a.xzijE  == -1) { printf("-xzijE "); count++;}
	if(a.yz     == +1) { printf("+yz    "); count++;}	if(a.yz     == -1) { printf("-yz    "); count++;}
	if(a.yzi    == +1) { printf("+yzi   "); count++;}	if(a.yzi    == -1) { printf("-yzi   "); count++;}
	if(a.yzj    == +1) { printf("+yzj   "); count++;}	if(a.yzj    == -1) { printf("-yzj   "); count++;}
	if(a.yzij   == +1) { printf("+yzij  "); count++;}	if(a.yzij   == -1) { printf("-yzij  "); count++;}
	if(a.yzE    == +1) { printf("+yzE   "); count++;}	if(a.yzE    == -1) { printf("-yzE   "); count++;}
	if(a.yziE   == +1) { printf("+yziE  "); count++;}	if(a.yziE   == -1) { printf("-yziE  "); count++;}
	if(a.yzjE   == +1) { printf("+yzjE  "); count++;}	if(a.yzjE   == -1) { printf("-yzjE  "); count++;}
	if(a.yzijE  == +1) { printf("+yzijE "); count++;}	if(a.yzijE  == -1) { printf("-yzijE "); count++;}
	if(a.xyz    == +1) { printf("+xyz   "); count++;}	if(a.xyz    == -1) { printf("-xyz   "); count++;}
	if(a.xyzi   == +1) { printf("+xyzi  "); count++;}	if(a.xyzi   == -1) { printf("-xyzi  "); count++;}
	if(a.xyzj   == +1) { printf("+xyzj  "); count++;}	if(a.xyzj   == -1) { printf("-xyzj  "); count++;}
	if(a.xyzij  == +1) { printf("+xyzij "); count++;}	if(a.xyzij  == -1) { printf("-xyzij "); count++;}
	if(a.xyzE   == +1) { printf("+xyzE  "); count++;}	if(a.xyzE   == -1) { printf("-xyzE  "); count++;}
	if(a.xyziE  == +1) { printf("+xyziE "); count++;}	if(a.xyziE  == -1) { printf("-xyziE "); count++;}
	if(a.xyzjE  == +1) { printf("+xyzjE "); count++;}	if(a.xyzjE  == -1) { printf("-xyzjE "); count++;}
	if(a.xyzijE == +1) { printf("+xyzijE"); count++;}	if(a.xyzijE == -1) { printf("-xyzijE"); count++;}

	
	if(count == 0) printf(" 0     ");
	if((count<0) || (count>1)) Program_Failed("PrintSimpleCHO - too many non-zero terms\n");
}

//////////////////////////////////////////////////////

void Print_C(C a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,    C_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyz,   C_Suffix[i++]);	printf("\n");
}

//////////////////////////////////////////////////////

void Print_H(H a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,  H_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xy, H_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yz, H_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xz, H_Suffix[i++]);	printf("\n");

}

//////////////////////////////////////////////////////

void Print_O(O a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,  O_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.i,   O_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.j,  O_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ij,  O_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.E,  O_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.iE,  O_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.jE, O_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ijE, O_Suffix[i++]);	printf("\n");

}

//////////////////////////////////////////////////////

void Print_CH(CH a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,  CH_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.x,   CH_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.y,  CH_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xy,  CH_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.z,  CH_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xz,  CH_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yz, CH_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyz, CH_Suffix[i++]);	printf("\n");

}

//////////////////////////////////////////////////////

void Print_CO(CO a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,    CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.i,     CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.j,    CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ij,    CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.E,    CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.iE,    CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.jE,   CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ijE,   CO_Suffix[i++]);	printf("\n");

	printf(" + (%10.3e)%-6s",a.xyz,   CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzi,   CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzj,  CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzij,  CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzE,  CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyziE,  CO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzjE, CO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzijE, CO_Suffix[i++]);	printf("\n");
}

//////////////////////////////////////////////////////

void Print_HO(HO a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,    HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.i,     HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.j,    HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ij,    HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.E,    HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.iE,    HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.jE,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ijE,   HO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.xy,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyi,   HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyj,  HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyij,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyE,  HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyiE,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyjE, HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyijE, HO_Suffix[i++]);	printf("\n");

	printf(" + (%10.3e)%-6s",a.yz,    HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzi,   HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzj,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzij,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzE,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yziE,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzjE,  HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzijE, HO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.xz,    HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzi,   HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzj,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzij,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzE,   HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xziE,  HO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzjE,  HO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzijE, HO_Suffix[i++]);	printf("\n");

}


//////////////////////////////////////////////////////

void Print_CHO(CHO a)
{
	int i = 0;
	printf(" + (%10.3e)%-6s",a.q,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.i,     CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.j,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ij,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.E,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.iE,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.jE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ijE,   CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.x,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xi,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xj,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xij,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xiE,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xjE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xijE,  CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.y,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yi,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yj,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yij,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yiE,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yjE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yijE,  CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.xy,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyi,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyj,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyij,  CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyiE,  CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyjE, CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyijE, CHO_Suffix[i++]);	printf("\n");

	printf(" + (%10.3e)%-6s",a.z,     CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.zi,     CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.zj,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.zij,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.zE,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.ziE,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.zjE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.zijE,   CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.xz,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzi,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzj,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzij,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xziE,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xzjE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xzijE,  CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.yz,    CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzi,    CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzj,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzij,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzE,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yziE,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.yzjE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.yzijE,  CHO_Suffix[i++]);	printf("\n");
	printf(" + (%10.3e)%-6s",a.xyz,   CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzi,   CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzj,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzij,  CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzE,  CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyziE,  CHO_Suffix[i++]);
	printf(" + (%10.3e)%-6s",a.xyzjE, CHO_Suffix[i++]);	printf(" + (%10.3e)%-6s",a.xyzijE, CHO_Suffix[i++]);	printf("\n");

}

//////////////////////////////////////////////////////

C C_Basis[2]; 

void Fill_C_Basis(void) {

	int i;

	for(i=0;i<2;i++) C_Basis[i] = C_Zero();

	C_Basis[ 0].q      = +1;	C_Basis[ 1].xyz    = +1;
}	

//////////////////////////////////////////////////////

H H_Basis[4]; 

void Fill_H_Basis(void) {

	int i;

	for(i=0;i<4;i++) H_Basis[i] = H_Zero();

	H_Basis[ 0].q = +1;	H_Basis[ 1].xy = +1;	H_Basis[ 2].yz = +1;	H_Basis[ 3].xz = +1;	
}	

//////////////////////////////////////////////////////

O O_Basis[8]; 

void Fill_O_Basis(void) {

	int i;

	for(i=0;i<8;i++) O_Basis[i] = O_Zero();

	O_Basis[ 0].q = +1;	O_Basis[ 1].i  = +1;	O_Basis[ 2].j  = +1;	O_Basis[ 3].ij  = +1;	
	O_Basis[ 4].E = +1;	O_Basis[ 5].iE = +1;	O_Basis[ 6].jE = +1;	O_Basis[ 7].ijE = +1;	
}	

//////////////////////////////////////////////////////

CH CH_Basis[8]; 

void Fill_CH_Basis(void) {

	int i;

	for(i=0;i<8;i++) CH_Basis[i] = CH_Zero();

	CH_Basis[ 0].q = +1;	CH_Basis[ 1].x  = +1;	CH_Basis[ 2].y  = +1;	CH_Basis[ 3].xy  = +1;	
	CH_Basis[ 4].z = +1;	CH_Basis[ 5].xz = +1;	CH_Basis[ 6].yz = +1;	CH_Basis[ 7].xyz = +1;	
}	

//////////////////////////////////////////////////////

CO CO_Basis[16]; 

void Fill_CO_Basis(void) {

	int i;

	for(i=0;i<16;i++) CO_Basis[i] = CO_Zero();

	CO_Basis[ 0].q      = +1;	CO_Basis[ 1].i      = +1;	CO_Basis[ 2].j      = +1;	CO_Basis[ 3].ij     = +1;	
	CO_Basis[ 4].E      = +1;	CO_Basis[ 5].iE     = +1;	CO_Basis[ 6].jE     = +1;	CO_Basis[ 7].ijE    = +1;	
	CO_Basis[ 8].xyz    = +1;	CO_Basis[ 9].xyzi   = +1;	CO_Basis[10].xyzj   = +1;	CO_Basis[11].xyzij  = +1;	
	CO_Basis[12].xyzE   = +1;	CO_Basis[13].xyziE  = +1;	CO_Basis[14].xyzjE  = +1;	CO_Basis[15].xyzijE = +1;	
}	

//////////////////////////////////////////////////////

HO HO_Basis[32]; 

void Fill_HO_Basis(void) {

	int i;

	for(i=0;i<32;i++) HO_Basis[i] = HO_Zero();

	HO_Basis[ 0].q      = +1;	HO_Basis[ 1].i      = +1;	HO_Basis[ 2].j      = +1;	HO_Basis[ 3].ij     = +1;	
	HO_Basis[ 4].E      = +1;	HO_Basis[ 5].iE     = +1;	HO_Basis[ 6].jE     = +1;	HO_Basis[ 7].ijE    = +1;	
	HO_Basis[ 8].xy     = +1;	HO_Basis[ 9].xyi    = +1;	HO_Basis[10].xyj    = +1;	HO_Basis[11].xyij   = +1;	
	HO_Basis[12].xyE    = +1;	HO_Basis[13].xyiE   = +1;	HO_Basis[14].xyjE   = +1;	HO_Basis[15].xyijE  = +1;	
	HO_Basis[16].yz     = +1;	HO_Basis[17].yzi    = +1;	HO_Basis[18].yzj    = +1;	HO_Basis[19].yzij   = +1;	
	HO_Basis[20].yzE    = +1;	HO_Basis[21].yziE   = +1;	HO_Basis[22].yzjE   = +1;	HO_Basis[23].yzijE  = +1;	
	HO_Basis[24].xz     = +1;	HO_Basis[25].xzi    = +1;	HO_Basis[26].xzj    = +1;	HO_Basis[27].xzij   = +1;	
	HO_Basis[28].xzE    = +1;	HO_Basis[29].xziE   = +1;	HO_Basis[30].xzjE   = +1;	HO_Basis[31].xzijE  = +1;	
}	
	

//////////////////////////////////////////////////////

CHO CHO_Basis[64]; 

void Fill_CHO_Basis(void) {

	int i;

	for(i=0;i<64;i++) CHO_Basis[i] = CHO_Zero();

	CHO_Basis[ 0].q      = +1;	CHO_Basis[ 1].i      = +1;	CHO_Basis[ 2].j      = +1;	CHO_Basis[ 3].ij     = +1;	
	CHO_Basis[ 4].E      = +1;	CHO_Basis[ 5].iE     = +1;	CHO_Basis[ 6].jE     = +1;	CHO_Basis[ 7].ijE    = +1;	
	CHO_Basis[ 8].x      = +1;	CHO_Basis[ 9].xi     = +1;	CHO_Basis[10].xj     = +1;	CHO_Basis[11].xij    = +1;	
	CHO_Basis[12].xE     = +1;	CHO_Basis[13].xiE    = +1;	CHO_Basis[14].xjE    = +1;	CHO_Basis[15].xijE   = +1;	
	CHO_Basis[16].y      = +1;	CHO_Basis[17].yi     = +1;	CHO_Basis[18].yj     = +1;	CHO_Basis[19].yij    = +1;	
	CHO_Basis[20].yE     = +1;	CHO_Basis[21].yiE    = +1;	CHO_Basis[22].yjE    = +1;	CHO_Basis[23].yijE   = +1;	
	CHO_Basis[24].xy     = +1;	CHO_Basis[25].xyi    = +1;	CHO_Basis[26].xyj    = +1;	CHO_Basis[27].xyij   = +1;	
	CHO_Basis[28].xyE    = +1;	CHO_Basis[29].xyiE   = +1;	CHO_Basis[30].xyjE   = +1;	CHO_Basis[31].xyijE  = +1;	
	CHO_Basis[32].z      = +1;	CHO_Basis[33].zi     = +1;	CHO_Basis[34].zj     = +1;	CHO_Basis[35].zij    = +1;	
	CHO_Basis[36].zE     = +1;	CHO_Basis[37].ziE    = +1;	CHO_Basis[38].zjE    = +1;	CHO_Basis[39].zijE   = +1;	
	CHO_Basis[40].xz     = +1;	CHO_Basis[41].xzi    = +1;	CHO_Basis[42].xzj    = +1;	CHO_Basis[43].xzij   = +1;	
	CHO_Basis[44].xzE    = +1;	CHO_Basis[45].xziE   = +1;	CHO_Basis[46].xzjE   = +1;	CHO_Basis[47].xzijE  = +1;	
	CHO_Basis[48].yz     = +1;	CHO_Basis[49].yzi    = +1;	CHO_Basis[50].yzj    = +1;	CHO_Basis[51].yzij   = +1;	
	CHO_Basis[52].yzE    = +1;	CHO_Basis[53].yziE   = +1;	CHO_Basis[54].yzjE   = +1;	CHO_Basis[55].yzijE  = +1;	
	CHO_Basis[56].xyz    = +1;	CHO_Basis[57].xyzi   = +1;	CHO_Basis[58].xyzj   = +1;	CHO_Basis[59].xyzij  = +1;	
	CHO_Basis[60].xyzE   = +1;	CHO_Basis[61].xyziE  = +1;	CHO_Basis[62].xyzjE  = +1;	CHO_Basis[63].xyzijE = +1;	
}	

//////////////////////////////////////////////////////////

char C_Sign[2][3] = {
	"++",
	"+-"
};
char H_Sign[4][5] = {
	"++++",
	"+-+-",
	"+--+",
	"++--"
};
char O_Sign[8][9] = {
	"++++++++", 
	"+-+-+--+", 
	"+--+++--", 
	"++--+-+-", 
	"+----+++", 
	"++-+---+", 
	"+++--+--", 
	"+-++--+-"
};
char CH_Sign[8][9] = {
	"++++++++",
	"++++++++",
	"+-+-+-+-",
	"+-+-+-+-",
	"+--++--+",
	"+--++--+",
	"++--++--",
	"++--++--"
};

char CO_Sign[16][17] = {
	"++++++++++++++++",
	"+-+-+--++-+-+--+",
	"+--+++--+--+++--",
	"++--+-+-++--+-+-",
	"+----++++----+++",
	"++-+---+++-+---+",
	"+++--+--+++--+--",
	"+-++--+-+-++--+-",
	"++++++++--------",
	"+-+-+--+-+-+-++-",
	"+--+++---++---++",
	"++--+-+---++-+-+",
	"+----+++-++++---",
	"++-+---+--+-+++-",
	"+++--+-----++-++",
	"+-++--+--+--++-+"
};

char HO_Sign[32][33] = {
	"++++++++++++++++++++++++++++++++",
	"+-+-+--++-+-+--++-+-+--++-+-+--+",
	"+--+++--+--+++--+--+++--+--+++--",
	"++--+-+-++--+-+-++--+-+-++--+-+-",
	"+----++++----++++----++++----+++",
	"++-+---+++-+---+++-+---+++-+---+",
	"+++--+--+++--+--+++--+--+++--+--",
	"+-++--+-+-++--+-+-++--+-+-++--+-",
	"++++++++--------++++++++--------",
	"+-+-+--+-+-+-++-+-+-+--+-+-+-++-",
	"+--+++---++---+++--+++---++---++",
	"++--+-+---++-+-+++--+-+---++-+-+",
	"+----+++-++++---+----+++-++++---",
	"++-+---+--+-+++-++-+---+--+-+++-",
	"+++--+-----++-+++++--+-----++-++",
	"+-++--+--+--++-++-++--+--+--++-+",
	"++++++++----------------++++++++",
	"+-+-+--+-+-+-++--+-+-++-+-+-+--+",
	"+--+++---++---++-++---+++--+++--",
	"++--+-+---++-+-+--++-+-+++--+-+-",
	"+----+++-++++----++++---+----+++",
	"++-+---+--+-+++---+-+++-++-+---+",
	"+++--+-----++-++---++-+++++--+--",
	"+-++--+--+--++-+-+--++-++-++--+-",
	"++++++++++++++++----------------",
	"+-+-+--++-+-+--+-+-+-++--+-+-++-",
	"+--+++--+--+++---++---++-++---++",
	"++--+-+-++--+-+---++-+-+--++-+-+",
	"+----++++----+++-++++----++++---",
	"++-+---+++-+---+--+-+++---+-+++-",
	"+++--+--+++--+-----++-++---++-++",
	"+-++--+-+-++--+--+--++-+-+--++-+"
};

char CHO_Sign[64][65] = {
	"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
	"+-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--+",
	"+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--",
	"++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-",
	"+----++++----++++----++++----++++----++++----++++----++++----+++",
	"++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+",
	"+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--",
	"+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-",
	"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
	"+-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--++-+-+--+",
	"+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--",
	"++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-",
	"+----++++----++++----++++----++++----++++----++++----++++----+++",
	"++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+++-+---+",
	"+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--+++--+--",
	"+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-+-++--+-",
	"++++++++--------++++++++--------++++++++--------++++++++--------",
	"+-+-+--+-+-+-++-+-+-+--+-+-+-++-+-+-+--+-+-+-++-+-+-+--+-+-+-++-",
	"+--+++---++---+++--+++---++---+++--+++---++---+++--+++---++---++",
	"++--+-+---++-+-+++--+-+---++-+-+++--+-+---++-+-+++--+-+---++-+-+",
	"+----+++-++++---+----+++-++++---+----+++-++++---+----+++-++++---",
	"++-+---+--+-+++-++-+---+--+-+++-++-+---+--+-+++-++-+---+--+-+++-",
	"+++--+-----++-+++++--+-----++-+++++--+-----++-+++++--+-----++-++",
	"+-++--+--+--++-++-++--+--+--++-++-++--+--+--++-++-++--+--+--++-+",
	"++++++++--------++++++++--------++++++++--------++++++++--------",
	"+-+-+--+-+-+-++-+-+-+--+-+-+-++-+-+-+--+-+-+-++-+-+-+--+-+-+-++-",
	"+--+++---++---+++--+++---++---+++--+++---++---+++--+++---++---++",
	"++--+-+---++-+-+++--+-+---++-+-+++--+-+---++-+-+++--+-+---++-+-+",
	"+----+++-++++---+----+++-++++---+----+++-++++---+----+++-++++---",
	"++-+---+--+-+++-++-+---+--+-+++-++-+---+--+-+++-++-+---+--+-+++-",
	"+++--+-----++-+++++--+-----++-+++++--+-----++-+++++--+-----++-++",
	"+-++--+--+--++-++-++--+--+--++-++-++--+--+--++-++-++--+--+--++-+",
	"++++++++----------------++++++++++++++++----------------++++++++",
	"+-+-+--+-+-+-++--+-+-++-+-+-+--++-+-+--+-+-+-++--+-+-++-+-+-+--+",
	"+--+++---++---++-++---+++--+++--+--+++---++---++-++---+++--+++--",
	"++--+-+---++-+-+--++-+-+++--+-+-++--+-+---++-+-+--++-+-+++--+-+-",
	"+----+++-++++----++++---+----++++----+++-++++----++++---+----+++",
	"++-+---+--+-+++---+-+++-++-+---+++-+---+--+-+++---+-+++-++-+---+",
	"+++--+-----++-++---++-+++++--+--+++--+-----++-++---++-+++++--+--",
	"+-++--+--+--++-+-+--++-++-++--+-+-++--+--+--++-+-+--++-++-++--+-",
	"++++++++----------------++++++++++++++++----------------++++++++",
	"+-+-+--+-+-+-++--+-+-++-+-+-+--++-+-+--+-+-+-++--+-+-++-+-+-+--+",
	"+--+++---++---++-++---+++--+++--+--+++---++---++-++---+++--+++--",
	"++--+-+---++-+-+--++-+-+++--+-+-++--+-+---++-+-+--++-+-+++--+-+-",
	"+----+++-++++----++++---+----++++----+++-++++----++++---+----+++",
	"++-+---+--+-+++---+-+++-++-+---+++-+---+--+-+++---+-+++-++-+---+",
	"+++--+-----++-++---++-+++++--+--+++--+-----++-++---++-+++++--+--",
	"+-++--+--+--++-+-+--++-++-++--+-+-++--+--+--++-+-+--++-++-++--+-",
	"++++++++++++++++----------------++++++++++++++++----------------",
	"+-+-+--++-+-+--+-+-+-++--+-+-++-+-+-+--++-+-+--+-+-+-++--+-+-++-",
	"+--+++--+--+++---++---++-++---+++--+++--+--+++---++---++-++---++",
	"++--+-+-++--+-+---++-+-+--++-+-+++--+-+-++--+-+---++-+-+--++-+-+",
	"+----++++----+++-++++----++++---+----++++----+++-++++----++++---",
	"++-+---+++-+---+--+-+++---+-+++-++-+---+++-+---+--+-+++---+-+++-",
	"+++--+--+++--+-----++-++---++-+++++--+--+++--+-----++-++---++-++",
	"+-++--+-+-++--+--+--++-+-+--++-++-++--+-+-++--+--+--++-+-+--++-+",
	"++++++++++++++++----------------++++++++++++++++----------------",
	"+-+-+--++-+-+--+-+-+-++--+-+-++-+-+-+--++-+-+--+-+-+-++--+-+-++-",
	"+--+++--+--+++---++---++-++---+++--+++--+--+++---++---++-++---++",
	"++--+-+-++--+-+---++-+-+--++-+-+++--+-+-++--+-+---++-+-+--++-+-+",
	"+----++++----+++-++++----++++---+----++++----+++-++++----++++---",
	"++-+---+++-+---+--+-+++---+-+++-++-+---+++-+---+--+-+++---+-+++-",
	"+++--+--+++--+-----++-++---++-+++++--+--+++--+-----++-++---++-++",
	"+-++--+-+-++--+--+--++-+-+--++-++-++--+-+-++--+--+--++-+-+--++-+"  
};

//////////////////////// C Product //////////////////////////////

C C_Product(C a, C b) {

	C c;

	c.q   = + a.q*b.q - a.xyz*b.xyz;
	c.xyz = + a.q*b.xyz + a.xyz*b.q;

	return c;
}

//////////////////////// H Product //////////////////////////////

H H_Product(H a, H b) {

	H c;

	c.q  =  + a.q*b.q  - a.xy*b.xy - a.yz*b.yz - a.xz*b.xz ;
	c.xy =  + a.q*b.xy + a.xy*b.q  + a.yz*b.xz - a.xz*b.yz ;
	c.yz =  + a.q*b.yz - a.xy*b.xz + a.yz*b.q  + a.xz*b.xy ;
	c.xz =  + a.q*b.xz + a.xy*b.yz - a.yz*b.xy + a.xz*b.q  ;

	return c;
}

//////////////////////// O Product //////////////////////////////

O O_Product(O a, O b) {

	O c;

	c.q   = + a.q*b.q   - a.i*b.i   - a.j*b.j   - a.ij*b.ij  - a.E*b.E   - a.iE*b.iE  - a.jE*b.jE  - a.ijE*b.ijE ; 
	c.i   = + a.q*b.i   + a.i*b.q   + a.j*b.ij  - a.ij*b.j   + a.E*b.iE  - a.iE*b.E   - a.jE*b.ijE + a.ijE*b.jE  ;
	c.j   = + a.q*b.j   - a.i*b.ij  + a.j*b.q   + a.ij*b.i   + a.E*b.jE  + a.iE*b.ijE - a.jE*b.E   - a.ijE*b.iE  ;
	c.ij  = + a.q*b.ij  + a.i*b.j   - a.j*b.i   + a.ij*b.q   + a.E*b.ijE - a.iE*b.jE  + a.jE*b.iE  - a.ijE*b.E   ;

	c.E   = + a.q*b.E   - a.i*b.iE  - a.j*b.jE  - a.ij*b.ijE + a.E*b.q   + a.iE*b.i   + a.jE*b.j   + a.ijE*b.ij  ;
	c.iE  = + a.q*b.iE  + a.i*b.E   - a.j*b.ijE + a.ij*b.jE  - a.E*b.i   + a.iE*b.q   - a.jE*b.ij  + a.ijE*b.j   ;
	c.jE  = + a.q*b.jE  + a.i*b.ijE + a.j*b.E   - a.ij*b.iE  - a.E*b.j   + a.iE*b.ij  + a.jE*b.q   - a.ijE*b.i   ;
	c.ijE = + a.q*b.ijE - a.i*b.jE  + a.j*b.iE  + a.ij*b.E   - a.E*b.ij  - a.iE*b.j   + a.jE*b.i   + a.ijE*b.q   ;

	return c;
}

//////////////////////// CH Product //////////////////////////////

CH CH_Product(CH a, CH b) {  // I am using the GA3E equations

	CH c;

	c.q   =  + a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.xz*b.xz  - a.yz*b.yz  - a.xyz*b.xyz ;
	c.x   =  + a.q*b.x   + a.x*b.q   - a.y*b.xy  - a.z*b.xz  + a.xy*b.y   + a.xz*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   =  + a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   + a.xz*b.xyz + a.yz*b.z   + a.xyz*b.xz ;
	c.xy  =  + a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   - a.xz*b.yz  + a.yz*b.xz  + a.xyz*b.z  ;

	c.z   =  + a.q*b.z   + a.x*b.xz  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz - a.xz*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xz  =  + a.q*b.xz  + a.x*b.z   - a.y*b.xyz - a.z*b.x   + a.xy*b.yz  + a.xz*b.q   - a.yz*b.xy  - a.xyz*b.y  ;
	c.yz  =  + a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   - a.xy*b.xz  + a.xz*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz =  + a.q*b.xyz + a.x*b.yz  - a.y*b.xz  + a.z*b.xy  + a.xy*b.z   - a.xz*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	return c;
}

//////////////////////// CO Product //////////////////////////////

CO CO_Product(CO a, CO b) {

	CO c;

c.q      =  + a.q*b.q      - a.i*b.i      - a.j*b.j      - a.ij*b.ij     - a.E*b.E      - a.iE*b.iE     - a.jE*b.jE     - a.ijE*b.ijE   
	 - a.xyz*b.xyz    + a.xyzi*b.xyzi   + a.xyzj*b.xyzj   + a.xyzij*b.xyzij  + a.xyzE*b.xyzE   + a.xyziE*b.xyziE  + a.xyzjE*b.xyzjE  + a.xyzijE*b.xyzijE ;
c.i      =  + a.q*b.i      + a.i*b.q      + a.j*b.ij     - a.ij*b.j      + a.E*b.iE     - a.iE*b.E      - a.jE*b.ijE    + a.ijE*b.jE    
	 - a.xyz*b.xyzi   - a.xyzi*b.xyz    - a.xyzj*b.xyzij  + a.xyzij*b.xyzj   - a.xyzE*b.xyziE  + a.xyziE*b.xyzE   + a.xyzjE*b.xyzijE - a.xyzijE*b.xyzjE  ;
c.j      =  + a.q*b.j      - a.i*b.ij     + a.j*b.q      + a.ij*b.i      + a.E*b.jE     + a.iE*b.ijE    - a.jE*b.E      - a.ijE*b.iE    
	 - a.xyz*b.xyzj   + a.xyzi*b.xyzij  - a.xyzj*b.xyz    - a.xyzij*b.xyzi   - a.xyzE*b.xyzjE  - a.xyziE*b.xyzijE + a.xyzjE*b.xyzE   + a.xyzijE*b.xyziE  ;
c.ij     =  + a.q*b.ij     + a.i*b.j      - a.j*b.i      + a.ij*b.q      + a.E*b.ijE    - a.iE*b.jE     + a.jE*b.iE     - a.ijE*b.E     
	 - a.xyz*b.xyzij  - a.xyzi*b.xyzj   + a.xyzj*b.xyzi   - a.xyzij*b.xyz    - a.xyzE*b.xyzijE + a.xyziE*b.xyzjE  - a.xyzjE*b.xyziE  + a.xyzijE*b.xyzE   ;
c.E      =  + a.q*b.E      - a.i*b.iE     - a.j*b.jE     - a.ij*b.ijE    + a.E*b.q      + a.iE*b.i      + a.jE*b.j      + a.ijE*b.ij    
	 - a.xyz*b.xyzE   + a.xyzi*b.xyziE  + a.xyzj*b.xyzjE  + a.xyzij*b.xyzijE - a.xyzE*b.xyz    - a.xyziE*b.xyzi   - a.xyzjE*b.xyzj   - a.xyzijE*b.xyzij  ;
c.iE     =  + a.q*b.iE     + a.i*b.E      - a.j*b.ijE    + a.ij*b.jE     - a.E*b.i      + a.iE*b.q      - a.jE*b.ij     + a.ijE*b.j     
	 - a.xyz*b.xyziE  - a.xyzi*b.xyzE   + a.xyzj*b.xyzijE - a.xyzij*b.xyzjE  + a.xyzE*b.xyzi   - a.xyziE*b.xyz    + a.xyzjE*b.xyzij  - a.xyzijE*b.xyzj   ;
c.jE     =  + a.q*b.jE     + a.i*b.ijE    + a.j*b.E      - a.ij*b.iE     - a.E*b.j      + a.iE*b.ij     + a.jE*b.q      - a.ijE*b.i     
	 - a.xyz*b.xyzjE  - a.xyzi*b.xyzijE - a.xyzj*b.xyzE   + a.xyzij*b.xyziE  + a.xyzE*b.xyzj   - a.xyziE*b.xyzij  - a.xyzjE*b.xyz    + a.xyzijE*b.xyzi   ;
c.ijE    =  + a.q*b.ijE    - a.i*b.jE     + a.j*b.iE     + a.ij*b.E      - a.E*b.ij     - a.iE*b.j      + a.jE*b.i      + a.ijE*b.q     
	 - a.xyz*b.xyzijE + a.xyzi*b.xyzjE  - a.xyzj*b.xyziE  - a.xyzij*b.xyzE   + a.xyzE*b.xyzij  + a.xyziE*b.xyzj   - a.xyzjE*b.xyzi   - a.xyzijE*b.xyz    ;
c.xyz    =  + a.q*b.xyz    - a.i*b.xyzi   - a.j*b.xyzj   - a.ij*b.xyzij  - a.E*b.xyzE   - a.iE*b.xyziE  - a.jE*b.xyzjE  - a.ijE*b.xyzijE
	 + a.xyz*b.q      - a.xyzi*b.i      - a.xyzj*b.j      - a.xyzij*b.ij     - a.xyzE*b.E      - a.xyziE*b.iE     - a.xyzjE*b.jE     - a.xyzijE*b.ijE    ;
c.xyzi   =  + a.q*b.xyzi   + a.i*b.xyz    + a.j*b.xyzij  - a.ij*b.xyzj   + a.E*b.xyziE  - a.iE*b.xyzE   - a.jE*b.xyzijE + a.ijE*b.xyzjE 
	 + a.xyz*b.i      + a.xyzi*b.q      + a.xyzj*b.ij     - a.xyzij*b.j      + a.xyzE*b.iE     - a.xyziE*b.E      - a.xyzjE*b.ijE    + a.xyzijE*b.jE     ;
c.xyzj   =  + a.q*b.xyzj   - a.i*b.xyzij  + a.j*b.xyz    + a.ij*b.xyzi   + a.E*b.xyzjE  + a.iE*b.xyzijE - a.jE*b.xyzE   - a.ijE*b.xyziE 
	 + a.xyz*b.j      - a.xyzi*b.ij     + a.xyzj*b.q      + a.xyzij*b.i      + a.xyzE*b.jE     + a.xyziE*b.ijE    - a.xyzjE*b.E      - a.xyzijE*b.iE     ;
c.xyzij  =  + a.q*b.xyzij  + a.i*b.xyzj   - a.j*b.xyzi   + a.ij*b.xyz    + a.E*b.xyzijE - a.iE*b.xyzjE  + a.jE*b.xyziE  - a.ijE*b.xyzE  
	 + a.xyz*b.ij     + a.xyzi*b.j      - a.xyzj*b.i      + a.xyzij*b.q      + a.xyzE*b.ijE    - a.xyziE*b.jE     + a.xyzjE*b.iE     - a.xyzijE*b.E      ;
c.xyzE   =  + a.q*b.xyzE   - a.i*b.xyziE  - a.j*b.xyzjE  - a.ij*b.xyzijE + a.E*b.xyz    + a.iE*b.xyzi   + a.jE*b.xyzj   + a.ijE*b.xyzij 
	 + a.xyz*b.E      - a.xyzi*b.iE     - a.xyzj*b.jE     - a.xyzij*b.ijE    + a.xyzE*b.q      + a.xyziE*b.i      + a.xyzjE*b.j      + a.xyzijE*b.ij     ;
c.xyziE  =  + a.q*b.xyziE  + a.i*b.xyzE   - a.j*b.xyzijE + a.ij*b.xyzjE  - a.E*b.xyzi   + a.iE*b.xyz    - a.jE*b.xyzij  + a.ijE*b.xyzj  
	 + a.xyz*b.iE     + a.xyzi*b.E      - a.xyzj*b.ijE    + a.xyzij*b.jE     - a.xyzE*b.i      + a.xyziE*b.q      - a.xyzjE*b.ij     + a.xyzijE*b.j      ;
c.xyzjE  =  + a.q*b.xyzjE  + a.i*b.xyzijE + a.j*b.xyzE   - a.ij*b.xyziE  - a.E*b.xyzj   + a.iE*b.xyzij  + a.jE*b.xyz    - a.ijE*b.xyzi  
	 + a.xyz*b.jE     + a.xyzi*b.ijE    + a.xyzj*b.E      - a.xyzij*b.iE     - a.xyzE*b.j      + a.xyziE*b.ij     + a.xyzjE*b.q      - a.xyzijE*b.i      ;
c.xyzijE =  + a.q*b.xyzijE - a.i*b.xyzjE  + a.j*b.xyziE  + a.ij*b.xyzE   - a.E*b.xyzij  - a.iE*b.xyzj   + a.jE*b.xyzi   + a.ijE*b.xyz   
	 + a.xyz*b.ijE    - a.xyzi*b.jE     + a.xyzj*b.iE     + a.xyzij*b.E      - a.xyzE*b.ij     - a.xyziE*b.j      + a.xyzjE*b.i      + a.xyzijE*b.q      ;

	return c;
}

//////////////////////// HO //////////////////////////////

HO HO_Product(HO a, HO b) {

	HO c;

c.q      =  + a.q*b.q      - a.i*b.i      - a.j*b.j      - a.ij*b.ij     - a.E*b.E      - a.iE*b.iE     - a.jE*b.jE     - a.ijE*b.ijE   
	 - a.xy*b.xy     + a.xyi*b.xyi    + a.xyj*b.xyj    + a.xyij*b.xyij   + a.xyE*b.xyE    + a.xyiE*b.xyiE   + a.xyjE*b.xyjE   + a.xyijE*b.xyijE 
	 - a.yz*b.yz     + a.yzi*b.yzi    + a.yzj*b.yzj    + a.yzij*b.yzij   + a.yzE*b.yzE    + a.yziE*b.yziE   + a.yzjE*b.yzjE   + a.yzijE*b.yzijE 
	 - a.xz*b.xz     + a.xzi*b.xzi    + a.xzj*b.xzj    + a.xzij*b.xzij   + a.xzE*b.xzE    + a.xziE*b.xziE   + a.xzjE*b.xzjE   + a.xzijE*b.xzijE  ;
c.i      =  + a.q*b.i      + a.i*b.q      + a.j*b.ij     - a.ij*b.j      + a.E*b.iE     - a.iE*b.E      - a.jE*b.ijE    + a.ijE*b.jE    
	 - a.xy*b.xyi    - a.xyi*b.xy     - a.xyj*b.xyij   + a.xyij*b.xyj    - a.xyE*b.xyiE   + a.xyiE*b.xyE    + a.xyjE*b.xyijE  - a.xyijE*b.xyjE  
	 - a.yz*b.yzi    - a.yzi*b.yz     - a.yzj*b.yzij   + a.yzij*b.yzj    - a.yzE*b.yziE   + a.yziE*b.yzE    + a.yzjE*b.yzijE  - a.yzijE*b.yzjE  
	 - a.xz*b.xzi    - a.xzi*b.xz     - a.xzj*b.xzij   + a.xzij*b.xzj    - a.xzE*b.xziE   + a.xziE*b.xzE    + a.xzjE*b.xzijE  - a.xzijE*b.xzjE   ;
c.j      =  + a.q*b.j      - a.i*b.ij     + a.j*b.q      + a.ij*b.i      + a.E*b.jE     + a.iE*b.ijE    - a.jE*b.E      - a.ijE*b.iE    
	 - a.xy*b.xyj    + a.xyi*b.xyij   - a.xyj*b.xy     - a.xyij*b.xyi    - a.xyE*b.xyjE   - a.xyiE*b.xyijE  + a.xyjE*b.xyE    + a.xyijE*b.xyiE  
	 - a.yz*b.yzj    + a.yzi*b.yzij   - a.yzj*b.yz     - a.yzij*b.yzi    - a.yzE*b.yzjE   - a.yziE*b.yzijE  + a.yzjE*b.yzE    + a.yzijE*b.yziE  
	 - a.xz*b.xzj    + a.xzi*b.xzij   - a.xzj*b.xz     - a.xzij*b.xzi    - a.xzE*b.xzjE   - a.xziE*b.xzijE  + a.xzjE*b.xzE    + a.xzijE*b.xziE   ;
c.ij     =  + a.q*b.ij     + a.i*b.j      - a.j*b.i      + a.ij*b.q      + a.E*b.ijE    - a.iE*b.jE     + a.jE*b.iE     - a.ijE*b.E     
	 - a.xy*b.xyij   - a.xyi*b.xyj    + a.xyj*b.xyi    - a.xyij*b.xy     - a.xyE*b.xyijE  + a.xyiE*b.xyjE   - a.xyjE*b.xyiE   + a.xyijE*b.xyE   
	 - a.yz*b.yzij   - a.yzi*b.yzj    + a.yzj*b.yzi    - a.yzij*b.yz     - a.yzE*b.yzijE  + a.yziE*b.yzjE   - a.yzjE*b.yziE   + a.yzijE*b.yzE   
	 - a.xz*b.xzij   - a.xzi*b.xzj    + a.xzj*b.xzi    - a.xzij*b.xz     - a.xzE*b.xzijE  + a.xziE*b.xzjE   - a.xzjE*b.xziE   + a.xzijE*b.xzE    ;
c.E      =  + a.q*b.E      - a.i*b.iE     - a.j*b.jE     - a.ij*b.ijE    + a.E*b.q      + a.iE*b.i      + a.jE*b.j      + a.ijE*b.ij    
	 - a.xy*b.xyE    + a.xyi*b.xyiE   + a.xyj*b.xyjE   + a.xyij*b.xyijE  - a.xyE*b.xy     - a.xyiE*b.xyi    - a.xyjE*b.xyj    - a.xyijE*b.xyij  
	 - a.yz*b.yzE    + a.yzi*b.yziE   + a.yzj*b.yzjE   + a.yzij*b.yzijE  - a.yzE*b.yz     - a.yziE*b.yzi    - a.yzjE*b.yzj    - a.yzijE*b.yzij  
	 - a.xz*b.xzE    + a.xzi*b.xziE   + a.xzj*b.xzjE   + a.xzij*b.xzijE  - a.xzE*b.xz     - a.xziE*b.xzi    - a.xzjE*b.xzj    - a.xzijE*b.xzij   ;
c.iE     =  + a.q*b.iE     + a.i*b.E      - a.j*b.ijE    + a.ij*b.jE     - a.E*b.i      + a.iE*b.q      - a.jE*b.ij     + a.ijE*b.j     
	 - a.xy*b.xyiE   - a.xyi*b.xyE    + a.xyj*b.xyijE  - a.xyij*b.xyjE   + a.xyE*b.xyi    - a.xyiE*b.xy     + a.xyjE*b.xyij   - a.xyijE*b.xyj   
	 - a.yz*b.yziE   - a.yzi*b.yzE    + a.yzj*b.yzijE  - a.yzij*b.yzjE   + a.yzE*b.yzi    - a.yziE*b.yz     + a.yzjE*b.yzij   - a.yzijE*b.yzj   
	 - a.xz*b.xziE   - a.xzi*b.xzE    + a.xzj*b.xzijE  - a.xzij*b.xzjE   + a.xzE*b.xzi    - a.xziE*b.xz     + a.xzjE*b.xzij   - a.xzijE*b.xzj    ;
c.jE     =  + a.q*b.jE     + a.i*b.ijE    + a.j*b.E      - a.ij*b.iE     - a.E*b.j      + a.iE*b.ij     + a.jE*b.q      - a.ijE*b.i     
	 - a.xy*b.xyjE   - a.xyi*b.xyijE  - a.xyj*b.xyE    + a.xyij*b.xyiE   + a.xyE*b.xyj    - a.xyiE*b.xyij   - a.xyjE*b.xy     + a.xyijE*b.xyi   
	 - a.yz*b.yzjE   - a.yzi*b.yzijE  - a.yzj*b.yzE    + a.yzij*b.yziE   + a.yzE*b.yzj    - a.yziE*b.yzij   - a.yzjE*b.yz     + a.yzijE*b.yzi   
	 - a.xz*b.xzjE   - a.xzi*b.xzijE  - a.xzj*b.xzE    + a.xzij*b.xziE   + a.xzE*b.xzj    - a.xziE*b.xzij   - a.xzjE*b.xz     + a.xzijE*b.xzi    ;
c.ijE    =  + a.q*b.ijE    - a.i*b.jE     + a.j*b.iE     + a.ij*b.E      - a.E*b.ij     - a.iE*b.j      + a.jE*b.i      + a.ijE*b.q     
	 - a.xy*b.xyijE  + a.xyi*b.xyjE   - a.xyj*b.xyiE   - a.xyij*b.xyE    + a.xyE*b.xyij   + a.xyiE*b.xyj    - a.xyjE*b.xyi    - a.xyijE*b.xy    
	 - a.yz*b.yzijE  + a.yzi*b.yzjE   - a.yzj*b.yziE   - a.yzij*b.yzE    + a.yzE*b.yzij   + a.yziE*b.yzj    - a.yzjE*b.yzi    - a.yzijE*b.yz    
	 - a.xz*b.xzijE  + a.xzi*b.xzjE   - a.xzj*b.xziE   - a.xzij*b.xzE    + a.xzE*b.xzij   + a.xziE*b.xzj    - a.xzjE*b.xzi    - a.xzijE*b.xz     ;
c.xy     =  + a.q*b.xy     - a.i*b.xyi    - a.j*b.xyj    - a.ij*b.xyij   - a.E*b.xyE    - a.iE*b.xyiE   - a.jE*b.xyjE   - a.ijE*b.xyijE 
	 + a.xy*b.q      - a.xyi*b.i      - a.xyj*b.j      - a.xyij*b.ij     - a.xyE*b.E      - a.xyiE*b.iE     - a.xyjE*b.jE     - a.xyijE*b.ijE   
	 + a.yz*b.xz     - a.yzi*b.xzi    - a.yzj*b.xzj    - a.yzij*b.xzij   - a.yzE*b.xzE    - a.yziE*b.xziE   - a.yzjE*b.xzjE   - a.yzijE*b.xzijE 
	 - a.xz*b.yz     + a.xzi*b.yzi    + a.xzj*b.yzj    + a.xzij*b.yzij   + a.xzE*b.yzE    + a.xziE*b.yziE   + a.xzjE*b.yzjE   + a.xzijE*b.yzijE  ;
c.xyi    =  + a.q*b.xyi    + a.i*b.xy     + a.j*b.xyij   - a.ij*b.xyj    + a.E*b.xyiE   - a.iE*b.xyE    - a.jE*b.xyijE  + a.ijE*b.xyjE  
	 + a.xy*b.i      + a.xyi*b.q      + a.xyj*b.ij     - a.xyij*b.j      + a.xyE*b.iE     - a.xyiE*b.E      - a.xyjE*b.ijE    + a.xyijE*b.jE    
	 + a.yz*b.xzi    + a.yzi*b.xz     + a.yzj*b.xzij   - a.yzij*b.xzj    + a.yzE*b.xziE   - a.yziE*b.xzE    - a.yzjE*b.xzijE  + a.yzijE*b.xzjE  
	 - a.xz*b.yzi    - a.xzi*b.yz     - a.xzj*b.yzij   + a.xzij*b.yzj    - a.xzE*b.yziE   + a.xziE*b.yzE    + a.xzjE*b.yzijE  - a.xzijE*b.yzjE   ;
c.xyj    =  + a.q*b.xyj    - a.i*b.xyij   + a.j*b.xy     + a.ij*b.xyi    + a.E*b.xyjE   + a.iE*b.xyijE  - a.jE*b.xyE    - a.ijE*b.xyiE  
	 + a.xy*b.j      - a.xyi*b.ij     + a.xyj*b.q      + a.xyij*b.i      + a.xyE*b.jE     + a.xyiE*b.ijE    - a.xyjE*b.E      - a.xyijE*b.iE    
	 + a.yz*b.xzj    - a.yzi*b.xzij   + a.yzj*b.xz     + a.yzij*b.xzi    + a.yzE*b.xzjE   + a.yziE*b.xzijE  - a.yzjE*b.xzE    - a.yzijE*b.xziE  
	 - a.xz*b.yzj    + a.xzi*b.yzij   - a.xzj*b.yz     - a.xzij*b.yzi    - a.xzE*b.yzjE   - a.xziE*b.yzijE  + a.xzjE*b.yzE    + a.xzijE*b.yziE   ;
c.xyij   =  + a.q*b.xyij   + a.i*b.xyj    - a.j*b.xyi    + a.ij*b.xy     + a.E*b.xyijE  - a.iE*b.xyjE   + a.jE*b.xyiE   - a.ijE*b.xyE   
	 + a.xy*b.ij     + a.xyi*b.j      - a.xyj*b.i      + a.xyij*b.q      + a.xyE*b.ijE    - a.xyiE*b.jE     + a.xyjE*b.iE     - a.xyijE*b.E     
	 + a.yz*b.xzij   + a.yzi*b.xzj    - a.yzj*b.xzi    + a.yzij*b.xz     + a.yzE*b.xzijE  - a.yziE*b.xzjE   + a.yzjE*b.xziE   - a.yzijE*b.xzE   
	 - a.xz*b.yzij   - a.xzi*b.yzj    + a.xzj*b.yzi    - a.xzij*b.yz     - a.xzE*b.yzijE  + a.xziE*b.yzjE   - a.xzjE*b.yziE   + a.xzijE*b.yzE    ;
c.xyE    =  + a.q*b.xyE    - a.i*b.xyiE   - a.j*b.xyjE   - a.ij*b.xyijE  + a.E*b.xy     + a.iE*b.xyi    + a.jE*b.xyj    + a.ijE*b.xyij  
	 + a.xy*b.E      - a.xyi*b.iE     - a.xyj*b.jE     - a.xyij*b.ijE    + a.xyE*b.q      + a.xyiE*b.i      + a.xyjE*b.j      + a.xyijE*b.ij    
	 + a.yz*b.xzE    - a.yzi*b.xziE   - a.yzj*b.xzjE   - a.yzij*b.xzijE  + a.yzE*b.xz     + a.yziE*b.xzi    + a.yzjE*b.xzj    + a.yzijE*b.xzij  
	 - a.xz*b.yzE    + a.xzi*b.yziE   + a.xzj*b.yzjE   + a.xzij*b.yzijE  - a.xzE*b.yz     - a.xziE*b.yzi    - a.xzjE*b.yzj    - a.xzijE*b.yzij   ;
c.xyiE   =  + a.q*b.xyiE   + a.i*b.xyE    - a.j*b.xyijE  + a.ij*b.xyjE   - a.E*b.xyi    + a.iE*b.xy     - a.jE*b.xyij   + a.ijE*b.xyj   
	 + a.xy*b.iE     + a.xyi*b.E      - a.xyj*b.ijE    + a.xyij*b.jE     - a.xyE*b.i      + a.xyiE*b.q      - a.xyjE*b.ij     + a.xyijE*b.j     
	 + a.yz*b.xziE   + a.yzi*b.xzE    - a.yzj*b.xzijE  + a.yzij*b.xzjE   - a.yzE*b.xzi    + a.yziE*b.xz     - a.yzjE*b.xzij   + a.yzijE*b.xzj   
	 - a.xz*b.yziE   - a.xzi*b.yzE    + a.xzj*b.yzijE  - a.xzij*b.yzjE   + a.xzE*b.yzi    - a.xziE*b.yz     + a.xzjE*b.yzij   - a.xzijE*b.yzj    ;
c.xyjE   =  + a.q*b.xyjE   + a.i*b.xyijE  + a.j*b.xyE    - a.ij*b.xyiE   - a.E*b.xyj    + a.iE*b.xyij   + a.jE*b.xy     - a.ijE*b.xyi   
	 + a.xy*b.jE     + a.xyi*b.ijE    + a.xyj*b.E      - a.xyij*b.iE     - a.xyE*b.j      + a.xyiE*b.ij     + a.xyjE*b.q      - a.xyijE*b.i     
	 + a.yz*b.xzjE   + a.yzi*b.xzijE  + a.yzj*b.xzE    - a.yzij*b.xziE   - a.yzE*b.xzj    + a.yziE*b.xzij   + a.yzjE*b.xz     - a.yzijE*b.xzi   
	 - a.xz*b.yzjE   - a.xzi*b.yzijE  - a.xzj*b.yzE    + a.xzij*b.yziE   + a.xzE*b.yzj    - a.xziE*b.yzij   - a.xzjE*b.yz     + a.xzijE*b.yzi    ;
c.xyijE  =  + a.q*b.xyijE  - a.i*b.xyjE   + a.j*b.xyiE   + a.ij*b.xyE    - a.E*b.xyij   - a.iE*b.xyj    + a.jE*b.xyi    + a.ijE*b.xy    
	 + a.xy*b.ijE    - a.xyi*b.jE     + a.xyj*b.iE     + a.xyij*b.E      - a.xyE*b.ij     - a.xyiE*b.j      + a.xyjE*b.i      + a.xyijE*b.q     
	 + a.yz*b.xzijE  - a.yzi*b.xzjE   + a.yzj*b.xziE   + a.yzij*b.xzE    - a.yzE*b.xzij   - a.yziE*b.xzj    + a.yzjE*b.xzi    + a.yzijE*b.xz    
	 - a.xz*b.yzijE  + a.xzi*b.yzjE   - a.xzj*b.yziE   - a.xzij*b.yzE    + a.xzE*b.yzij   + a.xziE*b.yzj    - a.xzjE*b.yzi    - a.xzijE*b.yz     ;
c.yz     =  + a.q*b.yz     - a.i*b.yzi    - a.j*b.yzj    - a.ij*b.yzij   - a.E*b.yzE    - a.iE*b.yziE   - a.jE*b.yzjE   - a.ijE*b.yzijE 
	 - a.xy*b.xz     + a.xyi*b.xzi    + a.xyj*b.xzj    + a.xyij*b.xzij   + a.xyE*b.xzE    + a.xyiE*b.xziE   + a.xyjE*b.xzjE   + a.xyijE*b.xzijE 
	 + a.yz*b.q      - a.yzi*b.i      - a.yzj*b.j      - a.yzij*b.ij     - a.yzE*b.E      - a.yziE*b.iE     - a.yzjE*b.jE     - a.yzijE*b.ijE   
	 + a.xz*b.xy     - a.xzi*b.xyi    - a.xzj*b.xyj    - a.xzij*b.xyij   - a.xzE*b.xyE    - a.xziE*b.xyiE   - a.xzjE*b.xyjE   - a.xzijE*b.xyijE  ;
c.yzi    =  + a.q*b.yzi    + a.i*b.yz     + a.j*b.yzij   - a.ij*b.yzj    + a.E*b.yziE   - a.iE*b.yzE    - a.jE*b.yzijE  + a.ijE*b.yzjE  
	 - a.xy*b.xzi    - a.xyi*b.xz     - a.xyj*b.xzij   + a.xyij*b.xzj    - a.xyE*b.xziE   + a.xyiE*b.xzE    + a.xyjE*b.xzijE  - a.xyijE*b.xzjE  
	 + a.yz*b.i      + a.yzi*b.q      + a.yzj*b.ij     - a.yzij*b.j      + a.yzE*b.iE     - a.yziE*b.E      - a.yzjE*b.ijE    + a.yzijE*b.jE    
	 + a.xz*b.xyi    + a.xzi*b.xy     + a.xzj*b.xyij   - a.xzij*b.xyj    + a.xzE*b.xyiE   - a.xziE*b.xyE    - a.xzjE*b.xyijE  + a.xzijE*b.xyjE   ;
c.yzj    =  + a.q*b.yzj    - a.i*b.yzij   + a.j*b.yz     + a.ij*b.yzi    + a.E*b.yzjE   + a.iE*b.yzijE  - a.jE*b.yzE    - a.ijE*b.yziE  
	 - a.xy*b.xzj    + a.xyi*b.xzij   - a.xyj*b.xz     - a.xyij*b.xzi    - a.xyE*b.xzjE   - a.xyiE*b.xzijE  + a.xyjE*b.xzE    + a.xyijE*b.xziE  
	 + a.yz*b.j      - a.yzi*b.ij     + a.yzj*b.q      + a.yzij*b.i      + a.yzE*b.jE     + a.yziE*b.ijE    - a.yzjE*b.E      - a.yzijE*b.iE    
	 + a.xz*b.xyj    - a.xzi*b.xyij   + a.xzj*b.xy     + a.xzij*b.xyi    + a.xzE*b.xyjE   + a.xziE*b.xyijE  - a.xzjE*b.xyE    - a.xzijE*b.xyiE   ;
c.yzij   =  + a.q*b.yzij   + a.i*b.yzj    - a.j*b.yzi    + a.ij*b.yz     + a.E*b.yzijE  - a.iE*b.yzjE   + a.jE*b.yziE   - a.ijE*b.yzE   
	 - a.xy*b.xzij   - a.xyi*b.xzj    + a.xyj*b.xzi    - a.xyij*b.xz     - a.xyE*b.xzijE  + a.xyiE*b.xzjE   - a.xyjE*b.xziE   + a.xyijE*b.xzE   
	 + a.yz*b.ij     + a.yzi*b.j      - a.yzj*b.i      + a.yzij*b.q      + a.yzE*b.ijE    - a.yziE*b.jE     + a.yzjE*b.iE     - a.yzijE*b.E     
	 + a.xz*b.xyij   + a.xzi*b.xyj    - a.xzj*b.xyi    + a.xzij*b.xy     + a.xzE*b.xyijE  - a.xziE*b.xyjE   + a.xzjE*b.xyiE   - a.xzijE*b.xyE    ;
c.yzE    =  + a.q*b.yzE    - a.i*b.yziE   - a.j*b.yzjE   - a.ij*b.yzijE  + a.E*b.yz     + a.iE*b.yzi    + a.jE*b.yzj    + a.ijE*b.yzij  
	 - a.xy*b.xzE    + a.xyi*b.xziE   + a.xyj*b.xzjE   + a.xyij*b.xzijE  - a.xyE*b.xz     - a.xyiE*b.xzi    - a.xyjE*b.xzj    - a.xyijE*b.xzij  
	 + a.yz*b.E      - a.yzi*b.iE     - a.yzj*b.jE     - a.yzij*b.ijE    + a.yzE*b.q      + a.yziE*b.i      + a.yzjE*b.j      + a.yzijE*b.ij    
	 + a.xz*b.xyE    - a.xzi*b.xyiE   - a.xzj*b.xyjE   - a.xzij*b.xyijE  + a.xzE*b.xy     + a.xziE*b.xyi    + a.xzjE*b.xyj    + a.xzijE*b.xyij   ;
c.yziE   =  + a.q*b.yziE   + a.i*b.yzE    - a.j*b.yzijE  + a.ij*b.yzjE   - a.E*b.yzi    + a.iE*b.yz     - a.jE*b.yzij   + a.ijE*b.yzj   
	 - a.xy*b.xziE   - a.xyi*b.xzE    + a.xyj*b.xzijE  - a.xyij*b.xzjE   + a.xyE*b.xzi    - a.xyiE*b.xz     + a.xyjE*b.xzij   - a.xyijE*b.xzj   
	 + a.yz*b.iE     + a.yzi*b.E      - a.yzj*b.ijE    + a.yzij*b.jE     - a.yzE*b.i      + a.yziE*b.q      - a.yzjE*b.ij     + a.yzijE*b.j     
	 + a.xz*b.xyiE   + a.xzi*b.xyE    - a.xzj*b.xyijE  + a.xzij*b.xyjE   - a.xzE*b.xyi    + a.xziE*b.xy     - a.xzjE*b.xyij   + a.xzijE*b.xyj    ;
c.yzjE   =  + a.q*b.yzjE   + a.i*b.yzijE  + a.j*b.yzE    - a.ij*b.yziE   - a.E*b.yzj    + a.iE*b.yzij   + a.jE*b.yz     - a.ijE*b.yzi   
	 - a.xy*b.xzjE   - a.xyi*b.xzijE  - a.xyj*b.xzE    + a.xyij*b.xziE   + a.xyE*b.xzj    - a.xyiE*b.xzij   - a.xyjE*b.xz     + a.xyijE*b.xzi   
	 + a.yz*b.jE     + a.yzi*b.ijE    + a.yzj*b.E      - a.yzij*b.iE     - a.yzE*b.j      + a.yziE*b.ij     + a.yzjE*b.q      - a.yzijE*b.i     
	 + a.xz*b.xyjE   + a.xzi*b.xyijE  + a.xzj*b.xyE    - a.xzij*b.xyiE   - a.xzE*b.xyj    + a.xziE*b.xyij   + a.xzjE*b.xy     - a.xzijE*b.xyi    ;
c.yzijE  =  + a.q*b.yzijE  - a.i*b.yzjE   + a.j*b.yziE   + a.ij*b.yzE    - a.E*b.yzij   - a.iE*b.yzj    + a.jE*b.yzi    + a.ijE*b.yz    
	 - a.xy*b.xzijE  + a.xyi*b.xzjE   - a.xyj*b.xziE   - a.xyij*b.xzE    + a.xyE*b.xzij   + a.xyiE*b.xzj    - a.xyjE*b.xzi    - a.xyijE*b.xz    
	 + a.yz*b.ijE    - a.yzi*b.jE     + a.yzj*b.iE     + a.yzij*b.E      - a.yzE*b.ij     - a.yziE*b.j      + a.yzjE*b.i      + a.yzijE*b.q     
	 + a.xz*b.xyijE  - a.xzi*b.xyjE   + a.xzj*b.xyiE   + a.xzij*b.xyE    - a.xzE*b.xyij   - a.xziE*b.xyj    + a.xzjE*b.xyi    + a.xzijE*b.xy     ;
c.xz     =  + a.q*b.xz     - a.i*b.xzi    - a.j*b.xzj    - a.ij*b.xzij   - a.E*b.xzE    - a.iE*b.xziE   - a.jE*b.xzjE   - a.ijE*b.xzijE 
	 + a.xy*b.yz     - a.xyi*b.yzi    - a.xyj*b.yzj    - a.xyij*b.yzij   - a.xyE*b.yzE    - a.xyiE*b.yziE   - a.xyjE*b.yzjE   - a.xyijE*b.yzijE 
	 - a.yz*b.xy     + a.yzi*b.xyi    + a.yzj*b.xyj    + a.yzij*b.xyij   + a.yzE*b.xyE    + a.yziE*b.xyiE   + a.yzjE*b.xyjE   + a.yzijE*b.xyijE 
	 + a.xz*b.q      - a.xzi*b.i      - a.xzj*b.j      - a.xzij*b.ij     - a.xzE*b.E      - a.xziE*b.iE     - a.xzjE*b.jE     - a.xzijE*b.ijE    ;
c.xzi    =  + a.q*b.xzi    + a.i*b.xz     + a.j*b.xzij   - a.ij*b.xzj    + a.E*b.xziE   - a.iE*b.xzE    - a.jE*b.xzijE  + a.ijE*b.xzjE  
	 + a.xy*b.yzi    + a.xyi*b.yz     + a.xyj*b.yzij   - a.xyij*b.yzj    + a.xyE*b.yziE   - a.xyiE*b.yzE    - a.xyjE*b.yzijE  + a.xyijE*b.yzjE  
	 - a.yz*b.xyi    - a.yzi*b.xy     - a.yzj*b.xyij   + a.yzij*b.xyj    - a.yzE*b.xyiE   + a.yziE*b.xyE    + a.yzjE*b.xyijE  - a.yzijE*b.xyjE  
	 + a.xz*b.i      + a.xzi*b.q      + a.xzj*b.ij     - a.xzij*b.j      + a.xzE*b.iE     - a.xziE*b.E      - a.xzjE*b.ijE    + a.xzijE*b.jE     ;
c.xzj    =  + a.q*b.xzj    - a.i*b.xzij   + a.j*b.xz     + a.ij*b.xzi    + a.E*b.xzjE   + a.iE*b.xzijE  - a.jE*b.xzE    - a.ijE*b.xziE  
	 + a.xy*b.yzj    - a.xyi*b.yzij   + a.xyj*b.yz     + a.xyij*b.yzi    + a.xyE*b.yzjE   + a.xyiE*b.yzijE  - a.xyjE*b.yzE    - a.xyijE*b.yziE  
	 - a.yz*b.xyj    + a.yzi*b.xyij   - a.yzj*b.xy     - a.yzij*b.xyi    - a.yzE*b.xyjE   - a.yziE*b.xyijE  + a.yzjE*b.xyE    + a.yzijE*b.xyiE  
	 + a.xz*b.j      - a.xzi*b.ij     + a.xzj*b.q      + a.xzij*b.i      + a.xzE*b.jE     + a.xziE*b.ijE    - a.xzjE*b.E      - a.xzijE*b.iE     ;
c.xzij   =  + a.q*b.xzij   + a.i*b.xzj    - a.j*b.xzi    + a.ij*b.xz     + a.E*b.xzijE  - a.iE*b.xzjE   + a.jE*b.xziE   - a.ijE*b.xzE   
	 + a.xy*b.yzij   + a.xyi*b.yzj    - a.xyj*b.yzi    + a.xyij*b.yz     + a.xyE*b.yzijE  - a.xyiE*b.yzjE   + a.xyjE*b.yziE   - a.xyijE*b.yzE   
	 - a.yz*b.xyij   - a.yzi*b.xyj    + a.yzj*b.xyi    - a.yzij*b.xy     - a.yzE*b.xyijE  + a.yziE*b.xyjE   - a.yzjE*b.xyiE   + a.yzijE*b.xyE   
	 + a.xz*b.ij     + a.xzi*b.j      - a.xzj*b.i      + a.xzij*b.q      + a.xzE*b.ijE    - a.xziE*b.jE     + a.xzjE*b.iE     - a.xzijE*b.E      ;
c.xzE    =  + a.q*b.xzE    - a.i*b.xziE   - a.j*b.xzjE   - a.ij*b.xzijE  + a.E*b.xz     + a.iE*b.xzi    + a.jE*b.xzj    + a.ijE*b.xzij  
	 + a.xy*b.yzE    - a.xyi*b.yziE   - a.xyj*b.yzjE   - a.xyij*b.yzijE  + a.xyE*b.yz     + a.xyiE*b.yzi    + a.xyjE*b.yzj    + a.xyijE*b.yzij  
	 - a.yz*b.xyE    + a.yzi*b.xyiE   + a.yzj*b.xyjE   + a.yzij*b.xyijE  - a.yzE*b.xy     - a.yziE*b.xyi    - a.yzjE*b.xyj    - a.yzijE*b.xyij  
	 + a.xz*b.E      - a.xzi*b.iE     - a.xzj*b.jE     - a.xzij*b.ijE    + a.xzE*b.q      + a.xziE*b.i      + a.xzjE*b.j      + a.xzijE*b.ij     ;
c.xziE   =  + a.q*b.xziE   + a.i*b.xzE    - a.j*b.xzijE  + a.ij*b.xzjE   - a.E*b.xzi    + a.iE*b.xz     - a.jE*b.xzij   + a.ijE*b.xzj   
	 + a.xy*b.yziE   + a.xyi*b.yzE    - a.xyj*b.yzijE  + a.xyij*b.yzjE   - a.xyE*b.yzi    + a.xyiE*b.yz     - a.xyjE*b.yzij   + a.xyijE*b.yzj   
	 - a.yz*b.xyiE   - a.yzi*b.xyE    + a.yzj*b.xyijE  - a.yzij*b.xyjE   + a.yzE*b.xyi    - a.yziE*b.xy     + a.yzjE*b.xyij   - a.yzijE*b.xyj   
	 + a.xz*b.iE     + a.xzi*b.E      - a.xzj*b.ijE    + a.xzij*b.jE     - a.xzE*b.i      + a.xziE*b.q      - a.xzjE*b.ij     + a.xzijE*b.j      ;
c.xzjE   =  + a.q*b.xzjE   + a.i*b.xzijE  + a.j*b.xzE    - a.ij*b.xziE   - a.E*b.xzj    + a.iE*b.xzij   + a.jE*b.xz     - a.ijE*b.xzi   
	 + a.xy*b.yzjE   + a.xyi*b.yzijE  + a.xyj*b.yzE    - a.xyij*b.yziE   - a.xyE*b.yzj    + a.xyiE*b.yzij   + a.xyjE*b.yz     - a.xyijE*b.yzi   
	 - a.yz*b.xyjE   - a.yzi*b.xyijE  - a.yzj*b.xyE    + a.yzij*b.xyiE   + a.yzE*b.xyj    - a.yziE*b.xyij   - a.yzjE*b.xy     + a.yzijE*b.xyi   
	 + a.xz*b.jE     + a.xzi*b.ijE    + a.xzj*b.E      - a.xzij*b.iE     - a.xzE*b.j      + a.xziE*b.ij     + a.xzjE*b.q      - a.xzijE*b.i      ;
c.xzijE  =  + a.q*b.xzijE  - a.i*b.xzjE   + a.j*b.xziE   + a.ij*b.xzE    - a.E*b.xzij   - a.iE*b.xzj    + a.jE*b.xzi    + a.ijE*b.xz    
	 + a.xy*b.yzijE  - a.xyi*b.yzjE   + a.xyj*b.yziE   + a.xyij*b.yzE    - a.xyE*b.yzij   - a.xyiE*b.yzj    + a.xyjE*b.yzi    + a.xyijE*b.yz    
	 - a.yz*b.xyijE  + a.yzi*b.xyjE   - a.yzj*b.xyiE   - a.yzij*b.xyE    + a.yzE*b.xyij   + a.yziE*b.xyj    - a.yzjE*b.xyi    - a.yzijE*b.xy    
	 + a.xz*b.ijE    - a.xzi*b.jE     + a.xzj*b.iE     + a.xzij*b.E      - a.xzE*b.ij     - a.xziE*b.j      + a.xzjE*b.i      + a.xzijE*b.q      ;

	return c;
}

//////////////////////////////////////////////////////

CHO CHO_Product(CHO a, CHO b) {  // I am now using the GA3EO product

	CHO c;

c.q =  + a.q*b.q - a.i*b.i - a.j*b.j - a.ij*b.ij - a.E*b.E - a.iE*b.iE - a.jE*b.jE - a.ijE*b.ijE
	 + a.x*b.x - a.xi*b.xi - a.xj*b.xj - a.xij*b.xij - a.xE*b.xE - a.xiE*b.xiE - a.xjE*b.xjE - a.xijE*b.xijE
	 + a.y*b.y - a.yi*b.yi - a.yj*b.yj - a.yij*b.yij - a.yE*b.yE - a.yiE*b.yiE - a.yjE*b.yjE - a.yijE*b.yijE
	 - a.xy*b.xy + a.xyi*b.xyi + a.xyj*b.xyj + a.xyij*b.xyij + a.xyE*b.xyE + a.xyiE*b.xyiE + a.xyjE*b.xyjE + a.xyijE*b.xyijE
	 + a.z*b.z - a.zi*b.zi - a.zj*b.zj - a.zij*b.zij - a.zE*b.zE - a.ziE*b.ziE - a.zjE*b.zjE - a.zijE*b.zijE
	 - a.xz*b.xz + a.xzi*b.xzi + a.xzj*b.xzj + a.xzij*b.xzij + a.xzE*b.xzE + a.xziE*b.xziE + a.xzjE*b.xzjE + a.xzijE*b.xzijE
	 - a.yz*b.yz + a.yzi*b.yzi + a.yzj*b.yzj + a.yzij*b.yzij + a.yzE*b.yzE + a.yziE*b.yziE + a.yzjE*b.yzjE + a.yzijE*b.yzijE
	 - a.xyz*b.xyz + a.xyzi*b.xyzi + a.xyzj*b.xyzj + a.xyzij*b.xyzij + a.xyzE*b.xyzE + a.xyziE*b.xyziE + a.xyzjE*b.xyzjE + a.xyzijE*b.xyzijE ;

c.i =  + a.q*b.i + a.i*b.q + a.j*b.ij - a.ij*b.j + a.E*b.iE - a.iE*b.E - a.jE*b.ijE + a.ijE*b.jE
	 + a.x*b.xi + a.xi*b.x + a.xj*b.xij - a.xij*b.xj + a.xE*b.xiE - a.xiE*b.xE - a.xjE*b.xijE + a.xijE*b.xjE
	 + a.y*b.yi + a.yi*b.y + a.yj*b.yij - a.yij*b.yj + a.yE*b.yiE - a.yiE*b.yE - a.yjE*b.yijE + a.yijE*b.yjE
	 - a.xy*b.xyi - a.xyi*b.xy - a.xyj*b.xyij + a.xyij*b.xyj - a.xyE*b.xyiE + a.xyiE*b.xyE + a.xyjE*b.xyijE - a.xyijE*b.xyjE
	 + a.z*b.zi + a.zi*b.z + a.zj*b.zij - a.zij*b.zj + a.zE*b.ziE - a.ziE*b.zE - a.zjE*b.zijE + a.zijE*b.zjE
	 - a.xz*b.xzi - a.xzi*b.xz - a.xzj*b.xzij + a.xzij*b.xzj - a.xzE*b.xziE + a.xziE*b.xzE + a.xzjE*b.xzijE - a.xzijE*b.xzjE
	 - a.yz*b.yzi - a.yzi*b.yz - a.yzj*b.yzij + a.yzij*b.yzj - a.yzE*b.yziE + a.yziE*b.yzE + a.yzjE*b.yzijE - a.yzijE*b.yzjE
	 - a.xyz*b.xyzi - a.xyzi*b.xyz - a.xyzj*b.xyzij + a.xyzij*b.xyzj - a.xyzE*b.xyziE + a.xyziE*b.xyzE + a.xyzjE*b.xyzijE - a.xyzijE*b.xyzjE ;

c.j =  + a.q*b.j - a.i*b.ij + a.j*b.q + a.ij*b.i + a.E*b.jE + a.iE*b.ijE - a.jE*b.E - a.ijE*b.iE
	 + a.x*b.xj - a.xi*b.xij + a.xj*b.x + a.xij*b.xi + a.xE*b.xjE + a.xiE*b.xijE - a.xjE*b.xE - a.xijE*b.xiE
	 + a.y*b.yj - a.yi*b.yij + a.yj*b.y + a.yij*b.yi + a.yE*b.yjE + a.yiE*b.yijE - a.yjE*b.yE - a.yijE*b.yiE
	 - a.xy*b.xyj + a.xyi*b.xyij - a.xyj*b.xy - a.xyij*b.xyi - a.xyE*b.xyjE - a.xyiE*b.xyijE + a.xyjE*b.xyE + a.xyijE*b.xyiE
	 + a.z*b.zj - a.zi*b.zij + a.zj*b.z + a.zij*b.zi + a.zE*b.zjE + a.ziE*b.zijE - a.zjE*b.zE - a.zijE*b.ziE
	 - a.xz*b.xzj + a.xzi*b.xzij - a.xzj*b.xz - a.xzij*b.xzi - a.xzE*b.xzjE - a.xziE*b.xzijE + a.xzjE*b.xzE + a.xzijE*b.xziE
	 - a.yz*b.yzj + a.yzi*b.yzij - a.yzj*b.yz - a.yzij*b.yzi - a.yzE*b.yzjE - a.yziE*b.yzijE + a.yzjE*b.yzE + a.yzijE*b.yziE
	 - a.xyz*b.xyzj + a.xyzi*b.xyzij - a.xyzj*b.xyz - a.xyzij*b.xyzi - a.xyzE*b.xyzjE - a.xyziE*b.xyzijE + a.xyzjE*b.xyzE + a.xyzijE*b.xyziE ;

c.ij =  + a.q*b.ij + a.i*b.j - a.j*b.i + a.ij*b.q + a.E*b.ijE - a.iE*b.jE + a.jE*b.iE - a.ijE*b.E
	 + a.x*b.xij + a.xi*b.xj - a.xj*b.xi + a.xij*b.x + a.xE*b.xijE - a.xiE*b.xjE + a.xjE*b.xiE - a.xijE*b.xE
	 + a.y*b.yij + a.yi*b.yj - a.yj*b.yi + a.yij*b.y + a.yE*b.yijE - a.yiE*b.yjE + a.yjE*b.yiE - a.yijE*b.yE
	 - a.xy*b.xyij - a.xyi*b.xyj + a.xyj*b.xyi - a.xyij*b.xy - a.xyE*b.xyijE + a.xyiE*b.xyjE - a.xyjE*b.xyiE + a.xyijE*b.xyE
	 + a.z*b.zij + a.zi*b.zj - a.zj*b.zi + a.zij*b.z + a.zE*b.zijE - a.ziE*b.zjE + a.zjE*b.ziE - a.zijE*b.zE
	 - a.xz*b.xzij - a.xzi*b.xzj + a.xzj*b.xzi - a.xzij*b.xz - a.xzE*b.xzijE + a.xziE*b.xzjE - a.xzjE*b.xziE + a.xzijE*b.xzE
	 - a.yz*b.yzij - a.yzi*b.yzj + a.yzj*b.yzi - a.yzij*b.yz - a.yzE*b.yzijE + a.yziE*b.yzjE - a.yzjE*b.yziE + a.yzijE*b.yzE
	 - a.xyz*b.xyzij - a.xyzi*b.xyzj + a.xyzj*b.xyzi - a.xyzij*b.xyz - a.xyzE*b.xyzijE + a.xyziE*b.xyzjE - a.xyzjE*b.xyziE + a.xyzijE*b.xyzE ;

c.E =  + a.q*b.E - a.i*b.iE - a.j*b.jE - a.ij*b.ijE + a.E*b.q + a.iE*b.i + a.jE*b.j + a.ijE*b.ij
	 + a.x*b.xE - a.xi*b.xiE - a.xj*b.xjE - a.xij*b.xijE + a.xE*b.x + a.xiE*b.xi + a.xjE*b.xj + a.xijE*b.xij
	 + a.y*b.yE - a.yi*b.yiE - a.yj*b.yjE - a.yij*b.yijE + a.yE*b.y + a.yiE*b.yi + a.yjE*b.yj + a.yijE*b.yij
	 - a.xy*b.xyE + a.xyi*b.xyiE + a.xyj*b.xyjE + a.xyij*b.xyijE - a.xyE*b.xy - a.xyiE*b.xyi - a.xyjE*b.xyj - a.xyijE*b.xyij
	 + a.z*b.zE - a.zi*b.ziE - a.zj*b.zjE - a.zij*b.zijE + a.zE*b.z + a.ziE*b.zi + a.zjE*b.zj + a.zijE*b.zij
	 - a.xz*b.xzE + a.xzi*b.xziE + a.xzj*b.xzjE + a.xzij*b.xzijE - a.xzE*b.xz - a.xziE*b.xzi - a.xzjE*b.xzj - a.xzijE*b.xzij
	 - a.yz*b.yzE + a.yzi*b.yziE + a.yzj*b.yzjE + a.yzij*b.yzijE - a.yzE*b.yz - a.yziE*b.yzi - a.yzjE*b.yzj - a.yzijE*b.yzij
	 - a.xyz*b.xyzE + a.xyzi*b.xyziE + a.xyzj*b.xyzjE + a.xyzij*b.xyzijE - a.xyzE*b.xyz - a.xyziE*b.xyzi - a.xyzjE*b.xyzj - a.xyzijE*b.xyzij ;

c.iE =  + a.q*b.iE + a.i*b.E - a.j*b.ijE + a.ij*b.jE - a.E*b.i + a.iE*b.q - a.jE*b.ij + a.ijE*b.j
	 + a.x*b.xiE + a.xi*b.xE - a.xj*b.xijE + a.xij*b.xjE - a.xE*b.xi + a.xiE*b.x - a.xjE*b.xij + a.xijE*b.xj
	 + a.y*b.yiE + a.yi*b.yE - a.yj*b.yijE + a.yij*b.yjE - a.yE*b.yi + a.yiE*b.y - a.yjE*b.yij + a.yijE*b.yj
	 - a.xy*b.xyiE - a.xyi*b.xyE + a.xyj*b.xyijE - a.xyij*b.xyjE + a.xyE*b.xyi - a.xyiE*b.xy + a.xyjE*b.xyij - a.xyijE*b.xyj
	 + a.z*b.ziE + a.zi*b.zE - a.zj*b.zijE + a.zij*b.zjE - a.zE*b.zi + a.ziE*b.z - a.zjE*b.zij + a.zijE*b.zj
	 - a.xz*b.xziE - a.xzi*b.xzE + a.xzj*b.xzijE - a.xzij*b.xzjE + a.xzE*b.xzi - a.xziE*b.xz + a.xzjE*b.xzij - a.xzijE*b.xzj
	 - a.yz*b.yziE - a.yzi*b.yzE + a.yzj*b.yzijE - a.yzij*b.yzjE + a.yzE*b.yzi - a.yziE*b.yz + a.yzjE*b.yzij - a.yzijE*b.yzj
	 - a.xyz*b.xyziE - a.xyzi*b.xyzE + a.xyzj*b.xyzijE - a.xyzij*b.xyzjE + a.xyzE*b.xyzi - a.xyziE*b.xyz + a.xyzjE*b.xyzij - a.xyzijE*b.xyzj ;

c.jE =  + a.q*b.jE + a.i*b.ijE + a.j*b.E - a.ij*b.iE - a.E*b.j + a.iE*b.ij + a.jE*b.q - a.ijE*b.i
	 + a.x*b.xjE + a.xi*b.xijE + a.xj*b.xE - a.xij*b.xiE - a.xE*b.xj + a.xiE*b.xij + a.xjE*b.x - a.xijE*b.xi
	 + a.y*b.yjE + a.yi*b.yijE + a.yj*b.yE - a.yij*b.yiE - a.yE*b.yj + a.yiE*b.yij + a.yjE*b.y - a.yijE*b.yi
	 - a.xy*b.xyjE - a.xyi*b.xyijE - a.xyj*b.xyE + a.xyij*b.xyiE + a.xyE*b.xyj - a.xyiE*b.xyij - a.xyjE*b.xy + a.xyijE*b.xyi
	 + a.z*b.zjE + a.zi*b.zijE + a.zj*b.zE - a.zij*b.ziE - a.zE*b.zj + a.ziE*b.zij + a.zjE*b.z - a.zijE*b.zi
	 - a.xz*b.xzjE - a.xzi*b.xzijE - a.xzj*b.xzE + a.xzij*b.xziE + a.xzE*b.xzj - a.xziE*b.xzij - a.xzjE*b.xz + a.xzijE*b.xzi
	 - a.yz*b.yzjE - a.yzi*b.yzijE - a.yzj*b.yzE + a.yzij*b.yziE + a.yzE*b.yzj - a.yziE*b.yzij - a.yzjE*b.yz + a.yzijE*b.yzi
	 - a.xyz*b.xyzjE - a.xyzi*b.xyzijE - a.xyzj*b.xyzE + a.xyzij*b.xyziE + a.xyzE*b.xyzj - a.xyziE*b.xyzij - a.xyzjE*b.xyz + a.xyzijE*b.xyzi ;

c.ijE =  + a.q*b.ijE - a.i*b.jE + a.j*b.iE + a.ij*b.E - a.E*b.ij - a.iE*b.j + a.jE*b.i + a.ijE*b.q
	 + a.x*b.xijE - a.xi*b.xjE + a.xj*b.xiE + a.xij*b.xE - a.xE*b.xij - a.xiE*b.xj + a.xjE*b.xi + a.xijE*b.x
	 + a.y*b.yijE - a.yi*b.yjE + a.yj*b.yiE + a.yij*b.yE - a.yE*b.yij - a.yiE*b.yj + a.yjE*b.yi + a.yijE*b.y
	 - a.xy*b.xyijE + a.xyi*b.xyjE - a.xyj*b.xyiE - a.xyij*b.xyE + a.xyE*b.xyij + a.xyiE*b.xyj - a.xyjE*b.xyi - a.xyijE*b.xy
	 + a.z*b.zijE - a.zi*b.zjE + a.zj*b.ziE + a.zij*b.zE - a.zE*b.zij - a.ziE*b.zj + a.zjE*b.zi + a.zijE*b.z
	 - a.xz*b.xzijE + a.xzi*b.xzjE - a.xzj*b.xziE - a.xzij*b.xzE + a.xzE*b.xzij + a.xziE*b.xzj - a.xzjE*b.xzi - a.xzijE*b.xz
	 - a.yz*b.yzijE + a.yzi*b.yzjE - a.yzj*b.yziE - a.yzij*b.yzE + a.yzE*b.yzij + a.yziE*b.yzj - a.yzjE*b.yzi - a.yzijE*b.yz
	 - a.xyz*b.xyzijE + a.xyzi*b.xyzjE - a.xyzj*b.xyziE - a.xyzij*b.xyzE + a.xyzE*b.xyzij + a.xyziE*b.xyzj - a.xyzjE*b.xyzi - a.xyzijE*b.xyz ;

c.x =  + a.q*b.x - a.i*b.xi - a.j*b.xj - a.ij*b.xij - a.E*b.xE - a.iE*b.xiE - a.jE*b.xjE - a.ijE*b.xijE
	 + a.x*b.q - a.xi*b.i - a.xj*b.j - a.xij*b.ij - a.xE*b.E - a.xiE*b.iE - a.xjE*b.jE - a.xijE*b.ijE
	 - a.y*b.xy + a.yi*b.xyi + a.yj*b.xyj + a.yij*b.xyij + a.yE*b.xyE + a.yiE*b.xyiE + a.yjE*b.xyjE + a.yijE*b.xyijE
	 + a.xy*b.y - a.xyi*b.yi - a.xyj*b.yj - a.xyij*b.yij - a.xyE*b.yE - a.xyiE*b.yiE - a.xyjE*b.yjE - a.xyijE*b.yijE
	 - a.z*b.xz + a.zi*b.xzi + a.zj*b.xzj + a.zij*b.xzij + a.zE*b.xzE + a.ziE*b.xziE + a.zjE*b.xzjE + a.zijE*b.xzijE
	 + a.xz*b.z - a.xzi*b.zi - a.xzj*b.zj - a.xzij*b.zij - a.xzE*b.zE - a.xziE*b.ziE - a.xzjE*b.zjE - a.xzijE*b.zijE
	 - a.yz*b.xyz + a.yzi*b.xyzi + a.yzj*b.xyzj + a.yzij*b.xyzij + a.yzE*b.xyzE + a.yziE*b.xyziE + a.yzjE*b.xyzjE + a.yzijE*b.xyzijE
	 - a.xyz*b.yz + a.xyzi*b.yzi + a.xyzj*b.yzj + a.xyzij*b.yzij + a.xyzE*b.yzE + a.xyziE*b.yziE + a.xyzjE*b.yzjE + a.xyzijE*b.yzijE ;

c.xi =  + a.q*b.xi + a.i*b.x + a.j*b.xij - a.ij*b.xj + a.E*b.xiE - a.iE*b.xE - a.jE*b.xijE + a.ijE*b.xjE
	 + a.x*b.i + a.xi*b.q + a.xj*b.ij - a.xij*b.j + a.xE*b.iE - a.xiE*b.E - a.xjE*b.ijE + a.xijE*b.jE
	 - a.y*b.xyi - a.yi*b.xy - a.yj*b.xyij + a.yij*b.xyj - a.yE*b.xyiE + a.yiE*b.xyE + a.yjE*b.xyijE - a.yijE*b.xyjE
	 + a.xy*b.yi + a.xyi*b.y + a.xyj*b.yij - a.xyij*b.yj + a.xyE*b.yiE - a.xyiE*b.yE - a.xyjE*b.yijE + a.xyijE*b.yjE
	 - a.z*b.xzi - a.zi*b.xz - a.zj*b.xzij + a.zij*b.xzj - a.zE*b.xziE + a.ziE*b.xzE + a.zjE*b.xzijE - a.zijE*b.xzjE
	 + a.xz*b.zi + a.xzi*b.z + a.xzj*b.zij - a.xzij*b.zj + a.xzE*b.ziE - a.xziE*b.zE - a.xzjE*b.zijE + a.xzijE*b.zjE
	 - a.yz*b.xyzi - a.yzi*b.xyz - a.yzj*b.xyzij + a.yzij*b.xyzj - a.yzE*b.xyziE + a.yziE*b.xyzE + a.yzjE*b.xyzijE - a.yzijE*b.xyzjE
	 - a.xyz*b.yzi - a.xyzi*b.yz - a.xyzj*b.yzij + a.xyzij*b.yzj - a.xyzE*b.yziE + a.xyziE*b.yzE + a.xyzjE*b.yzijE - a.xyzijE*b.yzjE ;

c.xj =  + a.q*b.xj - a.i*b.xij + a.j*b.x + a.ij*b.xi + a.E*b.xjE + a.iE*b.xijE - a.jE*b.xE - a.ijE*b.xiE
	 + a.x*b.j - a.xi*b.ij + a.xj*b.q + a.xij*b.i + a.xE*b.jE + a.xiE*b.ijE - a.xjE*b.E - a.xijE*b.iE
	 - a.y*b.xyj + a.yi*b.xyij - a.yj*b.xy - a.yij*b.xyi - a.yE*b.xyjE - a.yiE*b.xyijE + a.yjE*b.xyE + a.yijE*b.xyiE
	 + a.xy*b.yj - a.xyi*b.yij + a.xyj*b.y + a.xyij*b.yi + a.xyE*b.yjE + a.xyiE*b.yijE - a.xyjE*b.yE - a.xyijE*b.yiE
	 - a.z*b.xzj + a.zi*b.xzij - a.zj*b.xz - a.zij*b.xzi - a.zE*b.xzjE - a.ziE*b.xzijE + a.zjE*b.xzE + a.zijE*b.xziE
	 + a.xz*b.zj - a.xzi*b.zij + a.xzj*b.z + a.xzij*b.zi + a.xzE*b.zjE + a.xziE*b.zijE - a.xzjE*b.zE - a.xzijE*b.ziE
	 - a.yz*b.xyzj + a.yzi*b.xyzij - a.yzj*b.xyz - a.yzij*b.xyzi - a.yzE*b.xyzjE - a.yziE*b.xyzijE + a.yzjE*b.xyzE + a.yzijE*b.xyziE
	 - a.xyz*b.yzj + a.xyzi*b.yzij - a.xyzj*b.yz - a.xyzij*b.yzi - a.xyzE*b.yzjE - a.xyziE*b.yzijE + a.xyzjE*b.yzE + a.xyzijE*b.yziE ;

c.xij =  + a.q*b.xij + a.i*b.xj - a.j*b.xi + a.ij*b.x + a.E*b.xijE - a.iE*b.xjE + a.jE*b.xiE - a.ijE*b.xE
	 + a.x*b.ij + a.xi*b.j - a.xj*b.i + a.xij*b.q + a.xE*b.ijE - a.xiE*b.jE + a.xjE*b.iE - a.xijE*b.E
	 - a.y*b.xyij - a.yi*b.xyj + a.yj*b.xyi - a.yij*b.xy - a.yE*b.xyijE + a.yiE*b.xyjE - a.yjE*b.xyiE + a.yijE*b.xyE
	 + a.xy*b.yij + a.xyi*b.yj - a.xyj*b.yi + a.xyij*b.y + a.xyE*b.yijE - a.xyiE*b.yjE + a.xyjE*b.yiE - a.xyijE*b.yE
	 - a.z*b.xzij - a.zi*b.xzj + a.zj*b.xzi - a.zij*b.xz - a.zE*b.xzijE + a.ziE*b.xzjE - a.zjE*b.xziE + a.zijE*b.xzE
	 + a.xz*b.zij + a.xzi*b.zj - a.xzj*b.zi + a.xzij*b.z + a.xzE*b.zijE - a.xziE*b.zjE + a.xzjE*b.ziE - a.xzijE*b.zE
	 - a.yz*b.xyzij - a.yzi*b.xyzj + a.yzj*b.xyzi - a.yzij*b.xyz - a.yzE*b.xyzijE + a.yziE*b.xyzjE - a.yzjE*b.xyziE + a.yzijE*b.xyzE
	 - a.xyz*b.yzij - a.xyzi*b.yzj + a.xyzj*b.yzi - a.xyzij*b.yz - a.xyzE*b.yzijE + a.xyziE*b.yzjE - a.xyzjE*b.yziE + a.xyzijE*b.yzE ;

c.xE =  + a.q*b.xE - a.i*b.xiE - a.j*b.xjE - a.ij*b.xijE + a.E*b.x + a.iE*b.xi + a.jE*b.xj + a.ijE*b.xij
	 + a.x*b.E - a.xi*b.iE - a.xj*b.jE - a.xij*b.ijE + a.xE*b.q + a.xiE*b.i + a.xjE*b.j + a.xijE*b.ij
	 - a.y*b.xyE + a.yi*b.xyiE + a.yj*b.xyjE + a.yij*b.xyijE - a.yE*b.xy - a.yiE*b.xyi - a.yjE*b.xyj - a.yijE*b.xyij
	 + a.xy*b.yE - a.xyi*b.yiE - a.xyj*b.yjE - a.xyij*b.yijE + a.xyE*b.y + a.xyiE*b.yi + a.xyjE*b.yj + a.xyijE*b.yij
	 - a.z*b.xzE + a.zi*b.xziE + a.zj*b.xzjE + a.zij*b.xzijE - a.zE*b.xz - a.ziE*b.xzi - a.zjE*b.xzj - a.zijE*b.xzij
	 + a.xz*b.zE - a.xzi*b.ziE - a.xzj*b.zjE - a.xzij*b.zijE + a.xzE*b.z + a.xziE*b.zi + a.xzjE*b.zj + a.xzijE*b.zij
	 - a.yz*b.xyzE + a.yzi*b.xyziE + a.yzj*b.xyzjE + a.yzij*b.xyzijE - a.yzE*b.xyz - a.yziE*b.xyzi - a.yzjE*b.xyzj - a.yzijE*b.xyzij
	 - a.xyz*b.yzE + a.xyzi*b.yziE + a.xyzj*b.yzjE + a.xyzij*b.yzijE - a.xyzE*b.yz - a.xyziE*b.yzi - a.xyzjE*b.yzj - a.xyzijE*b.yzij ;

c.xiE =  + a.q*b.xiE + a.i*b.xE - a.j*b.xijE + a.ij*b.xjE - a.E*b.xi + a.iE*b.x - a.jE*b.xij + a.ijE*b.xj
	 + a.x*b.iE + a.xi*b.E - a.xj*b.ijE + a.xij*b.jE - a.xE*b.i + a.xiE*b.q - a.xjE*b.ij + a.xijE*b.j
	 - a.y*b.xyiE - a.yi*b.xyE + a.yj*b.xyijE - a.yij*b.xyjE + a.yE*b.xyi - a.yiE*b.xy + a.yjE*b.xyij - a.yijE*b.xyj
	 + a.xy*b.yiE + a.xyi*b.yE - a.xyj*b.yijE + a.xyij*b.yjE - a.xyE*b.yi + a.xyiE*b.y - a.xyjE*b.yij + a.xyijE*b.yj
	 - a.z*b.xziE - a.zi*b.xzE + a.zj*b.xzijE - a.zij*b.xzjE + a.zE*b.xzi - a.ziE*b.xz + a.zjE*b.xzij - a.zijE*b.xzj
	 + a.xz*b.ziE + a.xzi*b.zE - a.xzj*b.zijE + a.xzij*b.zjE - a.xzE*b.zi + a.xziE*b.z - a.xzjE*b.zij + a.xzijE*b.zj
	 - a.yz*b.xyziE - a.yzi*b.xyzE + a.yzj*b.xyzijE - a.yzij*b.xyzjE + a.yzE*b.xyzi - a.yziE*b.xyz + a.yzjE*b.xyzij - a.yzijE*b.xyzj
	 - a.xyz*b.yziE - a.xyzi*b.yzE + a.xyzj*b.yzijE - a.xyzij*b.yzjE + a.xyzE*b.yzi - a.xyziE*b.yz + a.xyzjE*b.yzij - a.xyzijE*b.yzj ;

c.xjE =  + a.q*b.xjE + a.i*b.xijE + a.j*b.xE - a.ij*b.xiE - a.E*b.xj + a.iE*b.xij + a.jE*b.x - a.ijE*b.xi
	 + a.x*b.jE + a.xi*b.ijE + a.xj*b.E - a.xij*b.iE - a.xE*b.j + a.xiE*b.ij + a.xjE*b.q - a.xijE*b.i
	 - a.y*b.xyjE - a.yi*b.xyijE - a.yj*b.xyE + a.yij*b.xyiE + a.yE*b.xyj - a.yiE*b.xyij - a.yjE*b.xy + a.yijE*b.xyi
	 + a.xy*b.yjE + a.xyi*b.yijE + a.xyj*b.yE - a.xyij*b.yiE - a.xyE*b.yj + a.xyiE*b.yij + a.xyjE*b.y - a.xyijE*b.yi
	 - a.z*b.xzjE - a.zi*b.xzijE - a.zj*b.xzE + a.zij*b.xziE + a.zE*b.xzj - a.ziE*b.xzij - a.zjE*b.xz + a.zijE*b.xzi
	 + a.xz*b.zjE + a.xzi*b.zijE + a.xzj*b.zE - a.xzij*b.ziE - a.xzE*b.zj + a.xziE*b.zij + a.xzjE*b.z - a.xzijE*b.zi
	 - a.yz*b.xyzjE - a.yzi*b.xyzijE - a.yzj*b.xyzE + a.yzij*b.xyziE + a.yzE*b.xyzj - a.yziE*b.xyzij - a.yzjE*b.xyz + a.yzijE*b.xyzi
	 - a.xyz*b.yzjE - a.xyzi*b.yzijE - a.xyzj*b.yzE + a.xyzij*b.yziE + a.xyzE*b.yzj - a.xyziE*b.yzij - a.xyzjE*b.yz + a.xyzijE*b.yzi ;

c.xijE =  + a.q*b.xijE - a.i*b.xjE + a.j*b.xiE + a.ij*b.xE - a.E*b.xij - a.iE*b.xj + a.jE*b.xi + a.ijE*b.x
	 + a.x*b.ijE - a.xi*b.jE + a.xj*b.iE + a.xij*b.E - a.xE*b.ij - a.xiE*b.j + a.xjE*b.i + a.xijE*b.q
	 - a.y*b.xyijE + a.yi*b.xyjE - a.yj*b.xyiE - a.yij*b.xyE + a.yE*b.xyij + a.yiE*b.xyj - a.yjE*b.xyi - a.yijE*b.xy
	 + a.xy*b.yijE - a.xyi*b.yjE + a.xyj*b.yiE + a.xyij*b.yE - a.xyE*b.yij - a.xyiE*b.yj + a.xyjE*b.yi + a.xyijE*b.y
	 - a.z*b.xzijE + a.zi*b.xzjE - a.zj*b.xziE - a.zij*b.xzE + a.zE*b.xzij + a.ziE*b.xzj - a.zjE*b.xzi - a.zijE*b.xz
	 + a.xz*b.zijE - a.xzi*b.zjE + a.xzj*b.ziE + a.xzij*b.zE - a.xzE*b.zij - a.xziE*b.zj + a.xzjE*b.zi + a.xzijE*b.z
	 - a.yz*b.xyzijE + a.yzi*b.xyzjE - a.yzj*b.xyziE - a.yzij*b.xyzE + a.yzE*b.xyzij + a.yziE*b.xyzj - a.yzjE*b.xyzi - a.yzijE*b.xyz
	 - a.xyz*b.yzijE + a.xyzi*b.yzjE - a.xyzj*b.yziE - a.xyzij*b.yzE + a.xyzE*b.yzij + a.xyziE*b.yzj - a.xyzjE*b.yzi - a.xyzijE*b.yz ;

c.y =  + a.q*b.y - a.i*b.yi - a.j*b.yj - a.ij*b.yij - a.E*b.yE - a.iE*b.yiE - a.jE*b.yjE - a.ijE*b.yijE
	 + a.x*b.xy - a.xi*b.xyi - a.xj*b.xyj - a.xij*b.xyij - a.xE*b.xyE - a.xiE*b.xyiE - a.xjE*b.xyjE - a.xijE*b.xyijE
	 + a.y*b.q - a.yi*b.i - a.yj*b.j - a.yij*b.ij - a.yE*b.E - a.yiE*b.iE - a.yjE*b.jE - a.yijE*b.ijE
	 - a.xy*b.x + a.xyi*b.xi + a.xyj*b.xj + a.xyij*b.xij + a.xyE*b.xE + a.xyiE*b.xiE + a.xyjE*b.xjE + a.xyijE*b.xijE
	 - a.z*b.yz + a.zi*b.yzi + a.zj*b.yzj + a.zij*b.yzij + a.zE*b.yzE + a.ziE*b.yziE + a.zjE*b.yzjE + a.zijE*b.yzijE
	 + a.xz*b.xyz - a.xzi*b.xyzi - a.xzj*b.xyzj - a.xzij*b.xyzij - a.xzE*b.xyzE - a.xziE*b.xyziE - a.xzjE*b.xyzjE - a.xzijE*b.xyzijE
	 + a.yz*b.z - a.yzi*b.zi - a.yzj*b.zj - a.yzij*b.zij - a.yzE*b.zE - a.yziE*b.ziE - a.yzjE*b.zjE - a.yzijE*b.zijE
	 + a.xyz*b.xz - a.xyzi*b.xzi - a.xyzj*b.xzj - a.xyzij*b.xzij - a.xyzE*b.xzE - a.xyziE*b.xziE - a.xyzjE*b.xzjE - a.xyzijE*b.xzijE ;

c.yi =  + a.q*b.yi + a.i*b.y + a.j*b.yij - a.ij*b.yj + a.E*b.yiE - a.iE*b.yE - a.jE*b.yijE + a.ijE*b.yjE
	 + a.x*b.xyi + a.xi*b.xy + a.xj*b.xyij - a.xij*b.xyj + a.xE*b.xyiE - a.xiE*b.xyE - a.xjE*b.xyijE + a.xijE*b.xyjE
	 + a.y*b.i + a.yi*b.q + a.yj*b.ij - a.yij*b.j + a.yE*b.iE - a.yiE*b.E - a.yjE*b.ijE + a.yijE*b.jE
	 - a.xy*b.xi - a.xyi*b.x - a.xyj*b.xij + a.xyij*b.xj - a.xyE*b.xiE + a.xyiE*b.xE + a.xyjE*b.xijE - a.xyijE*b.xjE
	 - a.z*b.yzi - a.zi*b.yz - a.zj*b.yzij + a.zij*b.yzj - a.zE*b.yziE + a.ziE*b.yzE + a.zjE*b.yzijE - a.zijE*b.yzjE
	 + a.xz*b.xyzi + a.xzi*b.xyz + a.xzj*b.xyzij - a.xzij*b.xyzj + a.xzE*b.xyziE - a.xziE*b.xyzE - a.xzjE*b.xyzijE + a.xzijE*b.xyzjE
	 + a.yz*b.zi + a.yzi*b.z + a.yzj*b.zij - a.yzij*b.zj + a.yzE*b.ziE - a.yziE*b.zE - a.yzjE*b.zijE + a.yzijE*b.zjE
	 + a.xyz*b.xzi + a.xyzi*b.xz + a.xyzj*b.xzij - a.xyzij*b.xzj + a.xyzE*b.xziE - a.xyziE*b.xzE - a.xyzjE*b.xzijE + a.xyzijE*b.xzjE ;

c.yj =  + a.q*b.yj - a.i*b.yij + a.j*b.y + a.ij*b.yi + a.E*b.yjE + a.iE*b.yijE - a.jE*b.yE - a.ijE*b.yiE
	 + a.x*b.xyj - a.xi*b.xyij + a.xj*b.xy + a.xij*b.xyi + a.xE*b.xyjE + a.xiE*b.xyijE - a.xjE*b.xyE - a.xijE*b.xyiE
	 + a.y*b.j - a.yi*b.ij + a.yj*b.q + a.yij*b.i + a.yE*b.jE + a.yiE*b.ijE - a.yjE*b.E - a.yijE*b.iE
	 - a.xy*b.xj + a.xyi*b.xij - a.xyj*b.x - a.xyij*b.xi - a.xyE*b.xjE - a.xyiE*b.xijE + a.xyjE*b.xE + a.xyijE*b.xiE
	 - a.z*b.yzj + a.zi*b.yzij - a.zj*b.yz - a.zij*b.yzi - a.zE*b.yzjE - a.ziE*b.yzijE + a.zjE*b.yzE + a.zijE*b.yziE
	 + a.xz*b.xyzj - a.xzi*b.xyzij + a.xzj*b.xyz + a.xzij*b.xyzi + a.xzE*b.xyzjE + a.xziE*b.xyzijE - a.xzjE*b.xyzE - a.xzijE*b.xyziE
	 + a.yz*b.zj - a.yzi*b.zij + a.yzj*b.z + a.yzij*b.zi + a.yzE*b.zjE + a.yziE*b.zijE - a.yzjE*b.zE - a.yzijE*b.ziE
	 + a.xyz*b.xzj - a.xyzi*b.xzij + a.xyzj*b.xz + a.xyzij*b.xzi + a.xyzE*b.xzjE + a.xyziE*b.xzijE - a.xyzjE*b.xzE - a.xyzijE*b.xziE ;

c.yij =  + a.q*b.yij + a.i*b.yj - a.j*b.yi + a.ij*b.y + a.E*b.yijE - a.iE*b.yjE + a.jE*b.yiE - a.ijE*b.yE
	 + a.x*b.xyij + a.xi*b.xyj - a.xj*b.xyi + a.xij*b.xy + a.xE*b.xyijE - a.xiE*b.xyjE + a.xjE*b.xyiE - a.xijE*b.xyE
	 + a.y*b.ij + a.yi*b.j - a.yj*b.i + a.yij*b.q + a.yE*b.ijE - a.yiE*b.jE + a.yjE*b.iE - a.yijE*b.E
	 - a.xy*b.xij - a.xyi*b.xj + a.xyj*b.xi - a.xyij*b.x - a.xyE*b.xijE + a.xyiE*b.xjE - a.xyjE*b.xiE + a.xyijE*b.xE
	 - a.z*b.yzij - a.zi*b.yzj + a.zj*b.yzi - a.zij*b.yz - a.zE*b.yzijE + a.ziE*b.yzjE - a.zjE*b.yziE + a.zijE*b.yzE
	 + a.xz*b.xyzij + a.xzi*b.xyzj - a.xzj*b.xyzi + a.xzij*b.xyz + a.xzE*b.xyzijE - a.xziE*b.xyzjE + a.xzjE*b.xyziE - a.xzijE*b.xyzE
	 + a.yz*b.zij + a.yzi*b.zj - a.yzj*b.zi + a.yzij*b.z + a.yzE*b.zijE - a.yziE*b.zjE + a.yzjE*b.ziE - a.yzijE*b.zE
	 + a.xyz*b.xzij + a.xyzi*b.xzj - a.xyzj*b.xzi + a.xyzij*b.xz + a.xyzE*b.xzijE - a.xyziE*b.xzjE + a.xyzjE*b.xziE - a.xyzijE*b.xzE ;

c.yE =  + a.q*b.yE - a.i*b.yiE - a.j*b.yjE - a.ij*b.yijE + a.E*b.y + a.iE*b.yi + a.jE*b.yj + a.ijE*b.yij
	 + a.x*b.xyE - a.xi*b.xyiE - a.xj*b.xyjE - a.xij*b.xyijE + a.xE*b.xy + a.xiE*b.xyi + a.xjE*b.xyj + a.xijE*b.xyij
	 + a.y*b.E - a.yi*b.iE - a.yj*b.jE - a.yij*b.ijE + a.yE*b.q + a.yiE*b.i + a.yjE*b.j + a.yijE*b.ij
	 - a.xy*b.xE + a.xyi*b.xiE + a.xyj*b.xjE + a.xyij*b.xijE - a.xyE*b.x - a.xyiE*b.xi - a.xyjE*b.xj - a.xyijE*b.xij
	 - a.z*b.yzE + a.zi*b.yziE + a.zj*b.yzjE + a.zij*b.yzijE - a.zE*b.yz - a.ziE*b.yzi - a.zjE*b.yzj - a.zijE*b.yzij
	 + a.xz*b.xyzE - a.xzi*b.xyziE - a.xzj*b.xyzjE - a.xzij*b.xyzijE + a.xzE*b.xyz + a.xziE*b.xyzi + a.xzjE*b.xyzj + a.xzijE*b.xyzij
	 + a.yz*b.zE - a.yzi*b.ziE - a.yzj*b.zjE - a.yzij*b.zijE + a.yzE*b.z + a.yziE*b.zi + a.yzjE*b.zj + a.yzijE*b.zij
	 + a.xyz*b.xzE - a.xyzi*b.xziE - a.xyzj*b.xzjE - a.xyzij*b.xzijE + a.xyzE*b.xz + a.xyziE*b.xzi + a.xyzjE*b.xzj + a.xyzijE*b.xzij ;

c.yiE =  + a.q*b.yiE + a.i*b.yE - a.j*b.yijE + a.ij*b.yjE - a.E*b.yi + a.iE*b.y - a.jE*b.yij + a.ijE*b.yj
	 + a.x*b.xyiE + a.xi*b.xyE - a.xj*b.xyijE + a.xij*b.xyjE - a.xE*b.xyi + a.xiE*b.xy - a.xjE*b.xyij + a.xijE*b.xyj
	 + a.y*b.iE + a.yi*b.E - a.yj*b.ijE + a.yij*b.jE - a.yE*b.i + a.yiE*b.q - a.yjE*b.ij + a.yijE*b.j
	 - a.xy*b.xiE - a.xyi*b.xE + a.xyj*b.xijE - a.xyij*b.xjE + a.xyE*b.xi - a.xyiE*b.x + a.xyjE*b.xij - a.xyijE*b.xj
	 - a.z*b.yziE - a.zi*b.yzE + a.zj*b.yzijE - a.zij*b.yzjE + a.zE*b.yzi - a.ziE*b.yz + a.zjE*b.yzij - a.zijE*b.yzj
	 + a.xz*b.xyziE + a.xzi*b.xyzE - a.xzj*b.xyzijE + a.xzij*b.xyzjE - a.xzE*b.xyzi + a.xziE*b.xyz - a.xzjE*b.xyzij + a.xzijE*b.xyzj
	 + a.yz*b.ziE + a.yzi*b.zE - a.yzj*b.zijE + a.yzij*b.zjE - a.yzE*b.zi + a.yziE*b.z - a.yzjE*b.zij + a.yzijE*b.zj
	 + a.xyz*b.xziE + a.xyzi*b.xzE - a.xyzj*b.xzijE + a.xyzij*b.xzjE - a.xyzE*b.xzi + a.xyziE*b.xz - a.xyzjE*b.xzij + a.xyzijE*b.xzj ;

c.yjE =  + a.q*b.yjE + a.i*b.yijE + a.j*b.yE - a.ij*b.yiE - a.E*b.yj + a.iE*b.yij + a.jE*b.y - a.ijE*b.yi
	 + a.x*b.xyjE + a.xi*b.xyijE + a.xj*b.xyE - a.xij*b.xyiE - a.xE*b.xyj + a.xiE*b.xyij + a.xjE*b.xy - a.xijE*b.xyi
	 + a.y*b.jE + a.yi*b.ijE + a.yj*b.E - a.yij*b.iE - a.yE*b.j + a.yiE*b.ij + a.yjE*b.q - a.yijE*b.i
	 - a.xy*b.xjE - a.xyi*b.xijE - a.xyj*b.xE + a.xyij*b.xiE + a.xyE*b.xj - a.xyiE*b.xij - a.xyjE*b.x + a.xyijE*b.xi
	 - a.z*b.yzjE - a.zi*b.yzijE - a.zj*b.yzE + a.zij*b.yziE + a.zE*b.yzj - a.ziE*b.yzij - a.zjE*b.yz + a.zijE*b.yzi
	 + a.xz*b.xyzjE + a.xzi*b.xyzijE + a.xzj*b.xyzE - a.xzij*b.xyziE - a.xzE*b.xyzj + a.xziE*b.xyzij + a.xzjE*b.xyz - a.xzijE*b.xyzi
	 + a.yz*b.zjE + a.yzi*b.zijE + a.yzj*b.zE - a.yzij*b.ziE - a.yzE*b.zj + a.yziE*b.zij + a.yzjE*b.z - a.yzijE*b.zi
	 + a.xyz*b.xzjE + a.xyzi*b.xzijE + a.xyzj*b.xzE - a.xyzij*b.xziE - a.xyzE*b.xzj + a.xyziE*b.xzij + a.xyzjE*b.xz - a.xyzijE*b.xzi ;

c.yijE =  + a.q*b.yijE - a.i*b.yjE + a.j*b.yiE + a.ij*b.yE - a.E*b.yij - a.iE*b.yj + a.jE*b.yi + a.ijE*b.y
	 + a.x*b.xyijE - a.xi*b.xyjE + a.xj*b.xyiE + a.xij*b.xyE - a.xE*b.xyij - a.xiE*b.xyj + a.xjE*b.xyi + a.xijE*b.xy
	 + a.y*b.ijE - a.yi*b.jE + a.yj*b.iE + a.yij*b.E - a.yE*b.ij - a.yiE*b.j + a.yjE*b.i + a.yijE*b.q
	 - a.xy*b.xijE + a.xyi*b.xjE - a.xyj*b.xiE - a.xyij*b.xE + a.xyE*b.xij + a.xyiE*b.xj - a.xyjE*b.xi - a.xyijE*b.x
	 - a.z*b.yzijE + a.zi*b.yzjE - a.zj*b.yziE - a.zij*b.yzE + a.zE*b.yzij + a.ziE*b.yzj - a.zjE*b.yzi - a.zijE*b.yz
	 + a.xz*b.xyzijE - a.xzi*b.xyzjE + a.xzj*b.xyziE + a.xzij*b.xyzE - a.xzE*b.xyzij - a.xziE*b.xyzj + a.xzjE*b.xyzi + a.xzijE*b.xyz
	 + a.yz*b.zijE - a.yzi*b.zjE + a.yzj*b.ziE + a.yzij*b.zE - a.yzE*b.zij - a.yziE*b.zj + a.yzjE*b.zi + a.yzijE*b.z
	 + a.xyz*b.xzijE - a.xyzi*b.xzjE + a.xyzj*b.xziE + a.xyzij*b.xzE - a.xyzE*b.xzij - a.xyziE*b.xzj + a.xyzjE*b.xzi + a.xyzijE*b.xz ;

c.xy =  + a.q*b.xy - a.i*b.xyi - a.j*b.xyj - a.ij*b.xyij - a.E*b.xyE - a.iE*b.xyiE - a.jE*b.xyjE - a.ijE*b.xyijE
	 + a.x*b.y - a.xi*b.yi - a.xj*b.yj - a.xij*b.yij - a.xE*b.yE - a.xiE*b.yiE - a.xjE*b.yjE - a.xijE*b.yijE
	 - a.y*b.x + a.yi*b.xi + a.yj*b.xj + a.yij*b.xij + a.yE*b.xE + a.yiE*b.xiE + a.yjE*b.xjE + a.yijE*b.xijE
	 + a.xy*b.q - a.xyi*b.i - a.xyj*b.j - a.xyij*b.ij - a.xyE*b.E - a.xyiE*b.iE - a.xyjE*b.jE - a.xyijE*b.ijE
	 + a.z*b.xyz - a.zi*b.xyzi - a.zj*b.xyzj - a.zij*b.xyzij - a.zE*b.xyzE - a.ziE*b.xyziE - a.zjE*b.xyzjE - a.zijE*b.xyzijE
	 - a.xz*b.yz + a.xzi*b.yzi + a.xzj*b.yzj + a.xzij*b.yzij + a.xzE*b.yzE + a.xziE*b.yziE + a.xzjE*b.yzjE + a.xzijE*b.yzijE
	 + a.yz*b.xz - a.yzi*b.xzi - a.yzj*b.xzj - a.yzij*b.xzij - a.yzE*b.xzE - a.yziE*b.xziE - a.yzjE*b.xzjE - a.yzijE*b.xzijE
	 + a.xyz*b.z - a.xyzi*b.zi - a.xyzj*b.zj - a.xyzij*b.zij - a.xyzE*b.zE - a.xyziE*b.ziE - a.xyzjE*b.zjE - a.xyzijE*b.zijE ;

c.xyi =  + a.q*b.xyi + a.i*b.xy + a.j*b.xyij - a.ij*b.xyj + a.E*b.xyiE - a.iE*b.xyE - a.jE*b.xyijE + a.ijE*b.xyjE
	 + a.x*b.yi + a.xi*b.y + a.xj*b.yij - a.xij*b.yj + a.xE*b.yiE - a.xiE*b.yE - a.xjE*b.yijE + a.xijE*b.yjE
	 - a.y*b.xi - a.yi*b.x - a.yj*b.xij + a.yij*b.xj - a.yE*b.xiE + a.yiE*b.xE + a.yjE*b.xijE - a.yijE*b.xjE
	 + a.xy*b.i + a.xyi*b.q + a.xyj*b.ij - a.xyij*b.j + a.xyE*b.iE - a.xyiE*b.E - a.xyjE*b.ijE + a.xyijE*b.jE
	 + a.z*b.xyzi + a.zi*b.xyz + a.zj*b.xyzij - a.zij*b.xyzj + a.zE*b.xyziE - a.ziE*b.xyzE - a.zjE*b.xyzijE + a.zijE*b.xyzjE
	 - a.xz*b.yzi - a.xzi*b.yz - a.xzj*b.yzij + a.xzij*b.yzj - a.xzE*b.yziE + a.xziE*b.yzE + a.xzjE*b.yzijE - a.xzijE*b.yzjE
	 + a.yz*b.xzi + a.yzi*b.xz + a.yzj*b.xzij - a.yzij*b.xzj + a.yzE*b.xziE - a.yziE*b.xzE - a.yzjE*b.xzijE + a.yzijE*b.xzjE
	 + a.xyz*b.zi + a.xyzi*b.z + a.xyzj*b.zij - a.xyzij*b.zj + a.xyzE*b.ziE - a.xyziE*b.zE - a.xyzjE*b.zijE + a.xyzijE*b.zjE ;

c.xyj =  + a.q*b.xyj - a.i*b.xyij + a.j*b.xy + a.ij*b.xyi + a.E*b.xyjE + a.iE*b.xyijE - a.jE*b.xyE - a.ijE*b.xyiE
	 + a.x*b.yj - a.xi*b.yij + a.xj*b.y + a.xij*b.yi + a.xE*b.yjE + a.xiE*b.yijE - a.xjE*b.yE - a.xijE*b.yiE
	 - a.y*b.xj + a.yi*b.xij - a.yj*b.x - a.yij*b.xi - a.yE*b.xjE - a.yiE*b.xijE + a.yjE*b.xE + a.yijE*b.xiE
	 + a.xy*b.j - a.xyi*b.ij + a.xyj*b.q + a.xyij*b.i + a.xyE*b.jE + a.xyiE*b.ijE - a.xyjE*b.E - a.xyijE*b.iE
	 + a.z*b.xyzj - a.zi*b.xyzij + a.zj*b.xyz + a.zij*b.xyzi + a.zE*b.xyzjE + a.ziE*b.xyzijE - a.zjE*b.xyzE - a.zijE*b.xyziE
	 - a.xz*b.yzj + a.xzi*b.yzij - a.xzj*b.yz - a.xzij*b.yzi - a.xzE*b.yzjE - a.xziE*b.yzijE + a.xzjE*b.yzE + a.xzijE*b.yziE
	 + a.yz*b.xzj - a.yzi*b.xzij + a.yzj*b.xz + a.yzij*b.xzi + a.yzE*b.xzjE + a.yziE*b.xzijE - a.yzjE*b.xzE - a.yzijE*b.xziE
	 + a.xyz*b.zj - a.xyzi*b.zij + a.xyzj*b.z + a.xyzij*b.zi + a.xyzE*b.zjE + a.xyziE*b.zijE - a.xyzjE*b.zE - a.xyzijE*b.ziE ;

c.xyij =  + a.q*b.xyij + a.i*b.xyj - a.j*b.xyi + a.ij*b.xy + a.E*b.xyijE - a.iE*b.xyjE + a.jE*b.xyiE - a.ijE*b.xyE
	 + a.x*b.yij + a.xi*b.yj - a.xj*b.yi + a.xij*b.y + a.xE*b.yijE - a.xiE*b.yjE + a.xjE*b.yiE - a.xijE*b.yE
	 - a.y*b.xij - a.yi*b.xj + a.yj*b.xi - a.yij*b.x - a.yE*b.xijE + a.yiE*b.xjE - a.yjE*b.xiE + a.yijE*b.xE
	 + a.xy*b.ij + a.xyi*b.j - a.xyj*b.i + a.xyij*b.q + a.xyE*b.ijE - a.xyiE*b.jE + a.xyjE*b.iE - a.xyijE*b.E
	 + a.z*b.xyzij + a.zi*b.xyzj - a.zj*b.xyzi + a.zij*b.xyz + a.zE*b.xyzijE - a.ziE*b.xyzjE + a.zjE*b.xyziE - a.zijE*b.xyzE
	 - a.xz*b.yzij - a.xzi*b.yzj + a.xzj*b.yzi - a.xzij*b.yz - a.xzE*b.yzijE + a.xziE*b.yzjE - a.xzjE*b.yziE + a.xzijE*b.yzE
	 + a.yz*b.xzij + a.yzi*b.xzj - a.yzj*b.xzi + a.yzij*b.xz + a.yzE*b.xzijE - a.yziE*b.xzjE + a.yzjE*b.xziE - a.yzijE*b.xzE
	 + a.xyz*b.zij + a.xyzi*b.zj - a.xyzj*b.zi + a.xyzij*b.z + a.xyzE*b.zijE - a.xyziE*b.zjE + a.xyzjE*b.ziE - a.xyzijE*b.zE ;

c.xyE =  + a.q*b.xyE - a.i*b.xyiE - a.j*b.xyjE - a.ij*b.xyijE + a.E*b.xy + a.iE*b.xyi + a.jE*b.xyj + a.ijE*b.xyij
	 + a.x*b.yE - a.xi*b.yiE - a.xj*b.yjE - a.xij*b.yijE + a.xE*b.y + a.xiE*b.yi + a.xjE*b.yj + a.xijE*b.yij
	 - a.y*b.xE + a.yi*b.xiE + a.yj*b.xjE + a.yij*b.xijE - a.yE*b.x - a.yiE*b.xi - a.yjE*b.xj - a.yijE*b.xij
	 + a.xy*b.E - a.xyi*b.iE - a.xyj*b.jE - a.xyij*b.ijE + a.xyE*b.q + a.xyiE*b.i + a.xyjE*b.j + a.xyijE*b.ij
	 + a.z*b.xyzE - a.zi*b.xyziE - a.zj*b.xyzjE - a.zij*b.xyzijE + a.zE*b.xyz + a.ziE*b.xyzi + a.zjE*b.xyzj + a.zijE*b.xyzij
	 - a.xz*b.yzE + a.xzi*b.yziE + a.xzj*b.yzjE + a.xzij*b.yzijE - a.xzE*b.yz - a.xziE*b.yzi - a.xzjE*b.yzj - a.xzijE*b.yzij
	 + a.yz*b.xzE - a.yzi*b.xziE - a.yzj*b.xzjE - a.yzij*b.xzijE + a.yzE*b.xz + a.yziE*b.xzi + a.yzjE*b.xzj + a.yzijE*b.xzij
	 + a.xyz*b.zE - a.xyzi*b.ziE - a.xyzj*b.zjE - a.xyzij*b.zijE + a.xyzE*b.z + a.xyziE*b.zi + a.xyzjE*b.zj + a.xyzijE*b.zij ;

c.xyiE =  + a.q*b.xyiE + a.i*b.xyE - a.j*b.xyijE + a.ij*b.xyjE - a.E*b.xyi + a.iE*b.xy - a.jE*b.xyij + a.ijE*b.xyj
	 + a.x*b.yiE + a.xi*b.yE - a.xj*b.yijE + a.xij*b.yjE - a.xE*b.yi + a.xiE*b.y - a.xjE*b.yij + a.xijE*b.yj
	 - a.y*b.xiE - a.yi*b.xE + a.yj*b.xijE - a.yij*b.xjE + a.yE*b.xi - a.yiE*b.x + a.yjE*b.xij - a.yijE*b.xj
	 + a.xy*b.iE + a.xyi*b.E - a.xyj*b.ijE + a.xyij*b.jE - a.xyE*b.i + a.xyiE*b.q - a.xyjE*b.ij + a.xyijE*b.j
	 + a.z*b.xyziE + a.zi*b.xyzE - a.zj*b.xyzijE + a.zij*b.xyzjE - a.zE*b.xyzi + a.ziE*b.xyz - a.zjE*b.xyzij + a.zijE*b.xyzj
	 - a.xz*b.yziE - a.xzi*b.yzE + a.xzj*b.yzijE - a.xzij*b.yzjE + a.xzE*b.yzi - a.xziE*b.yz + a.xzjE*b.yzij - a.xzijE*b.yzj
	 + a.yz*b.xziE + a.yzi*b.xzE - a.yzj*b.xzijE + a.yzij*b.xzjE - a.yzE*b.xzi + a.yziE*b.xz - a.yzjE*b.xzij + a.yzijE*b.xzj
	 + a.xyz*b.ziE + a.xyzi*b.zE - a.xyzj*b.zijE + a.xyzij*b.zjE - a.xyzE*b.zi + a.xyziE*b.z - a.xyzjE*b.zij + a.xyzijE*b.zj ;

c.xyjE =  + a.q*b.xyjE + a.i*b.xyijE + a.j*b.xyE - a.ij*b.xyiE - a.E*b.xyj + a.iE*b.xyij + a.jE*b.xy - a.ijE*b.xyi
	 + a.x*b.yjE + a.xi*b.yijE + a.xj*b.yE - a.xij*b.yiE - a.xE*b.yj + a.xiE*b.yij + a.xjE*b.y - a.xijE*b.yi
	 - a.y*b.xjE - a.yi*b.xijE - a.yj*b.xE + a.yij*b.xiE + a.yE*b.xj - a.yiE*b.xij - a.yjE*b.x + a.yijE*b.xi
	 + a.xy*b.jE + a.xyi*b.ijE + a.xyj*b.E - a.xyij*b.iE - a.xyE*b.j + a.xyiE*b.ij + a.xyjE*b.q - a.xyijE*b.i
	 + a.z*b.xyzjE + a.zi*b.xyzijE + a.zj*b.xyzE - a.zij*b.xyziE - a.zE*b.xyzj + a.ziE*b.xyzij + a.zjE*b.xyz - a.zijE*b.xyzi
	 - a.xz*b.yzjE - a.xzi*b.yzijE - a.xzj*b.yzE + a.xzij*b.yziE + a.xzE*b.yzj - a.xziE*b.yzij - a.xzjE*b.yz + a.xzijE*b.yzi
	 + a.yz*b.xzjE + a.yzi*b.xzijE + a.yzj*b.xzE - a.yzij*b.xziE - a.yzE*b.xzj + a.yziE*b.xzij + a.yzjE*b.xz - a.yzijE*b.xzi
	 + a.xyz*b.zjE + a.xyzi*b.zijE + a.xyzj*b.zE - a.xyzij*b.ziE - a.xyzE*b.zj + a.xyziE*b.zij + a.xyzjE*b.z - a.xyzijE*b.zi ;

c.xyijE =  + a.q*b.xyijE - a.i*b.xyjE + a.j*b.xyiE + a.ij*b.xyE - a.E*b.xyij - a.iE*b.xyj + a.jE*b.xyi + a.ijE*b.xy
	 + a.x*b.yijE - a.xi*b.yjE + a.xj*b.yiE + a.xij*b.yE - a.xE*b.yij - a.xiE*b.yj + a.xjE*b.yi + a.xijE*b.y
	 - a.y*b.xijE + a.yi*b.xjE - a.yj*b.xiE - a.yij*b.xE + a.yE*b.xij + a.yiE*b.xj - a.yjE*b.xi - a.yijE*b.x
	 + a.xy*b.ijE - a.xyi*b.jE + a.xyj*b.iE + a.xyij*b.E - a.xyE*b.ij - a.xyiE*b.j + a.xyjE*b.i + a.xyijE*b.q
	 + a.z*b.xyzijE - a.zi*b.xyzjE + a.zj*b.xyziE + a.zij*b.xyzE - a.zE*b.xyzij - a.ziE*b.xyzj + a.zjE*b.xyzi + a.zijE*b.xyz
	 - a.xz*b.yzijE + a.xzi*b.yzjE - a.xzj*b.yziE - a.xzij*b.yzE + a.xzE*b.yzij + a.xziE*b.yzj - a.xzjE*b.yzi - a.xzijE*b.yz
	 + a.yz*b.xzijE - a.yzi*b.xzjE + a.yzj*b.xziE + a.yzij*b.xzE - a.yzE*b.xzij - a.yziE*b.xzj + a.yzjE*b.xzi + a.yzijE*b.xz
	 + a.xyz*b.zijE - a.xyzi*b.zjE + a.xyzj*b.ziE + a.xyzij*b.zE - a.xyzE*b.zij - a.xyziE*b.zj + a.xyzjE*b.zi + a.xyzijE*b.z ;

c.z =  + a.q*b.z - a.i*b.zi - a.j*b.zj - a.ij*b.zij - a.E*b.zE - a.iE*b.ziE - a.jE*b.zjE - a.ijE*b.zijE
	 + a.x*b.xz - a.xi*b.xzi - a.xj*b.xzj - a.xij*b.xzij - a.xE*b.xzE - a.xiE*b.xziE - a.xjE*b.xzjE - a.xijE*b.xzijE
	 + a.y*b.yz - a.yi*b.yzi - a.yj*b.yzj - a.yij*b.yzij - a.yE*b.yzE - a.yiE*b.yziE - a.yjE*b.yzjE - a.yijE*b.yzijE
	 - a.xy*b.xyz + a.xyi*b.xyzi + a.xyj*b.xyzj + a.xyij*b.xyzij + a.xyE*b.xyzE + a.xyiE*b.xyziE + a.xyjE*b.xyzjE + a.xyijE*b.xyzijE
	 + a.z*b.q - a.zi*b.i - a.zj*b.j - a.zij*b.ij - a.zE*b.E - a.ziE*b.iE - a.zjE*b.jE - a.zijE*b.ijE
	 - a.xz*b.x + a.xzi*b.xi + a.xzj*b.xj + a.xzij*b.xij + a.xzE*b.xE + a.xziE*b.xiE + a.xzjE*b.xjE + a.xzijE*b.xijE
	 - a.yz*b.y + a.yzi*b.yi + a.yzj*b.yj + a.yzij*b.yij + a.yzE*b.yE + a.yziE*b.yiE + a.yzjE*b.yjE + a.yzijE*b.yijE
	 - a.xyz*b.xy + a.xyzi*b.xyi + a.xyzj*b.xyj + a.xyzij*b.xyij + a.xyzE*b.xyE + a.xyziE*b.xyiE + a.xyzjE*b.xyjE + a.xyzijE*b.xyijE ;

c.zi =  + a.q*b.zi + a.i*b.z + a.j*b.zij - a.ij*b.zj + a.E*b.ziE - a.iE*b.zE - a.jE*b.zijE + a.ijE*b.zjE
	 + a.x*b.xzi + a.xi*b.xz + a.xj*b.xzij - a.xij*b.xzj + a.xE*b.xziE - a.xiE*b.xzE - a.xjE*b.xzijE + a.xijE*b.xzjE
	 + a.y*b.yzi + a.yi*b.yz + a.yj*b.yzij - a.yij*b.yzj + a.yE*b.yziE - a.yiE*b.yzE - a.yjE*b.yzijE + a.yijE*b.yzjE
	 - a.xy*b.xyzi - a.xyi*b.xyz - a.xyj*b.xyzij + a.xyij*b.xyzj - a.xyE*b.xyziE + a.xyiE*b.xyzE + a.xyjE*b.xyzijE - a.xyijE*b.xyzjE
	 + a.z*b.i + a.zi*b.q + a.zj*b.ij - a.zij*b.j + a.zE*b.iE - a.ziE*b.E - a.zjE*b.ijE + a.zijE*b.jE
	 - a.xz*b.xi - a.xzi*b.x - a.xzj*b.xij + a.xzij*b.xj - a.xzE*b.xiE + a.xziE*b.xE + a.xzjE*b.xijE - a.xzijE*b.xjE
	 - a.yz*b.yi - a.yzi*b.y - a.yzj*b.yij + a.yzij*b.yj - a.yzE*b.yiE + a.yziE*b.yE + a.yzjE*b.yijE - a.yzijE*b.yjE
	 - a.xyz*b.xyi - a.xyzi*b.xy - a.xyzj*b.xyij + a.xyzij*b.xyj - a.xyzE*b.xyiE + a.xyziE*b.xyE + a.xyzjE*b.xyijE - a.xyzijE*b.xyjE ;

c.zj =  + a.q*b.zj - a.i*b.zij + a.j*b.z + a.ij*b.zi + a.E*b.zjE + a.iE*b.zijE - a.jE*b.zE - a.ijE*b.ziE
	 + a.x*b.xzj - a.xi*b.xzij + a.xj*b.xz + a.xij*b.xzi + a.xE*b.xzjE + a.xiE*b.xzijE - a.xjE*b.xzE - a.xijE*b.xziE
	 + a.y*b.yzj - a.yi*b.yzij + a.yj*b.yz + a.yij*b.yzi + a.yE*b.yzjE + a.yiE*b.yzijE - a.yjE*b.yzE - a.yijE*b.yziE
	 - a.xy*b.xyzj + a.xyi*b.xyzij - a.xyj*b.xyz - a.xyij*b.xyzi - a.xyE*b.xyzjE - a.xyiE*b.xyzijE + a.xyjE*b.xyzE + a.xyijE*b.xyziE
	 + a.z*b.j - a.zi*b.ij + a.zj*b.q + a.zij*b.i + a.zE*b.jE + a.ziE*b.ijE - a.zjE*b.E - a.zijE*b.iE
	 - a.xz*b.xj + a.xzi*b.xij - a.xzj*b.x - a.xzij*b.xi - a.xzE*b.xjE - a.xziE*b.xijE + a.xzjE*b.xE + a.xzijE*b.xiE
	 - a.yz*b.yj + a.yzi*b.yij - a.yzj*b.y - a.yzij*b.yi - a.yzE*b.yjE - a.yziE*b.yijE + a.yzjE*b.yE + a.yzijE*b.yiE
	 - a.xyz*b.xyj + a.xyzi*b.xyij - a.xyzj*b.xy - a.xyzij*b.xyi - a.xyzE*b.xyjE - a.xyziE*b.xyijE + a.xyzjE*b.xyE + a.xyzijE*b.xyiE ;

c.zij =  + a.q*b.zij + a.i*b.zj - a.j*b.zi + a.ij*b.z + a.E*b.zijE - a.iE*b.zjE + a.jE*b.ziE - a.ijE*b.zE
	 + a.x*b.xzij + a.xi*b.xzj - a.xj*b.xzi + a.xij*b.xz + a.xE*b.xzijE - a.xiE*b.xzjE + a.xjE*b.xziE - a.xijE*b.xzE
	 + a.y*b.yzij + a.yi*b.yzj - a.yj*b.yzi + a.yij*b.yz + a.yE*b.yzijE - a.yiE*b.yzjE + a.yjE*b.yziE - a.yijE*b.yzE
	 - a.xy*b.xyzij - a.xyi*b.xyzj + a.xyj*b.xyzi - a.xyij*b.xyz - a.xyE*b.xyzijE + a.xyiE*b.xyzjE - a.xyjE*b.xyziE + a.xyijE*b.xyzE
	 + a.z*b.ij + a.zi*b.j - a.zj*b.i + a.zij*b.q + a.zE*b.ijE - a.ziE*b.jE + a.zjE*b.iE - a.zijE*b.E
	 - a.xz*b.xij - a.xzi*b.xj + a.xzj*b.xi - a.xzij*b.x - a.xzE*b.xijE + a.xziE*b.xjE - a.xzjE*b.xiE + a.xzijE*b.xE
	 - a.yz*b.yij - a.yzi*b.yj + a.yzj*b.yi - a.yzij*b.y - a.yzE*b.yijE + a.yziE*b.yjE - a.yzjE*b.yiE + a.yzijE*b.yE
	 - a.xyz*b.xyij - a.xyzi*b.xyj + a.xyzj*b.xyi - a.xyzij*b.xy - a.xyzE*b.xyijE + a.xyziE*b.xyjE - a.xyzjE*b.xyiE + a.xyzijE*b.xyE ;

c.zE =  + a.q*b.zE - a.i*b.ziE - a.j*b.zjE - a.ij*b.zijE + a.E*b.z + a.iE*b.zi + a.jE*b.zj + a.ijE*b.zij
	 + a.x*b.xzE - a.xi*b.xziE - a.xj*b.xzjE - a.xij*b.xzijE + a.xE*b.xz + a.xiE*b.xzi + a.xjE*b.xzj + a.xijE*b.xzij
	 + a.y*b.yzE - a.yi*b.yziE - a.yj*b.yzjE - a.yij*b.yzijE + a.yE*b.yz + a.yiE*b.yzi + a.yjE*b.yzj + a.yijE*b.yzij
	 - a.xy*b.xyzE + a.xyi*b.xyziE + a.xyj*b.xyzjE + a.xyij*b.xyzijE - a.xyE*b.xyz - a.xyiE*b.xyzi - a.xyjE*b.xyzj - a.xyijE*b.xyzij
	 + a.z*b.E - a.zi*b.iE - a.zj*b.jE - a.zij*b.ijE + a.zE*b.q + a.ziE*b.i + a.zjE*b.j + a.zijE*b.ij
	 - a.xz*b.xE + a.xzi*b.xiE + a.xzj*b.xjE + a.xzij*b.xijE - a.xzE*b.x - a.xziE*b.xi - a.xzjE*b.xj - a.xzijE*b.xij
	 - a.yz*b.yE + a.yzi*b.yiE + a.yzj*b.yjE + a.yzij*b.yijE - a.yzE*b.y - a.yziE*b.yi - a.yzjE*b.yj - a.yzijE*b.yij
	 - a.xyz*b.xyE + a.xyzi*b.xyiE + a.xyzj*b.xyjE + a.xyzij*b.xyijE - a.xyzE*b.xy - a.xyziE*b.xyi - a.xyzjE*b.xyj - a.xyzijE*b.xyij ;

c.ziE =  + a.q*b.ziE + a.i*b.zE - a.j*b.zijE + a.ij*b.zjE - a.E*b.zi + a.iE*b.z - a.jE*b.zij + a.ijE*b.zj
	 + a.x*b.xziE + a.xi*b.xzE - a.xj*b.xzijE + a.xij*b.xzjE - a.xE*b.xzi + a.xiE*b.xz - a.xjE*b.xzij + a.xijE*b.xzj
	 + a.y*b.yziE + a.yi*b.yzE - a.yj*b.yzijE + a.yij*b.yzjE - a.yE*b.yzi + a.yiE*b.yz - a.yjE*b.yzij + a.yijE*b.yzj
	 - a.xy*b.xyziE - a.xyi*b.xyzE + a.xyj*b.xyzijE - a.xyij*b.xyzjE + a.xyE*b.xyzi - a.xyiE*b.xyz + a.xyjE*b.xyzij - a.xyijE*b.xyzj
	 + a.z*b.iE + a.zi*b.E - a.zj*b.ijE + a.zij*b.jE - a.zE*b.i + a.ziE*b.q - a.zjE*b.ij + a.zijE*b.j
	 - a.xz*b.xiE - a.xzi*b.xE + a.xzj*b.xijE - a.xzij*b.xjE + a.xzE*b.xi - a.xziE*b.x + a.xzjE*b.xij - a.xzijE*b.xj
	 - a.yz*b.yiE - a.yzi*b.yE + a.yzj*b.yijE - a.yzij*b.yjE + a.yzE*b.yi - a.yziE*b.y + a.yzjE*b.yij - a.yzijE*b.yj
	 - a.xyz*b.xyiE - a.xyzi*b.xyE + a.xyzj*b.xyijE - a.xyzij*b.xyjE + a.xyzE*b.xyi - a.xyziE*b.xy + a.xyzjE*b.xyij - a.xyzijE*b.xyj ;

c.zjE =  + a.q*b.zjE + a.i*b.zijE + a.j*b.zE - a.ij*b.ziE - a.E*b.zj + a.iE*b.zij + a.jE*b.z - a.ijE*b.zi
	 + a.x*b.xzjE + a.xi*b.xzijE + a.xj*b.xzE - a.xij*b.xziE - a.xE*b.xzj + a.xiE*b.xzij + a.xjE*b.xz - a.xijE*b.xzi
	 + a.y*b.yzjE + a.yi*b.yzijE + a.yj*b.yzE - a.yij*b.yziE - a.yE*b.yzj + a.yiE*b.yzij + a.yjE*b.yz - a.yijE*b.yzi
	 - a.xy*b.xyzjE - a.xyi*b.xyzijE - a.xyj*b.xyzE + a.xyij*b.xyziE + a.xyE*b.xyzj - a.xyiE*b.xyzij - a.xyjE*b.xyz + a.xyijE*b.xyzi
	 + a.z*b.jE + a.zi*b.ijE + a.zj*b.E - a.zij*b.iE - a.zE*b.j + a.ziE*b.ij + a.zjE*b.q - a.zijE*b.i
	 - a.xz*b.xjE - a.xzi*b.xijE - a.xzj*b.xE + a.xzij*b.xiE + a.xzE*b.xj - a.xziE*b.xij - a.xzjE*b.x + a.xzijE*b.xi
	 - a.yz*b.yjE - a.yzi*b.yijE - a.yzj*b.yE + a.yzij*b.yiE + a.yzE*b.yj - a.yziE*b.yij - a.yzjE*b.y + a.yzijE*b.yi
	 - a.xyz*b.xyjE - a.xyzi*b.xyijE - a.xyzj*b.xyE + a.xyzij*b.xyiE + a.xyzE*b.xyj - a.xyziE*b.xyij - a.xyzjE*b.xy + a.xyzijE*b.xyi ;

c.zijE =  + a.q*b.zijE - a.i*b.zjE + a.j*b.ziE + a.ij*b.zE - a.E*b.zij - a.iE*b.zj + a.jE*b.zi + a.ijE*b.z
	 + a.x*b.xzijE - a.xi*b.xzjE + a.xj*b.xziE + a.xij*b.xzE - a.xE*b.xzij - a.xiE*b.xzj + a.xjE*b.xzi + a.xijE*b.xz
	 + a.y*b.yzijE - a.yi*b.yzjE + a.yj*b.yziE + a.yij*b.yzE - a.yE*b.yzij - a.yiE*b.yzj + a.yjE*b.yzi + a.yijE*b.yz
	 - a.xy*b.xyzijE + a.xyi*b.xyzjE - a.xyj*b.xyziE - a.xyij*b.xyzE + a.xyE*b.xyzij + a.xyiE*b.xyzj - a.xyjE*b.xyzi - a.xyijE*b.xyz
	 + a.z*b.ijE - a.zi*b.jE + a.zj*b.iE + a.zij*b.E - a.zE*b.ij - a.ziE*b.j + a.zjE*b.i + a.zijE*b.q
	 - a.xz*b.xijE + a.xzi*b.xjE - a.xzj*b.xiE - a.xzij*b.xE + a.xzE*b.xij + a.xziE*b.xj - a.xzjE*b.xi - a.xzijE*b.x
	 - a.yz*b.yijE + a.yzi*b.yjE - a.yzj*b.yiE - a.yzij*b.yE + a.yzE*b.yij + a.yziE*b.yj - a.yzjE*b.yi - a.yzijE*b.y
	 - a.xyz*b.xyijE + a.xyzi*b.xyjE - a.xyzj*b.xyiE - a.xyzij*b.xyE + a.xyzE*b.xyij + a.xyziE*b.xyj - a.xyzjE*b.xyi - a.xyzijE*b.xy ;

c.xz =  + a.q*b.xz - a.i*b.xzi - a.j*b.xzj - a.ij*b.xzij - a.E*b.xzE - a.iE*b.xziE - a.jE*b.xzjE - a.ijE*b.xzijE
	 + a.x*b.z - a.xi*b.zi - a.xj*b.zj - a.xij*b.zij - a.xE*b.zE - a.xiE*b.ziE - a.xjE*b.zjE - a.xijE*b.zijE
	 - a.y*b.xyz + a.yi*b.xyzi + a.yj*b.xyzj + a.yij*b.xyzij + a.yE*b.xyzE + a.yiE*b.xyziE + a.yjE*b.xyzjE + a.yijE*b.xyzijE
	 + a.xy*b.yz - a.xyi*b.yzi - a.xyj*b.yzj - a.xyij*b.yzij - a.xyE*b.yzE - a.xyiE*b.yziE - a.xyjE*b.yzjE - a.xyijE*b.yzijE
	 - a.z*b.x + a.zi*b.xi + a.zj*b.xj + a.zij*b.xij + a.zE*b.xE + a.ziE*b.xiE + a.zjE*b.xjE + a.zijE*b.xijE
	 + a.xz*b.q - a.xzi*b.i - a.xzj*b.j - a.xzij*b.ij - a.xzE*b.E - a.xziE*b.iE - a.xzjE*b.jE - a.xzijE*b.ijE
	 - a.yz*b.xy + a.yzi*b.xyi + a.yzj*b.xyj + a.yzij*b.xyij + a.yzE*b.xyE + a.yziE*b.xyiE + a.yzjE*b.xyjE + a.yzijE*b.xyijE
	 - a.xyz*b.y + a.xyzi*b.yi + a.xyzj*b.yj + a.xyzij*b.yij + a.xyzE*b.yE + a.xyziE*b.yiE + a.xyzjE*b.yjE + a.xyzijE*b.yijE ;

c.xzi =  + a.q*b.xzi + a.i*b.xz + a.j*b.xzij - a.ij*b.xzj + a.E*b.xziE - a.iE*b.xzE - a.jE*b.xzijE + a.ijE*b.xzjE
	 + a.x*b.zi + a.xi*b.z + a.xj*b.zij - a.xij*b.zj + a.xE*b.ziE - a.xiE*b.zE - a.xjE*b.zijE + a.xijE*b.zjE
	 - a.y*b.xyzi - a.yi*b.xyz - a.yj*b.xyzij + a.yij*b.xyzj - a.yE*b.xyziE + a.yiE*b.xyzE + a.yjE*b.xyzijE - a.yijE*b.xyzjE
	 + a.xy*b.yzi + a.xyi*b.yz + a.xyj*b.yzij - a.xyij*b.yzj + a.xyE*b.yziE - a.xyiE*b.yzE - a.xyjE*b.yzijE + a.xyijE*b.yzjE
	 - a.z*b.xi - a.zi*b.x - a.zj*b.xij + a.zij*b.xj - a.zE*b.xiE + a.ziE*b.xE + a.zjE*b.xijE - a.zijE*b.xjE
	 + a.xz*b.i + a.xzi*b.q + a.xzj*b.ij - a.xzij*b.j + a.xzE*b.iE - a.xziE*b.E - a.xzjE*b.ijE + a.xzijE*b.jE
	 - a.yz*b.xyi - a.yzi*b.xy - a.yzj*b.xyij + a.yzij*b.xyj - a.yzE*b.xyiE + a.yziE*b.xyE + a.yzjE*b.xyijE - a.yzijE*b.xyjE
	 - a.xyz*b.yi - a.xyzi*b.y - a.xyzj*b.yij + a.xyzij*b.yj - a.xyzE*b.yiE + a.xyziE*b.yE + a.xyzjE*b.yijE - a.xyzijE*b.yjE ;

c.xzj =  + a.q*b.xzj - a.i*b.xzij + a.j*b.xz + a.ij*b.xzi + a.E*b.xzjE + a.iE*b.xzijE - a.jE*b.xzE - a.ijE*b.xziE
	 + a.x*b.zj - a.xi*b.zij + a.xj*b.z + a.xij*b.zi + a.xE*b.zjE + a.xiE*b.zijE - a.xjE*b.zE - a.xijE*b.ziE
	 - a.y*b.xyzj + a.yi*b.xyzij - a.yj*b.xyz - a.yij*b.xyzi - a.yE*b.xyzjE - a.yiE*b.xyzijE + a.yjE*b.xyzE + a.yijE*b.xyziE
	 + a.xy*b.yzj - a.xyi*b.yzij + a.xyj*b.yz + a.xyij*b.yzi + a.xyE*b.yzjE + a.xyiE*b.yzijE - a.xyjE*b.yzE - a.xyijE*b.yziE
	 - a.z*b.xj + a.zi*b.xij - a.zj*b.x - a.zij*b.xi - a.zE*b.xjE - a.ziE*b.xijE + a.zjE*b.xE + a.zijE*b.xiE
	 + a.xz*b.j - a.xzi*b.ij + a.xzj*b.q + a.xzij*b.i + a.xzE*b.jE + a.xziE*b.ijE - a.xzjE*b.E - a.xzijE*b.iE
	 - a.yz*b.xyj + a.yzi*b.xyij - a.yzj*b.xy - a.yzij*b.xyi - a.yzE*b.xyjE - a.yziE*b.xyijE + a.yzjE*b.xyE + a.yzijE*b.xyiE
	 - a.xyz*b.yj + a.xyzi*b.yij - a.xyzj*b.y - a.xyzij*b.yi - a.xyzE*b.yjE - a.xyziE*b.yijE + a.xyzjE*b.yE + a.xyzijE*b.yiE ;

c.xzij =  + a.q*b.xzij + a.i*b.xzj - a.j*b.xzi + a.ij*b.xz + a.E*b.xzijE - a.iE*b.xzjE + a.jE*b.xziE - a.ijE*b.xzE
	 + a.x*b.zij + a.xi*b.zj - a.xj*b.zi + a.xij*b.z + a.xE*b.zijE - a.xiE*b.zjE + a.xjE*b.ziE - a.xijE*b.zE
	 - a.y*b.xyzij - a.yi*b.xyzj + a.yj*b.xyzi - a.yij*b.xyz - a.yE*b.xyzijE + a.yiE*b.xyzjE - a.yjE*b.xyziE + a.yijE*b.xyzE
	 + a.xy*b.yzij + a.xyi*b.yzj - a.xyj*b.yzi + a.xyij*b.yz + a.xyE*b.yzijE - a.xyiE*b.yzjE + a.xyjE*b.yziE - a.xyijE*b.yzE
	 - a.z*b.xij - a.zi*b.xj + a.zj*b.xi - a.zij*b.x - a.zE*b.xijE + a.ziE*b.xjE - a.zjE*b.xiE + a.zijE*b.xE
	 + a.xz*b.ij + a.xzi*b.j - a.xzj*b.i + a.xzij*b.q + a.xzE*b.ijE - a.xziE*b.jE + a.xzjE*b.iE - a.xzijE*b.E
	 - a.yz*b.xyij - a.yzi*b.xyj + a.yzj*b.xyi - a.yzij*b.xy - a.yzE*b.xyijE + a.yziE*b.xyjE - a.yzjE*b.xyiE + a.yzijE*b.xyE
	 - a.xyz*b.yij - a.xyzi*b.yj + a.xyzj*b.yi - a.xyzij*b.y - a.xyzE*b.yijE + a.xyziE*b.yjE - a.xyzjE*b.yiE + a.xyzijE*b.yE ;

c.xzE =  + a.q*b.xzE - a.i*b.xziE - a.j*b.xzjE - a.ij*b.xzijE + a.E*b.xz + a.iE*b.xzi + a.jE*b.xzj + a.ijE*b.xzij
	 + a.x*b.zE - a.xi*b.ziE - a.xj*b.zjE - a.xij*b.zijE + a.xE*b.z + a.xiE*b.zi + a.xjE*b.zj + a.xijE*b.zij
	 - a.y*b.xyzE + a.yi*b.xyziE + a.yj*b.xyzjE + a.yij*b.xyzijE - a.yE*b.xyz - a.yiE*b.xyzi - a.yjE*b.xyzj - a.yijE*b.xyzij
	 + a.xy*b.yzE - a.xyi*b.yziE - a.xyj*b.yzjE - a.xyij*b.yzijE + a.xyE*b.yz + a.xyiE*b.yzi + a.xyjE*b.yzj + a.xyijE*b.yzij
	 - a.z*b.xE + a.zi*b.xiE + a.zj*b.xjE + a.zij*b.xijE - a.zE*b.x - a.ziE*b.xi - a.zjE*b.xj - a.zijE*b.xij
	 + a.xz*b.E - a.xzi*b.iE - a.xzj*b.jE - a.xzij*b.ijE + a.xzE*b.q + a.xziE*b.i + a.xzjE*b.j + a.xzijE*b.ij
	 - a.yz*b.xyE + a.yzi*b.xyiE + a.yzj*b.xyjE + a.yzij*b.xyijE - a.yzE*b.xy - a.yziE*b.xyi - a.yzjE*b.xyj - a.yzijE*b.xyij
	 - a.xyz*b.yE + a.xyzi*b.yiE + a.xyzj*b.yjE + a.xyzij*b.yijE - a.xyzE*b.y - a.xyziE*b.yi - a.xyzjE*b.yj - a.xyzijE*b.yij ;

c.xziE =  + a.q*b.xziE + a.i*b.xzE - a.j*b.xzijE + a.ij*b.xzjE - a.E*b.xzi + a.iE*b.xz - a.jE*b.xzij + a.ijE*b.xzj
	 + a.x*b.ziE + a.xi*b.zE - a.xj*b.zijE + a.xij*b.zjE - a.xE*b.zi + a.xiE*b.z - a.xjE*b.zij + a.xijE*b.zj
	 - a.y*b.xyziE - a.yi*b.xyzE + a.yj*b.xyzijE - a.yij*b.xyzjE + a.yE*b.xyzi - a.yiE*b.xyz + a.yjE*b.xyzij - a.yijE*b.xyzj
	 + a.xy*b.yziE + a.xyi*b.yzE - a.xyj*b.yzijE + a.xyij*b.yzjE - a.xyE*b.yzi + a.xyiE*b.yz - a.xyjE*b.yzij + a.xyijE*b.yzj
	 - a.z*b.xiE - a.zi*b.xE + a.zj*b.xijE - a.zij*b.xjE + a.zE*b.xi - a.ziE*b.x + a.zjE*b.xij - a.zijE*b.xj
	 + a.xz*b.iE + a.xzi*b.E - a.xzj*b.ijE + a.xzij*b.jE - a.xzE*b.i + a.xziE*b.q - a.xzjE*b.ij + a.xzijE*b.j
	 - a.yz*b.xyiE - a.yzi*b.xyE + a.yzj*b.xyijE - a.yzij*b.xyjE + a.yzE*b.xyi - a.yziE*b.xy + a.yzjE*b.xyij - a.yzijE*b.xyj
	 - a.xyz*b.yiE - a.xyzi*b.yE + a.xyzj*b.yijE - a.xyzij*b.yjE + a.xyzE*b.yi - a.xyziE*b.y + a.xyzjE*b.yij - a.xyzijE*b.yj ;

c.xzjE =  + a.q*b.xzjE + a.i*b.xzijE + a.j*b.xzE - a.ij*b.xziE - a.E*b.xzj + a.iE*b.xzij + a.jE*b.xz - a.ijE*b.xzi
	 + a.x*b.zjE + a.xi*b.zijE + a.xj*b.zE - a.xij*b.ziE - a.xE*b.zj + a.xiE*b.zij + a.xjE*b.z - a.xijE*b.zi
	 - a.y*b.xyzjE - a.yi*b.xyzijE - a.yj*b.xyzE + a.yij*b.xyziE + a.yE*b.xyzj - a.yiE*b.xyzij - a.yjE*b.xyz + a.yijE*b.xyzi
	 + a.xy*b.yzjE + a.xyi*b.yzijE + a.xyj*b.yzE - a.xyij*b.yziE - a.xyE*b.yzj + a.xyiE*b.yzij + a.xyjE*b.yz - a.xyijE*b.yzi
	 - a.z*b.xjE - a.zi*b.xijE - a.zj*b.xE + a.zij*b.xiE + a.zE*b.xj - a.ziE*b.xij - a.zjE*b.x + a.zijE*b.xi
	 + a.xz*b.jE + a.xzi*b.ijE + a.xzj*b.E - a.xzij*b.iE - a.xzE*b.j + a.xziE*b.ij + a.xzjE*b.q - a.xzijE*b.i
	 - a.yz*b.xyjE - a.yzi*b.xyijE - a.yzj*b.xyE + a.yzij*b.xyiE + a.yzE*b.xyj - a.yziE*b.xyij - a.yzjE*b.xy + a.yzijE*b.xyi
	 - a.xyz*b.yjE - a.xyzi*b.yijE - a.xyzj*b.yE + a.xyzij*b.yiE + a.xyzE*b.yj - a.xyziE*b.yij - a.xyzjE*b.y + a.xyzijE*b.yi ;

c.xzijE =  + a.q*b.xzijE - a.i*b.xzjE + a.j*b.xziE + a.ij*b.xzE - a.E*b.xzij - a.iE*b.xzj + a.jE*b.xzi + a.ijE*b.xz
	 + a.x*b.zijE - a.xi*b.zjE + a.xj*b.ziE + a.xij*b.zE - a.xE*b.zij - a.xiE*b.zj + a.xjE*b.zi + a.xijE*b.z
	 - a.y*b.xyzijE + a.yi*b.xyzjE - a.yj*b.xyziE - a.yij*b.xyzE + a.yE*b.xyzij + a.yiE*b.xyzj - a.yjE*b.xyzi - a.yijE*b.xyz
	 + a.xy*b.yzijE - a.xyi*b.yzjE + a.xyj*b.yziE + a.xyij*b.yzE - a.xyE*b.yzij - a.xyiE*b.yzj + a.xyjE*b.yzi + a.xyijE*b.yz
	 - a.z*b.xijE + a.zi*b.xjE - a.zj*b.xiE - a.zij*b.xE + a.zE*b.xij + a.ziE*b.xj - a.zjE*b.xi - a.zijE*b.x
	 + a.xz*b.ijE - a.xzi*b.jE + a.xzj*b.iE + a.xzij*b.E - a.xzE*b.ij - a.xziE*b.j + a.xzjE*b.i + a.xzijE*b.q
	 - a.yz*b.xyijE + a.yzi*b.xyjE - a.yzj*b.xyiE - a.yzij*b.xyE + a.yzE*b.xyij + a.yziE*b.xyj - a.yzjE*b.xyi - a.yzijE*b.xy
	 - a.xyz*b.yijE + a.xyzi*b.yjE - a.xyzj*b.yiE - a.xyzij*b.yE + a.xyzE*b.yij + a.xyziE*b.yj - a.xyzjE*b.yi - a.xyzijE*b.y ;

c.yz =  + a.q*b.yz - a.i*b.yzi - a.j*b.yzj - a.ij*b.yzij - a.E*b.yzE - a.iE*b.yziE - a.jE*b.yzjE - a.ijE*b.yzijE
	 + a.x*b.xyz - a.xi*b.xyzi - a.xj*b.xyzj - a.xij*b.xyzij - a.xE*b.xyzE - a.xiE*b.xyziE - a.xjE*b.xyzjE - a.xijE*b.xyzijE
	 + a.y*b.z - a.yi*b.zi - a.yj*b.zj - a.yij*b.zij - a.yE*b.zE - a.yiE*b.ziE - a.yjE*b.zjE - a.yijE*b.zijE
	 - a.xy*b.xz + a.xyi*b.xzi + a.xyj*b.xzj + a.xyij*b.xzij + a.xyE*b.xzE + a.xyiE*b.xziE + a.xyjE*b.xzjE + a.xyijE*b.xzijE
	 - a.z*b.y + a.zi*b.yi + a.zj*b.yj + a.zij*b.yij + a.zE*b.yE + a.ziE*b.yiE + a.zjE*b.yjE + a.zijE*b.yijE
	 + a.xz*b.xy - a.xzi*b.xyi - a.xzj*b.xyj - a.xzij*b.xyij - a.xzE*b.xyE - a.xziE*b.xyiE - a.xzjE*b.xyjE - a.xzijE*b.xyijE
	 + a.yz*b.q - a.yzi*b.i - a.yzj*b.j - a.yzij*b.ij - a.yzE*b.E - a.yziE*b.iE - a.yzjE*b.jE - a.yzijE*b.ijE
	 + a.xyz*b.x - a.xyzi*b.xi - a.xyzj*b.xj - a.xyzij*b.xij - a.xyzE*b.xE - a.xyziE*b.xiE - a.xyzjE*b.xjE - a.xyzijE*b.xijE ;

c.yzi =  + a.q*b.yzi + a.i*b.yz + a.j*b.yzij - a.ij*b.yzj + a.E*b.yziE - a.iE*b.yzE - a.jE*b.yzijE + a.ijE*b.yzjE
	 + a.x*b.xyzi + a.xi*b.xyz + a.xj*b.xyzij - a.xij*b.xyzj + a.xE*b.xyziE - a.xiE*b.xyzE - a.xjE*b.xyzijE + a.xijE*b.xyzjE
	 + a.y*b.zi + a.yi*b.z + a.yj*b.zij - a.yij*b.zj + a.yE*b.ziE - a.yiE*b.zE - a.yjE*b.zijE + a.yijE*b.zjE
	 - a.xy*b.xzi - a.xyi*b.xz - a.xyj*b.xzij + a.xyij*b.xzj - a.xyE*b.xziE + a.xyiE*b.xzE + a.xyjE*b.xzijE - a.xyijE*b.xzjE
	 - a.z*b.yi - a.zi*b.y - a.zj*b.yij + a.zij*b.yj - a.zE*b.yiE + a.ziE*b.yE + a.zjE*b.yijE - a.zijE*b.yjE
	 + a.xz*b.xyi + a.xzi*b.xy + a.xzj*b.xyij - a.xzij*b.xyj + a.xzE*b.xyiE - a.xziE*b.xyE - a.xzjE*b.xyijE + a.xzijE*b.xyjE
	 + a.yz*b.i + a.yzi*b.q + a.yzj*b.ij - a.yzij*b.j + a.yzE*b.iE - a.yziE*b.E - a.yzjE*b.ijE + a.yzijE*b.jE
	 + a.xyz*b.xi + a.xyzi*b.x + a.xyzj*b.xij - a.xyzij*b.xj + a.xyzE*b.xiE - a.xyziE*b.xE - a.xyzjE*b.xijE + a.xyzijE*b.xjE ;

c.yzj =  + a.q*b.yzj - a.i*b.yzij + a.j*b.yz + a.ij*b.yzi + a.E*b.yzjE + a.iE*b.yzijE - a.jE*b.yzE - a.ijE*b.yziE
	 + a.x*b.xyzj - a.xi*b.xyzij + a.xj*b.xyz + a.xij*b.xyzi + a.xE*b.xyzjE + a.xiE*b.xyzijE - a.xjE*b.xyzE - a.xijE*b.xyziE
	 + a.y*b.zj - a.yi*b.zij + a.yj*b.z + a.yij*b.zi + a.yE*b.zjE + a.yiE*b.zijE - a.yjE*b.zE - a.yijE*b.ziE
	 - a.xy*b.xzj + a.xyi*b.xzij - a.xyj*b.xz - a.xyij*b.xzi - a.xyE*b.xzjE - a.xyiE*b.xzijE + a.xyjE*b.xzE + a.xyijE*b.xziE
	 - a.z*b.yj + a.zi*b.yij - a.zj*b.y - a.zij*b.yi - a.zE*b.yjE - a.ziE*b.yijE + a.zjE*b.yE + a.zijE*b.yiE
	 + a.xz*b.xyj - a.xzi*b.xyij + a.xzj*b.xy + a.xzij*b.xyi + a.xzE*b.xyjE + a.xziE*b.xyijE - a.xzjE*b.xyE - a.xzijE*b.xyiE
	 + a.yz*b.j - a.yzi*b.ij + a.yzj*b.q + a.yzij*b.i + a.yzE*b.jE + a.yziE*b.ijE - a.yzjE*b.E - a.yzijE*b.iE
	 + a.xyz*b.xj - a.xyzi*b.xij + a.xyzj*b.x + a.xyzij*b.xi + a.xyzE*b.xjE + a.xyziE*b.xijE - a.xyzjE*b.xE - a.xyzijE*b.xiE ;

c.yzij =  + a.q*b.yzij + a.i*b.yzj - a.j*b.yzi + a.ij*b.yz + a.E*b.yzijE - a.iE*b.yzjE + a.jE*b.yziE - a.ijE*b.yzE
	 + a.x*b.xyzij + a.xi*b.xyzj - a.xj*b.xyzi + a.xij*b.xyz + a.xE*b.xyzijE - a.xiE*b.xyzjE + a.xjE*b.xyziE - a.xijE*b.xyzE
	 + a.y*b.zij + a.yi*b.zj - a.yj*b.zi + a.yij*b.z + a.yE*b.zijE - a.yiE*b.zjE + a.yjE*b.ziE - a.yijE*b.zE
	 - a.xy*b.xzij - a.xyi*b.xzj + a.xyj*b.xzi - a.xyij*b.xz - a.xyE*b.xzijE + a.xyiE*b.xzjE - a.xyjE*b.xziE + a.xyijE*b.xzE
	 - a.z*b.yij - a.zi*b.yj + a.zj*b.yi - a.zij*b.y - a.zE*b.yijE + a.ziE*b.yjE - a.zjE*b.yiE + a.zijE*b.yE
	 + a.xz*b.xyij + a.xzi*b.xyj - a.xzj*b.xyi + a.xzij*b.xy + a.xzE*b.xyijE - a.xziE*b.xyjE + a.xzjE*b.xyiE - a.xzijE*b.xyE
	 + a.yz*b.ij + a.yzi*b.j - a.yzj*b.i + a.yzij*b.q + a.yzE*b.ijE - a.yziE*b.jE + a.yzjE*b.iE - a.yzijE*b.E
	 + a.xyz*b.xij + a.xyzi*b.xj - a.xyzj*b.xi + a.xyzij*b.x + a.xyzE*b.xijE - a.xyziE*b.xjE + a.xyzjE*b.xiE - a.xyzijE*b.xE ;

c.yzE =  + a.q*b.yzE - a.i*b.yziE - a.j*b.yzjE - a.ij*b.yzijE + a.E*b.yz + a.iE*b.yzi + a.jE*b.yzj + a.ijE*b.yzij
	 + a.x*b.xyzE - a.xi*b.xyziE - a.xj*b.xyzjE - a.xij*b.xyzijE + a.xE*b.xyz + a.xiE*b.xyzi + a.xjE*b.xyzj + a.xijE*b.xyzij
	 + a.y*b.zE - a.yi*b.ziE - a.yj*b.zjE - a.yij*b.zijE + a.yE*b.z + a.yiE*b.zi + a.yjE*b.zj + a.yijE*b.zij
	 - a.xy*b.xzE + a.xyi*b.xziE + a.xyj*b.xzjE + a.xyij*b.xzijE - a.xyE*b.xz - a.xyiE*b.xzi - a.xyjE*b.xzj - a.xyijE*b.xzij
	 - a.z*b.yE + a.zi*b.yiE + a.zj*b.yjE + a.zij*b.yijE - a.zE*b.y - a.ziE*b.yi - a.zjE*b.yj - a.zijE*b.yij
	 + a.xz*b.xyE - a.xzi*b.xyiE - a.xzj*b.xyjE - a.xzij*b.xyijE + a.xzE*b.xy + a.xziE*b.xyi + a.xzjE*b.xyj + a.xzijE*b.xyij
	 + a.yz*b.E - a.yzi*b.iE - a.yzj*b.jE - a.yzij*b.ijE + a.yzE*b.q + a.yziE*b.i + a.yzjE*b.j + a.yzijE*b.ij
	 + a.xyz*b.xE - a.xyzi*b.xiE - a.xyzj*b.xjE - a.xyzij*b.xijE + a.xyzE*b.x + a.xyziE*b.xi + a.xyzjE*b.xj + a.xyzijE*b.xij ;

c.yziE =  + a.q*b.yziE + a.i*b.yzE - a.j*b.yzijE + a.ij*b.yzjE - a.E*b.yzi + a.iE*b.yz - a.jE*b.yzij + a.ijE*b.yzj
	 + a.x*b.xyziE + a.xi*b.xyzE - a.xj*b.xyzijE + a.xij*b.xyzjE - a.xE*b.xyzi + a.xiE*b.xyz - a.xjE*b.xyzij + a.xijE*b.xyzj
	 + a.y*b.ziE + a.yi*b.zE - a.yj*b.zijE + a.yij*b.zjE - a.yE*b.zi + a.yiE*b.z - a.yjE*b.zij + a.yijE*b.zj
	 - a.xy*b.xziE - a.xyi*b.xzE + a.xyj*b.xzijE - a.xyij*b.xzjE + a.xyE*b.xzi - a.xyiE*b.xz + a.xyjE*b.xzij - a.xyijE*b.xzj
	 - a.z*b.yiE - a.zi*b.yE + a.zj*b.yijE - a.zij*b.yjE + a.zE*b.yi - a.ziE*b.y + a.zjE*b.yij - a.zijE*b.yj
	 + a.xz*b.xyiE + a.xzi*b.xyE - a.xzj*b.xyijE + a.xzij*b.xyjE - a.xzE*b.xyi + a.xziE*b.xy - a.xzjE*b.xyij + a.xzijE*b.xyj
	 + a.yz*b.iE + a.yzi*b.E - a.yzj*b.ijE + a.yzij*b.jE - a.yzE*b.i + a.yziE*b.q - a.yzjE*b.ij + a.yzijE*b.j
	 + a.xyz*b.xiE + a.xyzi*b.xE - a.xyzj*b.xijE + a.xyzij*b.xjE - a.xyzE*b.xi + a.xyziE*b.x - a.xyzjE*b.xij + a.xyzijE*b.xj ;

c.yzjE =  + a.q*b.yzjE + a.i*b.yzijE + a.j*b.yzE - a.ij*b.yziE - a.E*b.yzj + a.iE*b.yzij + a.jE*b.yz - a.ijE*b.yzi
	 + a.x*b.xyzjE + a.xi*b.xyzijE + a.xj*b.xyzE - a.xij*b.xyziE - a.xE*b.xyzj + a.xiE*b.xyzij + a.xjE*b.xyz - a.xijE*b.xyzi
	 + a.y*b.zjE + a.yi*b.zijE + a.yj*b.zE - a.yij*b.ziE - a.yE*b.zj + a.yiE*b.zij + a.yjE*b.z - a.yijE*b.zi
	 - a.xy*b.xzjE - a.xyi*b.xzijE - a.xyj*b.xzE + a.xyij*b.xziE + a.xyE*b.xzj - a.xyiE*b.xzij - a.xyjE*b.xz + a.xyijE*b.xzi
	 - a.z*b.yjE - a.zi*b.yijE - a.zj*b.yE + a.zij*b.yiE + a.zE*b.yj - a.ziE*b.yij - a.zjE*b.y + a.zijE*b.yi
	 + a.xz*b.xyjE + a.xzi*b.xyijE + a.xzj*b.xyE - a.xzij*b.xyiE - a.xzE*b.xyj + a.xziE*b.xyij + a.xzjE*b.xy - a.xzijE*b.xyi
	 + a.yz*b.jE + a.yzi*b.ijE + a.yzj*b.E - a.yzij*b.iE - a.yzE*b.j + a.yziE*b.ij + a.yzjE*b.q - a.yzijE*b.i
	 + a.xyz*b.xjE + a.xyzi*b.xijE + a.xyzj*b.xE - a.xyzij*b.xiE - a.xyzE*b.xj + a.xyziE*b.xij + a.xyzjE*b.x - a.xyzijE*b.xi ;

c.yzijE =  + a.q*b.yzijE - a.i*b.yzjE + a.j*b.yziE + a.ij*b.yzE - a.E*b.yzij - a.iE*b.yzj + a.jE*b.yzi + a.ijE*b.yz
	 + a.x*b.xyzijE - a.xi*b.xyzjE + a.xj*b.xyziE + a.xij*b.xyzE - a.xE*b.xyzij - a.xiE*b.xyzj + a.xjE*b.xyzi + a.xijE*b.xyz
	 + a.y*b.zijE - a.yi*b.zjE + a.yj*b.ziE + a.yij*b.zE - a.yE*b.zij - a.yiE*b.zj + a.yjE*b.zi + a.yijE*b.z
	 - a.xy*b.xzijE + a.xyi*b.xzjE - a.xyj*b.xziE - a.xyij*b.xzE + a.xyE*b.xzij + a.xyiE*b.xzj - a.xyjE*b.xzi - a.xyijE*b.xz
	 - a.z*b.yijE + a.zi*b.yjE - a.zj*b.yiE - a.zij*b.yE + a.zE*b.yij + a.ziE*b.yj - a.zjE*b.yi - a.zijE*b.y
	 + a.xz*b.xyijE - a.xzi*b.xyjE + a.xzj*b.xyiE + a.xzij*b.xyE - a.xzE*b.xyij - a.xziE*b.xyj + a.xzjE*b.xyi + a.xzijE*b.xy
	 + a.yz*b.ijE - a.yzi*b.jE + a.yzj*b.iE + a.yzij*b.E - a.yzE*b.ij - a.yziE*b.j + a.yzjE*b.i + a.yzijE*b.q
	 + a.xyz*b.xijE - a.xyzi*b.xjE + a.xyzj*b.xiE + a.xyzij*b.xE - a.xyzE*b.xij - a.xyziE*b.xj + a.xyzjE*b.xi + a.xyzijE*b.x ;

c.xyz =  + a.q*b.xyz - a.i*b.xyzi - a.j*b.xyzj - a.ij*b.xyzij - a.E*b.xyzE - a.iE*b.xyziE - a.jE*b.xyzjE - a.ijE*b.xyzijE
	 + a.x*b.yz - a.xi*b.yzi - a.xj*b.yzj - a.xij*b.yzij - a.xE*b.yzE - a.xiE*b.yziE - a.xjE*b.yzjE - a.xijE*b.yzijE
	 - a.y*b.xz + a.yi*b.xzi + a.yj*b.xzj + a.yij*b.xzij + a.yE*b.xzE + a.yiE*b.xziE + a.yjE*b.xzjE + a.yijE*b.xzijE
	 + a.xy*b.z - a.xyi*b.zi - a.xyj*b.zj - a.xyij*b.zij - a.xyE*b.zE - a.xyiE*b.ziE - a.xyjE*b.zjE - a.xyijE*b.zijE
	 + a.z*b.xy - a.zi*b.xyi - a.zj*b.xyj - a.zij*b.xyij - a.zE*b.xyE - a.ziE*b.xyiE - a.zjE*b.xyjE - a.zijE*b.xyijE
	 - a.xz*b.y + a.xzi*b.yi + a.xzj*b.yj + a.xzij*b.yij + a.xzE*b.yE + a.xziE*b.yiE + a.xzjE*b.yjE + a.xzijE*b.yijE
	 + a.yz*b.x - a.yzi*b.xi - a.yzj*b.xj - a.yzij*b.xij - a.yzE*b.xE - a.yziE*b.xiE - a.yzjE*b.xjE - a.yzijE*b.xijE
	 + a.xyz*b.q - a.xyzi*b.i - a.xyzj*b.j - a.xyzij*b.ij - a.xyzE*b.E - a.xyziE*b.iE - a.xyzjE*b.jE - a.xyzijE*b.ijE ;

c.xyzi =  + a.q*b.xyzi + a.i*b.xyz + a.j*b.xyzij - a.ij*b.xyzj + a.E*b.xyziE - a.iE*b.xyzE - a.jE*b.xyzijE + a.ijE*b.xyzjE
	 + a.x*b.yzi + a.xi*b.yz + a.xj*b.yzij - a.xij*b.yzj + a.xE*b.yziE - a.xiE*b.yzE - a.xjE*b.yzijE + a.xijE*b.yzjE
	 - a.y*b.xzi - a.yi*b.xz - a.yj*b.xzij + a.yij*b.xzj - a.yE*b.xziE + a.yiE*b.xzE + a.yjE*b.xzijE - a.yijE*b.xzjE
	 + a.xy*b.zi + a.xyi*b.z + a.xyj*b.zij - a.xyij*b.zj + a.xyE*b.ziE - a.xyiE*b.zE - a.xyjE*b.zijE + a.xyijE*b.zjE
	 + a.z*b.xyi + a.zi*b.xy + a.zj*b.xyij - a.zij*b.xyj + a.zE*b.xyiE - a.ziE*b.xyE - a.zjE*b.xyijE + a.zijE*b.xyjE
	 - a.xz*b.yi - a.xzi*b.y - a.xzj*b.yij + a.xzij*b.yj - a.xzE*b.yiE + a.xziE*b.yE + a.xzjE*b.yijE - a.xzijE*b.yjE
	 + a.yz*b.xi + a.yzi*b.x + a.yzj*b.xij - a.yzij*b.xj + a.yzE*b.xiE - a.yziE*b.xE - a.yzjE*b.xijE + a.yzijE*b.xjE
	 + a.xyz*b.i + a.xyzi*b.q + a.xyzj*b.ij - a.xyzij*b.j + a.xyzE*b.iE - a.xyziE*b.E - a.xyzjE*b.ijE + a.xyzijE*b.jE ;

c.xyzj =  + a.q*b.xyzj - a.i*b.xyzij + a.j*b.xyz + a.ij*b.xyzi + a.E*b.xyzjE + a.iE*b.xyzijE - a.jE*b.xyzE - a.ijE*b.xyziE
	 + a.x*b.yzj - a.xi*b.yzij + a.xj*b.yz + a.xij*b.yzi + a.xE*b.yzjE + a.xiE*b.yzijE - a.xjE*b.yzE - a.xijE*b.yziE
	 - a.y*b.xzj + a.yi*b.xzij - a.yj*b.xz - a.yij*b.xzi - a.yE*b.xzjE - a.yiE*b.xzijE + a.yjE*b.xzE + a.yijE*b.xziE
	 + a.xy*b.zj - a.xyi*b.zij + a.xyj*b.z + a.xyij*b.zi + a.xyE*b.zjE + a.xyiE*b.zijE - a.xyjE*b.zE - a.xyijE*b.ziE
	 + a.z*b.xyj - a.zi*b.xyij + a.zj*b.xy + a.zij*b.xyi + a.zE*b.xyjE + a.ziE*b.xyijE - a.zjE*b.xyE - a.zijE*b.xyiE
	 - a.xz*b.yj + a.xzi*b.yij - a.xzj*b.y - a.xzij*b.yi - a.xzE*b.yjE - a.xziE*b.yijE + a.xzjE*b.yE + a.xzijE*b.yiE
	 + a.yz*b.xj - a.yzi*b.xij + a.yzj*b.x + a.yzij*b.xi + a.yzE*b.xjE + a.yziE*b.xijE - a.yzjE*b.xE - a.yzijE*b.xiE
	 + a.xyz*b.j - a.xyzi*b.ij + a.xyzj*b.q + a.xyzij*b.i + a.xyzE*b.jE + a.xyziE*b.ijE - a.xyzjE*b.E - a.xyzijE*b.iE ;

c.xyzij =  + a.q*b.xyzij + a.i*b.xyzj - a.j*b.xyzi + a.ij*b.xyz + a.E*b.xyzijE - a.iE*b.xyzjE + a.jE*b.xyziE - a.ijE*b.xyzE
	 + a.x*b.yzij + a.xi*b.yzj - a.xj*b.yzi + a.xij*b.yz + a.xE*b.yzijE - a.xiE*b.yzjE + a.xjE*b.yziE - a.xijE*b.yzE
	 - a.y*b.xzij - a.yi*b.xzj + a.yj*b.xzi - a.yij*b.xz - a.yE*b.xzijE + a.yiE*b.xzjE - a.yjE*b.xziE + a.yijE*b.xzE
	 + a.xy*b.zij + a.xyi*b.zj - a.xyj*b.zi + a.xyij*b.z + a.xyE*b.zijE - a.xyiE*b.zjE + a.xyjE*b.ziE - a.xyijE*b.zE
	 + a.z*b.xyij + a.zi*b.xyj - a.zj*b.xyi + a.zij*b.xy + a.zE*b.xyijE - a.ziE*b.xyjE + a.zjE*b.xyiE - a.zijE*b.xyE
	 - a.xz*b.yij - a.xzi*b.yj + a.xzj*b.yi - a.xzij*b.y - a.xzE*b.yijE + a.xziE*b.yjE - a.xzjE*b.yiE + a.xzijE*b.yE
	 + a.yz*b.xij + a.yzi*b.xj - a.yzj*b.xi + a.yzij*b.x + a.yzE*b.xijE - a.yziE*b.xjE + a.yzjE*b.xiE - a.yzijE*b.xE
	 + a.xyz*b.ij + a.xyzi*b.j - a.xyzj*b.i + a.xyzij*b.q + a.xyzE*b.ijE - a.xyziE*b.jE + a.xyzjE*b.iE - a.xyzijE*b.E ;

c.xyzE =  + a.q*b.xyzE - a.i*b.xyziE - a.j*b.xyzjE - a.ij*b.xyzijE + a.E*b.xyz + a.iE*b.xyzi + a.jE*b.xyzj + a.ijE*b.xyzij
	 + a.x*b.yzE - a.xi*b.yziE - a.xj*b.yzjE - a.xij*b.yzijE + a.xE*b.yz + a.xiE*b.yzi + a.xjE*b.yzj + a.xijE*b.yzij
	 - a.y*b.xzE + a.yi*b.xziE + a.yj*b.xzjE + a.yij*b.xzijE - a.yE*b.xz - a.yiE*b.xzi - a.yjE*b.xzj - a.yijE*b.xzij
	 + a.xy*b.zE - a.xyi*b.ziE - a.xyj*b.zjE - a.xyij*b.zijE + a.xyE*b.z + a.xyiE*b.zi + a.xyjE*b.zj + a.xyijE*b.zij
	 + a.z*b.xyE - a.zi*b.xyiE - a.zj*b.xyjE - a.zij*b.xyijE + a.zE*b.xy + a.ziE*b.xyi + a.zjE*b.xyj + a.zijE*b.xyij
	 - a.xz*b.yE + a.xzi*b.yiE + a.xzj*b.yjE + a.xzij*b.yijE - a.xzE*b.y - a.xziE*b.yi - a.xzjE*b.yj - a.xzijE*b.yij
	 + a.yz*b.xE - a.yzi*b.xiE - a.yzj*b.xjE - a.yzij*b.xijE + a.yzE*b.x + a.yziE*b.xi + a.yzjE*b.xj + a.yzijE*b.xij
	 + a.xyz*b.E - a.xyzi*b.iE - a.xyzj*b.jE - a.xyzij*b.ijE + a.xyzE*b.q + a.xyziE*b.i + a.xyzjE*b.j + a.xyzijE*b.ij ;

c.xyziE =  + a.q*b.xyziE + a.i*b.xyzE - a.j*b.xyzijE + a.ij*b.xyzjE - a.E*b.xyzi + a.iE*b.xyz - a.jE*b.xyzij + a.ijE*b.xyzj
	 + a.x*b.yziE + a.xi*b.yzE - a.xj*b.yzijE + a.xij*b.yzjE - a.xE*b.yzi + a.xiE*b.yz - a.xjE*b.yzij + a.xijE*b.yzj
	 - a.y*b.xziE - a.yi*b.xzE + a.yj*b.xzijE - a.yij*b.xzjE + a.yE*b.xzi - a.yiE*b.xz + a.yjE*b.xzij - a.yijE*b.xzj
	 + a.xy*b.ziE + a.xyi*b.zE - a.xyj*b.zijE + a.xyij*b.zjE - a.xyE*b.zi + a.xyiE*b.z - a.xyjE*b.zij + a.xyijE*b.zj
	 + a.z*b.xyiE + a.zi*b.xyE - a.zj*b.xyijE + a.zij*b.xyjE - a.zE*b.xyi + a.ziE*b.xy - a.zjE*b.xyij + a.zijE*b.xyj
	 - a.xz*b.yiE - a.xzi*b.yE + a.xzj*b.yijE - a.xzij*b.yjE + a.xzE*b.yi - a.xziE*b.y + a.xzjE*b.yij - a.xzijE*b.yj
	 + a.yz*b.xiE + a.yzi*b.xE - a.yzj*b.xijE + a.yzij*b.xjE - a.yzE*b.xi + a.yziE*b.x - a.yzjE*b.xij + a.yzijE*b.xj
	 + a.xyz*b.iE + a.xyzi*b.E - a.xyzj*b.ijE + a.xyzij*b.jE - a.xyzE*b.i + a.xyziE*b.q - a.xyzjE*b.ij + a.xyzijE*b.j ;

c.xyzjE =  + a.q*b.xyzjE + a.i*b.xyzijE + a.j*b.xyzE - a.ij*b.xyziE - a.E*b.xyzj + a.iE*b.xyzij + a.jE*b.xyz - a.ijE*b.xyzi
	 + a.x*b.yzjE + a.xi*b.yzijE + a.xj*b.yzE - a.xij*b.yziE - a.xE*b.yzj + a.xiE*b.yzij + a.xjE*b.yz - a.xijE*b.yzi
	 - a.y*b.xzjE - a.yi*b.xzijE - a.yj*b.xzE + a.yij*b.xziE + a.yE*b.xzj - a.yiE*b.xzij - a.yjE*b.xz + a.yijE*b.xzi
	 + a.xy*b.zjE + a.xyi*b.zijE + a.xyj*b.zE - a.xyij*b.ziE - a.xyE*b.zj + a.xyiE*b.zij + a.xyjE*b.z - a.xyijE*b.zi
	 + a.z*b.xyjE + a.zi*b.xyijE + a.zj*b.xyE - a.zij*b.xyiE - a.zE*b.xyj + a.ziE*b.xyij + a.zjE*b.xy - a.zijE*b.xyi
	 - a.xz*b.yjE - a.xzi*b.yijE - a.xzj*b.yE + a.xzij*b.yiE + a.xzE*b.yj - a.xziE*b.yij - a.xzjE*b.y + a.xzijE*b.yi
	 + a.yz*b.xjE + a.yzi*b.xijE + a.yzj*b.xE - a.yzij*b.xiE - a.yzE*b.xj + a.yziE*b.xij + a.yzjE*b.x - a.yzijE*b.xi
	 + a.xyz*b.jE + a.xyzi*b.ijE + a.xyzj*b.E - a.xyzij*b.iE - a.xyzE*b.j + a.xyziE*b.ij + a.xyzjE*b.q - a.xyzijE*b.i ;

c.xyzijE =  + a.q*b.xyzijE - a.i*b.xyzjE + a.j*b.xyziE + a.ij*b.xyzE - a.E*b.xyzij - a.iE*b.xyzj + a.jE*b.xyzi + a.ijE*b.xyz
	 + a.x*b.yzijE - a.xi*b.yzjE + a.xj*b.yziE + a.xij*b.yzE - a.xE*b.yzij - a.xiE*b.yzj + a.xjE*b.yzi + a.xijE*b.yz
	 - a.y*b.xzijE + a.yi*b.xzjE - a.yj*b.xziE - a.yij*b.xzE + a.yE*b.xzij + a.yiE*b.xzj - a.yjE*b.xzi - a.yijE*b.xz
	 + a.xy*b.zijE - a.xyi*b.zjE + a.xyj*b.ziE + a.xyij*b.zE - a.xyE*b.zij - a.xyiE*b.zj + a.xyjE*b.zi + a.xyijE*b.z
	 + a.z*b.xyijE - a.zi*b.xyjE + a.zj*b.xyiE + a.zij*b.xyE - a.zE*b.xyij - a.ziE*b.xyj + a.zjE*b.xyi + a.zijE*b.xy
	 - a.xz*b.yijE + a.xzi*b.yjE - a.xzj*b.yiE - a.xzij*b.yE + a.xzE*b.yij + a.xziE*b.yj - a.xzjE*b.yi - a.xzijE*b.y
	 + a.yz*b.xijE - a.yzi*b.xjE + a.yzj*b.xiE + a.yzij*b.xE - a.yzE*b.xij - a.yziE*b.xj + a.yzjE*b.xi + a.yzijE*b.x
	 + a.xyz*b.ijE - a.xyzi*b.jE + a.xyzj*b.iE + a.xyzij*b.E - a.xyzE*b.ij - a.xyziE*b.j + a.xyzjE*b.i + a.xyzijE*b.q ;

//////////////////////////////////////////////////////

	return c;
}




