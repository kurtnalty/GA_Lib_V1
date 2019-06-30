# GA_Lib_V1
Numeric Geometric Algebra routines in C, and GiNaC symbolic routines in C++

These are a collection of Geometric Algebra routines for 2D, 3D, 4D, Minkowski 3+1 spacetime and 4+1 fivespace.

The C routines are numeric. The C++ routines are symbolic, using the GiNaC library.

*****************

Prequisites for numerical C routines: libm and headers

Prequisites for symbolic C++ routines: We need libm, libgmp, libcln, and libginac and associated headers installed.

Note: GiNac has different initialization syntax for matrices between older and newer version. Should you get compile errors
in the demo programs, see the note at the end of this file for the easy, included fix.

*****************

Simple usage: 

    We include the actual source code in our main file in this implementation. We do not build a binary library.
    
For the C routines, simply include the desired geometric algebra C routine files. Example snippet from GA2E_Demo.c

    #include <stdio.h>
    #include <stdlib.h>

    #include "GA2E_Routines.c"

    /////////////////////////////////// Demonstrate GA2E Routines ////////////////////////////////

    int main(void)
    {

For the C++ routines, include the ginac releated headers ans the C++ routines. Example snippet from GA2E_Test_Harness.cp
    
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <time.h>
    #include <iostream>
    #include <fstream>
    #include <ginac/ginac.h>
    using namespace std;
    using namespace GiNaC;

    #include "GA2E_Routines.cp"

     /////////////////////////////////// Demonstrate GA2E Routines ////////////////////////////////

    int main(void)
    {
    
Compile instruction example, assuming gcc compiler for numeric C work

gcc GA2E_Demo.c -o GA2E_Demo GA2E_Demo.c -l m

Compile command example, assuming g++ compiler for symbolic C++ work

g++ GA2E_Test_Harness.cp -o GA2E_Test_Harness -l ginac -l cln -l gmp -l m

*****************

Supported Algebras:

    GA2E - Two Dimensional Euclidean Space Geometric Algebra < --- Simplest!
    
    GA3E - Three Dimensional Euclidean Space Geometric Algebra using xz bivector base convention
    
    GA3Ezx - Three Dimensional Euclidean Space Geometric Algebra using zx bivector base convention  < --- Good Starting point for beginners
    
    GA4E - Four Dimensional Euclidean Space Geometric Algebra
    
    Mink - Four Dimensional SpaceTime (Minkowski SpaceTime) using +++- metric for xyzt
    
    GA5_4_1 - Five Dimensional SpaceTime matching the full Dirac Gamma matrix set. 
    
*****************
Thank-you to David Hestenes, Alan MacDonald, Eric Lenyel, Leo Dorst and many others for championing this field!

Special thanks to Alan Bromborsky for the symbolic GA routines in python, which inspired and were used as a check for my routines.

*****************************************

GiNaC matrix initialization routines have changed between older and newest versions of GiNaC.

All libraries have two copies of the affected routines, with new version commented out, and older enabled by default.

Sample code snippet from GA3E_Routines.cp

    //////////////////////////////////////////////////////
    /*
    matrix GA3E_To_Matrix(GA3E MV) // newer ginac calling sequence
    {

Notice the /* block commenting the newer calling sequence routines? Should the demo fail to compile, use your text editor to comment out the older routine, and 
use the newer routine.



