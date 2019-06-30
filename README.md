# GA_Lib_V1
Numeric Geometric Algebra routines in C, and GiNaC symbolic routines in C++

These are a collection of Geometric Algebra routines for 2D, 3D, 4D, Minkowski 3+1 spacetime and 4+1 fivespace.

The C routines are numeric. The C++ routines are symbolic, using the GiNaC library.

*****************

Prequisites for numerical C routines: libm and headers

Prequisites for symbolic C++ routines: We need libm, libgmp, libcln, and libginac and associated headers installed.

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
    
Compile instruction example, assume gcc compiler for numeric C work

gcc GA2E_Demo.c -o GA2E_Demo GA2E_Demo.c -l m

Compile command example, assuming g++ compiler for symbolic C++ work

g++ GA2E_Test_Harness.cp -o GA2E_Test_Harness -l ginac -l cln -l gmp -l m

