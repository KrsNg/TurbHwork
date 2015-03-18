#include <stdio.h>
#include <iostream>
#include <stddef.h>
#include <string>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <stdlib.h>
#include <sstream>

using std::cout;
using std::endl;
using std::ios;


// Function prototypes====================================================================================================
double evalFunc(double x);

//=================================================================================================

// Main function
int main(){
  

  double tol = 1e-3; 
  double c = 0;
  double a= 0;
  double b = 10;  
  double f_a = 0, f_b = 0, f_c = 0, root = 0;
          
//   f_a = evalFunc(U_ave, Re_D, a);
//   f_b = evalFunc(U_ave, Re_D, b);

  /* Check that that neither end-point is a root and if f(a) and f(b) have the same sign, throw an exception. */
  if ( f_a == 0 ){
    root = a;
    std::cout<<"root is "<<root<<endl;
  } else if ( f_b == 0 ){
     root = b;
     std::cout<<"root is "<<root<<endl;
  } else if ( f_a * f_b > 0 ){
     std::cout<<"f(a) and f(b) do not have opposite signs"<<endl;
  }
  
  int iter = 0;

  // Begin the iterations-------------------------------------
	
//    while ( fabs(b - a) > tol ){
     for (int i = 0; i< 20; i++){
//         a = 700; b = 0;                                     //  just for testing
        c = 0.5*(a + b);
               
	f_a = evalFunc(a);
	f_b = evalFunc(b);	
	f_c = evalFunc(c);
	
	 std::cout<<" a is = "<<a<<endl;
	 std::cout<<"  b is = "<<b<<endl;
	 std::cout<<"   c is = "<<c<<endl;
	 std::cout<<" f(a) is = "<<f_a<<endl;
	 std::cout<<"  f(b) is = "<<f_b<<endl;
	 std::cout<<"   f(c) is = "<<f_c<<endl;
	 
	/* Check if we found a root or whether or not  we should continue with:
                   [a, c] if f(a) and f(c) have opposite signs, or
                   [c, b] if f(c) and f(b) have opposite signs. */
	 
	if (f_c == 0){
// 	  break;
	} else if (f_a * f_c < 0){
	  b = c;
	} else{
	  a = c;
	}
	
	/* If |b - a| < eps_step, check whether or not
               |f(a)| < |f(b)| and |f(a)| < eps_abs and return 'a', or
               |f(b)| < eps_abs and return 'b'. */

        if ( b - a < tol ){
            if ( fabs(f_a) < fabs(f_b) && fabs(f_a) < tol ){
                root = a;
		std::cout<<"root is "<<root<<endl;
		break;
	    }
            else if ( fabs(f_b) < tol ){
                root = b;
		std::cout<<"root is "<<root<<endl;
		break;
	    }
	}
	
	iter = iter + 1;
// 	std::cout<<"                  root = "<<root<<endl;

   }
   
   	std::cout<<"                  Final results root = "<<root<<endl;
	std::cout<<"Num iter = "<<iter<<endl;
   
return 0;

}

double evalFunc(double x){
  
        double func = 0;
	
	func = 3*(x*x*x*x) - 2*(x*x) - 7;
	
	return func;
  
}

