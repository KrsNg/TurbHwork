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
double evalFuncDeriv(double x);
//=================================================================================================

// Main function
int main(){
  
  double tol = 1e-5; 
  double c = 0;
  double a= 0;
  double b = 1000;  
  double f_a = 0, f_b = 0, f_c = 0, root = 0;

  f_a = evalFunc(a);
  f_b = evalFunc(b);

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
	
   while ( fabs(b - a) > tol ){

        c = 0.5*(a + b);     
	f_a = evalFunc(a);
	f_b = evalFunc(b);	
	f_c = evalFunc(c);

	if (f_c == 0){
// 	  break;
	} else if (f_a * f_c < 0){
	  b = c;
	} else{
	  a = c;
	}

	if ( b - a < tol ){
            if ( fabs(f_a) < fabs(f_b) && fabs(f_a) < tol ){
                root = a;
		std::cout<<"root is "<<root<<endl;
		break;
	    }   else if ( fabs(f_b) < tol ){
                root = b;
		std::cout<<"root is "<<root<<endl;
		break;
	    }
	}
	iter = iter + 1;
// 	std::cout<<"root is "<<c<<endl;
	root = c;
   }
    
   	std::cout<<"Final results Bisection method root = "<<root<<endl;
	std::cout<<"Num iter = "<<iter<<endl;
	
	iter = 0;
	double y = 0;
       	double delta = 300;
//  for (int i = 0; i<40; i++){
 while (delta > tol){

//    y = b - evalFunc(b)/evalFuncDeriv(b); // update to b
   
   double D = (evalFunc(b + evalFunc(b)) - evalFunc(b)) /evalFunc(b);
   y = b - evalFunc(b)/D;
   
   std::cout<<"root is "<<y<<endl;
   delta = fabs(y-b);
   if (fabs(y - b) < tol){
    break;
   }
   b = y;
   root = y;
   iter  = iter + 1;
   
 }
	
	std::cout<<"Final results Newton method root = "<<root<<endl;
	std::cout<<"Num iter = "<<iter<<endl;

return 0;

}

double evalFunc(double x){
  
        double func = 0;
	
	func = 3*(x*x*x*x) - 2*(x*x) - 7;
	
	return func;
  
}

double evalFuncDeriv(double x){
  
     double deriv = 0;
     
     deriv = 12*x*x*x - 4*x;
     
     return deriv;
  
}

