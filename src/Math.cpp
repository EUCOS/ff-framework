#include "Math.h"

int FF::GCD(int a, int b)
{
  return (b!=0 ? GCD(b,a%b):a);//Euclid's algorithm.
}//Close GCD.

int FF::LCM(int a, int b)
{
  return ((a*b)/GCD(a,b));
}//Close LCM.

bool FF::abseq(int a, int b)
{
	return (abs(a) == abs(b));
}//Close abseq

