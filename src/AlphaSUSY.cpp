#include "AlphaSUSY.h"

//PUBLIC.
AlphaSUSY::AlphaSUSY(const std::vector<int>& Numerator,
		       int Denominator,
		       const std::vector<int>& Coefficients):
  Alpha(Numerator, Denominator, Coefficients)
{
  ;
}//Close constructor.

AlphaSUSY::AlphaSUSY(const AlphaSUSY& New_AlphaSUSY):
  Alpha(New_AlphaSUSY)
{
  ;
}//Close copy constructor.

//INTERFACE.
char AlphaSUSY::Type() const
{
  //s for SUSY.
  return's';
}
