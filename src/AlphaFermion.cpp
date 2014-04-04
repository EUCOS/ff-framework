#include "AlphaFermion.h"

//PUBLIC.
AlphaFermion::AlphaFermion(const std::vector<int>& Numerator, 
			     int Denominator,
			     const std::vector<int>& Coefficients):
  Alpha(Numerator, Denominator, Coefficients)
{
  ;
}//Close constructor.

AlphaFermion::AlphaFermion(const AlphaFermion& New_AlphaFermion):
  Alpha(New_AlphaFermion)
{
  ;
}//Close copy constructor.

//INTERFACE.
char AlphaFermion::Type() const
{
  //f for Fermion.
  return 'f';
}
