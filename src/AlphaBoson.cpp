#include "AlphaBoson.h"

//PUBLIC.
AlphaBoson::AlphaBoson(const std::vector<int>& Numerator, 
			 int Denominator, 
			 const std::vector<int>& Coefficients):
  Alpha(Numerator, Denominator, Coefficients)
{
  ;
}//Close constructor.

AlphaBoson::AlphaBoson(const AlphaBoson& New_AlphaBoson):Alpha(New_AlphaBoson)
{
  ;
}//Close copy constructor.


//INTERFACE.
char AlphaBoson::Type() const
{
  //b for Boson.
  return 'b';
}//Close Type.
