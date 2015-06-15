#include "Alpha.h"

Alpha::Alpha(const std::vector<int>& Numerator, int Denominator, 
	     const std::vector<int>& Coefficients):
  BasisAlpha(Numerator, Denominator)
{
  Coefficients_ = Coefficients;
}//Close constructor.

Alpha::Alpha(const Alpha& New_Alpha):BasisAlpha(New_Alpha)
{
  Coefficients_ = New_Alpha.Coefficients();
}//Close copy constructor.

//INTERFACE.
//just a test of push
int Alpha::Mass_Left()
{
  int Mass_Left = 0;
  for(int a=0; a<LM_Size(); a++)
    Mass_Left+=Numerator().at(a)*Numerator().at(a);

  return Mass_Left;
}//Close Mass_Left.

int Alpha::Mass_Right()
{
  int Mass_Right = 0;

  for(int a=LM_Size(); a<static_cast<int>(Numerator().size()); a++)
    Mass_Right+=Numerator().at(a)*Numerator().at(a);

  return Mass_Right;
}//Close Mass_Right.

char Alpha::Type() const
{
  //n for Normal.
  return 'n';
}

//DEBUG.
void Alpha::Display_Coefficients() const
{
  std::cout<<"Coefficients: ";
  for(int a=0; a<static_cast<int>(Coefficients().size()); a++)
    std::cout<<Coefficients().at(a)<<" ";
}//Close Display_Coefficients.
