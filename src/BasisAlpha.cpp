#include "BasisAlpha.h"
//CONSTRUCTORS.
BasisAlpha::BasisAlpha(const std::vector<int>& Numerator, 
			 int Denominator)
{
  Numerator_ = Numerator;
  Denominator_ = Denominator;

  Calculate_LM_Size();
  Calculate_RM_Compact_Size();
}// Close constructor.

BasisAlpha::BasisAlpha(const std::vector<int>& Numerator, 
			 int Denominator, const BasisVector& BV)
{
  Numerator_ = Numerator;
  Denominator_ = Denominator;

  LM_Size_ = BV.LM_Size();
  RM_Compact_Size_ = BV.RM_Compact_Size();
}//Close second constructor.

BasisAlpha::BasisAlpha(const std::vector<int>& Numerator, 
			 int Denominator, int Large_ST_Dimensions)
{
  Numerator_ = Numerator;
  Denominator_ = Denominator;

  LM_Size_ = 28 - 2*Large_ST_Dimensions;
  RM_Compact_Size_ = 2*(10 - Large_ST_Dimensions);
}//Close third constructor.

BasisAlpha::BasisAlpha(const BasisAlpha& Old_BasisAlpha)
{
  Numerator_ = Old_BasisAlpha.Numerator();
  Denominator_ = Old_BasisAlpha.Denominator();

  LM_Size_ = Old_BasisAlpha.LM_Size();
  RM_Compact_Size_ = Old_BasisAlpha.RM_Compact_Size();
}//Close copy constructor.

//INTERFACE.
int BasisAlpha::Lorentz_Dot(const BasisAlpha& BasisAlpha_2) const
{
  int LM = 0;
  int RM = 0;

  for(int a=0; a<LM_Size(); a++)
    LM += Numerator().at(a)*BasisAlpha_2.Numerator().at(a);

  for(int a=LM_Size(); a<static_cast<int>(Numerator().size()); a++)
    RM += Numerator().at(a)*BasisAlpha_2.Numerator().at(a);

  return LM - RM;

}//Close Lorentz_Dot.

//DEBUG.
void BasisAlpha::Display() const
{

  std::cout<<Denominator()<<": ";
  for(int a=0; a<LM_Size(); a++)
    std::cout<<Numerator().at(a)<<" ";
  std::cout<<"|| ";
  for(int a=LM_Size(); a<static_cast<int>(Numerator().size()); a++)
    std::cout<<Numerator().at(a)<<" ";
  std::cout<<std::endl;
}// Close Display_Numerator.

//PRIVATE.
//HELPERS.
void BasisAlpha::Calculate_LM_Size()
{
  LM_Size_ = (Numerator().size() - 24)/2;
}//Close Calculate_LM_Size.

void BasisAlpha::Calculate_RM_Compact_Size()
{
  if(LM_Size() == 0)
    Calculate_LM_Size();

  if(Numerator().size() != 40)
    RM_Compact_Size_ = LM_Size() - 8;
  else
    RM_Compact_Size_ = 0;
}//Close Calculate_RM_Compact_Size.
