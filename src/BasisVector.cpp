#include "BasisVector.h"

//CONSTRUCTORS.
BasisVector::BasisVector(const std::vector<int>& BV, int Order)
{
  BV_ = BV;
  Order_ = Order;

  Calculate_LM_Size();
  Calculate_RM_Compact_Size();
}//Close constructor.

BasisVector::BasisVector(const std::vector<int>& BV, int Order,
			   int Large_ST_Dimensions)
{
  BV_ = BV;
  Order_ = Order;

  LM_Size_ = 28 - 2*Large_ST_Dimensions;
  RM_Compact_Size_ = 2*(10 - Large_ST_Dimensions);
}//Close second constructor.

BasisVector::BasisVector(const BasisVector& New_Basis_Vector)
{
  BV_ = New_Basis_Vector.BV();
  Order_ = New_Basis_Vector.Order();
  LM_Size_ = New_Basis_Vector.LM_Size();
  RM_Compact_Size_ = New_Basis_Vector.RM_Compact_Size();

}//Close copy constructor.

//DEBUG.
void BasisVector::Display() const
{
  std::cout<<Order()<<": ";
  for(int a=0; a<LM_Size(); a++)
    std::cout<<BV().at(a)<<" ";
  std::cout<<"|| ";
  for(int a=LM_Size(); a<static_cast<int>(BV().size()); a++)
    std::cout<<BV().at(a)<<" ";
  std::cout<<std::endl;
}//Close Display_BV.

//PRIVATE.
//HELPERS.
void BasisVector::Calculate_LM_Size()
{
  LM_Size_ = (BV().size() - 24)/2;
}//Close Calculate_Last_LM_Index.

void BasisVector::Calculate_RM_Compact_Size()
{
  if(LM_Size() == 0)
    Calculate_LM_Size();

  if(BV().size() != 40)
    RM_Compact_Size_ = LM_Size() - 8;
  else
    RM_Compact_Size_ = 0;//Automatically terminates
  //loops that use this as an upper limit.
}//Close Calculate_Last_RM_Compact_Index.

