#include "GSOCoefficientMatrix.h"

//PUBLIC.
GSOCoefficientMatrix::GSOCoefficientMatrix
(const std::vector<std::vector<int> >& Numerators, 
 const std::vector<int>& Denominators)
{
  Numerators_ = Numerators;
  Denominators_ = Denominators;
}//Close second constructor.

GSOCoefficientMatrix::GSOCoefficientMatrix(const GSOCoefficientMatrix& 
					       New_GSOCoefficientMatrix)
{
  Numerators_ = New_GSOCoefficientMatrix.Numerators();
  Denominators_ = New_GSOCoefficientMatrix.Denominators();
}//Close copy constructor.

//INTERFACE.
void GSOCoefficientMatrix::Load_GSOCoefficientMatrix_Row
(const std::vector<int>& New_Row)
{
  Numerators_.push_back(New_Row);
}//Close Load_GSOCoefficientMatrix_Row.

void GSOCoefficientMatrix::Load_GSOCoefficientMatrix_Order (int New_Order)
{
  Denominators_.push_back(New_Order);
}//Close Load_GSOCoefficientMatrix_Order.

//DEBUG.
void GSOCoefficientMatrix::Display() const
{
  //Put the denominators on top, since they correspond to the column.
  for(int a=0; a<static_cast<int>(Denominators().size()); ++a)
    std::cout<<Denominators().at(a)<<" ";
  std::cout<<std::endl;

  for(int a=0; a<static_cast<int>(Numerators().size()); ++a)
    std::cout<<"--";
  std::cout<<std::endl;

  //Now the numerators.
  for(int a=0; a<static_cast<int>(Numerators().size()); a++)
    {
      for(int b=0; b<static_cast<int>(Numerators().at(a).size()); b++)
	std::cout<<Numerators().at(a).at(b)<< " ";
      std::cout<<std::endl;
    }//Close for loop on Numerators.
  std::cout<<std::endl;
}//Close Display.

