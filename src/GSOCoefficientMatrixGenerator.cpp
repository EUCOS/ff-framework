#include "GSOCoefficientMatrixGenerator.h"
GSOCoefficientMatrixGenerator::GSOCoefficientMatrixGenerator
(std::vector<int> Orders, int Layer)
{
  Orders_ = Orders;
  Layer_ = Layer;
}//Close constructor.

GSOCoefficientMatrixGenerator::GSOCoefficientMatrixGenerator
(const GSOCoefficientMatrixGenerator& New_GSOCoefficientMatrixGenerator)
{
  Orders_ = New_GSOCoefficientMatrixGenerator.Orders();
  Layer_ = New_GSOCoefficientMatrixGenerator.Layer();
  GSO_Coefficient_Matrix_Extensions_ = New_GSOCoefficientMatrixGenerator.
    GSO_Coefficient_Matrix_Extensions();
}//Close copy constructor.

//INTERFACE.
void GSOCoefficientMatrixGenerator::Build_GSO_Coefficient_Extensions()
{
  std::vector<int> Row(Orders().size() - 1, 0);//For the 1 vector.
  Extend_GSO_Coefficient_Matrix(Row, 0);
}//Close Build_GSO_Coefficient_Extensions.

//DEBUG.
void GSOCoefficientMatrixGenerator::
Display_GSO_Coefficient_Matrix_Extensions() const
{
  std::cout<<"GSO Coefficient Matrix Extensions: "<<
    GSO_Coefficient_Matrix_Extensions().size()<<std::endl;

  std::list<std::vector<int> >::const_iterator itGSO_Extensions = 
    GSO_Coefficient_Matrix_Extensions_.begin();
  std::list<std::vector<int> >::const_iterator itGSO_Extensions_End = 
    GSO_Coefficient_Matrix_Extensions_.end();

  for(; itGSO_Extensions != itGSO_Extensions_End; ++itGSO_Extensions)
    {
      for(int a=0; a<static_cast<int>(itGSO_Extensions->size()); ++a)
	std::cout<<itGSO_Extensions->at(a)<<" ";
      std::cout<<std::endl;
    }
  std::cout<<std::endl;
}//Close Display_GSO_Coefficient_Matrix_Extensions.

//PRIVATE.
void GSOCoefficientMatrixGenerator::Extend_GSO_Coefficient_Matrix
(std::vector<int> Row, int Index)
{
  if(Index < static_cast<int>(Row.size()))
    {
      for(int a=0; a<Orders().at(Index); ++a)
	{
	  Row.at(Index) = a;
	  Extend_GSO_Coefficient_Matrix(Row, Index+1);
	}//Close for loop on matrix values.
    }else
    GSO_Coefficient_Matrix_Extensions_.push_back(Row);
}//Close Extend_GSO_Coefficient_Matrix.
