#include "GSOCoefficientMatrixBuilder.h"
//CONSTRUCTORS.
GSOCoefficientMatrixBuilder::GSOCoefficientMatrixBuilder
(const GSOCoefficientMatrix& Half_GSO_Matrix, 
 const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  Half_GSO_Matrix_ = Half_GSO_Matrix;
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  GSO_Coefficient_Orders_ = Half_GSO_Matrix.Denominators();
  Consistent_GSO_Matrix_ = true;
}//Close constructor.

GSOCoefficientMatrixBuilder::GSOCoefficientMatrixBuilder
(const Model& FFHS_Model, const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  Half_GSO_Matrix_ = FFHS_Model.k_ij();
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  GSO_Coefficient_Orders_ = FFHS_Model.k_ij().Denominators();
  Consistent_GSO_Matrix_ = true;
}//Close second constructor.

GSOCoefficientMatrixBuilder::GSOCoefficientMatrixBuilder
(const std::vector<std::vector<int> >& Half_Numerators,
 const std::vector<int>& Half_Denominators, 
 const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  Half_GSO_Matrix_ = GSOCoefficientMatrix(Half_Numerators, Half_Denominators);
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  GSO_Coefficient_Orders_ = Half_Denominators;
  Consistent_GSO_Matrix_ = true;
}//Close third constructor.

GSOCoefficientMatrixBuilder::GSOCoefficientMatrixBuilder
(const GSOCoefficientMatrixBuilder& New_GSOCoefficientMatrixBuilder)
{
  Half_GSO_Matrix_ = New_GSOCoefficientMatrixBuilder.Half_GSO_Matrix();
  Common_Basis_Alphas_ = New_GSOCoefficientMatrixBuilder.Common_Basis_Alphas();
  GSO_Coefficient_Orders_ = 
    New_GSOCoefficientMatrixBuilder.GSO_Coefficient_Orders();
  Complete_GSO_Matrix_ = New_GSOCoefficientMatrixBuilder.Complete_GSO_Matrix();
  Consistent_GSO_Matrix_ = 
    New_GSOCoefficientMatrixBuilder.Consistent_GSO_Matrix();
}//Close copy constructor.

//INTERFACE.
void GSOCoefficientMatrixBuilder::Build_Complete_GSO_Matrix()
{
  Convert_Half_GSO_Matrix();
  Complete_Half_GSO_Matrix();
}//Close Build_Complete_GSO_Matrix.

//DEBUG.
void GSOCoefficientMatrixBuilder::Display_Half_GSO_Matrix() const
{
  std::cout<<"Half GSO Coefficient Matrix: "<<std::endl;
  Half_GSO_Matrix().Display();
}//Close Display_Half_GSO_Matrix.

void GSOCoefficientMatrixBuilder::Display_Common_Basis_Alphas() const
{
  std::cout<<"Common Basis Alphas: "<<std::endl;
  for(int a=0; a<static_cast<int>(Common_Basis_Alphas().size()); a++)
    Common_Basis_Alphas().at(a).Display();
  std::cout<<std::endl;
}//Close Display_Common_Basis_Alphas.

void GSOCoefficientMatrixBuilder::Display_Complete_GSO_Matrix() const
{
  std::cout<<"Complete GSO Coefficient Matrix: "<<std::endl;
  Complete_GSO_Matrix().Display();
}//Close Display_Complete_GSO_Matrix.

//PRIVATE.
//HELPERS.
void GSOCoefficientMatrixBuilder::Convert_Half_GSO_Matrix()
{
  int CommonDenom = Common_Basis_Alphas().at(0).Denominator();
  std::vector<std::vector<int> > Half_GSO_Numerators = 
    Half_GSO_Matrix().Numerators();

  for(int a=0; a<static_cast<int>(GSO_Coefficient_Orders().size())-1; ++a)
    {
      int Column_Loop_Start = a;
      if(a > 0)
      	Column_Loop_Start++;
      int Convert = CommonDenom / GSO_Coefficient_Orders().at(a);
      for(int b=Column_Loop_Start; b<static_cast<int>(Half_GSO_Numerators.size()); ++b)
	{
	  Half_GSO_Numerators.at(b).at(a) =  Reset_GSO_Coefficient_Range
	    (2*Half_GSO_Numerators.at(b).at(a)*Convert, CommonDenom);
	  //GSO Coefficients are coded 2m/N.
	}

    }//Close for loop on GSO_Coefficient_Orders().size()-1.

  Half_GSO_Matrix_ = GSOCoefficientMatrix
    (Half_GSO_Numerators, std::vector<int>
     (Half_GSO_Numerators.size(), CommonDenom));

}//Close Convert_Half_GSO_Matrix. 

void GSOCoefficientMatrixBuilder::Complete_Half_GSO_Matrix()
{
  std::vector<std::vector<int> > Complete_GSO_Numerators;
  for(int a=0; a<static_cast<int>(Half_GSO_Matrix().Numerators().size()); a++)
    {
      std::vector<int> Complete_GSO_Numerator;
      for(int b=0; b<static_cast<int>(Half_GSO_Matrix().Numerators().size()); b++)
	{
	  if((b<a)||(b==0))
	    Complete_GSO_Numerator.push_back
	      (Half_GSO_Matrix().Numerators().at(a).at(b));
	  else if(b==a)//Build diagonal element.
	    {
	      int Converted_GSO_Coefficient = Compute_Diagonal_GSO_Element
		(Common_Basis_Alphas().at(a), 
		 Half_GSO_Matrix().Numerators().at(a).at(0));

	      Converted_GSO_Coefficient = Reset_GSO_Coefficient_Range
		(Converted_GSO_Coefficient, 
		 Common_Basis_Alphas().at(a).Denominator());
		
	      Complete_GSO_Numerator.push_back(Converted_GSO_Coefficient);

	      if(Consistent_GSO_Matrix())
		Consistent_GSO_Matrix_ = 
		  Check_GSO_Element_Consistency(Converted_GSO_Coefficient, b);
	    }else//Build upper triangular element.
	    {
	      int Converted_GSO_Coefficient = Compute_Off_Diagonal_GSO_Element
		(Common_Basis_Alphas().at(a), Common_Basis_Alphas().at(b),
		 Half_GSO_Matrix().Numerators().at(b).at(a));

	      Converted_GSO_Coefficient = Reset_GSO_Coefficient_Range
		(Converted_GSO_Coefficient, 
		 Common_Basis_Alphas().at(a).Denominator());
		
	      Complete_GSO_Numerator.push_back(Converted_GSO_Coefficient);
	      if(Consistent_GSO_Matrix())
		Consistent_GSO_Matrix_ = 
		  Check_GSO_Element_Consistency(Converted_GSO_Coefficient, b);
	    }
	}//Close for loop on Half_GSO_Matrix().Numerators().at(a).size().
      Complete_GSO_Numerators.push_back(Complete_GSO_Numerator);
      Complete_GSO_Numerator.clear();
    }//Close for loop on Half_GSO_Matrix().Numerators.size().
  Complete_GSO_Matrix_ = GSOCoefficientMatrix 
    (Complete_GSO_Numerators, std::vector<int> 
     (Complete_GSO_Numerators.size(), Common_Basis_Alphas().at(0).Denominator()));
}//Close Complete_Half_GSO_Matrix.

int GSOCoefficientMatrixBuilder::Compute_Diagonal_GSO_Element
(const BasisAlpha& The_Alpha, int First_Row_GSO_Coefficient)
{
  int Converted_GSO_Denominator = The_Alpha.Denominator();
  int Converted_GSO_Coefficient = The_Alpha.Lorentz_Dot(The_Alpha);

  Converted_GSO_Coefficient /= (8*Converted_GSO_Denominator);
  //Always an integer by modular invariance.

  Converted_GSO_Coefficient -= The_Alpha.Numerator().at(0);

  Converted_GSO_Coefficient -= First_Row_GSO_Coefficient;

  Converted_GSO_Coefficient %= (2*Converted_GSO_Denominator);

  return Converted_GSO_Coefficient;
}//Close Compute_Diagonal_GSO_Element.

int GSOCoefficientMatrixBuilder::Compute_Off_Diagonal_GSO_Element
(const BasisAlpha& Alpha_a, const BasisAlpha& Alpha_b, 
 int Half_GSO_Coefficient)
{
  int Converted_GSO_Denominator = Alpha_a.Denominator();

  int Converted_GSO_Coefficient = Alpha_a.Lorentz_Dot(Alpha_b);

  Converted_GSO_Coefficient /= (4*Converted_GSO_Denominator);
  //Always an integer by modular invariance.

  Converted_GSO_Coefficient -= Half_GSO_Coefficient;

  Converted_GSO_Coefficient %= (2*Converted_GSO_Denominator);

  return Converted_GSO_Coefficient;
}//Close Compute_Off_Diagonal_GSO_Element.

int GSOCoefficientMatrixBuilder::Reset_GSO_Coefficient_Range
(int Converted_GSO_Coefficient, int Converted_GSO_Denominator)
{
  
  if(Converted_GSO_Coefficient == -Converted_GSO_Denominator)
    Converted_GSO_Coefficient *= -1;
  else if(Converted_GSO_Coefficient < -Converted_GSO_Denominator)
    Converted_GSO_Coefficient += (2*Converted_GSO_Denominator);
  else if(Converted_GSO_Coefficient > Converted_GSO_Denominator)
    Converted_GSO_Coefficient -= (2*Converted_GSO_Denominator);

  return Converted_GSO_Coefficient;
}//Close Reset_GSO_Coefficient_Range.

bool GSOCoefficientMatrixBuilder::Check_GSO_Element_Consistency(int GSO_Element,
								   int column)
{
  if(Common_Basis_Alphas().at(column).Numerator().at(0) == 
     Common_Basis_Alphas().at(column).Denominator())//ST Fermion
    return (FF::LCM(2,GSO_Coefficient_Orders().at(column))*GSO_Element)%
      (2*Common_Basis_Alphas().at(0).Denominator())==0;
  else //ST Boson.
    return(GSO_Coefficient_Orders().at(column)*GSO_Element)%
      (2*Common_Basis_Alphas().at(0).Denominator())==0;
}//Close Check_GSO_Element_Consistency.
