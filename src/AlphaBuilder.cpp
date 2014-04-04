#include "AlphaBuilder.h"

//PUBLIC.
AlphaBuilder::AlphaBuilder(const std::vector<BasisAlpha>& Common_BasisAlphas,
			     const std::vector<BasisAlpha>& BasisAlphas)
{
  Common_BasisAlphas_ = Common_BasisAlphas;
  Coefficient_Limits_ = Get_Coefficient_Limits(BasisAlphas);

  Linearly_Independent_Alphas_ = true;
}//Close constructor.

AlphaBuilder::AlphaBuilder(const AlphaBuilder& New_AlphaBuilder)
{
  Common_BasisAlphas_ = New_AlphaBuilder.Common_BasisAlphas();
  Coefficient_Limits_ = New_AlphaBuilder.Coefficient_Limits();
  Alpha_Bosons_ = New_AlphaBuilder.Alpha_Bosons();
  Alpha_Fermions_ = New_AlphaBuilder.Alpha_Fermions();
  Alpha_SUSYs_ = New_AlphaBuilder.Alpha_SUSYs();
  Linearly_Independent_Alphas_ = New_AlphaBuilder.Linearly_Independent_Alphas();
}//Close copy constructor.

//INTERFACE.
void AlphaBuilder::Build_Alphas()
{
  std::vector<int> Starting_Coefficients
    (Common_BasisAlphas().size(),0);
  std::vector<int> Empty_Numerator
    (Common_BasisAlphas().at(0).Numerator().size(), 0);

  Alpha Empty_Alpha(Empty_Numerator, Common_BasisAlphas().at(0).Denominator(),
		    Starting_Coefficients);

  Add_Alphas(0, Empty_Alpha);
}//Close Build_Linear_Combinations.

//DEBUG.
void AlphaBuilder::Display_Common_BasisAlphas() const
{
  int Common_BasisAlphas_Size = Common_BasisAlphas().size();
  std::cout<<"Common basis alphas: "<<Common_BasisAlphas_Size<<std::endl;
  for(int a=0; a<Common_BasisAlphas_Size; ++a)
    Common_BasisAlphas().at(a).Display();

}//Close Display_Common_BasisAlpha_Set.

void AlphaBuilder::Display_Coefficient_Limits() const
{
  int Coefficient_Limits_Size = Coefficient_Limits().size();
  std::cout<<"Coefficient limits: "<<std::endl;
  for(int a=0; a<Coefficient_Limits_Size; ++a)
    std::cout<<Coefficient_Limits().at(a)<<" ";
  std::cout<<std::endl;
}//Close Display_Coefficient_Limits.

void AlphaBuilder::Display_Alpha_Bosons() const
{
  int Alpha_Bosons_Size = Alpha_Bosons().size();
  std::cout<<"Boson Alphas: "<<Alpha_Bosons_Size<<std::endl;
  std::set<AlphaBoson>::iterator itAlpha_Boson = Alpha_Bosons_.begin();
  std::set<AlphaBoson>::iterator Alpha_Bosons_End = Alpha_Bosons_.end();
  for(; itAlpha_Boson != Alpha_Bosons_End; ++itAlpha_Boson)
    itAlpha_Boson->Display();
}//Close Display_Alpha_Bosons.

void AlphaBuilder::Display_Alpha_Fermions() const
{
  int Alpha_Fermions_Size = Alpha_Fermions().size();
  std::cout<<"Fermion Alphas: "<<Alpha_Fermions_Size<<std::endl;
  std::set<AlphaFermion>::iterator itAlpha_Fermion = Alpha_Fermions_.begin();
  std::set<AlphaFermion>::iterator Alpha_Fermions_End = Alpha_Fermions_.end();
  for(; itAlpha_Fermion != Alpha_Fermions_End; ++itAlpha_Fermion)
    itAlpha_Fermion->Display();
}//Close Display_Alpha_Fermions.

void AlphaBuilder::Display_Alpha_SUSYs() const
{
  std::cout<<"SUSY Alphas: "<<Alpha_SUSYs().size()<<std::endl;
  std::set<AlphaSUSY>::iterator itAlpha_SUSY = Alpha_SUSYs_.begin();
  for(; itAlpha_SUSY != Alpha_SUSYs_.end(); itAlpha_SUSY++)
    itAlpha_SUSY->Display();
}//Close Display_Alpha_SUSYs.

void AlphaBuilder::Display_All_Alphas() const
{
  Display_Alpha_Bosons();
  std::cout<<std::endl;
  Display_Alpha_Fermions();
  std::cout<<std::endl;
  Display_Alpha_SUSYs();
  std::cout<<std::endl;
}//Close Display_All_Alphas.

//PRIVATE.
//HELPERS.

std::vector<int> AlphaBuilder::Get_Coefficient_Limits
(const std::vector<BasisAlpha>& BasisAlphas)
{
  std::vector<int> Coefficient_Limits;
  int BasisAlphas_Size = BasisAlphas.size();
  for(int a=0; a<BasisAlphas_Size; ++a)
    Coefficient_Limits.push_back(BasisAlphas.at(a).Denominator());
  return Coefficient_Limits;
}//Close Get_Coefficient_Limits.

void AlphaBuilder::Add_Alphas(int Layer, Alpha Last_Alpha)
{
  if(Layer < static_cast<int>(Common_BasisAlphas().size()))
    {
      int Coefficient_Limit = Coefficient_Limits().at(Layer);
      for(int a=0; a<Coefficient_Limit; ++a)
	{
	  std::vector<int> New_Numerator = Last_Alpha.Numerator();
	  std::vector<int> New_Coefficients = Last_Alpha.Coefficients();

	  //Set the coefficients.
	  New_Coefficients.at(Layer) = a;

	  //Now add the basis alpha multiplied by the coefficient.
	  int New_Numerator_Size = New_Numerator.size();
	  for(int b=0; b<New_Numerator_Size; ++b)
	    New_Numerator.at(b) += 
	      a*Common_BasisAlphas().at(Layer).Numerator().at(b);

	  //Adjust the range.
	  New_Numerator = Adjust_Alpha_Range(New_Numerator);

	  //Now add the linear combination.
	  Alpha New_Alpha(New_Numerator, 
			  Common_BasisAlphas().at(Layer).Denominator(), 
			  New_Coefficients);
	  Add_Alphas(Layer+1, New_Alpha);
	}//Close for loop on Coefficient_Limits().at(Layer).
    }else//Close if statement.
    {
      int Mass_Denominator = Last_Alpha.Denominator()*Last_Alpha.Denominator();

      //Check for linear independence.
      if(!Linearly_Independent_Alpha(Last_Alpha))
	Linearly_Independent_Alphas_ = false;

      //Get the Masses of the left and right movers.
      int Mass_Left = Last_Alpha.Mass_Left();
      int Mass_Right = Last_Alpha.Mass_Right();

      //Check if the alpha is a SUSY sector.
      if(Mass_Left == 8*Mass_Denominator && 
	 Mass_Right == 0 && Last_Alpha.Numerator().at(0) != 0)
	//SUSY sector, also automatically a fermion sector.
	{
	  Alpha_SUSYs_.insert(AlphaSUSY(Last_Alpha.Numerator(), 
					 Last_Alpha.Denominator(),
					 Last_Alpha.Coefficients()));
	  Alpha_Fermions_.insert(AlphaFermion(Last_Alpha.Numerator(),
					       Last_Alpha.Denominator(),
					       Last_Alpha.Coefficients()));
	} else if(Mass_Left == 0 && 
		  Mass_Right <= 16*Mass_Denominator)//Boson sector.
	Alpha_Bosons_.insert(AlphaBoson(Last_Alpha.Numerator(),
					 Last_Alpha.Denominator(),
					 Last_Alpha.Coefficients()));
      else if(Mass_Left <= 8*Mass_Denominator && 
	      Mass_Right <= 16*Mass_Denominator &&
	      Last_Alpha.Numerator().at(0) != 0)//Fermion sector.
	Alpha_Fermions_.insert(AlphaFermion(Last_Alpha.Numerator(),
					     Last_Alpha.Denominator(),
					     Last_Alpha.Coefficients()));

    }//Close else statement.
}//Close Add_Alphas.

std::vector<int> AlphaBuilder::Adjust_Alpha_Range
(std::vector<int> Alpha_Numerator)
{
  int Denominator = Common_BasisAlphas().at(0).Denominator();
  //Because they're all the same anyway.
  int Alpha_Numerator_Size = Alpha_Numerator.size();
  for(int a=0; a<Alpha_Numerator_Size; ++a)
    {
      Alpha_Numerator.at(a)%=(2*Denominator);

      if(Alpha_Numerator.at(a) > Denominator)
	Alpha_Numerator.at(a) -= (2*Denominator);
      else if(Alpha_Numerator.at(a) <= -Denominator)
	Alpha_Numerator.at(a) += (2*Denominator);
    }//Close for loop on New_Numerator.
  return Alpha_Numerator;
}//Close Adjust_Alpha_Range.

bool AlphaBuilder::Linearly_Independent_Alpha(const Alpha& Last_Alpha)
{
  if(Last_Alpha.Numerator() != 
     std::vector<int>(Last_Alpha.Numerator().size(),0))
    return true;
  else if(Last_Alpha.Coefficients() == 
	  std::vector<int>(Last_Alpha.Coefficients().size(), 0))
    return true;
  else
    return false;
}
