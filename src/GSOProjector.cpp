#include <cassert>
#include "GSOProjector.h"

//PUBLIC.
GSOProjector::GSOProjector(const std::vector<BasisAlpha>& Common_Basis_Alphas,
			     char Alpha_Type,
			     const GSOCoefficientMatrix& k_ij,
			     const std::vector<int>& Coefficients,
			     const std::map<int, int>& Fermion_Mode_Map)
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Alpha_Type_ = Alpha_Type;
  k_ij_ = k_ij;
  Coefficients_ = Coefficients;
  Fermion_Mode_Map_ = Fermion_Mode_Map;
}//Close constructor.

GSOProjector::GSOProjector(const std::vector<BasisAlpha>& Common_Basis_Alphas,
			     const Alpha& The_Alpha, 
			     const GSOCoefficientMatrix& k_ij,
			     const std::map<int, int>& Fermion_Mode_Map)
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Alpha_Type_ = The_Alpha.Type();
  Coefficients_ = The_Alpha.Coefficients();
  k_ij_ = k_ij;
  Fermion_Mode_Map_ = Fermion_Mode_Map;

}//Close second constructor.

GSOProjector::GSOProjector(const GSOProjector& New_GSOProjector)
{
  Common_Basis_Alphas_ = New_GSOProjector.Common_Basis_Alphas();
  Alpha_Type_ = New_GSOProjector.Alpha_Type();
  k_ij_ = New_GSOProjector.k_ij();
  Coefficients_ = New_GSOProjector.Coefficients();
  Fermion_Mode_Map_ = New_GSOProjector.Fermion_Mode_Map();
}//Close copy constructor.


//INTERFACE.
bool GSOProjector::GSOP(const State& The_State) const
{
  if(Alpha_Type() == 'b')
    return GSOP_Boson(The_State);
  else if((Alpha_Type() == 'f') || (Alpha_Type() == 's'))
    return GSOP_Fermion(The_State);
  std::cout << "*** Error: unexpected edge case found." << std::endl;
  assert(false);
}//Close GSOP.

//PRIVATE.
//HELPERS.

bool GSOProjector::GSOP_Boson(const State& The_State) const
{
  bool GSO_Ready = true;
  std::vector<int> State_Numerator = The_State.Numerator();
  int State_Denominator = The_State.Denominator();
  int LM_Size = Common_Basis_Alphas().at(0).LM_Size();
  int Alpha_Denominator = Common_Basis_Alphas().at(0).Denominator();

  for(int a=0; a<static_cast<int>(Common_Basis_Alphas().size()); a++)
    {
      if(GSO_Ready == true)
	{
	  int RHS_Numerator = 0;
	  std::vector<int> Alpha_Numerator = 
	    Common_Basis_Alphas().at(a).Numerator();
	  //Firstly, the summation over k_ij.
	  for(int b=0; b<static_cast<int>(Coefficients().size()); b++)
	    RHS_Numerator += (k_ij().Numerators().at(a).at(b)*Coefficients().at(b)*
			      2*State_Denominator);

	  //Then the dot product between Q and Alpha.
	  for(int b=LM_Size; b<static_cast<int>(State_Numerator.size()); b++)
	    RHS_Numerator += (Alpha_Numerator.at(b)*State_Numerator.at(b));

	  //Finally, the correction factor for real elements.
	  std::map<int, int>::const_iterator itMap = Fermion_Mode_Map_.begin();
	  for(; itMap != Fermion_Mode_Map_.end(); ++itMap)
	    { 
	      if((itMap->first) < LM_Size && (itMap->second) > (LM_Size-1))
		{
		  RHS_Numerator += (Alpha_Numerator.at(itMap->second)*
				    State_Numerator.at(itMap->second));
		}//Close correction for real elements.
	    }//Close for loop on itMap.
	  

	  GSO_Ready = ((RHS_Numerator%(2*Alpha_Denominator*2*State_Denominator)) 
		       == 0);
	}else
	return GSO_Ready;
    }//Close for loop on Common_Basis_Alphas.

  return GSO_Ready;
}//Close GSOP_Boson.

bool GSOProjector::GSOP_Fermion(const State& The_State) const
{
  bool GSO_Ready = true;
  int LM_Size = Common_Basis_Alphas().at(0).LM_Size();
  int Alpha_Denominator = Common_Basis_Alphas().at(0).Denominator();
  std::vector<int> State_Numerator = The_State.Numerator();
  int State_Denominator = The_State.Denominator();

  int Common_Basis_Alphas_Size = Common_Basis_Alphas().size();
  for(int a=0; a<Common_Basis_Alphas_Size; ++a)
    {
      if(GSO_Ready == true)
	{
	  std::vector<int> Alpha_Numerator = 
	    Common_Basis_Alphas().at(a).Numerator();
	  int RHS_Numerator = 0;
	  int State_Numerator_Size = State_Numerator.size();
	  //First the Lorentz dot product between the state and the alpha.
	  for(int b=0; b<LM_Size; ++b)
	    RHS_Numerator += (State_Numerator.at(b)*Alpha_Numerator.at(b));
	  for(int b=LM_Size; b<State_Numerator_Size;++b)
	    RHS_Numerator -= (State_Numerator.at(b)*Alpha_Numerator.at(b));

	  //Then the correction for the real part.
	  std::map<int, int>::const_iterator itMap = Fermion_Mode_Map_.begin();
	  std::map<int, int>::const_iterator itMap_End = Fermion_Mode_Map_.end();
	  for(; itMap != itMap_End; ++itMap)
	    {
	      if((itMap->first)<LM_Size && (itMap->second)>(LM_Size-1))
		{
		  RHS_Numerator += (State_Numerator.at(itMap->first)*
				    Alpha_Numerator.at(itMap->first));
		  RHS_Numerator -= (State_Numerator.at(itMap->second)*
				    Alpha_Numerator.at(itMap->second));
		}//Close if statement for correcting the real part.
	    }//Close for loop on the map.

	  //Now perform the summation over k_ij.
	  int Coefficients_Size = Coefficients().size();
	  for(int b=0; b<Coefficients_Size; ++b)
	    RHS_Numerator -= (k_ij().Numerators().at(a).at(b)*Coefficients().at(b)*
			      2*State_Denominator);
	  //Now subtract the ST part.
	  RHS_Numerator -= (Alpha_Numerator.at(0)*
			    2*The_State.Denominator());

	  GSO_Ready = (RHS_Numerator%(2*2*Alpha_Denominator*
				      The_State.Denominator()) == 0);
			      
	}else
	return GSO_Ready;
    }//Close for loop on Common_Basis_Alphas.
  return GSO_Ready;
}//Close GSOP_Fermion.

